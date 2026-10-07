#!/usr/bin/env python3
"""Pro Planet-Hunt v2 — mit Tytos Verbesserungen.

- Transit-maskiertes Flattening (Radius wird nicht kleingeglaettet)
- Echter Sternradius aus TIC
- Sekundaer-Eclipse-Check (Phase 0.5)
- V/U-Form-Heuristik
- Harmonischen-Check (P/2, 2P)
- Odd-Even-Test
- Automatischer Abgleich mit TOI/bestaetigten Katalogen
"""
import sys
import numpy as np
import lightkurve as lk
from astropy.stats import sigma_clipped_stats
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RSUN_IN_REARTH = 109.1
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

def tap_csv(query):
    r = requests.get(TAP, params={"query": query, "format": "csv"}, timeout=90)
    lines = r.text.strip().splitlines()
    if len(lines) < 2:
        return []
    hdr = lines[0].split(",")
    out = []
    for ln in lines[1:]:
        vals = [v.strip().strip('"') for v in ln.split(",")]
        out.append(dict(zip(hdr, vals)))
    return out

def get_star(target, tic_id=None):
    try:
        from astroquery.mast import Catalogs
        q = f"TIC {tic_id}" if tic_id else target
        cat = Catalogs.query_object(q, catalog="TIC", radius=0.002)
        if len(cat) > 0:
            rad = float(cat["rad"][0]); teff = float(cat["Teff"][0])
            tic = int(cat["ID"][0])
            if np.isfinite(rad):
                print(f"[*] TIC {tic}: R*={rad:.3f} Rsun, Teff={teff:.0f} K")
                return rad, tic
    except Exception as e:
        print(f"[!] TIC query failed: {e}")
    return 1.0, tic_id

def load_masked(target, cadence=120, max_sectors=12):
    print(f"[*] Suche Lightcurves fuer {target} ...")
    sr = lk.search_lightcurve(target, mission="TESS", author="SPOC", exptime=cadence)
    if len(sr) == 0:
        sr = lk.search_lightcurve(target, mission="TESS")
    if len(sr) > max_sectors:
        sr = sr[:max_sectors]
    print(f"[*] Lade {len(sr)} Files ({cadence}s)")
    lcc = sr.download_all()
    lc = lcc.stitch().remove_nans().remove_outliers(sigma=5)
    # Neu aufbauen aus reinen numpy-Arrays -> keine masked-array-Konflikte in flatten
    tt = np.asarray(lc.time.value, dtype=float)
    ff = np.asarray(lc.flux, dtype=float)
    good = np.isfinite(tt) & np.isfinite(ff)
    lc = lk.LightCurve(time=tt[good], flux=ff[good])

    # 1. Durchgang: grobes flatten nur um Transits zu FINDEN
    rough = lc.flatten(window_length=401).remove_nans()
    bls = rough.to_periodogram(method="bls",
                               period=np.linspace(0.5, 20, 6000),
                               frequency_factor=1000)
    P = bls.period_at_max_power.value
    t0 = bls.transit_time_at_max_power.value
    dur = bls.duration_at_max_power.value

    # 2. maskiertes flatten (eigene Implementierung, da lk.flatten(mask=) buggy):
    #    Transit-Punkte beim Trend ueberbruecken, savgol-Trend, dann teilen.
    from scipy.signal import savgol_filter
    tt = np.asarray(lc.time.value, dtype=float)
    ff = np.asarray(lc.flux, dtype=float)
    phase = ((tt - t0 + 0.5*P) % P) - 0.5*P
    in_tr = np.abs(phase) < 1.0*dur
    ff_fit = ff.copy()
    if in_tr.any() and (~in_tr).sum() > 10:
        ff_fit[in_tr] = np.interp(tt[in_tr], tt[~in_tr], ff[~in_tr])
    win = 401 if len(ff_fit) > 401 else (len(ff_fit)//2*2-1)
    trend = savgol_filter(ff_fit, win, 2)
    flat = lk.LightCurve(time=tt, flux=ff/trend).remove_nans()
    span = (flat.time.max() - flat.time.min()).value
    print(f"[*] {len(flat.time)} Punkte ueber {span:.1f} Tage "
          f"(maskiertes Flattening, grob P={P:.3f}d)")
    return flat

def diagnostics(flat, P, t0, dur, rstar):
    t = flat.time.value; f = flat.flux.value
    phase = ((t - t0 + 0.5*P) % P) - 0.5*P
    intr = np.abs(phase) < 0.5*dur
    oot  = np.abs(phase) > 1.5*dur
    _, med_oot, std_oot = sigma_clipped_stats(f[oot])
    depth = med_oot - np.median(f[intr])
    snr = depth/(std_oot/np.sqrt(max(intr.sum(),1))) if std_oot>0 else 0.0
    rprs = np.sqrt(abs(depth))
    rp = rprs * rstar * RSUN_IN_REARTH

    # Odd-Even
    tn = np.round((t-t0)/P)
    odd = intr & (tn%2==1); even = intr & (tn%2==0)
    d_odd = med_oot-np.median(f[odd]) if odd.sum()>2 else np.nan
    d_even= med_oot-np.median(f[even]) if even.sum()>2 else np.nan
    oe = abs(d_odd-d_even)/((d_odd+d_even)/2) if (np.isfinite(d_odd) and np.isfinite(d_even) and d_odd+d_even>0) else np.nan

    # Sekundaer-Eclipse bei Phase 0.5
    sec = np.abs(np.abs(phase)-0.5*P) < 0.5*dur
    if sec.sum()>3:
        d_sec = med_oot-np.median(f[sec])
        snr_sec = d_sec/(std_oot/np.sqrt(sec.sum())) if std_oot>0 else 0.0
    else:
        d_sec, snr_sec = np.nan, np.nan

    # V/U-Form: Verhaeltnis Tiefe im Kern (25%) vs. ganze Transitbreite
    core = np.abs(phase) < 0.25*dur
    wing = (np.abs(phase)>=0.25*dur)&(np.abs(phase)<0.5*dur)
    if core.sum()>2 and wing.sum()>2:
        d_core = med_oot-np.median(f[core])
        d_wing = med_oot-np.median(f[wing])
        vshape = d_core/d_wing if d_wing>0 else np.nan  # ~1=U (flach), >>1=V spitz
    else:
        vshape = np.nan

    # Harmonischen-Check
    def depth_at(per):
        ph=((t-t0+0.5*per)%per)-0.5*per
        ii=np.abs(ph)<0.5*dur; oo=np.abs(ph)>1.5*dur
        if ii.sum()<3 or oo.sum()<10: return 0.0
        _,m,s=sigma_clipped_stats(f[oo])
        return (m-np.median(f[ii]))/(s/np.sqrt(ii.sum())) if s>0 else 0.0
    snr_half, snr_double = depth_at(P/2), depth_at(P*2)

    return dict(depth=depth, snr=snr, rp=rp, rprs=rprs, oe=oe,
                d_sec=d_sec, snr_sec=snr_sec, vshape=vshape,
                snr_half=snr_half, snr_double=snr_double,
                d_odd=d_odd, d_even=d_even)

def catalog_match(tic, P):
    """Pruefe ob Periode schon als TOI oder bestaetigt bekannt ist."""
    hits = []
    for row in tap_csv(f"select toi,tfopwg_disp,pl_orbper from toi where tid={tic}"):
        try:
            p = float(row["pl_orbper"])
            if p>0 and abs(p-P)/P < 0.02:
                hits.append(f"TOI-{row['toi']} ({row['tfopwg_disp']}, P={p:.4f}d)")
        except: pass
    return hits

def verdict(d):
    if d["snr"] < 7 or abs(d["snr"]) < 7:
        return "RAUSCHEN (SNR<7)"
    if np.isfinite(d["snr_sec"]) and d["snr_sec"] > 5:
        return "VERDACHT: Doppelstern (Sekundaer-Eclipse)"
    if np.isfinite(d["oe"]) and d["oe"] > 0.5:
        return "VERDACHT: Doppelstern (Odd-Even)"
    if np.isfinite(d["vshape"]) and d["vshape"] > 1.8:
        return "VERDACHT: V-foermig (evtl. Doppelstern)"
    if d["snr_double"] > d["snr"]*1.1:
        return "UNKLAR: 2P staerker (Harmonische? echte Periode evtl. 2x)"
    return "KANDIDAT (planetenartig)"

def run(target, max_p=20.0, tic_id=None):
    rstar, tic = get_star(target, tic_id)
    flat = load_masked(target)
    bls = flat.to_periodogram(method="bls",
                              period=np.linspace(0.5, max_p, 6000),
                              frequency_factor=1000)
    P = bls.period_at_max_power.value
    t0 = bls.transit_time_at_max_power.value
    dur = bls.duration_at_max_power.value
    d = diagnostics(flat, P, t0, dur, rstar)
    hits = catalog_match(tic, P) if tic else []

    print(f"\n==== ERGEBNIS ====")
    print(f"Periode    : {P:.5f} d")
    print(f"Tiefe      : {d['depth']*1e6:.0f} ppm")
    print(f"Dauer      : {dur*24:.2f} h")
    print(f"Rp         : ~{d['rp']:.2f} Rerde (R*={rstar:.2f} Rsun)")
    print(f"SNR        : {d['snr']:.1f}")
    print(f"Odd-Even   : {d['oe']:.2f}" if np.isfinite(d['oe']) else "Odd-Even   : n/a")
    print(f"Sek.Eclipse: SNR={d['snr_sec']:.1f}" if np.isfinite(d['snr_sec']) else "Sek.Eclipse: n/a")
    print(f"V/U-Form   : {d['vshape']:.2f} (~1=U/Planet, >1.8=V/Binary)" if np.isfinite(d['vshape']) else "V/U-Form   : n/a")
    print(f"Harmonische: SNR(P/2)={d['snr_half']:.1f}  SNR(2P)={d['snr_double']:.1f}")
    print(f"Katalog    : {', '.join(hits) if hits else 'KEIN bekannter TOI bei dieser Periode'}")
    print(f"URTEIL     : {verdict(d)}")

    fig, ax = plt.subplots(figsize=(9,5))
    folded = flat.fold(period=P, epoch_time=t0)
    folded.scatter(ax=ax, s=2, c="0.6")
    try: folded.bin(time_bin_size=dur/4).plot(ax=ax, c="red", lw=2)
    except: pass
    ax.set_xlim(-0.2, 0.2)
    ax.set_title(f"{target}: P={P:.4f}d ~{d['rp']:.1f}Re SNR={d['snr']:.0f} | {verdict(d)}")
    out=f"v2_{target.replace(' ','_')}.png"
    plt.tight_layout(); plt.savefig(out, dpi=110)
    print(f"[*] Plot: {out}")

if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv)>1 else "Pi Mensae"
    maxp = float(sys.argv[2]) if len(sys.argv)>2 else 20.0
    tic = sys.argv[3] if len(sys.argv)>3 else None
    run(target, max_p=maxp, tic_id=tic)
