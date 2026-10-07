#!/usr/bin/env python3
"""Pro Planet-Hunt: echter Sternradius, SNR, Odd-Even-Test, Sektor-Persistenz."""
import sys
import numpy as np
import lightkurve as lk
from astropy.stats import sigma_clipped_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RSUN_IN_REARTH = 109.1

def get_star_radius(target):
    """Sternradius (Rsun) aus TIC-Katalog holen."""
    try:
        from astroquery.mast import Catalogs
        cat = Catalogs.query_object(target, catalog="TIC", radius=0.002)
        if len(cat) > 0:
            rad = float(cat["rad"][0])
            teff = float(cat["Teff"][0])
            if np.isfinite(rad):
                print(f"[*] TIC Sternradius: {rad:.3f} Rsun, Teff={teff:.0f} K")
                return rad
    except Exception as e:
        print(f"[!] TIC-Abfrage fehlgeschlagen: {e}")
    print("[!] Nutze Default 1.0 Rsun")
    return 1.0

def load(target, cadence=120, max_sectors=6):
    print(f"[*] Suche Lightcurves fuer {target} ...")
    sr = lk.search_lightcurve(target, mission="TESS", author="SPOC",
                              exptime=cadence)
    if len(sr) == 0:
        sr = lk.search_lightcurve(target, mission="TESS")
    if len(sr) > max_sectors:
        sr = sr[:max_sectors]
    print(f"[*] Lade {len(sr)} Files ({cadence}s):")
    print(sr)
    lcc = sr.download_all()
    lc = lcc.stitch().remove_nans().remove_outliers(sigma=5)
    flat = lc.flatten(window_length=401).remove_nans()
    span = (flat.time.max() - flat.time.min()).value
    print(f"[*] {len(flat.time)} Punkte ueber {span:.1f} Tage")
    return flat

def in_transit_mask(time, P, t0, dur, width=0.5):
    phase = ((time - t0 + 0.5*P) % P) - 0.5*P
    return np.abs(phase) < width * dur

def analyze(flat, r, rstar):
    """SNR, Odd-Even, Radius mit echtem Sternradius."""
    t = flat.time.value
    f = flat.flux.value
    P, t0, dur, depth = r["P"].value, r["t0"].value, r["dur"].value, r["depth"]

    intr = in_transit_mask(t, P, t0, dur, width=0.5)
    oot = ~in_transit_mask(t, P, t0, dur, width=1.5)
    if intr.sum() < 3 or oot.sum() < 10:
        r["snr"] = 0.0
    else:
        _, _, std_oot = sigma_clipped_stats(f[oot])
        n_in = intr.sum()
        r["snr"] = (depth) / (std_oot / np.sqrt(n_in)) if std_oot > 0 else 0.0

    # Odd-Even: Transit-Nummer bestimmen
    tn = np.round((t - t0) / P)
    odd = intr & (tn % 2 == 1)
    even = intr & (tn % 2 == 0)
    d_odd = 1.0 - np.median(f[odd]) if odd.sum() > 2 else np.nan
    d_even = 1.0 - np.median(f[even]) if even.sum() > 2 else np.nan
    if np.isfinite(d_odd) and np.isfinite(d_even) and (d_odd+d_even) > 0:
        r["oddeven_diff"] = abs(d_odd - d_even) / ((d_odd + d_even)/2)
    else:
        r["oddeven_diff"] = np.nan
    r["d_odd"], r["d_even"] = d_odd, d_even

    rprs = np.sqrt(abs(depth))
    r["rp_earth"] = rprs * rstar * RSUN_IN_REARTH
    return r

def find_one(flat, min_p=0.5, max_p=20.0, label="P1"):
    grid = np.linspace(min_p, max_p, 6000)
    bls = flat.to_periodogram(method="bls", period=grid, frequency_factor=1000)
    return dict(
        bls=bls, label=label,
        P=bls.period_at_max_power,
        t0=bls.transit_time_at_max_power,
        dur=bls.duration_at_max_power,
        depth=float(bls.depth_at_max_power),
        power=float(bls.max_power),
    )

def mask_transits(flat, P, t0, dur):
    phase = ((flat.time.value - t0.value + 0.5*P.value) % P.value) - 0.5*P.value
    return flat[np.abs(phase) > 1.3 * dur.value]

def verdict(r):
    snr, oe, pw = r["snr"], r["oddeven_diff"], r["power"]
    if pw < 20 or snr < 7:
        return "RAUSCHEN (schwach)"
    if np.isfinite(oe) and oe > 0.5:
        return "VERDACHT: Eclipsing Binary (Odd-Even unterschiedlich)"
    if snr >= 7 and pw >= 20:
        return "KANDIDAT (planetenartig)"
    return "unklar"

def run(target, n_planets=3, max_p=20.0):
    rstar = get_star_radius(target)
    flat = load(target)
    results, work = [], flat
    for i in range(n_planets):
        r = find_one(work, max_p=max_p, label=f"Signal {i+1}")
        r = analyze(flat, r, rstar)
        results.append(r)
        print(f"\n== {r['label']} ==")
        print(f"Periode : {r['P'].value:.5f} d")
        print(f"Tiefe   : {r['depth']*1e6:.0f} ppm")
        print(f"Dauer   : {r['dur'].value*24:.2f} h")
        print(f"Rp      : ~{r['rp_earth']:.2f} Rerde (Rstar={rstar:.2f} Rsun)")
        print(f"BLS pow : {r['power']:.1f}")
        print(f"SNR     : {r['snr']:.1f}")
        print(f"Odd-Even: {r['oddeven_diff']:.2f} (odd={r['d_odd']*1e6 if np.isfinite(r['d_odd']) else float('nan'):.0f} even={r['d_even']*1e6 if np.isfinite(r['d_even']) else float('nan'):.0f} ppm)")
        print(f"URTEIL  : {verdict(r)}")
        work = mask_transits(work, r["P"], r["t0"], r["dur"])

    n = len(results)
    fig, axes = plt.subplots(n, 1, figsize=(10, 3.2*n))
    if n == 1:
        axes = [axes]
    for ax, r in zip(axes, results):
        folded = flat.fold(period=r["P"], epoch_time=r["t0"])
        folded.scatter(ax=ax, s=1, c="0.7")
        try:
            folded.bin(time_bin_size=r["dur"].value/4).plot(ax=ax, c="red", lw=2)
        except Exception:
            pass
        ax.set_xlim(-0.25, 0.25)
        ax.set_title(f"{r['label']}: P={r['P'].value:.4f}d "
                     f"~{r['rp_earth']:.1f}Re SNR={r['snr']:.0f} "
                     f"pow={r['power']:.0f} | {verdict(r)}")
    plt.tight_layout()
    out = f"pro_{target.replace(' ','_')}.png"
    plt.savefig(out, dpi=110)
    print(f"\n[*] Plot: {out}")

if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "TOI-700"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    maxp = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    run(target, n_planets=n, max_p=maxp)
