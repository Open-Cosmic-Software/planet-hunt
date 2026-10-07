#!/usr/bin/env python3
"""Multi-Planet-Suche: BLS -> maskieren -> wiederholen."""
import sys
import numpy as np
import lightkurve as lk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

def find_one(flat, min_p=0.5, max_p=20.0, label="P1"):
    grid = np.linspace(min_p, max_p, 6000)
    bls = flat.to_periodogram(method="bls", period=grid, frequency_factor=1000)
    P = bls.period_at_max_power
    t0 = bls.transit_time_at_max_power
    dur = bls.duration_at_max_power
    depth = float(bls.depth_at_max_power)
    power = float(bls.max_power)
    rprs = np.sqrt(abs(depth))
    rp = rprs * 109.0
    print(f"\n== {label} ==")
    print(f"Periode : {P.value:.5f} d")
    print(f"t0      : {t0.value:.4f}")
    print(f"Dauer   : {dur.value*24:.2f} h")
    print(f"Tiefe   : {depth*1e6:.0f} ppm")
    print(f"Rp/Rs   : {rprs:.4f}  -> ~{rp:.2f} Rerde (bei 1 Rsun)")
    print(f"BLS power: {power:.1f}")
    return dict(bls=bls, P=P, t0=t0, dur=dur, depth=depth, power=power,
                rp=rp, label=label)

def mask_transits(flat, P, t0, dur):
    # Punkte innerhalb 1.3x Transitdauer um jeden erwarteten Transit raus
    phase = ((flat.time.value - t0.value + 0.5*P.value) % P.value) - 0.5*P.value
    keep = np.abs(phase) > 1.3 * dur.value
    return flat[keep]

def run(target, n_planets=3, max_p=20.0):
    flat = load(target)
    results = []
    work = flat
    for i in range(n_planets):
        r = find_one(work, max_p=max_p, label=f"Signal {i+1}")
        results.append(r)
        work = mask_transits(work, r["P"], r["t0"], r["dur"])
    # Plot: fuer jedes Signal ein Phase-Fold
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
        ax.set_title(f"{r['label']}: P={r['P'].value:.4f}d, "
                     f"{r['depth']*1e6:.0f}ppm, ~{r['rp']:.1f}Re, "
                     f"pow={r['power']:.0f}")
    plt.tight_layout()
    out = f"multi_{target.replace(' ','_')}.png"
    plt.savefig(out, dpi=110)
    print(f"\n[*] Plot: {out}")
    print("\n==== ZUSAMMENFASSUNG ====")
    for r in results:
        print(f"{r['label']}: P={r['P'].value:.4f}d  "
              f"depth={r['depth']*1e6:.0f}ppm  ~{r['rp']:.1f}Re  "
              f"power={r['power']:.0f}")

if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "TOI-700"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    maxp = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    run(target, n_planets=n, max_p=maxp)
