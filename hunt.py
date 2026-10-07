#!/usr/bin/env python3
"""Planet-Hunt: TESS Lightcurve -> Detrend -> BLS -> Transit-Kandidat."""
import sys
import numpy as np
import lightkurve as lk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def hunt(target, max_period=15.0, min_period=0.5, outprefix="out",
         cadence=120, max_sectors=4):
    print(f"[*] Suche Lightcurves fuer {target} ...")
    sr = lk.search_lightcurve(target, mission="TESS", author="SPOC",
                              exptime=cadence)
    if len(sr) == 0:
        print("[!] Keine SPOC-Lightcurves gefunden, versuche QLP ...")
        sr = lk.search_lightcurve(target, mission="TESS")
    # Nur die ersten max_sectors nehmen (ueber Jahre verteilt reicht)
    if len(sr) > max_sectors:
        sr = sr[:max_sectors]
    print(f"[*] Lade {len(sr)} Lightcurve-Files ({cadence}s cadence):")
    print(sr)

    lcc = sr.download_all()
    # Stitch + bereinigen + normalisieren
    lc = lcc.stitch().remove_nans().remove_outliers(sigma=5)
    # Detrending (flatten) gegen Sternrauschen/Trends
    flat = lc.flatten(window_length=401)

    print(f"[*] {len(flat.time)} Datenpunkte ueber "
          f"{(flat.time.max()-flat.time.min()).value:.1f} Tage")

    # Box-Least-Squares Periodensuche
    period_grid = np.linspace(min_period, max_period, 5000)
    bls = flat.to_periodogram(method="bls", period=period_grid,
                              frequency_factor=1000)
    best_period = bls.period_at_max_power
    best_t0 = bls.transit_time_at_max_power
    best_dur = bls.duration_at_max_power
    depth = bls.depth_at_max_power

    print("\n==== BLS ERGEBNIS ====")
    print(f"Beste Periode : {best_period.value:.5f} d")
    print(f"Transit t0    : {best_t0.value:.5f}")
    print(f"Dauer         : {best_dur.value*24:.3f} h")
    print(f"Tiefe         : {float(depth)*1e6:.0f} ppm "
          f"({float(depth)*100:.4f} %)")

    # Radius-Schaetzung (relativ zum Stern): Rp/Rs = sqrt(depth)
    rprs = np.sqrt(float(depth))
    # Grobe Annahme Sternradius ~ 1 Rsun -> Rp in Erdradien (Rsun=109 Rerde)
    rp_earth = rprs * 109.0
    print(f"Rp/Rs         : {rprs:.4f}")
    print(f"Rp (bei 1 Rsun): ~{rp_earth:.2f} Erdradien")

    # Plots
    fig, axes = plt.subplots(3, 1, figsize=(11, 12))
    bls.plot(ax=axes[0])
    axes[0].set_title(f"{target} - BLS Periodogram (P={best_period:.4f})")
    flat.scatter(ax=axes[1], s=1)
    axes[1].set_title("Detrended Lightcurve")
    folded = flat.fold(period=best_period, epoch_time=best_t0)
    folded.scatter(ax=axes[2], s=1)
    folded.bin(time_bin_size=0.005).plot(ax=axes[2], color="red", lw=2)
    axes[2].set_title(f"Phase-Folded @ {best_period:.4f}")
    axes[2].set_xlim(-0.3, 0.3)
    plt.tight_layout()
    out = f"{outprefix}_{target.replace(' ','_')}.png"
    plt.savefig(out, dpi=110)
    print(f"[*] Plot gespeichert: {out}")
    return bls

if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "Pi Mensae"
    maxp = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0
    hunt(target, max_period=maxp)
