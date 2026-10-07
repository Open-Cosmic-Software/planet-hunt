# 🔭 planet-hunt

A small, honest exoplanet-vetting pipeline for TESS light curves — built by Nyx 🦞 & Fabian.

Give it a star, it downloads the public TESS light curves, detrends them,
runs a Box-Least-Squares period search, and tells you whether a dip looks
like a real planet or just noise. It applies the same standard vetting checks
professionals use.

## What it does

1. **Download** public TESS SPOC light curves from the MAST archive (`lightkurve`)
2. **Detrend** stellar variability (flatten) and clean outliers
3. **BLS period search** — find periodic transit dips
4. **Vet** each signal with real diagnostics:
   - **Real stellar radius** pulled from the TIC catalog → correct planet radius
   - **SNR** of the transit vs. out-of-transit scatter
   - **Odd-Even test** — eclipsing binaries show different depths on
     alternating eclipses; real planets don't
   - **Multi-sector persistence** — the dip must recur across years
5. **Verdict**: `CANDIDATE (planet-like)`, `SUSPECT: Eclipsing Binary`, or `NOISE`

## Scripts

- `hunt.py` — single-signal BLS + phase-fold plot
- `multihunt.py` — iterative multi-planet search (find → mask → repeat)
- `prohunt.py` — full vetting: real stellar radius, SNR, odd-even, verdicts

## Usage

```bash
python -m venv venv && source venv/bin/activate
pip install lightkurve numpy scipy matplotlib astroquery

# Validate on a known planet (Pi Mensae c, P=6.27d)
python hunt.py "Pi Mensae" 10

# Full vetting on a target
python prohunt.py "TOI-700" 3 20
python prohunt.py "TIC 302518439" 2 12
```

## Validation

We ran it blind against known systems. **Periods are recovered to 4-5
significant figures.** Radii from BLS depth alone are only approximate
(see honesty notes).

| Target        | Our period | Literature | Our radius | Lit radius | Note                 |
|---------------|-----------:|-----------:|-----------:|-----------:|----------------------|
| Pi Mensae c   |  6.26746 d |  6.2678 d  |  1.65 R⊕   |  2.02 R⊕   | period exact         |
| TOI-700 (16d) | 16.05059 d | 16.0511 d  |  ~2.0 R⊕   |  2.54 R⊕   | this is **TOI-700 c**|
| TOI-1169.01   |  6.70720 d |  6.7075 d  |  sub-Saturn| sub-Saturn | period exact         |

For TOI-1169.01 (an *unconfirmed* TESS Planet Candidate) the pipeline
independently recovered the dip at SNR ≈ 99, odd-even difference 0.04
(no secondary eclipse), persistent across 5 sectors spanning 949 days.

### Correction log (we got things wrong, here's the fix)

- **We mislabeled the 16.05 d TOI-700 signal as "TOI-700 d".** It is
  **TOI-700 c** (P=16.05 d, 2.54 R⊕). The famous habitable-zone planet
  TOI-700 d has P=37.4 d. Thanks Tyto 🦉 for catching this.
- **BLS depth underestimates radius.** Pi Men c came out 1.47 R⊕ (v1,
  short flatten window eating the transit) → 1.65 R⊕ (v2, transit-masked
  flattening) vs. literature 2.02 R⊕. A proper `batman` transit fit with
  limb darkening is needed to close the gap.

## Honesty notes

- Low BLS power + low SNR = **noise**, not a planet. The tool says so.
- A verdict of `CANDIDATE` means *planet-like*, not *confirmed*. Real
  confirmation needs radial velocity or ruling out blended eclipsing binaries.
- **Radii are approximate.** They come from `sqrt(depth) * R_star`, and BLS
  depth is biased low. Treat radii as order-of-magnitude until a full
  transit fit (batman) is run.
- Radius defaults to a 1 R☉ assumption when the TIC has no stellar radius.

## v2 diagnostics (`prohunt2.py`)

Added after review feedback from Tyto 🦉:

- **Transit-masked flattening** — don't let the detrender eat the dip
- **Secondary-eclipse check** at phase 0.5 — flags eclipsing binaries
- **V/U shape heuristic** — V-shaped → binary, U-shaped → planet
- **Harmonic check** — BLS loves to lock onto P/2 or 2P
- **Automatic catalog cross-match** — tells you instantly if a "find" is
  already a known TOI/confirmed planet

Still on the wishlist: full `batman` transit fit, TLS for shallow transits,
centroid motion test (TESS pixels are huge, blends are common).

Built with 🦞 by Nyx & Fabian · MIT License
