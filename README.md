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

We ran it blind against known systems to prove it works:

| Target        | Our period | Literature | Our radius | Verdict     |
|---------------|-----------:|-----------:|-----------:|-------------|
| Pi Mensae c   |  6.26765 d |  6.2679 d  |  1.47 R⊕   | matches     |
| TOI-700 d     | 16.05059 d | 16.0510 d  |  1.97 R⊕   | CANDIDATE   |
| TOI-1169.01   |  6.70720 d |  6.7075 d  |  sub-Saturn| CANDIDATE   |

For TOI-1169.01 (an *unconfirmed* TESS Planet Candidate) the pipeline
independently recovered the dip at SNR ≈ 99, odd-even difference 0.04
(no binary), persistent across 5 sectors spanning 949 days (2021–2024).

## Honesty notes

- Low BLS power + low SNR = **noise**, not a planet. The tool says so.
- A verdict of `CANDIDATE` means *planet-like*, not *confirmed*. Real
  confirmation needs radial velocity or ruling out blended eclipsing binaries.
- Radius defaults to a 1 R☉ assumption when the TIC has no stellar radius.

Built with 🦞 by Nyx & Fabian · MIT License
