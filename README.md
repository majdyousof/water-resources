# Catchment water balance and simulation

Water Resources Engineering coursework by Majd Yousof, Imperial College London, December 2024.

The analysis follows catchment 39008 near Oxford from rainfall and crop water demand to river discharge, groundwater conditions and reservoir storage. It uses the supplied 1970–2015 hydrometeorological record, soil-moisture observations and groundwater data.

Open [the notebook](MajdYousofCWSubmission.ipynb) for the calculations and figures, or [the report](report.pdf) for the discussion.

## Run

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run these commands from this directory:

```bash
uv sync --locked
uv run jupyter lab MajdYousofCWSubmission.ipynb
```

In Jupyter, choose **Restart Kernel and Run All Cells**. Run from the project root so `data/` and `model.py` are available. The first simulation compiles the model with Numba; the full Nelder–Mead calibration can take several minutes. Figures are regenerated in `outputs/plots/`.

To execute without opening Jupyter:

```bash
uv run jupyter nbconvert --execute --to notebook MajdYousofCWSubmission.ipynb \
  --output MajdYousofCWSubmission.executed.ipynb --output-dir outputs \
  --ExecutePreprocessor.timeout=600
```

Python 3.12 and all dependencies are recorded in `.python-version` and `uv.lock`.

## Analysis

1. Estimate crop coefficients and monthly water balance with and without soil-water stress.
2. Calculate the annual precipitation–crop-demand balance.
3. Calibrate the supplied lumped hydrological model using Nelder–Mead and a chronological 80:20 split.
4. Compare observed groundwater with simulated storage using standardisation and the Standardised Groundwater Index (SGI).
5. Calculate reservoir storage using the Waitt curve for 150,000 people, demand of 142 litres/person/day and a 2% river abstraction allowance.

## Results

| Metric | Value |
| --- | ---: |
| Uncalibrated discharge sum of squared errors | 13.193142 |
| Calibrated discharge sum of squared errors | 2.600952 |
| Test-record discharge sum of squared errors | 0.855774 |
| Groundwater SGI sum of squared errors | 286.775993 |
| Reservoir capacity | 11,264,480.64 m³ |
| Critical-period index | 55 months |

## Limitations

- Monthly soil-water stress factors cycle every 12 daily rows rather than matching calendar months. The 365-day crop coefficient sequence also starts in October and does not account for leap days.
- The annual balance is precipitation minus crop PET. Positive values represent a surplus rather than irrigation demand. The first and last calendar years contain only part of a year.
- Calibration permits porosities above one and changes the time step despite daily observations. Spin-up treatment differs between the objective and reported errors; the test simulation restarts from the original initial states. Unnormalised sums of squares over unequal record lengths are not directly comparable.
- The Waitt calculation has a one-month mismatch between inflow-window lengths and cumulative demand. Its reported 55 months is a critical-period index, not a return period. The report also labels the daily demand of 21,300 m³ as monthly; the code aggregates daily demand correctly.

These limitations affect the interpretation of irrigation requirements and reservoir capacity.

## Files and sources

- `MajdYousofCWSubmission.ipynb`: coursework analysis, figures and results.
- `model.py`: the lumped hydrological model supplied for the coursework (`Model_B.py`).
- `data/`: the supplied CSVs and shapefiles.
- `report.pdf`: coursework report.
- `outputs/`: generated execution copies and figures, excluded from Git.

Hydrometeorological and catchment data are from **CAMELS-GB**, described by [Coxon et al. (2020)](https://doi.org/10.5194/essd-12-2459-2020). The coursework also supplied COSMOS soil-moisture data, groundwater observations and boundary files. Further methodological references are in the report and notebook. Third-party data and supplied code retain their original attribution and terms.

Code checks:

```bash
uv run ruff check .
uv run ruff format --check .
uv run pyright
```
