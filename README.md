# Water Resources Engineering

Coursework for CIVE70020 Water Resources Engineering at Imperial College London, analysing catchment 39008 near Oxford. The notebook covers catchment water balance, crop water demand, rainfall–runoff model calibration, groundwater comparison and reservoir storage estimation using the Waitt curve.

The accompanying [report](report/report.pdf) presents the methods, results and discussion.

```text
water-resources/
├── MajdYousofCWSubmission.ipynb  # Analysis, figures and results
├── model.py                     # Hydrological model and calibration functions
├── report/
│   └── report.pdf               # Coursework report
├── data/
│   ├── Catchment/               # CAMELS-GB catchment attributes
│   │   └── shapefiles/
│   │       ├── 39008_boundary/  # Catchment boundary
│   │       └── gb/              # Great Britain boundaries
│   └── Hydrology/               # Weather, river flow, soil moisture and groundwater
├── outputs/                     # Generated notebooks and plots (not tracked)
├── pyproject.toml               # Dependencies and code-check configuration
├── uv.lock                      # Locked dependency versions
├── .python-version              # Python version
├── .gitignore                   # Files excluded from Git
├── .gitattributes               # File-handling rules
└── README.md
```
