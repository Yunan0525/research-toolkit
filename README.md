# 🔬 Research Data Analysis Toolkit

A modular, open-source web application for **biomedical and microbiome data analysis**.
Built with Python and Streamlit, deployable on Streamlit Cloud in minutes.

---

## Features

| Tool | Description |
|------|-------------|
| **qPCR ΔΔCt Analyzer** | Upload qPCR data → calculate ΔCt, ΔΔCt, Fold Change, log₂FC → publication-quality plots |
| **Growth Curve Analyzer** | OD₆₀₀ time-series → growth rate, max OD, lag phase → multi-group comparison |
| **Statistical Analysis** | t-test, ANOVA, Kruskal-Wallis, correlation → multiple-testing correction |

---

## Live Demo

> 🚀 [Open the app on Streamlit Cloud](https://research-toolkit-yunanhu.streamlit.app/)


---

## Project Structure

```
research-toolkit/
├── app.py                        # Homepage and entry point
├── requirements.txt              # Python dependencies
├── .streamlit/
│   └── config.toml               # Theme and layout settings
│
├── pages/                        # One file = one tool page (auto-discovered)
│   ├── 01_qPCR_Analyzer.py
│   ├── 02_Growth_Curve_Analyzer.py
│   └── 03_Statistical_Analysis.py
│
├── analysis/                     # Pure analysis logic (no Streamlit code)
│   ├── qpcr.py
│   ├── growth_curve.py
│   └── statistics.py
│
├── visualization/                # Plotly figure builders
│   ├── qpcr_plots.py
│   ├── growth_plots.py
│   └── stats_plots.py
│
├── utils/                        # Shared helpers
│   ├── file_io.py
│   ├── validators.py
│   └── export.py
│
├── components/                   # Reusable Streamlit UI widgets
│   ├── sidebar.py
│   ├── upload_widget.py
│   └── result_display.py
│
├── data/                         # Sample datasets for testing
│   ├── sample_qpcr.csv
│   ├── sample_growth_curve.csv
│   └── sample_stats_data.csv
│
└── tests/                        # pytest unit tests
    ├── test_qpcr.py
    ├── test_growth_curve.py
    └── test_statistics.py
```

---

## Local Installation

### Prerequisites
- Python 3.10 or newer ([download](https://www.python.org/downloads/))
- Git ([download](https://git-scm.com/downloads))

### Step-by-step

```bash
# 1. Clone the repository
git clone https://github.com/Yunan0525/research-toolkit.git
cd research-toolkit

# 2. Create a virtual environment (keeps your system Python clean)
python -m venv venv

# 3. Activate the virtual environment
#    On macOS / Linux:
source venv/bin/activate
#    On Windows:
venv\Scripts\activate

# 4. Install all dependencies
pip install -r requirements.txt

# 5. Launch the app
streamlit run app.py
```

The app opens automatically at **http://localhost:8501**.

---

## Deployment on Streamlit Cloud

1. Push this repository to GitHub (must be **public** for the free tier).
2. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with GitHub.
3. Click **New app**.
4. Select your repository, branch (`main`), and main file (`app.py`).
5. Click **Deploy** — your app is live in ~2 minutes.

Every `git push` to `main` triggers an automatic redeploy.

---

## Running Tests

```bash
# From the project root with the virtual environment active:
pytest tests/ -v
```

---

## Input File Formats

### qPCR Analyzer (`sample_qpcr.csv`)
| Column | Description | Example |
|--------|-------------|---------|
| `Sample` | Sample identifier | `WT_1` |
| `Gene` | Target or reference gene name | `IL6`, `GAPDH` |
| `Ct` | Ct (threshold cycle) value | `24.5` |
| `Group` | Experimental group | `Control`, `Treatment` |

### Growth Curve Analyzer (`sample_growth_curve.csv`)
| Column | Description | Example |
|--------|-------------|---------|
| `Time_h` | Time in hours | `0`, `1`, `2` |
| `OD600` | Optical density at 600 nm | `0.05` |
| `Replicate` | Replicate identifier | `R1`, `R2` |
| `Group` | Treatment group | `Control`, `Drug_A` |

### Statistical Analysis (`sample_stats_data.csv`)
| Column | Description | Example |
|--------|-------------|---------|
| `Group` | Group label | `Control`, `Treatment` |
| `Value` | Numeric measurement | `1.23` |
| `Subject` | Subject ID (for paired tests) | `P01` |

---

## Contributing

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/new-module`.
3. Commit your changes: `git commit -m "Add new-module"`.
4. Push and open a Pull Request.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgements

Built with [Streamlit](https://streamlit.io), [Plotly](https://plotly.com),
[SciPy](https://scipy.org), [pandas](https://pandas.pydata.org),
and [pingouin](https://pingouin-stats.org).
