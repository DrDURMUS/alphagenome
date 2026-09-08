# AlphaGenome Studio & Atlas CLI 🧬

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-green.svg)](https://opensource.org/licenses/Apache-2.0)
[![AlphaGenome](https://img.shields.io/badge/Powered%20By-Google%20DeepMind%20AlphaGenome-6366f1.svg)](https://github.com/google-deepmind/alphagenome)
[![Status](https://img.shields.io/badge/status-active-success.svg)]()

**AlphaGenome Studio** is an integrated bioinformatics platform and web interface for analyzing human genome variant effects using **Google DeepMind's AlphaGenome** and **AlphaGenome Atlas** (9 billion precomputed SNVs across 1 Petabyte of functional genomics data).

The platform provides an intuitive **Web Studio (Drag & Drop)**, a high-performance **CLI tool**, local caching, and one-click desktop launchers for rapid clinical and molecular genetic research.

---

## 🌟 Key Features

- ⚡ **AlphaGenome Atlas Engine (9 Billion SNVs)**:
  - Instant query access to over 9 billion precomputed single nucleotide variants (SNVs) across the entire human genome.
  - Covers 30+ regulatory scorers including gene expression log-fold change (`GeneMaskLFCScorer`), gene activation (`GeneMaskActiveScorer`), splicing alterations (`SpliceJunctionScorer`), and chromatin accessibility (`CenterMaskScorer`).
  - Automatic persistent local disk caching (`~/.cache/alphagenome`) for millisecond repeat queries.
- 🧬 **Direct Model Inference**:
  - Live deep learning sequence-to-function predictions for custom intervals, novel variants, and multi-track genomic profiles (RNA-seq, CAGE, DNase, ChIP-seq, Splice sites).
- 🖥️ **AlphaGenome Web Studio**:
  - **Drag & Drop**: Analyze `.vcf`, `.vcf.gz`, `.tsv`, `.csv`, and `.txt` files effortlessly with automatic preview and validation.
  - **Decision Support**: Automated variant effect tiering (🔴 Distinct Effect, 🟠 Moderate, 🟡 Weak, 🟢 Negligible) and mechanism breakdowns.
  - **Interactive Table**: Search, filter by biological impact, sort, and export to TSV/CSV.
  - **Scorer Catalog**: Searchable documentation of all 30+ Atlas scorers.
  - **Zero External UI Dependencies**: Powered by Python's built-in `http.server.ThreadingHTTPServer`.
- 🖱️ **Double-Click Launchers (macOS)**:
  - Launch the web application instantly without touching the command line via `AlphaGenome Studio.command` or the native desktop application generator (`create_desktop_app.sh`).

---

## 📁 Repository Structure

```text
alphagenome/
├── alphagenome_cli.py          # Unified CLI for Atlas queries & model inference
├── web_app.py                  # Web Studio backend server (standard library HTTP server)
├── ui/                         # Modern Glassmorphic Web UI
│   ├── index.html              # Responsive single-page application
│   ├── style.css               # Design system & dark theme
│   └── app.js                  # Frontend state & interactive data tables
├── start_studio.sh             # Portable launcher shell script
├── create_desktop_app.sh       # macOS Desktop shortcut/applet generator
├── AlphaGenome Studio.command  # Double-clickable macOS launcher script
├── test_atlas.py               # Test suite for Atlas client, parser, & web API
├── example_variants.tsv        # Sample genomic input data
├── KULLANIM_KILAVUZU.md        # Comprehensive user manual (Türkçe)
└── README.md                   # Project documentation
```

---

## 🚀 Quickstart & Installation

### 1. Clone the Repository
```bash
git clone https://github.com/DrDURMUS/alphagenome.git
cd alphagenome
```

### 2. Set Up Python Environment
Create and activate a Python 3.10+ Conda environment:
```bash
conda create -n alphagenome python=3.11 -y
conda activate alphagenome
```

Install the required packages:
```bash
pip install -r requirements.txt
```

### 3. AlphaGenome API Key
Obtain an API key from the [Google DeepMind AlphaGenome portal](https://alphagenome.deepmind.google/).

Set it in your terminal environment:
```bash
export ALPHAGENOME_API_KEY="your_api_key_here"
```
*(Alternatively, you can enter your API key directly through the Web Studio settings modal).*

---

## 💻 Usage

### 🌐 Option A: Web Studio (Recommended)

#### 1. Double-Click Launch (macOS)
- Simply double-click **`AlphaGenome Studio.command`** in the project folder.
- *Optional:* Run `./create_desktop_app.sh` to generate a native **`AlphaGenome Studio.app`** directly on your Desktop!

#### 2. Terminal Launch
```bash
conda activate alphagenome
python web_app.py --port 8000
```
Or via the CLI:
```bash
python alphagenome_cli.py ui --port 8000
```
Your browser will automatically open to `http://127.0.0.1:8000`.

---

### ⌨️ Option B: Command Line Interface (CLI)

#### 1. Query Single Variant from Atlas (Instant)
```bash
# Query a single SNV across all scorers
python alphagenome_cli.py atlas variant chr22:36201698:A>C

# Filter by tissue or cell type
python alphagenome_cli.py atlas variant chr22:36201698:A>C --tissue colon

# Filter by gene name
python alphagenome_cli.py atlas variant chr22:36201698:A>C --gene APOL4
```

#### 2. Query Genomic Interval from Atlas
```bash
python alphagenome_cli.py atlas interval chr22:35677410-35678410 --tissue liver
```

#### 3. Batch Variant Analysis (VCF or TSV)
```bash
# Analyze a VCF file with Atlas
python alphagenome_cli.py atlas vcf input.vcf --output results.tsv

# Analyze a TSV variant list
python alphagenome_cli.py atlas list variants.tsv --output results.tsv
```

#### 4. List All Available Scorers
```bash
python alphagenome_cli.py atlas scorers
```

#### 5. Live Deep Learning Inference (Model Mode)
```bash
# Direct sequence forward prediction
python alphagenome_cli.py predict chr22:35677410-36725986 --tissue colon

# In-silico variant scoring with local model
python alphagenome_cli.py variant chr22:36201698:A>C --score
```

---

## 🧪 Running Tests

To run the complete automated test suite:
```bash
python test_atlas.py
```
This tests:
- `tidy_atlas_scores`: AnnData conversion to clean pandas DataFrames
- VCF and TSV variant file parsers
- Local disk caching (`CachedAtlasClient`)
- Variant and interval coordinate parsers
- Web Studio API endpoints (`/api/status`, `/api/scorers`, etc.)

---

## 📖 Clinical & Biological Interpretation Guide

When analyzing regulatory variants, AlphaGenome Studio organizes predictions into functional mechanism tiers:

| Impact Tier | Description | Typical Metrics |
| :--- | :--- | :--- |
| 🔴 **Distinct (Belirgin)** | Strong regulatory or splicing perturbation | High quantile score (≥ 0.95), absolute LFC ≥ 0.5, or SpliceSite delta ≥ 0.2 |
| 🟠 **Moderate (Orta)** | Detectable regulatory effect | Quantile score 0.85 - 0.95 |
| 🟡 **Weak (Zayıf)** | Borderline/minor regulatory change | Quantile score 0.70 - 0.85 |
| 🟢 **Negligible (İhmal)** | Benign / baseline background variation | Quantile score < 0.70 |

> **Note**: In-silico predictions are intended for research and decision-support purposes. Clinical findings should be corroborated with established guidelines (ACMG/AMP/ClinGen) and targeted experimental validation (e.g., SpliceAI, RNA-seq).

---

## 📄 License & Disclaimer

This project is licensed under the Apache 2.0 License.

**Disclaimer**: AlphaGenome and AlphaGenome Atlas are research tools developed by Google DeepMind. This tool is provided for academic, biomedical, and clinical research purposes.
