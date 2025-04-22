# 🔬 Sepsis Diagnostics Project

This repository hosts a computational pipeline for identifying and prioritizing genetic biomarkers for sepsis using supervised machine learning, statistical analysis, and biological network interpretation. The project aims to support the development of robust diagnostic tools for sepsis through multi-dataset integration and rigorous biological validation.

---

## 📚 Project Overview

Sepsis is a life-threatening condition caused by a dysregulated host response to infection. Early diagnosis remains a major challenge due to its heterogeneous nature. This project explores gene expression data from multiple public datasets to discover consistent transcriptomic signals associated with sepsis, and to define a compact, interpretable diagnostic signature.

Key objectives of this project:
- Detect robust sepsis-associated genes across heterogeneous transcriptomic datasets.
- Develop a compact, high-performing gene panel using ML-based feature selection.
- Validate the biological coherence of the gene signature using protein–protein interaction networks and pathway enrichment.

---

## 🧠 Methodological Pipeline

The analysis pipeline involves the following stages:

### 🧮 1. Statistical Testing
- **Mann–Whitney U Test**: Evaluate differential expression between sepsis and control samples.
- **Effect Size Calculation**: Quantify biological relevance using Cohen’s d.

### 🌲 2. Machine Learning Feature Selection
- **Random Forest**: Compute variable importance across multiple datasets.
- **Stability Ranking**: Aggregate importance across resampled runs for robustness.

### 🧬 3. Signature Construction
- Ensemble ranking across statistical and ML-based metrics.
- Final 15-gene panel constructed by integrating consistent top-ranking genes.

### 🔗 4. Biological Validation
- **STRING PPI Network**: Explore interactions and functional modules.
- **Enrichment Analysis**: Validate pathways using Gene Ontology and KEGG.
- Group genes into meaningful modules:
  - *Innate immune activation*
  - *Immunosuppression*
  - *Neutrophil-driven inflammation*

---

📊 Datasets
This project uses publicly available gene expression datasets from the NCBI Gene Expression Omnibus (GEO). Eleven datasets were used to train the diagnostic model, and three datasets were used for validation.

🔧 Training Datasets

GSE185263 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185263  
GSE65682  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682  
GSE236713 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236713  
GSE131761 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131761  
GSE154918 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154918  
GSE69063  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69063  
GSE57065  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57065  
GSE100159 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100159  
GSE243217 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243217  
GSE28750  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28750  
GSE54514  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54514  



✅ Validation Datasets

GSE69528  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69528  
GSE60424  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424  
GSE63311  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63311  




## 📁 Repository Structure

The project directory is organized into logical folders for datasets, analysis outputs, and notebooks.  
Each training dataset folder contains both the dataset and the corresponding analysis code and results.

```text
├── Frequency Results/          # Frequency of gene selection across datasets by different methods
├── Functional Analysis/        # Enrichment analysis results (STRING, g:Profiler)
├── Genes/                      # Sepsis-related gene panels and annotation files
├── Tests outputs/              # Outputs from different tests (MWU, correlation)
├── Top Genes Validation/       # Evaluation of top-ranked gene sets across 11 investigated datasets
├── Training Datasets/          # Dataset-specific directories with code and outputs
│   ├── GSE185263/
│   │   ├── GSE185263code.R             # Code for preprocessing and analysis
│   │   ├── sepsis_dataGSE185263.csv    # Final filtered gene expression matrix
│   │   ├── RF-results/                 # Random Forest ranking code + results
│   │   ├── Mann-W-U test/             # Mann–Whitney U test code + results
│   │   ├── correlation test/          # Correlation filter code + results
│   │   └── top25percentile/           # Top 25% selection logic + results
│   └── ... (other GSE folders follow same format)
├── notebook/                  # Jupyter or RMarkdown notebooks for exploratory work
├── validation Datasets/      # Datasets used for external validation of diagnostic clusters
├── Project on genetic diagnostic signature.pdf # Project overview or report (PDF)
├── LICENSE.txt               # GPL-3.0 License file
└── README.md                 # Project overview and usage guide


