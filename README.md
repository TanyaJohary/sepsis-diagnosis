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
