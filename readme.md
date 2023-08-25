# The coenzyme A precursor pantethine restrains sarcoma growth through promotion of type 1 immunity

## Article information

**Title:** The coenzyme A precursor pantethine restrains sarcoma growth through promotion of type 1 immunity

**Authors:** Richard MIALLOT 1, Virginie MILLET 1, Anais ROGER 1, Romain FENOUIL 1, Catherine TARDIVEL 2, Jean-Charles MARTIN 2, Laetitia SHINTU 3, Paul BERCHARD 4, Juliane SOUSA LANZA 1, Bernard MALISSEN 1 5, Sandrine HENRI 1, Sophie UGOLINI 1, Aurélie DUTOUR 4, Pascal FINETTI 6, François BERTUCCI 6 7, Jean-Yves BLAY 4 8, Franck GALLAND 1, Philippe NAQUET 1

1. Aix-Marseille Université, INSERM, CNRS, Centre d'Immunologie de Marseille-Luminy, Marseille (CIML), France. 
2. Aix Marseille Université, INRAE, INSERM, C2VN, Marseille, France.
3. ISM2, Aix Marseille Université, CNRS, Centrale Marseille, Marseille, France.
4. INSERM 1052, CNRS 5286, Cancer Research Center of Lyon (CRCL), Childhood Cancers and Cell Death Laboratory, Lyon, France.
5. Aix Marseille Université, INSERM, CNRS, Centre d’Immunophénomique (CIPHE), Marseille, France 
6. Aix-Marseille Université, INSERM, CNRS, Centre de Recherche en Cancérologie de Marseille (CRCM), Institut Paoli‑Calmettes (IPC), Laboratory of Predictive Oncology, Marseille, France.
7. Institut Paoli‑Calmettes, Department of Medical Oncology, Marseille, France.
8. Université Lyon I, UNICANCER Centre Léon Bérard, Department of Medicine, Lyon, France.

**Summary:**
The tumor microenvironment is a dynamic network of stromal, cancer and immune cells that interact and compete for resources. We have previously identified the Vanin1 pathway as a tumor suppressor of sarcoma development via vitamin B5 and coenzyme A regeneration. Using an aggressive sarcoma cell line that lacks Vnn1 expression, we showed that the administration of pantethine, a vitamin B5 precursor, attenuates tumor growth in immunocompetent mice. Pantethine boosts anti-tumor immunity, including the polarization of myeloid and dendritic cells towards enhanced IFNγ-driven antigen presentation pathways and improved development of hypermetabolic effector CD8+ T cells endowed with potential anti-tumor activity. At later stages of treatment, the effect of pantethine was limited by the development of immune cell exhaustion. Nevertheless, its activity was comparable to that of anti-PD1 treatment in sensitive tumors. In humans, VNN1 expression correlates with improved survival and immune cell infiltration in soft tissue sarcomas, but not in osteosarcomas. Pantethine could be a potential therapeutic immunoadjuvant for the development of anti-tumor immunity. 

---

## Goal of the github
This github project contains Source code (scripts and dockerfiles) used to produce analyses reported in the article (and intermediate/additional analyses).
Required data is published on gene expression omnibus (GEO). Docker/Singularity images and results of analyses are available for download on Zenodo. 


---

## Description of the datasets

There are 4 datasets in this study, each being a replicate containing the 4 experimental conditions (split by HTO).

Indvidual samples were aligned with CellRanger and analyzed individually to check for quality and identify cell populations.
All samples were then unified in a single (merged) analysis for comparison of experimental conditions.

Resulting data, html reports, and files were uploaded to Zenodo with following DOIs.
RUN 1, cDNA 1 (align & QC): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8277073.svg)](https://doi.org/10.5281/zenodo.8277073)
RUN 1, cDNA 2 (align & QC): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8279739.svg)](https://doi.org/10.5281/zenodo.8279739)
RUN 2, cDNA 1 (align & QC): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8279752.svg)](https://doi.org/10.5281/zenodo.8279752)
RUN 2, cDNA 2 (align & QC): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8279754.svg)](https://doi.org/10.5281/zenodo.8279754)
Merged dataset (full analysis): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8279759.svg)](https://doi.org/10.5281/zenodo.8279759)


## Recommendations for reproducibility

Software environment were created with Docker/Singularity container technologies.
Dockerfiles are available in this repository and corresponding container as binary files on Zenodo archives.

In order to replicate existing results, one should download binary container images as rebuilding from source does not guarantee to result in the same execution environment.

## Download, adapt paths, compile reports

For each sample, analyses steps are organized in separate folders, numbered in the order of execution.
Each step is compiled individually.

The folders hierarchy containing code in this repository should replicate the results folder downloadable from Zenodo.
Code and data/results can (should) be stored in separate foder structures.

The file 'globalParams.R' defines path to the root folders of the data/analyses and should be adapted to the execution environment.

Each analysis step contains a script named 'launch_reports_compilation.R' that should be executed within the correpsonding container.
It loads the parameters/path from globalParams, and starts the compilation of HTML from corresponding Rmd file with Knitr.


