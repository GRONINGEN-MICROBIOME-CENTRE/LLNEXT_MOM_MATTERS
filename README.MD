# Mom Matters: Maternal influences on infant gut microbiome and health in 714 mother–infant pairs in the Lifelines NEXT cohort

Code accompanying:

> **Sinha T. et al. (2026).** *Maternal influences on infant gut microbiome and health in 714 mother–infant pairs in the Lifelines NEXT cohort.*

---

# Table of Contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
  - [General scripts](#general-scripts)
  - [Mother–infant gut microbiome in Lifelines NEXT](#motherinfant-gut-microbiome-in-lifelines-next)
  - [Delivery mode shapes the infant gut microbiome](#delivery-mode-shapes-the-infant-gut-microbiome)
  - [Birth environment and the infant microbiome](#birth-environment-and-the-infant-microbiome)
  - [Infant feeding and neuroactive potential](#infant-feeding-and-neuroactive-potential)
  - [Infant carbohydrate-active enzymes (CAZymes)](#infant-carbohydrate-active-enzymes-cazymes)
  - [Maternal microbiome predicts eczema](#maternal-microbiome-predicts-eczema)
  - [Mother–infant strain sharing](#motherinfant-strain-sharing)
  - [Infant intraspecies phylogenies](#infant-intraspecies-phylogenies)
- [Data availability](#data-availability)
- [Contact](#contact)

---

# Overview

This repository contains the scripts used for the analyses presented in:

**Sinha T. et al. (2026).** *Maternal influences on infant gut microbiome and health in 714 mother–infant pairs in the Lifelines NEXT cohort.*

The repository includes code for data preprocessing, statistical analyses, visualization, and figure generation used throughout the manuscript.

---

# Repository structure

## General scripts

General utility functions and helper scripts used throughout the manuscript, including:

- Shared functions
- Plotting functions
- Demo datasets

---

## Mother–infant gut microbiome in Lifelines NEXT

Scripts for the overall characterization of the maternal and infant gut microbiome.

Includes:

- Cleaning and preprocessing of metadata
- Creation of cohort summary statistics
- Principal Coordinates Analysis
- Overall microbiome composition analyses
- Associaion of maternal and infant microbiome with time point
- Inter- and intra-individual distance analyses
- Community clustering analyses

**Figures**

- Figure 1
- Extended Data 1
- Extended Data 2
- Extended Data 3

---

## Delivery mode shapes the infant gut microbiome

Scripts investigating the impact of delivery mode on infant microbiome development.

Includes:

- MEFISTO analysis in relation to predictors
- TCAM analysis in relation to predictors
- PERMANOVA analyses in relation to predictors
- Associations between clinical predictors and:
  - Species-level genome bins (SGBs)
  - Alpha diversity
  - Functional pathways

**Figures**

- Figure 2
- Extended Data 4

---

## Birth environment and the infant microbiome

Scripts investigating birth-related environmental factors influencing the infant microbiome.

**Figures**

- Figure 3

---

## Infant feeding and neuroactive potential

Scripts investigating associations between infant feeding and microbial neuroactive potential.

Includes:

- Associations between predictors and human milk oligosaccharides (HMOs)
- Associations between predictors and gut–brain modules (GBMs)

**Figures**

- Extended Data 5

---

## Infant carbohydrate-active enzymes (CAZymes)

Scripts related to carbohydrate-active enzyme analyses.

Includes:

- Cleaning CAZyme annotations
- Associations between CAZymes and SGBs
- Longitudinal analyses of CAZymes across time
- Associations of CAZymes with infant predictors
- CAZyme substrate analyses

**Figures**

- Figure 3
- Extended Data 5
- Extended Data 6

---

## Maternal microbiome predicts eczema

Scripts investigating maternal microbiome characteristics associated with infant eczema.

Includes:

- Associations between maternal predictors and:
  - SGBs
  - Alpha diversity
  - Functional pathways
- Maternal microbiome prediction models for infant eczema

**Figures**

- Figure 4

---

## Mother–infant strain sharing

Scripts investigating microbial strain transmission between mothers and infants.

Includes:

- Maternal gut strain transmission
- Breast milk strain transmission
- Vaginal strain transmission
- Associations between strain sharing, time, microbial abundance and persistance and clinical predictors

**Figures**

- Figure 5
- Extended Data 7
- Extended Data 8
- Extended Data 9

---

## Infant intraspecies phylogenies

Scripts for phylogenetic analyses of infant microbial populations.

Includes:

- PGLMM analyses

**Figures**

- Extended Data 10

---

# Data availability

Quality-trimmed, human-decontaminated shotgun metagenomic sequencing data are publicly available through the NCBI Sequence Read Archive (SRA) under accession **PRJNA1468137**.

Metadata (including questionnaire and clinical information) together with quality-trimmed, human-decontaminated sequencing reads are available through the European Genome-Phenome Archive (EGA):

- **EGAD50000001187** — Fecal samples
- **EGAD50000002469** — Breast milk and vaginal samples
- **EGAD50000002677** — Metadata

Access can be requested through the Lifelines NEXT application form:

https://forms.gle/A4Jem2rMnjcygWRD6

Access is available to qualified researchers under the Lifelines NEXT Data Access Agreement:

https://groningenmicrobiome.org/?page_id=2598

This controlled-access procedure ensures that metadata are used exclusively for scientific research purposes and in accordance with the informed consent provided by Lifelines NEXT participants.

Detailed instructions describing the EGA application procedure are available at:

https://github.com/GRONINGEN-MICROBIOME-CENTRE/LLNEXT_pilot/blob/main/Data_Access_EGA.md

No fees are associated with access to the sequencing data or metadata.

Metagenome-assembled genomes (MAGs) recovered from fecal samples are publicly available through Zenodo:

https://zenodo.org/records/18965233


---

# Contact

For questions regarding the analysis workflow or datasets, please contact:

**Trishla Sinha, MD, PhD**  
University Medical Center Groningen  
Groningen Microbiome Centre

📧 t.sinha@rug.nl
