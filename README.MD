# Mom Matters: Maternal influences on infant gut microbiome and health in 714 mother–infant pairs in the Lifelines NEXT cohort

Code accompanying:

> **Sinha T. et al. (2026).** *Maternal influences on infant gut microbiome and health*

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

Here we present a comprehensive characterization of maternal and infant gut microbiome development during pregnancy and the first year of life, and identify the maternal, perinatal, and environmental factors shaping infant gut microbiome assembly and early health. In this study we metagenomically sequenced **4,526 longitudinal fecal samples** from **714 mother–infant pairs** in the Dutch **Lifelines NEXT (LLNEXT)** birth cohort, spanning **12 weeks of pregnancy to one year postpartum**. We integrated these data with **474 clinical and exposure variables**, together with **ultra-deep shotgun metagenomic sequencing of breast milk and vaginal microbiomes**, enabling a comprehensive analysis of maternal and infant microbiome dynamics.

🤰 **Maternal gut microbiome dynamics:** The maternal gut microbiome undergoes only subtle compositional changes throughout pregnancy and the postpartum period and is primarily influenced by diet, infections, and pre-pregnancy smoking.

🦠 **Maternal origin of the infant gut microbiome:** The maternal gut microbiome is the principal reservoir for infant gut bacterial strains, whereas transmission from the vaginal and breast milk microbiomes occurs only occasionally.

🔗 **Determinants of strain transmission:** Mother–infant gut strain sharing is highly time-dependent, and higher maternal abundance of a bacterial species increases the likelihood of successful strain transmission to the infant. Maternal strains persist longer in infants and are enriched in functions potentially useful for infant health. 

🩺 **Maternal microbiome and infant health:** The maternal gut microbiome is predictive of infant eczema, highlighting the importance of maternal microbial composition for early-life health outcomes.

🚼 **Drivers of infant gut microbiome development:** Delivery mode and feeding mode are the primary determinants of infant gut microbiome composition and functional maturation, while maternal exposures also contribute to microbiome development.

🏡 **Birth environment:** Among 585 vaginally delivered infants, 155 were born at home. Home birth was only moderately associated with infant gut microbiome composition, with effects comparable to other birth-related factors such as duration of pushing and rupture of membranes.

🌱 **Central role of the mother:** Together, these findings demonstrate the central role of the mother and her microbiome in shaping the infant gut ecosystem, microbial transmission, functional development, and early health outcomes.

---

# Cohort

This study is based on the **Lifelines NEXT (LLNEXT)** cohort, a large prospective, population-based birth cohort in the Netherlands designed to investigate early-life determinants of health and disease.

### Cohort overview and data catalogue

https://umcgresearchdatacatalogue.nl/all/cohorts/LIFELINES_NEXT

### Original cohort description

https://pubmed.ncbi.nlm.nih.gov/32100173/

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

## Infant feeding, HMOs and neuroactive potential

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
