# Community assembly and biodiversity–ecosystem functioning analyses

This repository contains the code used in:

de León et al. (in prep.)
*Biodiversity-function relationship determined by the niche-neutral gradient of metacommunity assembly

## Repository structure

The analyses are organized into independent scripts that should be run sequentially.

### 1. CATS.R

Calculates the annual position of the metacommunity along the niche–neutral gradient using the Community Assembly via Trait Selection (CATS) framework.


---

### 2. BEF.R

Estimates the annual biodiversity–ecosystem functioning (BEF) relationship.

This script includes:

- Primary BEF estimation
- Environmental covariate models
- Linear vs quadratic robustness analyses
- Richness-restricted analyses (1–9 species)
- Quantile regression analyses
- Generation of Figs.S4-S13

---

### 3. Functional_diversity.R

Calculates annual functional diversity metrics (FRic, FEve, FDiv, FDis and Rao's Q) and performs mediation analyses relating community assembly, richness and functional diversity.

---

### 4. NST.R

Calculates the Normalized Stochasticity Ratio (NST) under different null models and regional species pools, and compares NST with CATS estimates.

This script includes:

- Generation of Figs.S1

---

### 5. Main_figures.R

Generates all figures included in the main manuscript and Figs.S2-S3.

---

## Running order

The scripts should be executed in the following order:

```
CATS.R
      ↓
BEF.R
      ↓
Functional_diversity.R
      ↓
NST.R
      ↓
Main_figures.R
```

---

## Data requirements

The analyses require the following input datasets:

- `br.xlsx`
- `bm.xlsx`
- `Traits.xlsx`
- `traits_numeric.csv`
- `Environmental.xlsx`
- `Ambientales_resumen.xlsx`
-  `Climatic.xlsx`

Dryad link where to find this datasets were shared exclusively with reviewers (NOT FOR PUBLICATION). Database will be published soon. 

---

## Software

Analyses were performed in R.

Main packages:

- glmmTMB
- MuMIn
- FD
- cluster
- NST
- quantreg
- visreg
- bestglm
- ggplot2
- ape
- vegan

---

## Notes

The primary manuscript results are based on the annual CATS estimates (`nn.1`, shared here as "CATS.csv") and the BEF estimates obtained from the environmental covariate models (`pendientes.covariables`).

The remaining analyses are included as robustness tests and supplementary analyses.
