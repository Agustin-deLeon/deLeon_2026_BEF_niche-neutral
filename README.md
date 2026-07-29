# Community assembly and biodiversity–ecosystem functioning analyses

This repository contains the code used in:

de León et al. (in prep.)
*Biodiversity-function relationship determined by the niche-neutral gradient of metacommunity assembly

## Repository structure

The analyses are organized into independent scripts that should be run sequentially.

### 1. CATS.R

Calculates the annual position of the metacommunity along the niche–neutral gradient using the Community Assembly via Trait Selection (CATS) framework.

Main outputs:

- `nn.1`
- `bm.nuevo`
- `br.muestreos`
- `profundidad`

---

### 2. BEF.R

Estimates the annual biodiversity–ecosystem functioning (BEF) relationship.

This script includes:

- Primary BEF estimation
- Environmental covariate models
- Linear vs quadratic robustness analyses
- Richness-restricted analyses (1–9 species)
- Quantile regression analyses

Main outputs:

- `pendientes.lineal`
- `pendientes.covariables`
- `slopes_df`
- `slopes_df_sat`
- `slopes_df_sat_linear`
- `slopes_covariates_q90`

---

### 3. Functional_diversity.R

Calculates annual functional diversity metrics (FRic, FEve, FDiv, FDis and Rao's Q) and performs mediation analyses relating community assembly, richness and functional diversity.

Main outputs:

- `observado`
- `datos_sem`
- `tabla_sem`

---

### 4. NST.R

Calculates the Normalized Stochasticity Ratio (NST) under different null models and regional species pools, and compares NST with CATS estimates.

Main outputs:

- `NST_out`

---

### 5. Main_figures.R

Generates all figures included in the main manuscript and Supplementary Information.

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
- `traits_numerico_n.csv`
- `Environmental.xlsx`
- `Ambientales_resumen.xlsx`

Additional annual climatic variables are imported from:

- `df_plot2`

---

## Software

Analyses were performed in R (version X.X.X).

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

The primary manuscript results are based on the annual CATS estimates (`nn.1`) and the BEF estimates obtained from the environmental covariate models (`pendientes.covariables`).

The remaining analyses are included as robustness tests and supplementary analyses.
