###############################################################
###############################################################

# de León et al. 
# Biodiversity-function relationship determined
# by the niche-neutral gradient of community assembly

#=========================================================
### -------- FUNCTIONAL DIVERSITY ANALYSES ----------- ###
#=========================================================
#
# This script calculates annual functional diversity metrics
# from community composition and species trait data.
#
# Specifically, it:
#   (1) computes annual values of FRic, FEve, FDiv, FDis,
#       and Rao's Q using Gower distances and dbFD;
#   (2) evaluates the relationships between the annual
#       position along the niche–neutral gradient
#       (R²Traits), species richness, and functional
#       diversity using path (mediation) analyses;
#   (3) summarizes standardized direct, indirect, and total
#       effects for each functional diversity metric; and
#   (4) generates the path diagrams used in the manuscript.
#
# Required objects:
#   - M_UM: species-by-sampling-unit community matrix.
#   - traits: species trait matrix (species in rows).
#   - nn.1: annual CATS results, including R2.exp.
#   - riqueza: annual mean species richness.
#
# The object nn.1 must be generated previously by running
# CATS.R.
#=========================================================

library(FD)
library(cluster)
library(dplyr)

M_obs <- br.muestreos_df %>%
  group_by(año) %>%
  summarise(
    across(starts_with("Acmella_"):last_col(), sum),
    .groups = "drop"
  )

M_obs <- M_obs[,c(TRUE, colSums(M_obs[,-1]) > 0)]
comm <- as.matrix(M_obs[,-1])
rownames(traits)[rownames(traits) == "Symphyotrichum_squamatum"] <-
  "Symphyotrichum_graminifolium"
comm <- comm[, !colnames(comm) %in% c(
  "Pterocaulon_polystachyum",
  "Sisyrinchium_sp"
)]
traits_obs <- traits[colnames(comm), , drop = FALSE]
stopifnot(all(colnames(comm) == rownames(traits_obs)))


D_obs <- daisy(traits_obs,metric = "gower",stand = TRUE)
fd_obs <- dbFD(x = D_obs,a = comm,corr = "cailliez",m = 10)

observado <- data.frame(
  Año  = M_obs$año,
  FRic = fd_obs$FRic,
  FEve = fd_obs$FEve,
  FDiv = fd_obs$FDiv,
  FDis = fd_obs$FDis,
  FRao = fd_obs$RaoQ
)

observado


#=========================================================
# ---- Mediation analyses
#=========================================================

datos_sem <- data.frame(
  año       = observado$Año,
  R2_traits = nn.1$R2.exp,
  riqueza   = as.numeric(riqueza),
  
  FRic = observado$FRic,
  FEve = observado$FEve,
  FDiv = observado$FDiv,
  FDis = observado$FDis,
  RaoQ = observado$RaoQ
)

# Chequeos
stopifnot(
  nrow(datos_sem) == length(R2.exp),
  nrow(datos_sem) == length(riqueza)
)

summary(datos_sem)
cor(datos_sem[, -1], use = "pairwise.complete.obs")


#=========================================================
# PATH ANALYSIS
#=========================================================

library(lavaan)

metricas <- c("FRic", "FEve", "FDiv", "FDis", "FRao")

fits <- list()
resultados <- list()

for(m in metricas){
  
  datos <- data.frame(
    R2.exp  = scale(nn.1$R2.exp)[,1],
    riqueza = scale(nn.1$riqueza)[,1],
    y        = scale(observado[[m]])[,1]
  )
  
  summary(datos)
  print(var(datos$y))
  modelo <- '
  
    riqueza ~ a*R2.exp
    
    y ~ b*riqueza + c*R2.exp
    
    indirect := a*b
    direct   := c
    total    := c + (a*b)
    
  '
  
  fit <- sem(modelo, data = datos)
  
  fits[[m]] <- fit
  
  resultados[[m]] <-
    parameterEstimates(
      fit,
      standardized = TRUE
    )
  
}

func_path <- do.call(
  rbind,
  lapply(names(resultados), function(m){
    
    x <- resultados[[m]]
    
    x <- subset(
      x,
      label %in% c("a","b","c","indirect","total")
    )
    
    data.frame(
      metrica = m,
      efecto = x$label,
      beta = round(x$std.all,3),
      p = signif(x$pvalue,3)
    )
    
  })
)

func_path




library(DiagrammeR)

plot_mediation <- function(resumen,
                           metrica,
                           titulo = NULL,
                           digits = 2){
  
  if(is.null(titulo))
    titulo <- metrica
  
  r <- resumen[resumen$metrica == metrica, ]
  
  a <- r$beta[r$efecto == "a"]
  b <- r$beta[r$efecto == "b"]
  c <- r$beta[r$efecto == "c"]
  
  pa <- r$p[r$efecto == "a"]
  pb <- r$p[r$efecto == "b"]
  pc <- r$p[r$efecto == "c"]
  
  ptxt <- function(x){
    
    if(is.na(x)) return("NA")
    
    if(x < 0.001){
      "p < 0.001"
    }else{
      paste0("p = ", signif(x,2))
    }
    
  }
  
  grViz(sprintf("
digraph {

graph [layout=dot, rankdir=LR]

node[
shape=box,
style=rounded,
fontsize=20,
fontname=Helvetica,
margin=0.25]

A [label='R²Traits']
B [label='Metacommunity\\nspecies richness']
C [label='%s']

A -> B [label='a = %.2f\\n%s', fontsize=18]
B -> C [label='b = %.2f\\n%s', fontsize=18]
A -> C [label='c = %.2f\\n%s', fontsize=18]

labelloc='t'
label=<<B>%s</B>>
fontsize=28

}
",
titulo,
a, ptxt(pa),
b, ptxt(pb),
c, ptxt(pc),
titulo))
}

plot_mediation(func_path, "FDis", "Functional dispersion")

plot_mediation(func_path, "FEve", "Functional evenness")

plot_mediation(func_path, "FDiv", "Functional divergence")

plot_mediation(func_path, metrica="FRao", "Rao´s Q")

plot_mediation(func_path, "FRic", "Functional richness")

