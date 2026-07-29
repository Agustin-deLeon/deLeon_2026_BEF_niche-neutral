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

M_obs <- M_UM %>%
  group_by(año) %>%
  summarise(
    across(starts_with("Acmella_"):last_col(), sum),
    .groups = "drop"
  )

M_obs <- M_obs[,c(TRUE, colSums(M_obs[,-1]) > 0)]
comm <- as.matrix(M_obs[,-1])
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
# FUNCIÓN DE MEDIACIÓN / PATH ANALYSIS
#=========================================================

analisis_vias <- function(
    datos,
    respuesta,
    n_boot = 10000,
    seed = 1234
) {
  
  vars <- c("R2_traits", "riqueza", respuesta)
  
  dat <- datos[
    complete.cases(datos[, vars]),
    vars,
    drop = FALSE
  ]
  
  # Estandarizar para obtener coeficientes comparables
  dat_std <- as.data.frame(scale(dat))
  
  #-------------------------------------------------------
  # Modelos del path
  #-------------------------------------------------------
  
  # Camino a:
  # R2_traits -> riqueza
  mod_a <- lm(
    riqueza ~ R2_traits,
    data = dat_std
  )
  
  # Caminos b y c':
  # riqueza -> FD
  # R2_traits -> FD
  mod_b <- lm(
    reformulate(
      c("riqueza", "R2_traits"),
      response = respuesta
    ),
    data = dat_std
  )
  
  # Efecto total c:
  mod_total <- lm(
    reformulate(
      "R2_traits",
      response = respuesta
    ),
    data = dat_std
  )
  
  a <- unname(coef(mod_a)["R2_traits"])
  b <- unname(coef(mod_b)["riqueza"])
  c_directo <- unname(coef(mod_b)["R2_traits"])
  c_total <- unname(coef(mod_total)["R2_traits"])
  indirecto <- a * b
  
  # Bootstrap del efecto indirecto
  set.seed(seed)
  boot_indirecto <- replicate(
    n_boot,
    {
      idx <- sample(
        seq_len(nrow(dat)),
        size = nrow(dat),
        replace = TRUE
      )
      
      boot_dat <- dat[idx, , drop = FALSE]
      
      if (
        sd(boot_dat$R2_traits) == 0 ||
        sd(boot_dat$riqueza) == 0 ||
        sd(boot_dat[[respuesta]]) == 0
      ) {
        return(NA_real_)
      }
      
      boot_std <- as.data.frame(scale(boot_dat))
      
      ma <- try(
        lm(
          riqueza ~ R2_traits,
          data = boot_std
        ),
        silent = TRUE
      )
      
      mb <- try(
        lm(
          reformulate(
            c("riqueza", "R2_traits"),
            response = respuesta
          ),
          data = boot_std
        ),
        silent = TRUE
      )
      
      if (
        inherits(ma, "try-error") ||
        inherits(mb, "try-error")
      ) {
        return(NA_real_)
      }
      
      aa <- coef(ma)["R2_traits"]
      bb <- coef(mb)["riqueza"]
      
      as.numeric(aa * bb)
    }
  )
  
  boot_indirecto <- boot_indirecto[is.finite(boot_indirecto)]
  
  IC_indirecto <- quantile(
    boot_indirecto,
    probs = c(0.025, 0.975),
    na.rm = TRUE)
  
  p_boot <- 2 * min(
    mean(boot_indirecto <= 0),
    mean(boot_indirecto >= 0))
  
  p_boot <- min(p_boot, 1)
  
  # Resultados
  coef_b <- summary(mod_b)$coefficients
  tabla <- data.frame(
    metrica = respuesta,
    a_R2_a_riqueza = a,
    p_a = summary(mod_a)$coefficients[
      "R2_traits",
      "Pr(>|t|)"
    ],
    b_riqueza_a_FD = b,
    p_b = coef_b[
      "riqueza",
      "Pr(>|t|)"
    ],
    
    efecto_directo_R2 = c_directo,
    p_directo = coef_b[
      "R2_traits",
      "Pr(>|t|)"
    ],
    
    efecto_indirecto = indirecto,
    IC_indirecto_2.5 = IC_indirecto[1],
    IC_indirecto_97.5 = IC_indirecto[2],
    p_indirecto_boot = p_boot,
    
    efecto_total_R2 = c_total,
    p_total = summary(mod_total)$coefficients[
      "R2_traits",
      "Pr(>|t|)"
    ],
    
    n = nrow(dat)
  )
  
  cat("\n=================================================\n")
  cat("MÉTRICA:", respuesta, "\n")
  cat("=================================================\n")
  
  cat("\nModelo de riqueza:\n")
  print(summary(mod_a))
  
  cat("\nModelo de diversidad funcional:\n")
  print(summary(mod_b))
  
  cat("\nModelo del efecto total:\n")
  print(summary(mod_total))
  
  cat("\nEfectos estandarizados:\n")
  print(tabla)
  
  invisible(
    list(
      tabla = tabla,
      modelo_riqueza = mod_a,
      modelo_FD = mod_b,
      modelo_total = mod_total,
      bootstrap_indirecto = boot_indirecto
    )
  )
}

sem_FDis <- analisis_vias(
  datos = datos_sem,
  respuesta = "FDis",
  n_boot = 10000,
  seed = 1234
)

sem_FDis$tabla

metricas_sem <- c(
  "FRic",
  "FEve",
  "FDiv",
  "FDis",
  "RaoQ"
)

resultados_sem <- lapply(
  metricas_sem,
  function(m) {
    analisis_vias(
      datos = datos_sem,
      respuesta = m,
      n_boot = 10000,
      seed = 1234
    )
  }
)

names(resultados_sem) <- metricas_sem

tabla_sem <- do.call(
  rbind,
  lapply(
    resultados_sem,
    function(x) x$tabla
  )
)

rownames(tabla_sem) <- NULL

tabla_sem

tabla_sem_redondeada <- tabla_sem

cols_num <- sapply(
  tabla_sem_redondeada,
  is.numeric
)

tabla_sem_redondeada[cols_num] <-
  lapply(
    tabla_sem_redondeada[cols_num],
    round,
    digits = 4
  )

tabla_sem_redondeada



graficar_vias <- function(resultado, nombre_metrica) {
  
  tab <- resultado$tabla
  
  plot(
    NA,
    xlim = c(0, 10),
    ylim = c(0, 10),
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )
  
  # Cajas
  rect(0.5, 4, 3, 6, lwd = 2)
  rect(4, 7, 6.5, 9, lwd = 2)
  rect(7, 4, 9.5, 6, lwd = 2)
  
  text(1.75, 5, expression(R[Traits]^2), cex = 1.3)
  text(5.25, 8, "Species richness", cex = 1.2)
  text(8.25, 5, nombre_metrica, cex = 1.3)
  
  # R2 -> riqueza
  arrows(
    3, 5.5,
    4.2, 7.4,
    length = 0.12,
    lwd = 2
  )
  
  text(
    3.2,
    7,
    paste0(
      "a = ",
      round(tab$a_R2_a_riqueza, 2),
      "\np = ",
      signif(tab$p_a, 2)
    )
  )
  
  # riqueza -> FD
  arrows(
    6.3, 7.4,
    7.2, 5.6,
    length = 0.12,
    lwd = 2
  )
  
  text(
    7.1,
    7,
    paste0(
      "b = ",
      round(tab$b_riqueza_a_FD, 2),
      "\np = ",
      signif(tab$p_b, 2)
    )
  )
  
  # R2 -> FD
  arrows(
    3, 5,
    7, 5,
    length = 0.12,
    lwd = 2
  )
  
  text(
    5,
    4.4,
    paste0(
      "c' = ",
      round(tab$efecto_directo_R2, 2),
      "\np = ",
      signif(tab$p_directo, 2)
    )
  )
  
  title(
    main = paste(
      "Path analysis:",
      nombre_metrica
    )
  )
}


