###############################################################
###############################################################

# de León et al. 
# Biodiversity-function relationship determined
# by the niche-neutral gradient of community assembly

#=========================================================
### -------- BIODIVERSITY–ECOSYSTEM FUNCTIONING -------- ###
#=========================================================
#
# This script estimates the annual biodiversity–ecosystem
# functioning (BEF) relationship using several analytical
# approaches.
#
# Analyses should be run sequentially, as later sections
# depend on objects generated previously.
#
# Main analysis:
#   (1) Annual linear GLMMs (full richness range).
#       Output:
#         - pendientes.lineal
#
# Robustness analyses (Supplementary Information):
#   (2) Linear models accounting for environmental
#       covariates.
#         - pendientes.covariables
#
#   (3) Linear vs. quadratic GLMMs.
#         - slopes_df
#
#   (4) Linear vs. quadratic GLMMs after restricting the
#       richness range to 1–9 species.
#         - slopes_df_sat
#
#   (5) Linear GLMMs restricted to 1–9 species.
#         - slopes_df_sat_linear
#
#   (6) 90th quantile regression models.
#         - slopes_covariates_q90
#
# The final section reproduces the supplementary figures
# associated with the robustness analyses.
#
# Required objects:
#   - bm.nuevo (generated in CATS.R)
#   - nn.1 (generated in CATS.R)
#=========================================================




####################################
#GLMM
library(glmmTMB)
library(MuMIn)
library(viridisLite)

par(mfrow=c(6,4), mar=c(5.6,6,3,2), cex.lab = 2)

pendientes.lineal <- list()
pendientes.charcos <- list()

for (ii in 2005:2025) {
  
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:26)]
  colnames(df.um.ch) <- c(
    "um.año", "ch.id", "um.biom", "um.rich",
    "lluvia_muestreo", "lluvia_anual", "temp_muestreo", "temp_anual",
    "DM", "ddmm", "Shape", "Islands", "log.Area", "log.Volumen",
    "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree",
    "log.Betweenness", "Closenness"
  )
  
  # >>> NUEVO: eliminar riqueza 0 (o menor) <<<
  df.um.ch <- df.um.ch[df.um.ch$um.rich > 0, ]
  
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$ch.id <- as.factor(df.um.ch$ch.id)
  
  titulo <- paste(ii)
  
  # Si queda poco dato o 1 solo charco, saltar prolijo
  if (nrow(df.um.ch) < 5 || length(unique(df.um.ch$ch.id)) < 2) {
    plot.new()
    title(main = titulo)
    mtext("Datos insuficientes", side = 3, line = -2)
    next
  }
  
  m1 <- glmmTMB(
    um.biom ~ um.rich + (um.rich | ch.id),
    family = gaussian(),
    data = df.um.ch
  )
  
  fixef_m1 <- fixef(m1)$cond
  ranef_intercepts <- ranef(m1)$cond$ch.id[, "(Intercept)"]
  ranef_slopes <- ranef(m1)$cond$ch.id[, "um.rich"]
  
  pendientes.charcos[[as.character(ii)]] <- ranef_slopes
  
  unique_ids <- levels(df.um.ch$ch.id)
  
  # =====================================================
  # Predicción marginal (efectos fijos)
  # =====================================================
  
  pred <- as.data.frame(
    ggpredict(
      m1,
      terms = "um.rich [all]"
    )
  )
  
  # =====================================================
  # Gráfico base
  # =====================================================
  
  colors <- rep("black", length(unique_ids))
  names(colors) <- unique_ids
  
  plot(
    df.um.ch$um.rich,
    df.um.ch$um.biom,
    
    xlab = "Richness",
    
    ylab = expression(
      "Biomass (g / 0.04 m"^2*")"
    ),
    
    main = titulo,
    
    pch = 19,
    cex = 0.8,
    
    col = adjustcolor(
      "black",
      alpha.f = 0.45
    ),
    
    #bty = "l",
    
    cex.lab = 2.2,
    cex.axis = 2,
    cex.main = 1.7
  )
  
  # =====================================================
  # Pendientes por charco
  # =====================================================
  
  for(i in seq_along(unique_ids)){
    
    id_i <- unique_ids[i]
    
    current_data <-
      df.um.ch[
        df.um.ch$ch.id == id_i,
      ]
    
    intercept <-
      fixef_m1["(Intercept)"] +
      ranef_intercepts[i]
    
    slope <-
      fixef_m1["um.rich"] +
      ranef_slopes[i]
    
    x_vals <- seq(
      min(current_data$um.rich),
      max(current_data$um.rich),
      length.out = 50
    )
    
    y_vals <-
      intercept +
      slope * x_vals
    
    lines(
      
      x_vals,
      y_vals,
      
      col = adjustcolor(
        "darkgreen",
        alpha.f = 0.5
      ), 
      
      lty = 2,
      
      lwd = 1.2
    )
  }
  
  # =====================================================
  # Línea fija
  # =====================================================
  
  lines(
    pred$x,
    pred$predicted,
    lwd = 3,
    col = "black"
  )
  
  # ---- Guardar métricas ----
  s <- summary(m1)
  sigma.temp <- s$sigma^2
  stddev_intercept <- attr(s$varcor$cond$ch.id, "stddev")["(Intercept)"]
  stddev_umrich    <- attr(s$varcor$cond$ch.id, "stddev")["um.rich"]
  r2_vals <- r.squaredGLMM(m1)
  
  pend.temp <- cbind(
    t(as.matrix(fixef_m1)),
    sigma.temp,
    r2_vals[1],
    r2_vals[2],
    stddev_intercept,
    stddev_umrich,
    s$coefficients$cond["um.rich", "Pr(>|z|)"]
  )
  
  colnames(pend.temp) <- c(
    "intercepto", "pendiente", "sigma", "R2 m", "R2 c",
    "stddev int", "stddev pend", "pvalor"
  )
  rownames(pend.temp) <- ii
  
  pendientes.lineal[[as.character(ii)]] <- pend.temp
}


pendientes.lineal <- do.call(rbind, pendientes.lineal)

pendientes_df <- do.call(rbind, lapply(names(pendientes.charcos), function(año) {
  data.frame(año = as.numeric(año), pendiente = pendientes.charcos[[año]])
}))

pendientes_df$year <- as.numeric(pendientes_df$año)





#######################################################
#### ------------ Con covariables ---------------- ####
#######################################################
find.x <- function(m, y_en, x_en, KK = 3) {
  out <- NULL
  require(bestglm)
  require(mgcv)
  
  X <- m[, x_en]
  y <- m[, y_en]
  Xy <- as.data.frame(cbind(X, y))
  
  # Model selection by AIC
  model <- bestglm(Xy, family = gaussian, IC = "AIC")
  model$BestModels -> models_tbl
  names(model$BestModel$coefficients)[-1] -> vars_best
  
  a.t <- NULL
  for (i in 1:nrow(models_tbl)) {
    vars_temp <- colnames(models_tbl)[which(models_tbl[i, ] == TRUE)]
    form <- as.formula(paste("y", "~", paste(vars_temp, collapse = "+")))
    fit <- lm(form, data = Xy)
    
    aa <- summary(fit)$coefficients
    a <- if ("um.rich" %in% rownames(aa)) aa["um.rich", 1] else NA_real_
    a.t <- rbind(a.t, cbind(row = i, a))
  }
  
  print(a.t)
  cbind(models_tbl, a = a.t[, 2]) -> models_tbl
  models_tbl[which(models_tbl$Criterion == median(models_tbl$Criterion)), ] -> temp
  
  vars_model <- colnames(temp)[which(temp[1, ] == TRUE)]
  out[[1]] <- vars_model
  out[[2]] <- models_tbl
  out[[3]] <- vars_best
  return(out)
}

# =========================================================
# Environmental variables
# =========================================================

library(visreg)

slopes_covariates <- NULL

par(
  mfrow = c(6,4),
  mar = c(5.5,6,3,2),
  cex.lab = 2
)

for (ii in 2005:2025) {
  
  # =====================================================
  # TEMPLATE
  # =====================================================
  
  slopes_cov_temp <- as.data.frame(
    matrix(
      NA,
      ncol = 18,
      nrow = 1
    )
  )
  
  colnames(slopes_cov_temp) <- c(
    
    "Intercept",
    "um.rich",
    "p.value",
    "R2",
    
    "Islands",
    "log.Area",
    "Mean.Depth",
    "CV.Depth",
    "Hydroperiod",
    "Degree",
    "log.Betweenness",
    
    "um.rich_Islands",
    "um.rich_log.Area",
    "um.rich_Mean.Depth",
    "um.rich_CV.Depth",
    "um.rich_Hydroperiod",
    "um.rich_Degree",
    "um.rich_log.Betweenness"
  )
  
  # =====================================================
  # DATOS
  # =====================================================
  
  bm.temp <- bm.nuevo[
    bm.nuevo[,1] == ii,
  ]
  
  df.um.ch <- bm.temp[, c(
    
    1,
    3,
    
    6,
    5,
    
    7:10,
    
    14:16,
    
    29,
    34,
    
    30:32,
    
    23:25
  )]
  
  colnames(df.um.ch) <- c(
    
    "um.year",
    "pond.id",
    
    "um.biom",
    "um.rich",
    
    "rain_sampling",
    "rain_annual",
    "temp_sampling",
    "temp_annual",
    
    "DM",
    "ddmm",
    "Shape",
    
    "Islands",
    "log.Area",
    
    "Mean.Depth",
    "Sd.Depth",
    "CV.Depth",
    
    "Hydroperiod",
    "Degree",
    "log.Betweenness"
  )
  
  # =====================================================
  # FILTRO
  # =====================================================
  
  df.um.ch <- df.um.ch[
    df.um.ch$um.rich > 0,
  ]
  
  df.um.ch <- na.omit(df.um.ch)
  
  df.um.ch$pond.id <- as.factor(
    df.um.ch$pond.id
  )
  
  title_txt <- paste(ii)
  
  # =====================================================
  # CONTROL
  # =====================================================
  
  if (
    nrow(df.um.ch) < 5 ||
    length(unique(df.um.ch$pond.id)) < 2
  ) {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "Insufficient data",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # =====================================================
  # ESCALADO SOLO COVARIABLES
  # =====================================================
  
  vars_scale <- c(
    
    "Islands",
    "log.Area",
    "Mean.Depth",
    "CV.Depth",
    "Hydroperiod",
    "Degree",
    "log.Betweenness"
  )
  
  df.um.ch[, vars_scale] <-
    scale(
      df.um.ch[, vars_scale]
    )
  
  # =====================================================
  # SELECCION VARIABLES
  # =====================================================
  
  variables <- find.x(
    m = df.um.ch,
    y_en = 3,
    x_en = c(
      4,
      12:19
    )
  )
  
  variables_model2 <- variables[[1]]
  
  # Force richness to be included in all models
  if (!("um.rich" %in% variables_model2)) {
    variables_model2 <- c("um.rich", variables_model2)
  }
  
  # =====================================================
  # SI NO HAY MODELO
  # =====================================================
  
  if(length(variables_model2) == 0){
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "No model selected",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # =====================================================
  # FORMULA
  # =====================================================
  
  form <- as.formula(
    paste(
      "um.biom ~",
      paste(
        variables_model2,
        collapse = " + "
      )
    )
  )
  
  # =====================================================
  # MODELO
  # =====================================================
  
  fit <- lm(
    form,
    data = df.um.ch
  )
  
  coefs <- summary(fit)$coefficients
  
  # =====================================================
  # SI richness NO ENTRA
  # =====================================================
  
  if(
    !("um.rich" %in% rownames(coefs))
  ){
    
    resid <- residuals(fit)
    
    df.resid <- cbind(
      resid,
      df.um.ch
    )
    
    fit.resid <- lm(
      resid ~ um.rich,
      data = df.resid
    )
    
    visreg(
      
      fit.resid,
      
      "um.rich",
      
      gg = FALSE,
      
      xlab = "Richness",
      
      ylab = expression(
        "Residual biomass"
      ),
      
      main = paste0(
        title_txt,
        " (Residual)"
      ),
      
      line.par = list(
        col = "black",
        lwd = 3
      ),
      
      fill.par = list(
        col = adjustcolor(
          "darkgreen",
          alpha.f = 0.2
        )
      ),
      
      points.par = list(
        pch = 19,
        cex = 0.8,
        col = adjustcolor(
          "black",
          alpha.f = 0.5
        )
      )
    )
    
    est <- summary(fit.resid)$coefficients[
      "um.rich",
      "Estimate"
    ]
    
    p_val <- summary(fit.resid)$coefficients[
      "um.rich",
      "Pr(>|t|)"
    ]
    
    texto <- paste0(
      "Residual"#,
      #round(est, 3)
    )
    
    mtext(
      texto,
      side = 3,
      line = 0,
      cex = 0.8
    )
    
  } else {
    
    # =================================================
    # MODELO NORMAL
    # =================================================
    
    est <- coefs[
      "um.rich",
      "Estimate"
    ]
    
    p_val <- coefs[
      "um.rich",
      "Pr(>|t|)"
    ]
    
    visreg(
      
      fit,
      
      "um.rich",
      
      gg = FALSE,
      
      xlab = "Richness",
      
      ylab = expression(
        "Biomass (g / 0.04 m"^2*")"
      ),
      
      main = title_txt,
      
      line.par = list(
        col = "black",
        lwd = 3
      ),
      
      fill.par = list(
        col = adjustcolor(
          "darkgreen",
          alpha.f = 0.2
        )
      ),
      
      points.par = list(
        pch = 19,
        cex = 0.8,
        col = adjustcolor(
          "black",
          alpha.f = 0.5
        )
      )
    )
    
    texto <- paste0(
      
      "Lineal | Est. = ",
      
      round(est, 3)
    )
    
    mtext(
      texto,
      side = 3,
      line = 0,
      cex = 0.8
    )
  }
  
  # =====================================================
  # GUARDAR
  # =====================================================
  
  if(
    "(Intercept)" %in% rownames(coefs)
  ){
    
    slopes_cov_temp$Intercept <-
      coefs["(Intercept)",1]
  }
  
  slopes_cov_temp$um.rich <- est
  
  slopes_cov_temp$p.value <- p_val
  
  slopes_cov_temp$R2 <-
    summary(fit)$adj.r.squared
  
  # =====================================================
  # COVARIABLES
  # =====================================================
  
  covariates <- c(
    
    "Islands",
    "log.Area",
    "Mean.Depth",
    "CV.Depth",
    "Hydroperiod",
    "Degree",
    "log.Betweenness"
  )
  
  for(v in covariates){
    
    if(
      (v %in% variables_model2) &&
      (v %in% rownames(coefs))
    ){
      
      slopes_cov_temp[[v]] <-
        coefs[v,1]
    }
  }
  
  # =====================================================
  # APPEND
  # =====================================================
  
  slopes_covariates <- rbind(
    slopes_covariates,
    slopes_cov_temp
  )
}

pendientes.covariables <- slopes_covariates
m_cov<-(lm(pendientes.covariables[,2]~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_cov)








library(glmmTMB)
library(MuMIn)
library(viridisLite)
library(visreg)

par(mfrow=c(6,4), mar=c(5.5,6,3,2), cex.lab = 2)

slopes <- list()
sigma_vec  <- numeric()

for (ii in 2005:2025) {
  
  # ---- 1) Filter year ----
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  
  # ---- 2) Build df ----
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:26)]
  
  colnames(df.um.ch) <- c(
    "um.year", "pond.id", "um.biom", "um.rich",
    "rain_sampling", "rain_annual", "temp_sampling", "temp_annual",
    "DM", "ddmm", "Shape", "Islands", "log.Area", "log.Volume",
    "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree",
    "log.Betweenness", "Closenness"
  )
  
  df.um.ch <- na.omit(df.um.ch)
  
  df.um.ch$pond.id <- as.factor(df.um.ch$pond.id)
  
  title_txt <- paste(ii)
  
  # ---- 3) Control ----
  if (
    nrow(df.um.ch) < 5 ||
    length(unique(df.um.ch$pond.id)) < 2
  ) {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "Insufficient data",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # ---- 4) Candidate models ----
  
  m1 <- try(
    glmmTMB(
      um.biom ~ um.rich + (1 + um.rich | pond.id),
      family = gaussian(),
      data = df.um.ch
    ),
    silent = TRUE
  )
  
  m12 <- try(
    glmmTMB(
      um.biom ~ um.rich + I(um.rich^2) +
        (1 + um.rich | pond.id),
      family = gaussian(),
      data = df.um.ch
    ),
    silent = TRUE
  )
  
  # ---- 5) AIC ----
  
  AIC_m1 <- if (
    !inherits(m1, "try-error")
  ) {
    AIC(m1)
  } else {
    NA
  }
  
  AIC_m12 <- if (
    !inherits(m12, "try-error")
  ) {
    AIC(m12)
  } else {
    NA
  }
  
  deltaAIC <- AIC_m1 - AIC_m12
  
  # ---- 6) Select model ----
  
  if (
    !inherits(m1, "try-error") &&
    !inherits(m12, "try-error")
  ) {
    
    if (
      !is.na(deltaAIC) &&
      deltaAIC > 2
    ) {
      
      best_model <- m12
      best_label <- "Quadratic"
      
    } else {
      
      best_model <- m1
      best_label <- "Linear"
    }
    
  } else if (
    !inherits(m1, "try-error")
  ) {
    
    best_model <- m1
    best_label <- "Linear"
    
  } else if (
    !inherits(m12, "try-error")
  ) {
    
    best_model <- m12
    best_label <- "Quadratic"
    
  } else {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "No model converged",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # ---- 7) Summary ----
  
  s <- summary(best_model)
  
  fixef_best <- fixef(best_model)$cond
  
  g0 <- fixef_best["(Intercept)"]
  g1 <- fixef_best["um.rich"]
  
  g2 <- if (
    "I(um.rich^2)" %in% names(fixef_best)
  ) {
    fixef_best["I(um.rich^2)"]
  } else {
    0
  }
  
  # ---- 8) VISREG ----
  
  visreg(
    best_model,
    "um.rich",
    
    gg = FALSE,
    
    re.form = NA,
    
    #partial = FALSE,
    
    xlab = "Richness",
    
    ylab = expression(
      "Biomass (g / 0.04 m"^2*")"
    ),
    
    main = paste0(
      title_txt#,
      # " (",
      # best_label,
      # ")"
    ),
    
    line.par = list(
      col = "black",
      lwd = 3
    ),
    
    fill.par = list(
      col = adjustcolor(
        "darkgreen",
        alpha.f = 0.2
      )
    ),
    
    points.par = list(
      pch = 19,
      cex = 0.8,
      col = adjustcolor(
        "black",
        alpha.f = 0.5
      )
    )
  )
  
  # ---- 9) Pond lines ----
  
  ran_int <- ranef(best_model)$cond$pond.id[, "(Intercept)"]
  ran_slp <- ranef(best_model)$cond$pond.id[, "um.rich"]
  
  unique_ids <- levels(df.um.ch$pond.id)
  
  colors <- viridis(
    length(unique_ids),
    option = "D"
  )
  
  names(colors) <- unique_ids
  
  for (k in seq_along(unique_ids)) {
    
    idk <- unique_ids[k]
    
    current_data <- df.um.ch[
      df.um.ch$pond.id == idk,
    ]
    
    x_min <- min(
      current_data$um.rich,
      na.rm = TRUE
    )
    
    x_max <- max(
      current_data$um.rich,
      na.rm = TRUE
    )
    
    intercept <- g0 + ran_int[idk]
    
    slope1 <- g1 + ran_slp[idk]
    
    slope2 <- g2
    
    curve(
      intercept +
        slope1 * x +
        slope2 * x^2,
      
      from = x_min,
      to = x_max,
      
      add = TRUE,
      
      col = colors[idk],
      lwd = 1,
      lty = 2
    )
  }
  
  # ---- 10) Annotation ----
  
  texto <- ifelse(
    
    best_label == "Linear",
    
    paste0(
      "Linear | Est. = ",
      round(g1, 3)
    ),
    
    paste0("Quadratic | ΔAIC = ",
           round(deltaAIC, 2))
  )
  
  mtext(
    texto,
    side = 3,
    line = 0,
    cex = 0.8
  )
  
  # ---- 11) Metrics ----
  
  p_umrich <- NA_real_
  
  if (
    "um.rich" %in% rownames(
      s$coefficients$cond
    )
  ) {
    
    p_umrich <- s$coefficients$cond[
      "um.rich",
      "Pr(>|z|)"
    ]
  }
  
  sigma2 <- s$sigma^2
  
  sigma_vec <- c(
    sigma_vec,
    sigma2
  )
  
  r2m <- r2c <- NA_real_
  
  r2_try <- try(
    r.squaredGLMM(best_model),
    silent = TRUE
  )
  
  if (!inherits(r2_try, "try-error")) {
    
    r2m <- as.numeric(r2_try[1])
    r2c <- as.numeric(r2_try[2])
  }
  
  sd_int <- sd_slp <- NA_real_
  
  vc <- try(
    summary(best_model)$varcor$cond$pond.id,
    silent = TRUE
  )
  
  if (!inherits(vc, "try-error")) {
    
    sd_int <- attr(
      vc,
      "stddev"
    )["(Intercept)"]
    
    sd_slp <- attr(
      vc,
      "stddev"
    )["um.rich"]
  }
  
  # ---- 12) Save ----
  
  row_out <- data.frame(
    
    year = ii,
    
    best_label = best_label,
    
    AIC_m1 = AIC_m1,
    AIC_m12 = AIC_m12,
    
    deltaAIC = deltaAIC,
    
    intercept = unname(g0),
    
    slope = unname(g1),
    
    quadratic = unname(g2),
    
    deriv_S1 = unname(
      g1 + 2*g2*1
    ),
    
    deriv_S2 = unname(
      g1 + 2*g2*2
    ),
    
    sigma = sigma2,
    
    R2_m = r2m,
    R2_c = r2c,
    
    p_value = p_umrich,
    
    stddev_int = sd_int,
    
    stddev_slope = sd_slp
  )
  
  slopes[[as.character(ii)]] <- row_out
}

slopes_df <- do.call(rbind, slopes)

slopes_df$P_quad <- ifelse(slopes_df$best_label=="Quadratic",1,0)

m_quad<-(glm(slopes_df$P_quad~nn.1$R2.unexp))
summary(m_quad)
m_d1<-(lm(slopes_df$deriv_S1~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d1)
m_d2<-(lm(slopes_df$deriv_S2~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d2)























library(glmmTMB)
library(MuMIn)
library(viridisLite)
library(visreg)

par(mfrow=c(6,4), mar=c(5.5,6,3,2), cex.lab = 2)

slopes <- list()
sigma_vec  <- numeric()

for (ii in 2005:2025) {
  
  # ---- 1) Filter year ----
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  
  # ---- 2) Build df ----
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:26)]
  
  colnames(df.um.ch) <- c(
    "um.year", "pond.id", "um.biom", "um.rich",
    "rain_sampling", "rain_annual", "temp_sampling", "temp_annual",
    "DM", "ddmm", "Shape", "Islands", "log.Area", "log.Volume",
    "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree",
    "log.Betweenness", "Closenness"
  )
  
  # ---- FILTRO RIQUEZA 1–9 ----
  
  df.um.ch <- df.um.ch[
    df.um.ch$um.rich >= 1 &
      df.um.ch$um.rich <= 9,
  ]
  
  df.um.ch <- na.omit(df.um.ch)
  
  df.um.ch$pond.id <- as.factor(df.um.ch$pond.id)
  
  title_txt <- paste(ii)
  
  # ---- 3) Control ----
  if (
    nrow(df.um.ch) < 5 ||
    length(unique(df.um.ch$pond.id)) < 2
  ) {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "Insufficient data",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # ---- 4) Candidate models ----
  
  m1 <- try(
    glmmTMB(
      um.biom ~ um.rich + (1 + um.rich | pond.id),
      family = gaussian(),
      data = df.um.ch
    ),
    silent = TRUE
  )
  
  m12 <- try(
    glmmTMB(
      um.biom ~ um.rich + I(um.rich^2) +
        (1 + um.rich | pond.id),
      family = gaussian(),
      data = df.um.ch
    ),
    silent = TRUE
  )
  
  # ---- 5) AIC ----
  
  AIC_m1 <- if (
    !inherits(m1, "try-error")
  ) {
    AIC(m1)
  } else {
    NA
  }
  
  AIC_m12 <- if (
    !inherits(m12, "try-error")
  ) {
    AIC(m12)
  } else {
    NA
  }
  
  deltaAIC <- AIC_m1 - AIC_m12
  
  # ---- 6) Select model ----
  
  if (
    !inherits(m1, "try-error") &&
    !inherits(m12, "try-error")
  ) {
    
    if (
      !is.na(deltaAIC) &&
      deltaAIC > 2
    ) {
      
      best_model <- m12
      best_label <- "Quadratic"
      
    } else {
      
      best_model <- m1
      best_label <- "Linear"
    }
    
  } else if (
    !inherits(m1, "try-error")
  ) {
    
    best_model <- m1
    best_label <- "Linear"
    
  } else if (
    !inherits(m12, "try-error")
  ) {
    
    best_model <- m12
    best_label <- "Quadratic"
    
  } else {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "No model converged",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # ---- 7) Summary ----
  
  s <- summary(best_model)
  
  fixef_best <- fixef(best_model)$cond
  
  g0 <- fixef_best["(Intercept)"]
  g1 <- fixef_best["um.rich"]
  
  g2 <- if (
    "I(um.rich^2)" %in% names(fixef_best)
  ) {
    fixef_best["I(um.rich^2)"]
  } else {
    0
  }
  
  # ---- 8) VISREG ----
  
  visreg(
    best_model,
    "um.rich",
    
    gg = FALSE,
    re.form = NA,
    #partial = FALSE,
    
    xlab = "Richness",
    
    ylab = expression(
      "Biomass (g / 0.04 m"^2*")"
    ),
    
    main = paste0(
      title_txt#,
      # " (",
      # best_label,
      # ")"
    ),
    
    line.par = list(
      col = "black",
      lwd = 3
    ),
    
    fill.par = list(
      col = adjustcolor(
        "darkgreen",
        alpha.f = 0.2
      )
    ),
    
    points.par = list(
      pch = 19,
      cex = 0.8,
      col = adjustcolor(
        "black",
        alpha.f = 0.5
      )
    )
  )
  
  # ---- 9) Pond lines ----
  
  ran_int <- ranef(best_model)$cond$pond.id[, "(Intercept)"]
  ran_slp <- ranef(best_model)$cond$pond.id[, "um.rich"]
  
  unique_ids <- levels(df.um.ch$pond.id)
  
  colors <- viridis(
    length(unique_ids),
    option = "D"
  )
  
  names(colors) <- unique_ids
  
  for (k in seq_along(unique_ids)) {
    
    idk <- unique_ids[k]
    
    current_data <- df.um.ch[
      df.um.ch$pond.id == idk,
    ]
    
    x_min <- min(
      current_data$um.rich,
      na.rm = TRUE
    )
    
    x_max <- max(
      current_data$um.rich,
      na.rm = TRUE
    )
    
    intercept <- g0 + ran_int[idk]
    
    slope1 <- g1 + ran_slp[idk]
    
    slope2 <- g2
    
    curve(
      intercept +
        slope1 * x +
        slope2 * x^2,
      
      from = x_min,
      to = x_max,
      
      add = TRUE,
      
      col = colors[idk],
      lwd = 1,
      lty = 2
    )
  }
  
  # ---- 10) Annotation ----
  
  texto <- ifelse(
    
    best_label == "Linear",
    
    paste0(
      "Linear | Est. = ",
      round(g1, 3)
    ),
    
    paste0("Quadratic | ΔAIC = ",
           round(deltaAIC, 2))
  )
  
  mtext(
    texto,
    side = 3,
    line = 0,
    cex = 0.8
  )
  
  # ---- 11) Metrics ----
  
  p_umrich <- NA_real_
  
  if (
    "um.rich" %in% rownames(
      s$coefficients$cond
    )
  ) {
    
    p_umrich <- s$coefficients$cond[
      "um.rich",
      "Pr(>|z|)"
    ]
  }
  
  sigma2 <- s$sigma^2
  
  sigma_vec <- c(
    sigma_vec,
    sigma2
  )
  
  r2m <- r2c <- NA_real_
  
  r2_try <- try(
    r.squaredGLMM(best_model),
    silent = TRUE
  )
  
  if (!inherits(r2_try, "try-error")) {
    
    r2m <- as.numeric(r2_try[1])
    r2c <- as.numeric(r2_try[2])
  }
  
  sd_int <- sd_slp <- NA_real_
  
  vc <- try(
    summary(best_model)$varcor$cond$pond.id,
    silent = TRUE
  )
  
  if (!inherits(vc, "try-error")) {
    
    sd_int <- attr(
      vc,
      "stddev"
    )["(Intercept)"]
    
    sd_slp <- attr(
      vc,
      "stddev"
    )["um.rich"]
  }
  
  # ---- 12) Save ----
  
  row_out <- data.frame(
    
    year = ii,
    
    best_label = best_label,
    
    AIC_m1 = AIC_m1,
    AIC_m12 = AIC_m12,
    
    deltaAIC = deltaAIC,
    
    intercept = unname(g0),
    
    slope = unname(g1),
    
    quadratic = unname(g2),
    
    deriv_S1 = unname(
      g1 + 2*g2*1
    ),
    
    deriv_S2 = unname(
      g1 + 2*g2*2
    ),
    
    sigma = sigma2,
    
    R2_m = r2m,
    R2_c = r2c,
    
    p_value = p_umrich,
    
    stddev_int = sd_int,
    
    stddev_slope = sd_slp
  )
  
  slopes[[as.character(ii)]] <- row_out
}

slopes_df_sat <- do.call(rbind, slopes)

slopes_df_sat$P_quad <- ifelse(slopes_df_sat$best_label=="Quadratic",1,0)

m_quad_sat<-(glm(slopes_df_sat$P_quad~nn.1$R2.unexp))
summary(m_quad_sat)
m_d1_sat<-(lm(slopes_df_sat$deriv_S1~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d1_sat)
m_d2_sat<-(lm(slopes_df_sat$deriv_S2~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d2_sat)
















library(glmmTMB)
library(MuMIn)
library(viridisLite)
library(visreg)

par(
  mfrow = c(6,4),
  mar = c(5.5,6,3,2),
  cex.lab = 2
)

slopes <- list()
sigma_vec <- numeric()

for (ii in 2005:2025) {
  
  # ---- 1) Filter year ----
  
  bm.temp <- bm.nuevo[
    bm.nuevo[,1] == ii,
  ]
  
  # ---- 2) Build df ----
  
  df.um.ch <- bm.temp[, c(
    
    1,
    3,
    
    6,
    5,
    
    7:10,
    
    14:26
  )]
  
  colnames(df.um.ch) <- c(
    
    "um.year",
    "pond.id",
    
    "um.biom",
    "um.rich",
    
    "rain_sampling",
    "rain_annual",
    "temp_sampling",
    "temp_annual",
    
    "DM",
    "ddmm",
    "Shape",
    
    "Islands",
    "log.Area",
    "log.Volume",
    
    "Mean.Depth",
    "Sd.Depth",
    "CV.Depth",
    
    "Hydroperiod",
    "Degree",
    "log.Betweenness",
    "Closenness"
  )
  
  # ---- FILTRO RIQUEZA 1–9 ----
  
  df.um.ch <- df.um.ch[
    df.um.ch$um.rich >= 1 &
      df.um.ch$um.rich <= 9,
  ]
  
  df.um.ch <- na.omit(df.um.ch)
  
  df.um.ch$pond.id <- as.factor(
    df.um.ch$pond.id
  )
  
  title_txt <- paste(ii)
  
  # ---- 3) Control ----
  
  if (
    nrow(df.um.ch) < 5 ||
    length(unique(df.um.ch$pond.id)) < 2
  ) {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "Insufficient data",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # ---- 4) Linear model ----
  
  best_model <- try(
    
    glmmTMB(
      
      um.biom ~ um.rich +
        (1 + um.rich | pond.id),
      
      family = gaussian(),
      
      data = df.um.ch
    ),
    
    silent = TRUE
  )
  
  if (
    inherits(best_model, "try-error")
  ) {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "Model did not converge",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # ---- 5) Summary ----
  
  s <- summary(best_model)
  
  fixef_best <- fixef(best_model)$cond
  
  g0 <- fixef_best["(Intercept)"]
  
  g1 <- fixef_best["um.rich"]
  
  # ---- 6) VISREG ----
  
  visreg(
    
    best_model,
    
    "um.rich",
    
    gg = FALSE,
    
    re.form = NA,
    
    xlab = "Richness",
    
    ylab = expression(
      "Biomass (g / 0.04 m"^2*")"
    ),
    
    main = title_txt,
    
    line.par = list(
      col = "black",
      lwd = 3
    ),
    
    fill.par = list(
      col = adjustcolor(
        "darkgreen",
        alpha.f = 0.2
      )
    ),
    
    points.par = list(
      pch = 19,
      cex = 0.8,
      col = adjustcolor(
        "black",
        alpha.f = 0.5
      )
    )
  )
  
  # ---- 7) Pond lines ----
  
  ran_int <- ranef(best_model)$cond$pond.id[, "(Intercept)"]
  
  ran_slp <- ranef(best_model)$cond$pond.id[, "um.rich"]
  
  unique_ids <- levels(
    df.um.ch$pond.id
  )
  
  colors <- viridis(
    length(unique_ids),
    option = "D"
  )
  
  names(colors) <- unique_ids
  
  for (k in seq_along(unique_ids)) {
    
    idk <- unique_ids[k]
    
    current_data <- df.um.ch[
      df.um.ch$pond.id == idk,
    ]
    
    x_min <- min(
      current_data$um.rich,
      na.rm = TRUE
    )
    
    x_max <- max(
      current_data$um.rich,
      na.rm = TRUE
    )
    
    intercept <- g0 + ran_int[idk]
    
    slope1 <- g1 + ran_slp[idk]
    
    curve(
      intercept +
        slope1 * x,
      
      from = x_min,
      to = x_max,
      
      add = TRUE,
      
      col = colors[idk],
      lwd = 1,
      lty = 2
    )
  }
  
  # ---- 8) Annotation ----
  
  texto <- paste0(
    
    "Linear | Est. = ",
    
    round(g1, 3)
  )
  
  mtext(
    texto,
    side = 3,
    line = 0,
    cex = 0.8
  )
  
  # ---- 9) Metrics ----
  
  p_umrich <- s$coefficients$cond[
    "um.rich",
    "Pr(>|z|)"
  ]
  
  sigma2 <- s$sigma^2
  
  sigma_vec <- c(
    sigma_vec,
    sigma2
  )
  
  r2m <- r2c <- NA_real_
  
  r2_try <- try(
    r.squaredGLMM(best_model),
    silent = TRUE
  )
  
  if (
    !inherits(r2_try, "try-error")
  ) {
    
    r2m <- as.numeric(r2_try[1])
    
    r2c <- as.numeric(r2_try[2])
  }
  
  sd_int <- sd_slp <- NA_real_
  
  vc <- try(
    summary(best_model)$varcor$cond$pond.id,
    silent = TRUE
  )
  
  if (
    !inherits(vc, "try-error")
  ) {
    
    sd_int <- attr(
      vc,
      "stddev"
    )["(Intercept)"]
    
    sd_slp <- attr(
      vc,
      "stddev"
    )["um.rich"]
  }
  
  # ---- 10) Save ----
  
  row_out <- data.frame(
    
    year = ii,
    
    intercept = unname(g0),
    
    slope = unname(g1),
    
    deriv_S1 = unname(g1),
    
    deriv_S2 = unname(g1),
    
    sigma = sigma2,
    
    R2_m = r2m,
    
    R2_c = r2c,
    
    p_value = p_umrich,
    
    stddev_int = sd_int,
    
    stddev_slope = sd_slp
  )
  
  slopes[[as.character(ii)]] <- row_out
}

slopes_df_sat_linear <- do.call(
  rbind,
  slopes)

m_lineal_sat <- lm(slopes_df_sat_linear$deriv_S2 ~
                     nn.1$R2.exp +
                     nn.1$lluvia_anual)
summary(m_lineal_sat)



######################################################################
######################################################################
################### ------ quantile 90 ---------- ####################
######################################################################


library(quantreg)
library(visreg)

slopes_covariates_q90 <- NULL

par(
  mfrow = c(6,4),
  mar = c(5.5,6,3,2),
  cex.lab = 2
)

for (ii in 2005:2025) {
  
  # =====================================================
  # TEMPLATE
  # =====================================================
  
  slopes_cov_temp <- as.data.frame(
    matrix(
      NA,
      ncol = 22,
      nrow = 1
    )
  )
  
  colnames(slopes_cov_temp) <- c(
    
    "Intercept",
    "um.rich",
    "I.um.rich2",
    
    "p.value",
    "AIC",
    "modelo",
    
    "deriv_S1",
    "deriv_S2",
    
    "Islands",
    "log.Area",
    "Mean.Depth",
    "CV.Depth",
    "Hydroperiod",
    "Degree",
    "log.Betweenness",
    
    "um.rich_Islands",
    "um.rich_log.Area",
    "um.rich_Mean.Depth",
    "um.rich_CV.Depth",
    "um.rich_Hydroperiod",
    "um.rich_Degree",
    "um.rich_log.Betweenness"
  )
  
  # =====================================================
  # DATOS
  # =====================================================
  
  bm.temp <- bm.nuevo[
    bm.nuevo[,1] == ii,
  ]
  
  df.um.ch <- bm.temp[, c(
    
    1,
    3,
    
    6,
    5,
    
    7:10,
    
    14:16,
    
    29,
    34,
    
    30:32,
    
    23:25
  )]
  
  colnames(df.um.ch) <- c(
    
    "um.year",
    "pond.id",
    
    "um.biom",
    "um.rich",
    
    "rain_sampling",
    "rain_annual",
    "temp_sampling",
    "temp_annual",
    
    "DM",
    "ddmm",
    "Shape",
    
    "Islands",
    "log.Area",
    
    "Mean.Depth",
    "Sd.Depth",
    "CV.Depth",
    
    "Hydroperiod",
    "Degree",
    "log.Betweenness"
  )
  
  # =====================================================
  # FILTRO
  # =====================================================
  
  df.um.ch <- df.um.ch[
    df.um.ch$um.rich > 0,
  ]
  
  df.um.ch <- na.omit(df.um.ch)
  
  df.um.ch$pond.id <- as.factor(
    df.um.ch$pond.id
  )
  
  title_txt <- paste(ii)
  
  # =====================================================
  # CONTROL
  # =====================================================
  
  if (
    nrow(df.um.ch) < 5 ||
    length(unique(df.um.ch$pond.id)) < 2
  ) {
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "Insufficient data",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # =====================================================
  # ESCALADO SOLO COVARIABLES
  # =====================================================
  
  vars_scale <- c(
    
    "Islands",
    "log.Area",
    "Mean.Depth",
    "CV.Depth",
    "Hydroperiod",
    "Degree",
    "log.Betweenness"
  )
  
  df.um.ch[, vars_scale] <-
    scale(
      df.um.ch[, vars_scale]
    )
  
  # =====================================================
  # SELECCION VARIABLES
  # =====================================================
  
  variables <- find.x(
    m = df.um.ch,
    y_en = 3,
    x_en = c(
      4,
      12:18
    )
  )
  
  variables_model2 <- variables[[1]]
  
  # Force richness to be included in all models
  if (!("um.rich" %in% variables_model2)) {
    variables_model2 <- c("um.rich", variables_model2)
  }
  
  # =====================================================
  # SI NO HAY MODELO
  # =====================================================
  
  if(length(variables_model2) == 0){
    
    plot.new()
    
    title(main = title_txt)
    
    mtext(
      "No model selected",
      side = 3,
      line = -2
    )
    
    next
  }
  
  # =====================================================
  # FORMULAS
  # =====================================================
  
  vars_lin <- variables_model2
  
  vars_quad <- variables_model2
  
  if(!("I(um.rich^2)" %in% vars_quad)){
    
    if("um.rich" %in% vars_quad){
      
      idx <- which(vars_quad == "um.rich")
      
      vars_quad <- append(
        vars_quad,
        "I(um.rich^2)",
        after = idx
      )
    }
  }
  
  form_lin <- as.formula(
    paste(
      "um.biom ~",
      paste(vars_lin, collapse = " + ")
    )
  )
  
  form_quad <- as.formula(
    paste(
      "um.biom ~",
      paste(vars_quad, collapse = " + ")
    )
  )
  
  # =====================================================
  # MODELOS QUANTILE
  # =====================================================
  
  fit_lin <- rq(
    form_lin,
    tau = 0.9,
    data = df.um.ch
  )
  
  fit_quad <- rq(
    form_quad,
    tau = 0.9,
    data = df.um.ch
  )
  
  # =====================================================
  # AIC
  # =====================================================
  
  AIC_lin <- AIC(fit_lin)
  
  AIC_quad <- AIC(fit_quad)
  
  deltaAIC <- AIC_lin - AIC_quad
  
  # =====================================================
  # SELECCION
  # =====================================================
  
  if(
    "um.rich" %in% variables_model2 &&
    deltaAIC > 2
  ){
    
    fit <- fit_quad
    
    modelo <- "Quadratic"
    
  } else {
    
    fit <- fit_lin
    
    modelo <- "Linear"
  }
  
  # =====================================================
  # SI richness NO ENTRA
  # =====================================================
  
  if(
    !("um.rich" %in% names(coef(fit)))
  ){
    
    resid <- residuals(fit)
    
    df.resid <- cbind(
      resid,
      df.um.ch
    )
    
    fit.resid.lin <- rq(
      resid ~ um.rich,
      tau = 0.9,
      data = df.resid
    )
    
    fit.resid.quad <- rq(
      resid ~ um.rich + I(um.rich^2),
      tau = 0.9,
      data = df.resid
    )
    
    AIC.resid.lin <- AIC(fit.resid.lin)
    
    AIC.resid.quad <- AIC(fit.resid.quad)
    
    deltaAIC.resid <- AIC.resid.lin - AIC.resid.quad
    
    if(deltaAIC.resid > 2){
      
      fit.resid <- fit.resid.quad
      
      modelo <- "Residual Quadratic"
      
    } else {
      
      fit.resid <- fit.resid.lin
      
      modelo <- "Residual Linear"
    }
    
    est <- coef(fit.resid)["um.rich"]
    
    coef_quad <- ifelse(
      "I(um.rich^2)" %in% names(coef(fit.resid)),
      coef(fit.resid)["I(um.rich^2)"],
      0
    )
    
    p_val <- summary(
      fit.resid,
      se = "nid"
    )$coefficients[
      "um.rich",
      "Pr(>|t|)"
    ]
    
    deriv_S1 <- est + 2*coef_quad*1
    
    deriv_S2 <- est + 2*coef_quad*2
    
    visreg(
      
      fit.resid,
      
      "um.rich",
      
      gg = FALSE,
      
      xlab = "Richness",
      
      ylab = expression(
        "Residual biomass"
      ),
      
      main = paste0(
        title_txt#,
        #" (Residual)"
      ),
      
      line.par = list(
        col = "black",
        lwd = 3
      ),
      
      fill.par = list(
        col = adjustcolor(
          "darkgreen",
          alpha.f = 0.2
        )
      ),
      
      points.par = list(
        pch = 19,
        cex = 0.8,
        col = adjustcolor(
          "black",
          alpha.f = 0.5
        )
      )
    )
    
    texto <- ifelse(
      
      modelo == "Residual Linear",
      
      paste0(
        "Residual Linear | Est. = ",
        round(est, 3)
      ),
      
      paste0(
        "Residual Quadratic | S1 = ",
        round(deriv_S1, 3),
        " | S2 = ",
        round(deriv_S2, 3)
      )
    )
    
    mtext(
      texto,
      side = 3,
      line = 0,
      cex = 0.8
    )
    
  } else {
    
    # =================================================
    # MODELO NORMAL
    # =================================================
    
    est <- coef(fit)["um.rich"]
    
    coef_quad <- ifelse(
      "I(um.rich^2)" %in% names(coef(fit)),
      coef(fit)["I(um.rich^2)"],
      0
    )
    
    p_val <- summary(
      fit,
      se = "nid"
    )$coefficients[
      "um.rich",
      "Pr(>|t|)"
    ]
    
    deriv_S1 <- est + 2*coef_quad*1
    
    deriv_S2 <- est + 2*coef_quad*2
    
    visreg(
      
      fit,
      
      "um.rich",
      
      gg = FALSE,
      
      xlab = "Richness",
      
      ylab = expression(
        "Biomass (g / 0.04 m"^2*")"
      ),
      
      main = paste0(
        title_txt#,
        #" (",
        #modelo,
        #")"
      ),
      
      line.par = list(
        col = "black",
        lwd = 3
      ),
      
      fill.par = list(
        col = adjustcolor(
          "darkgreen",
          alpha.f = 0.2
        )
      ),
      
      points.par = list(
        pch = 19,
        cex = 0.8,
        col = adjustcolor(
          "black",
          alpha.f = 0.5
        )
      )
    )
    
    texto <- ifelse(
      
      modelo == "Linear",
      
      paste0(
        "Linear | Est. = ",
        round(est, 3)
      ),
      
      paste0(
        "Quadratic | S1 = ",
        round(deriv_S1, 3),
        " | S2 = ",
        round(deriv_S2, 3)
      )
    )
    
    mtext(
      texto,
      side = 3,
      line = 0,
      cex = 0.8
    )
  }
  
  # =====================================================
  # GUARDAR
  # =====================================================
  
  if(
    "(Intercept)" %in% names(coef(fit))
  ){
    
    slopes_cov_temp$Intercept <-
      coef(fit)["(Intercept)"]
  }
  
  slopes_cov_temp$um.rich <- est
  
  slopes_cov_temp$I.um.rich2 <- coef_quad
  
  slopes_cov_temp$p.value <- p_val
  
  slopes_cov_temp$AIC <- AIC(fit)
  
  slopes_cov_temp$modelo <- modelo
  
  slopes_cov_temp$deriv_S1 <- deriv_S1
  
  slopes_cov_temp$deriv_S2 <- deriv_S2
  
  # =====================================================
  # COVARIABLES
  # =====================================================
  
  covariates <- c(
    
    "Islands",
    "log.Area",
    "Mean.Depth",
    "CV.Depth",
    "Hydroperiod",
    "Degree",
    "log.Betweenness"
  )
  
  for(v in covariates){
    
    if(
      (v %in% variables_model2) &&
      (v %in% names(coef(fit)))
    ){
      
      slopes_cov_temp[[v]] <-
        coef(fit)[v]
    }
  }
  
  # =====================================================
  # APPEND
  # =====================================================
  
  slopes_covariates_q90 <- rbind(
    slopes_covariates_q90,
    slopes_cov_temp
  )
}

slopes_covariates_q90$P_quad <-
  as.numeric(
    slopes_covariates_q90$modelo %in%
      c("Quadratic", "Residual Quadratic"))
slopes_covariates_q90$P_quad-> nn.1$P_quad
m_quad_q90<-(glm(P_quad~R2.exp, data=nn.1, family=binomial()))
summary(m_quad_q90)


slopes_covariates_q90$deriv_S1-> nn.1$BEF
m_q90_d1<- lm(BEF~R2.exp, data=nn.1)
summary(m_q90_d1)

slopes_covariates_q90$deriv_S2-> nn.1$BEF
m_q90_d2<- lm(BEF~R2.exp, data=nn.1)
summary(m_q90_d2)



# ---------- MODEL VISUALIZATION ---------- #

m_cov<-(lm(pendientes.covariables[,2]~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_cov)

m_quad<-(glm(slopes_df$P_quad~nn.1$R2.exp))
summary(m_quad)
m_d1<-(lm(slopes_df$deriv_S1~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d1)
m_d2<-(lm(slopes_df$deriv_S2~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d2)

m_lineal_sat <- lm(slopes_df_sat_linear$deriv_S2 ~
                     nn.1$R2.exp + nn.1$lluvia_anual)
summary(m_lineal_sat)

m_quad_sat<-(glm(slopes_df_sat$P_quad~nn.1$R2.exp))
summary(m_quad_sat)
m_d1_sat<-(lm(slopes_df_sat$deriv_S1~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d1_sat)
m_d2_sat<-(lm(slopes_df_sat$deriv_S2~nn.1$R2.exp+nn.1$lluvia_anual))
summary(m_d2_sat)


m_quad_q90<-(glm(slopes_covariates_q90$P_quad~nn.1$R2.exp, family=binomial()))
summary(m_quad_q90)
slopes_covariates_q90$deriv_S1-> nn.1$BEF
m_q90_d1<- lm(slopes_covariates_q90$deriv_S1~nn.1$R2.exp)
summary(m_q90_d1)
m_q90_d2<- lm(slopes_covariates_q90$deriv_S2~nn.1$R2.exp)
summary(m_q90_d2)



# ============
#
# Fig. S6
#
# ============

ax_col <- "#2E2E2E"
col_bef <- "#3A8EC1"
pt_bef  <- adjustcolor(col_bef, alpha.f = 0.6)
double_curve <- function(fun, from, to, col, halo_col = ax_col,
                         halo_lwd = 6, col_lwd = 4.5, ...) {
  curve(fun, from = from, to = to, add = TRUE, col = halo_col, lwd = halo_lwd, ...)
  curve(fun, from = from, to = to, add = TRUE, col = col,      lwd = col_lwd,  ...)
}


# ---- PANEL A: Model with environmental cov. --- #
pendientes.covariables[,2]-> nn.1$BEF
m_cov<- lm(BEF ~ R2.exp + lluvia_anual, data = nn.1)
p <- coefficients(m_cov)
par(mfrow = c(2,2), mar = c(9, 9, 4, 3), mgp = c(3.5, 1, 0),
    cex.lab = 2, cex.axis = 2)

BEF.neutral <- with(nn.1,
                    BEF - p["lluvia_anual"] * lluvia_anual +
                      mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))
plot(BEF.neutral ~ nn.1$R2.exp,
     xlab = "",
     ylab = "Richness-Biomass relationship (BEF)",
     bty = "l",
     pch = 19, col = pt_bef, cex = 3,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * nn.1$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_bef)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")")
  )), side = 1, line = 7, cex = 2, col = ax_col)
text(x = 0.3, y = 1.2, "F(2,18): 8.2\np: 0.003\nR²: 0.47",
     cex = 2.0, pos = 1, col = ax_col)
mtext("A", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.5, col = ax_col)
x.seq <- seq(
  min(nn.1$R2.exp, na.rm = TRUE),
  max(nn.1$R2.exp, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  R2.exp = x.seq,
  lluvia_anual = mean(nn.1$lluvia_anual, na.rm = TRUE))
pred <- predict(
  m_cov,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["lluvia_anual"] * nn.1$lluvia_anual,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_bef, alpha.f = 0.15))
lines(x.seq, fit, col = col_bef, lwd = 3)
lines(x.seq, lwr, col = col_bef, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_bef, lty = 2, lwd = 1.5)

# ============================================================
BEF.precipitation <- with(nn.1,
                          BEF - p["R2.exp"] * R2.exp +
                            mean(p["R2.exp"] * R2.exp, na.rm = TRUE))
plot(BEF.precipitation ~ nn.1$lluvia_anual,
     xlab = "",
     ylab = "",
     bty = "l",
     pch = 19, col = pt_bef, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
mtext("Annual precipitation",
      side = 1, line = 3.75, cex = 2, col = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["lluvia_anual"] * x +
    mean(p["R2.exp"] * nn.1$R2.exp, na.rm = TRUE),
  from = min(nn.1$lluvia_anual, na.rm = TRUE),
  to   = max(nn.1$lluvia_anual, na.rm = TRUE),
  col  = col_bef
)
x.seq <- seq(
  min(nn.1$lluvia_anual, na.rm = TRUE),
  max(nn.1$lluvia_anual, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  lluvia_anual = x.seq,
  R2.exp = mean(nn.1$R2.exp, na.rm = TRUE))
pred <- predict(
  m_cov,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["R2.exp"] * nn.1$R2.exp,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_bef, alpha.f = 0.15))
lines(x.seq, fit, col = col_bef, lwd = 3)
lines(x.seq, lwr, col = col_bef, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_bef, lty = 2, lwd = 1.5)

######################
# Fig. S6 - PANEL B
slopes_df_sat_linear$deriv_S2-> nn.1$BEF
ax_col <- "#2E2E2E"
col_E <- "darkolivegreen3"  # azul violáceo
pt_E <- adjustcolor(col_E, alpha.f = 0.6)
double_curve <- function(fun, from, to, col, halo_col = ax_col,
                         halo_lwd = 6, col_lwd = 4.5, ...) {
  curve(fun, from = from, to = to, add = TRUE, col = halo_col, lwd = halo_lwd, ...)
  curve(fun, from = from, to = to, add = TRUE, col = col,      lwd = col_lwd,  ...)}
m_lineal_sat<- lm(BEF ~ R2.exp + lluvia_anual, data = nn.1)
p <- coefficients(m_lineal_sat)
BEF.neutral <- with(nn.1,
                    BEF - p["lluvia_anual"] * lluvia_anual +
                      mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))
plot(BEF.neutral ~ nn.1$R2.exp,
     xlab = "",
     ylab = "Richness-Biomass relationship (BEF)",
     bty = "l",
     pch = 19, col = pt_E, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * nn.1$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_E)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")")
  )), side = 1, line = 7, cex = 2, col = ax_col)
text(x = 0.3, y = 0.8, "F(2,18): 5.3\np: 0.015\nR²: 0.37",
     cex = 2.0, pos = 1, col = ax_col)
mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.5, col = ax_col)
x.seq <- seq(
  min(nn.1$R2.exp, na.rm = TRUE),
  max(nn.1$R2.exp, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  R2.exp = x.seq,
  lluvia_anual = mean(nn.1$lluvia_anual, na.rm = TRUE))
pred <- predict(
  m_lineal_sat,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["lluvia_anual"] * nn.1$lluvia_anual,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_E, alpha.f = 0.15))
lines(x.seq, fit, col = col_E, lwd = 3)
lines(x.seq, lwr, col = col_E, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_E, lty = 2, lwd = 1.5)

# ================================================
BEF.precipitation <- with(nn.1,
                          BEF - p["R2.exp"] * R2.exp +
                            mean(p["R2.exp"] * R2.exp, na.rm = TRUE))
plot(BEF.precipitation ~ nn.1$lluvia_anual,
     xlab = "",
     ylab = "",
     bty = "l",
     pch = 19, col = pt_E, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
mtext("Annual precipitation",
      side = 1, line = 3.75, cex = 2, col = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["lluvia_anual"] * x +
    mean(p["R2.exp"] * nn.1$R2.exp, na.rm = TRUE),
  from = min(nn.1$lluvia_anual, na.rm = TRUE),
  to   = max(nn.1$lluvia_anual, na.rm = TRUE),
  col  = col_E)
x.seq <- seq(
  min(nn.1$lluvia_anual, na.rm = TRUE),
  max(nn.1$lluvia_anual, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  lluvia_anual = x.seq,
  R2.exp = mean(nn.1$R2.exp, na.rm = TRUE))
pred <- predict(
  m_lineal_sat,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["R2.exp"] * nn.1$R2.exp,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_E, alpha.f = 0.15))
lines(x.seq, fit, col = col_E, lwd = 3)
lines(x.seq, lwr, col = col_E, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_E, lty = 2, lwd = 1.5)





plot_glm_visreg_R2 <- function(modelo, xvar = "R2.exp",
                               xlab = 
                                 expression(
                                   atop(
                                     "Variance in species occurrences",
                                     paste("explained by traits (", R[plain(Traits)]^2, ")")
                                   )
                                 ),
                               ylab = "P(quadratic)",
                               col_line = "navy",
                               col_fill = adjustcolor("lightblue", 0.25),
                               ax_col = "#2E2E2E",
                               line_xlab = 6.5) {
  
  visreg::visreg(
    modelo,
    xvar, scale = "response", ylim = c(0,1),
    yaxt = "n", xlab = " ", ylab = " ", cex.lab = 3,
    line.par = list(
      col = col_line,
      lwd = 3),
    fill.par = list(
      col = col_fill))
  axis(2, at = seq(0,1,0.2),
    labels = seq(0,1,0.2), las = 1)
  mtext(ylab, side = 2,
    line = 5,   # probar 5, 6, 7...
    cex = 1.8, col = ax_col)
  mtext(xlab, side = 1,
    line = 7,#5.75,
    cex = 2,#2,
    col = ax_col
  )
}

####################################################
###  ------------ Robustness ------------------   ##
####################################################

col_A <- "#2E5EAA"  
col_B <- "#B2472F"  
col_C <- "#1F8A70"  
col_D <- "#6A78A8"  

pt_A <- adjustcolor(col_A, alpha.f = 0.6)
pt_B <- adjustcolor(col_B, alpha.f = 0.6)
pt_C <- adjustcolor(col_C, alpha.f = 0.6)
pt_D <- adjustcolor(col_D, alpha.f = 0.6)

ax_col <- "#2E2E2E"

# ============= #
#               #
#    Fig. S8    #
#               #
# ============= #

slopes_df$P_quad-> nn.1$P_quad
m_quad<-(glm(P_quad~R2.exp, data=nn.1, family=binomial()))
summary(m_quad)

slopes_df$deriv_S1-> nn.1$deriv_S1
m_d1<-(lm(deriv_S1~R2.exp+lluvia_anual, data=nn.1))
summary(m_d1)

slopes_df$deriv_S2-> nn.1$deriv_S2
m_d2<-(lm(deriv_S2~R2.exp+lluvia_anual, data=nn.1))
summary(m_d2)

par(mfrow = c(3,2), mar = c(9,8,4,2), cex.lab = 2.2, cex.axis = 1.6)

# A
p<- coefficients(m_d1)
slopes_df$deriv_S1-> nn.1$BEF
BEF.neutral <- with(
  nn.1,
  BEF -
    p["lluvia_anual"] * lluvia_anual +
    p["lluvia_anual"] * mean(lluvia_anual))
plot(BEF.neutral ~ nn.1$R2.exp,
     xlab = "",
     ylab = "Richness-Biomass relationship (BEF)\n(derivative S=1)",
     bty = "l",
     pch = 19, col = pt_A, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * nn.1$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_A
)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")")
  )), side = 1, line = 7, cex = 1.5, col = ax_col)
text(x = 0.3, y = 0.75, "F(2,18): 4.7\np: 0.022\nR²: 0.34", cex = 2.5,
     pos = 3, col = ax_col)
mtext("A", side = 3, line = 1.0, adj = 0, font = 2, cex = 1.5, col = ax_col)
x.seq <- seq(
  min(nn.1$R2.exp, na.rm = TRUE),
  max(nn.1$R2.exp, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  R2.exp = x.seq,
  lluvia_anual = mean(nn.1$lluvia_anual, na.rm = TRUE))
pred <- predict(
  m_d1,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["lluvia_anual"] * nn.1$lluvia_anual,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_A, alpha.f = 0.15))
lines(x.seq, fit, col = col_A, lwd = 3)
lines(x.seq, lwr, col = col_A, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_A, lty = 2, lwd = 1.5)

BEF.precipitation <- with(nn.1,
                          BEF - p["R2.exp"] * R2.exp +
                            mean(p["R2.exp"] * R2.exp, na.rm = TRUE))
plot(BEF.precipitation ~ nn.1$lluvia_anual,
     xlab = "",
     ylab = "",
     bty = "l",
     pch = 19, col = pt_A, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
mtext("Annual precipitation",
      side = 1, line = 3.75, cex = 1.5, col = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["lluvia_anual"] * x +
    mean(p["R2.exp"] * nn.1$R2.exp, na.rm = TRUE),
  from = min(nn.1$lluvia_anual, na.rm = TRUE),
  to   = max(nn.1$lluvia_anual, na.rm = TRUE),
  col  = col_A)
x.seq <- seq(
  min(nn.1$lluvia_anual, na.rm = TRUE),
  max(nn.1$lluvia_anual, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  lluvia_anual = x.seq,
  R2.exp = mean(nn.1$R2.exp, na.rm = TRUE))
pred <- predict(
  m_d1,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["R2.exp"] * nn.1$R2.exp,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_A, alpha.f = 0.15))
lines(x.seq, fit, col = col_A, lwd = 3)
lines(x.seq, lwr, col = col_A, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_A, lty = 2, lwd = 1.5)

# B
p<- coefficients(m_d2)
slopes_df$deriv_S2-> nn.1$BEF
BEF.neutral <- with(
  nn.1,
  BEF -
    p["lluvia_anual"] * lluvia_anual +
    p["lluvia_anual"] * mean(lluvia_anual)
)
plot(BEF.neutral ~ nn.1$R2.exp,
     xlab = "",
     ylab = "Richness-Biomass relationship (BEF)\n(derivative S=2)",
     bty = "l",
     pch = 19, col = pt_B, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * nn.1$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_B
)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")")
  )), side = 1, line = 7, cex = 1.5, col = ax_col)
text(x = 0.3, y = 0.75, "F(2,18): 5.7\np: 0.012\nR²: 0.39", cex = 2.5,
     pos = 3, col = ax_col)
mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(
  min(nn.1$R2.exp, na.rm = TRUE),
  max(nn.1$R2.exp, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  R2.exp = x.seq,
  lluvia_anual = mean(nn.1$lluvia_anual, na.rm = TRUE))
pred <- predict(
  m_d2,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["lluvia_anual"] * nn.1$lluvia_anual,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["lluvia_anual"] *
  mean(nn.1$lluvia_anual, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_B, alpha.f = 0.15))
lines(x.seq, fit, col = col_B, lwd = 3)
lines(x.seq, lwr, col = col_B, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_B, lty = 2, lwd = 1.5)


BEF.precipitation <- with(nn.1,
                          BEF - p["R2.exp"] * R2.exp +
                            mean(p["R2.exp"] * R2.exp, na.rm = TRUE))
plot(BEF.precipitation ~ nn.1$lluvia_anual,
     xlab = "",
     ylab = "",
     bty = "l",
     pch = 19, col = pt_B, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
mtext("Annual precipitation",
      side = 1, line = 3.75, cex = 1.5, col = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["lluvia_anual"] * x +
    mean(p["R2.exp"] * nn.1$R2.exp, na.rm = TRUE),
  from = min(nn.1$lluvia_anual, na.rm = TRUE),
  to   = max(nn.1$lluvia_anual, na.rm = TRUE),
  col  = col_B)
x.seq <- seq(
  min(nn.1$lluvia_anual, na.rm = TRUE),
  max(nn.1$lluvia_anual, na.rm = TRUE),
  length.out = 200)
newdat <- data.frame(
  lluvia_anual = x.seq,
  R2.exp = mean(nn.1$R2.exp, na.rm = TRUE))
pred <- predict(
  m_d2,
  newdata = newdat,
  interval = "confidence")
adj <- mean(
  p["R2.exp"] * nn.1$R2.exp,
  na.rm = TRUE)
fit <- pred[, "fit"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
lwr <- pred[, "lwr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
upr <- pred[, "upr"] -
  p["R2.exp"] *
  mean(nn.1$R2.exp, na.rm = TRUE) +
  adj
polygon(
  c(x.seq, rev(x.seq)),
  c(lwr, rev(upr)),
  border = NA,
  col = adjustcolor(col_B, alpha.f = 0.15))
lines(x.seq, fit, col = col_B, lwd = 3)
lines(x.seq, lwr, col = col_B, lty = 2, lwd = 1.5)
lines(x.seq, upr, col = col_B, lty = 2, lwd = 1.5)


# C
plot_glm_visreg_R2(
  modelo = m_quad)
text(x = 0.72, y = 0.9, "p: 0.87", cex = 2.5,
     pos = 1, col = ax_col)
mtext("C", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)



# ============
#
# Fig. S10
#
# ============

slopes_df_sat$P_quad-> nn.1$P_quad
m_quad_sat<-(glm(P_quad~R2.exp, data=nn.1, family=binomial()))
summary(m_quad_sat)

slopes_df_sat$deriv_S1-> nn.1$deriv_S1
m_d1_sat<-(lm(deriv_S1~R2.exp #+lluvia_anual
              , data=nn.1))
summary(m_d1_sat)

slopes_df_sat$deriv_S2-> nn.1$deriv_S2
m_d2_sat<-(lm(deriv_S2~R2.exp #+lluvia_anual
              , data=nn.1))
summary(m_d2_sat)

par(mfrow = c(2,2), mar = c(9,8,4,2), cex.lab = 2.2, cex.axis = 1.6)

# A
p<- coefficients(m_d1_sat)
slopes_df_sat$deriv_S1-> nn.1$BEF
plot(BEF ~ R2.exp, data = nn.1,
     ylab = "Richness-Biomass relationship (BEF)\n(derivative S=1)",
     xlab = " ",
     pch = 19, col = pt_C, cex = 2.5,
     bty = "l",
     col.axis = ax_col, col.lab = ax_col)
mtext(expression(atop("Variance in species occurrences",
                      paste("explained by traits (", R[plain(Traits)]^2, ")"))),
      side = 1, line = 7, cex = 2, col = ax_col)
double_curve(
  fun = function(x)
    p[1] +
    p[2]*x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_C)
text(x = 0.3, y = 1,
     "F(1,19): 3.5\np: 0.077\nR²: 0.16",
     cex = 2, pos = 3, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),max(nn.1$R2.exp, na.rm = TRUE),length.out = 200)
pred <- predict(m_d1_sat,newdata = data.frame(R2.exp = x.seq),interval = "confidence")
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),
        border = NA,col = adjustcolor(col_C, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_C,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_C,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_C,lty = 2,lwd = 1.5)


# B
p<- coefficients(m_d2_sat)
slopes_df_sat$deriv_S2-> nn.1$BEF
plot(BEF ~ R2.exp, data = nn.1,
     ylab = "Richness-Biomass relationship (BEF)\n(derivative S=2)",
     xlab = " ",
     pch = 19, col = pt_D, cex = 2.5,
     bty = "l",
     col.axis = ax_col, col.lab = ax_col)
mtext(expression(atop("Variance in species occurrences",
                      paste("explained by traits (", R[plain(Traits)]^2, ")"))),
      side = 1, line = 7, cex = 2, col = ax_col)
double_curve(
  fun = function(x)
    p[1] +
    p[2]*x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_D)
text(x = 0.3, y = 1,
     "F(1,19): 3.7\np: 0.07\nR²: 0.16",
     cex = 2, pos = 3, col = ax_col)
mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),max(nn.1$R2.exp, na.rm = TRUE),length.out = 200)
pred <- predict(m_d2_sat,newdata = data.frame(R2.exp = x.seq),interval = "confidence")
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),
        border = NA,col = adjustcolor(col_D, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_D,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_D,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_D,lty = 2,lwd = 1.5)

# C
plot_glm_visreg_R2(
  modelo = m_quad_sat)
text(x = 0.72, y = 0.9, "p: 0.81", cex = 2,
     pos = 1, col = ax_col)
mtext("C", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)



################
# Fig. S12
# 
# Quantile model

col_F <- "#A05C7B"  # malva apagado
col_G <- "#C28B2C"  # mostaza oscura
pt_F <- adjustcolor(col_F, alpha.f = 0.6)
pt_G <- adjustcolor(col_G, alpha.f = 0.6)

slopes_covariates_q90$P_quad-> nn.1$P_quad
m_quad_q90<-(glm(P_quad~R2.exp, data=nn.1, family=binomial()))
summary(m_quad_q90)

slopes_covariates_q90$deriv_S1-> nn.1$BEF
m_q90_d1<- lm(BEF~R2.exp, data=nn.1)
summary(m_q90_d1)

slopes_covariates_q90$deriv_S2-> nn.1$BEF
m_q90_d2<- lm(BEF~R2.exp, data=nn.1)
summary(m_q90_d2)


par(mfrow = c(2,2), mar = c(9, 9, 4, 3), mgp = c(3.5, 1, 0),
    cex.lab = 2, cex.axis = 2)

slopes_covariates_q90$deriv_S1-> nn.1$BEF
p<- coefficients(m_q90_d1)
summary(m_q90_d1)
plot(BEF ~ R2.exp, data=nn.1,
     xlab = "",
     ylab = "Richness-Biomass relationship (BEF)\n(q90, derivative S=1)",
     bty = "l",
     pch = 19, col = pt_F, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_F)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")")
  )), side = 1, line = 7, cex = 2, col = ax_col)
text(x = 0.3, y = 5, "F(1,19): 8.09\np: 0.01\nR²: 0.3", cex = 2,
     pos = 3, col = ax_col)
mtext("A", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),max(nn.1$R2.exp, na.rm = TRUE),length.out = 200)
pred <- predict(m_q90_d1,newdata = data.frame(R2.exp = x.seq),interval = "confidence")
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),
        border = NA,col = adjustcolor(col_F, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_F,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_F,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_F,lty = 2,lwd = 1.5)

slopes_covariates_q90$deriv_S2-> nn.1$BEF
p<- coefficients(m_q90_d2)
summary(m_q90_d2)
plot(BEF ~ R2.exp, data=nn.1,
     xlab = "",
     ylab = "Richness-Biomass relationship (BEF)\n(q90, derivative S=2)",
     bty = "l",
     pch = 19, col = pt_G, cex = 3,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_G)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")")
  )), side = 1, line = 7, cex = 2, col = ax_col)
text(x = 0.3, y = 3, "F(1,19): 7.5\np: 0.013\nR²: 0.28", cex = 2,
     pos = 3, col = ax_col)
mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),max(nn.1$R2.exp, na.rm = TRUE),length.out = 200)
pred <- predict(m_q90_d2,newdata = data.frame(R2.exp = x.seq),interval = "confidence")
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),border = NA,col = adjustcolor(col_G, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_G,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_G,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_G,lty = 2,lwd = 1.5)

# C
plot_glm_visreg_R2(
  modelo = m_quad_q90)
text(x = 0.28, y = 0.9, "p: 0.11", cex = 1.8,
     pos = 1, col = ax_col)
mtext("C", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
