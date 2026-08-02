###############################################################
###############################################################

# de León et al. 
# Biodiversity-function relationship determined
# by the niche-neutral gradient of community assembly

###############################################################
# Main figures - Manuscript
#
# This script generates all main figures presented in the paper.
#
# Before running this script, the following scripts must be
# executed to generate the required intermediate datasets:
#
# 1. BEF.R
#    - Estimates annual BEF relationships using the different
#      analytical approaches (main GLMM analysis and robustness
#      analyses presented in SUPP).
#
# 2. CATS.R
#    - Estimates the annual position of the metacommunity along
#      the niche–neutral gradient (R²Traits and related metrics).
#
# 3. Functional_diversity.R
#    - Calculates annual functional diversity metrics
#      (FRic, FEve, FDis, FDiv, Rao's Q).
#
# 4. NST.R
#    - Calculates the Normalized Stochasticity Ratio (NST) and
#      compares it with CATS estimates (Supplementary Figure S1).
#
# After these scripts have been executed, this script can be run
# to reproduce all figures included in the main manuscript.
###############################################################


library(dplyr)
library(readxl)
library(ggplot2)


# Preparing environmental variables we will use for the analyses
Clima_dropbox <- read_excel("Clima_dropbox.xlsx")
lluvia_anual <- Clima_dropbox %>%
  group_by(year) %>%
  summarise(lluvia_prom = median(Lluvia_mm, na.rm = TRUE), lluvia_total = sum(Lluvia_mm, na.rm = TRUE)) %>%
  arrange(year) %>%
  mutate(delta_lluvia_anual = lluvia_prom - lag(lluvia_prom))

lluvia_abril_sep <- Clima_dropbox %>%
  filter(mes >= 4 & mes <= 9) %>%   # abril a septiembre
  group_by(year) %>%
  summarise(lluvia_abril_sep = median(Lluvia_mm, na.rm = TRUE))

head(lluvia_abril_sep)
head(lluvia_anual)
lluvia_anual<- lluvia_anual[-c(1:3),]
lluvia_abril_sep<- lluvia_abril_sep[-c(1:3),]


lluvia_resumen <- Clima_dropbox[-c(289,290),] %>%
  group_by(year) %>%
  summarise(
    media   = mean(Lluvia_mm, na.rm = TRUE),
    mediana = median(Lluvia_mm, na.rm = TRUE),
    total   = sum(Lluvia_mm, na.rm = TRUE))

calc_lluvia <- function(data, mes_ini, mes_fin){
  data %>%
    filter(mes >= mes_ini & mes <= mes_fin) %>%
    group_by(year) %>%
    summarise(lluvia_mediana = median(Lluvia_mm, na.rm = TRUE), .groups = "drop") %>%
    mutate(ventana = paste0(mes_ini, "-", mes_fin))
}

ventanas <- list(c(4,9)) # Months period used
par(mfrow = c(2,3))
for(v in ventanas){
  clim <- calc_lluvia(Clima_dropbox, v[1], v[2])
  datos <- data.frame(
    year = nn.1$año,
    R2.exp = nn.1$R2.exp
  ) %>%
    left_join(clim, by = "year")
  mod <- lm(
    R2.exp ~ lluvia_mediana + I(lluvia_mediana^2),
    data = datos)
  xseq <- seq(min(datos$lluvia_mediana),
              max(datos$lluvia_mediana),
              length.out = 200)
  pred <- predict(mod,newdata = data.frame(lluvia_mediana = xseq))
  lines(xseq, pred, lwd = 2)
  print(summary(mod))
}


calc_delta_lluvia <- function(data, mes_ini, mes_fin){
  
  data %>%
    filter(mes >= mes_ini & mes <= mes_fin) %>%
    group_by(year) %>%
    summarise(
      lluvia_mediana = median(Lluvia_mm, na.rm = TRUE),
      HR_mediana = median(HR, na.rm = TRUE),
      .groups = "drop") %>%
    arrange(year) %>%
    mutate(
      delta_lluvia = lluvia_mediana - lag(lluvia_mediana),
      delta_HR = HR_mediana - lag(HR_mediana),
      ventana = paste0(mes_ini, "-", mes_fin))
}

ventanas <- list(c(4,9))
par(mfrow = c(1,1))
for(v in ventanas){
  clim <- calc_delta_lluvia(Clima_dropbox, v[1], v[2])
  datos <- data.frame(
    year = nn.1$año,
    R2.exp = nn.1$R2.exp) %>%
    left_join(clim, by = "year")
  mod <- lm(R2.exp ~ delta_lluvia + I(delta_lluvia^2),data = datos)
  xseq <- seq(min(datos$delta_lluvia, na.rm = TRUE),
              max(datos$delta_lluvia, na.rm = TRUE),
              length.out = 200)
  pred <- predict(mod, newdata = data.frame(delta_lluvia = xseq))
  plot(datos$delta_lluvia,
       datos$R2.exp,
       pch = 19,
       xlab = "Delta lluvia",
       ylab = "R2.exp",
       main = paste0(v[1], "-", v[2]))
  lines(xseq, pred, lwd = 2)
  print(summary(mod))
}






########## ------------ FIGURE 1. ------------- #########
#--------           SYSTEM DYNAMIC         -------------#

layout(matrix(1:2, nrow = 1))
par(mar = c(6, 8, 4, 6), mgp = c(3.5, 1, 0), cex.lab = 2.2, cex.axis = 1.8)

ax_col <- "#2E2E2E"
col_drift <- "darkgreen"
col_bef   <- "darkblue"

pt_drift <- adjustcolor(col_drift, alpha.f = 0.6)
pt_bef   <- adjustcolor(col_bef,   alpha.f = 0.6)
double_line <- function(x, y, col, halo_col = ax_col, halo_lwd = 4, col_lwd = 2.5) {
  lines(x, y, col = halo_col, lwd = halo_lwd)
  lines(x, y, col = col,      lwd = col_lwd)
}

#========================#
# ---- PANEL A: ---------#
df_f1 <- data.frame(
  year = df_plot$year,
  R2.exp = df_plot$R2.exp,
  BEF = pendientes.lineal[,2])
r1 <- range(df_f1$R2.exp, na.rm = TRUE)
r2 <- range(df_f1$BEF,    na.rm = TRUE)
a <- diff(r1) / diff(r2)
b <- r1[1] - a * r2[1]
df_f1$BEF.scaled <- a * df_f1$BEF + b

col_assembly <- adjustcolor("darkgreen", alpha.f = 1)
col_bef      <- adjustcolor("darkblue",  alpha.f = 1)

p1.A <- ggplot(df_f1, aes(x = year)) +
  geom_line(
    aes(y = R2.exp, colour = "Assembly"),
    linewidth = 1.3
  ) +
  geom_point(
    aes(y = R2.exp, colour = "Assembly"),
    size = 5
  ) +
  geom_line(
    aes(y = BEF.scaled, colour = "BEF"),
    linewidth = 1.3,
    linetype = 2
  ) +
  geom_point(
    aes(y = BEF.scaled, colour = "BEF"),
    shape = 18,
    size = 5
  ) +
  scale_colour_manual(
    values = c(
      Assembly = col_assembly,
      BEF      = col_bef
    ),
    name = NULL,
    labels = c(
      expression(R[Traits]^2),
      "BEF"
    )
  ) +
  scale_y_continuous(
    name = expression(
      atop(
        "Variance in species occurrences",
        paste("explained by traits (", R[plain(Traits)]^2, ")")
      )
    ),
    sec.axis = sec_axis(
      trans = ~ (. - b) / a,
      name = "Richness-Biomass relationship (BEF)"
    )
  ) +
  scale_x_continuous(
    breaks = seq(2005, 2025, 2)
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.80, 0.90),
    legend.background = element_rect(
      fill = adjustcolor("white", alpha.f = 0.6),
      colour = NA
    ),
    legend.key = element_blank(),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y.left = element_text(
      size = 20,
      colour = "darkgreen"
    ),
    axis.title.y.right = element_text(
      size = 20,
      colour = "darkblue"
    ),
    axis.text.x = element_text(
      size = 18,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y.left = element_text(
      size = 18,
      colour = "darkgreen"
    ),
    axis.text.y.right = element_text(
      size = 18,
      colour = "darkblue"
    )) +
  labs(x = "Year")

#========================#
# ---- PANEL B: ---------#
# Panel B. R2.traits as a function of interannual variation in precipitation

col_A <- "#2E5EAA"
col_B <- "#B2472F"
pt_A <- adjustcolor(col_A, alpha.f = 0.6)
pt_B <- adjustcolor(col_B, alpha.f = 0.6)

m1_mediana <-lm(formula = R2.exp ~ delta_lluvia + I(delta_lluvia^2), data = datos)
summary(m1_mediana)

x.seq <- seq(
  min(datos$delta_lluvia, na.rm = TRUE),
  max(datos$delta_lluvia, na.rm = TRUE),
  length.out = 300)
pred <- predict(
  m1_mediana,
  newdata = data.frame(
    delta_lluvia = x.seq),
  interval = "confidence")
pred.df <- data.frame(
  delta_lluvia = x.seq,
  fit = pred[, "fit"],
  lwr = pred[, "lwr"],
  upr = pred[, "upr"])

p1.B <-
  ggplot(datos,
         aes(delta_lluvia, R2.exp)) +
  geom_ribbon(
    data = pred.df,
    aes(
      x = delta_lluvia,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    fill = col_A,
    alpha = 0.20
  ) +
  geom_line(
    data = pred.df,
    aes(y = lwr),
    colour = col_A,
    linewidth = 0.8,
    linetype = 2
  ) +
  geom_line(
    data = pred.df,
    aes(
      y = upr),
    colour = col_A,
    linewidth = 0.8,
    linetype = 2
  ) +
  geom_line(
    data = pred.df,
    aes(y = fit),
    colour = col_A,
    linewidth = 1.3
  ) +
  geom_point(
    colour = pt_A,
    size = 4
  ) +
  labs(
    x = "Interannual change in pre-sampling precipitation\n(Apr-Sep)",
    y = expression(
      atop(
        "Variance in species occurrences",
        paste(
          "explained by traits (",
          R[plain(Traits)]^2,
          ")"
        )
      )
    )
  ) +
  annotate(
    "text",
    x = 50,
    y = 0.175,
    label = "F(2,18): 3.8\np = 0.04\nR² = 0.3",
    size = 7
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

p1.A <- p1.A +
  theme(
    plot.margin = margin(
      t = 5,
      r = 10,
      b = 5,
      l = 10))
p1.B <- p1.B +
  theme(
    plot.margin = margin(
      t = 5,
      r = 10,
      b = 5,
      l = 35))
ggarrange(
  p1.A,
  p1.B,
  ncol = 2,
  labels = c("A", "B"),
  label.x = 0.02,
  label.y = 0.98,
  hjust = 0,
  vjust = 1,
  font.label = list(
    size = 18,
    face = "bold"))



########## ------------ FIGURE 2. ------------- #########
# R2.traits effect on taxonomic and functional diversity # 

#############################################

col_A <- "#B2472F"  # riqueza
col_B <- "#3B6FB6"  # FRic
col_C <- "#8A6BBE"  # FEve
col_D <- "#2D8B74"  # FDis

pt_A <- adjustcolor(col_A, alpha.f = 0.70)
pt_B <- adjustcolor(col_B, alpha.f = 0.70)
pt_C <- adjustcolor(col_C, alpha.f = 0.70)
pt_D <- adjustcolor(col_D, alpha.f = 0.70)

ax_col <- "black"


#=========================================================
# B-C: Metacommunity species richness
#=========================================================

par(mfrow=c(1,1), mar=c(9,9,3,3))
nn.1$riqueza<- riqueza
m_rich <-lm(formula = riqueza ~ R2.exp, data = nn.1)
summary(m_rich)
p<- coefficients(m_rich)
plot(riqueza ~ R2.exp, data = nn.1,
     ylab = "Metacommunity species richness",
     xlab = " ",
     pch = 19, col = pt_B, cex = 2.5,
     bty = "l",
     col.axis = ax_col, col.lab = ax_col)
mtext(expression(atop("Variance in species occurrences",
                      paste("explained by traits (", R[plain(Traits)]^2, ")"))),
      side = 1, line = 6, cex = 2, col = ax_col)
double_curve(
  fun = function(x)
    p[1] +
    p[2]*x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_B)
text(x = 0.8, y = 95,
     "F(1,19): 5.9\np: 0.025\nR²: 0.24",
     cex = 2, pos = 3, col = ax_col)

x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),max(nn.1$R2.exp, na.rm = TRUE),length.out = 200)
pred <- predict(m_rich,newdata = data.frame(R2.exp = x.seq),interval = "confidence")
# banda semitransparente
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),
        border = NA,col = adjustcolor(col_B, alpha.f = 0.15))
# línea central
lines(x.seq,pred[, "fit"],col = col_B,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_B,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_B,lty = 2,lwd = 1.5)


#=========================================================
# B-C: Path results
#=========================================================

# func_path has results from Path analyses, created in "Functional_diversity.R" 
# plot_mediation is a function also created in the other script 
# for visualizing path results

#B
plot_mediation(func_path, "FDis", "Functional dispersion")
#C
plot_mediation(func_path, "FEve", "Functional evenness")





######### ---- FIGURE 3. ------ ########
# R2.traits effect on biomass and BEF # 

# Data we will use:

nn.bef <- nn.1[, c( "año", "BEF", "R2.exp", "prof_media", "prof_CV")]

clima.year <- Clima_dropbox %>%
  group_by(year) %>%
  summarise(
    lluvia = median(Lluvia_mm, na.rm = TRUE),
    temp = median(Temp_Med_C, na.rm = TRUE),
    hr = median(HR, na.rm = TRUE),
    intv = median(IntV, na.rm = TRUE)
  ) %>%
  arrange(year)

datos <- left_join(
  nn.bef,
  clima.year,
  by = c("año" = "year")
)
datos.z <- as.data.frame(scale(datos[,-1]))

#####################
## --- Biomasa
um_por_año <- M %>%
  group_by(año) %>%
  summarise(n_UM = n())
bm.nuevo_df<- as.data.frame(bm.nuevo[,1:6])
na.omit(bm.nuevo_df)-> bm.nuevo_df
biomasa_total_anual <- bm.nuevo_df %>%
  group_by(año) %>%
  summarise(BiomasaTotal = sum(biom))
biomasa_estandarizada <- biomasa_total_anual %>%
  left_join(um_por_año, by = "año") %>%
  mutate(Biomasa_Estandarizada = BiomasaTotal / n_UM)
biomasa_estandarizada<- as.matrix(biomasa_estandarizada[,4])
datos$biomasa<- biomasa_estandarizada
bestglm(Xy = datos[,c("R2.exp","prof_media","prof_CV","lluvia","temp","hr","intv", "biomasa")], 
        family = gaussian, IC = "AIC", nvmax=2)$BestModels


# Panel A. R2.traits ~ biomass

layout(matrix(c(1,3,2,4), nrow = 2))
par(mar = c(9, 9, 3, 3), mgp = c(3.5, 1, 0), cex.lab = 2.2, cex.axis = 1.8)
col_bef  <- "#3A8EC1"
col_biomasa <- "#E69F00"
pt_bef   <- adjustcolor(col_bef,  alpha.f = 0.6)
pt_biomasa  <- adjustcolor(col_biomasa, alpha.f = 0.6)
ax_col   <- "#2E2E2E"


m_biom <- lm(biomasa ~ R2.exp + lluvia, data = datos)
summary(m_biom)
p2 <- coefficients(m_biom)
biomasa.neutral <- with(datos,
                        biomasa - p2["lluvia"] * lluvia +
                          mean(p2["lluvia"] * lluvia, na.rm = TRUE))
plot(biomasa.neutral ~ datos$R2.exp,
     xlab = "", ylab = expression("Mean biomass (g / 0.04 m"^2*")"),
     bty = "l", pch = 19, col = pt_biomasa, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p2["(Intercept)"] +
    mean(p2["lluvia"] * datos$lluvia, na.rm = TRUE) +
    p2["R2.exp"] * x,
  from = min(datos$R2.exp, na.rm = TRUE),
  to   = max(datos$R2.exp, na.rm = TRUE),
  col  = col_biomasa
)

mtext(expression(atop("Variance in species occurrences",
                      paste("explained by traits (", R[plain(Traits)]^2, ")"))),
      side = 1, line = 7, cex = 2, col = ax_col)
text(x = 0.25, y = 8, "F(2,18): 7.3\np: 0.0047\nR²: 0.45", cex = 2, pos = 1, col = ax_col)

x.seq <- seq(min(datos$R2.exp, na.rm = TRUE),max(datos$R2.exp, na.rm = TRUE),length.out = 200)
newdat <- data.frame(
  R2.exp = x.seq,
  lluvia = mean(datos$lluvia, na.rm = TRUE))
pred <- predict(m_biom,newdata = newdat,interval = "confidence")
adj <- mean(p2["lluvia"] * datos$lluvia, na.rm = TRUE)
lines(x.seq, pred[, "fit"] - p2["lluvia"] * mean(datos$lluvia) + adj,
      col = col_biomasa, lwd = 3)
lines(x.seq, pred[, "lwr"] - p2["lluvia"] * mean(datos$lluvia) + adj,
      col = col_biomasa, lty = 2)
lines(x.seq, pred[, "upr"] - p2["lluvia"] * mean(datos$lluvia) + adj,
      col = col_biomasa, lty = 2)
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),border = NA,
  col = adjustcolor(col_biomasa, alpha.f = 0.2))
mtext("A", side = 3, line = 1.2, adj = 0, font = 2, cex = 2, col = ax_col)
# ----------  biomasa ~ precip ----------
biomasa.precipitation <- with(datos, biomasa - p2["R2.exp"] * R2.exp +
                                mean(p2["R2.exp"] * R2.exp, na.rm = TRUE))
plot(biomasa.precipitation ~ datos$lluvia,
     xlab = " ", ylab = "",
     bty = "l", pch = 19, col = pt_biomasa, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)
mtext("Annual precipitation", side = 1, line = 3.75, cex = 2, col = ax_col)
double_curve(
  fun  = function(x) p2["(Intercept)"] + p2["lluvia"] * x +
    mean(p2["R2.exp"] * datos$R2.exp, na.rm = TRUE),
  from = min(datos$lluvia, na.rm = TRUE),
  to   = max(datos$lluvia, na.rm = TRUE),col  = col_biomasa)
x.seq <- seq(min(datos$lluvia, na.rm = TRUE),max(datos$lluvia, na.rm = TRUE),length.out = 200)
newdat <- data.frame(
  lluvia = x.seq,R2.exp = mean(datos$R2.exp, na.rm = TRUE))
pred <- predict(m_biom,newdata = newdat,interval = "confidence")
adj <- mean(p2["R2.exp"] * datos$R2.exp, na.rm = TRUE)
lines(x.seq, pred[, "fit"] - p2["R2.exp"] * mean(datos$R2.exp) + adj,
      col = col_biomasa, lwd = 3)
lines(x.seq, pred[, "lwr"] - p2["R2.exp"] * mean(datos$R2.exp) + adj,
      col = col_biomasa, lty = 2)
lines(x.seq, pred[, "upr"] - p2["R2.exp"] * mean(datos$R2.exp) + adj,
      col = col_biomasa, lty = 2)
polygon(c(x.seq, rev(x.seq)),c(pred[, "lwr"], rev(pred[, "upr"])),border = NA, 
  col = adjustcolor(col_biomasa, alpha.f = 0.2))



###############
# --- BEF 
datos$BEF <- pendientes.lineal[,2] # Here you can choose the BEF estimate you want to use for the analyses
bestglm(
  Xy = datos[, c(
    "R2.exp",
    "prof_media",
    "prof_CV",
    "lluvia",
    "temp",
    "hr",
    "intv",
    "BEF"
  )],
  IC = "AIC",
  nvmax = 2)$BestModels
r2(lm(BEF~R2.exp+lluvia,data=datos))
r2(lm(BEF~lluvia,data=datos))

m_bef <- lm(BEF ~ R2.exp + lluvia, data = datos)
summary(m_bef)
p2 <- coefficients(m_bef)
BEF.neutral <- with(datos,
                    BEF - p2["lluvia"] * lluvia +
                      mean(p2["lluvia"] * lluvia, na.rm = TRUE))
plot(BEF.neutral ~ datos$R2.exp,
     xlab = "", ylab = "Richness-biomass relationship (BEF)",
     bty = "l", pch = 19, col = pt_bef, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p2["(Intercept)"] +
    mean(p2["lluvia"] * datos$lluvia, na.rm = TRUE) +
    p2["R2.exp"] * x,
  from = min(datos$R2.exp, na.rm = TRUE),
  to   = max(datos$R2.exp, na.rm = TRUE),
  col  = col_bef
)
mtext(expression(atop("Variance in species occurrences",
                      paste("explained by traits (", R[plain(Traits)]^2, ")"))),
      side = 1, line = 7, cex = 2, col = ax_col)
text(x = 0.25, y = 0.9, "F(2,18): 6.2\np: 0.009\nR²: 0.41", cex = 2, pos = 1, col = ax_col)

# grilla
x.seq <- seq(
  min(datos$R2.exp, na.rm = TRUE),
  max(datos$R2.exp, na.rm = TRUE),
  length.out = 200
)
# predicción manteniendo lluvia en su valor medio
newdat <- data.frame(
  R2.exp = x.seq,
  lluvia = mean(datos$lluvia, na.rm = TRUE)
)
pred <- predict(
  m_bef,
  newdata = newdat,
  interval = "confidence"
)
# transformar al mismo espacio que BEF.neutral
adj <- mean(p2["lluvia"] * datos$lluvia, na.rm = TRUE)
lines(x.seq, pred[, "fit"] - p2["lluvia"] * mean(datos$lluvia) + adj,
      col = col_bef, lwd = 3)
lines(x.seq, pred[, "lwr"] - p2["lluvia"] * mean(datos$lluvia) + adj,
      col = col_bef, lty = 2)
lines(x.seq, pred[, "upr"] - p2["lluvia"] * mean(datos$lluvia) + adj,
      col = col_bef, lty = 2)
polygon(
  c(x.seq, rev(x.seq)),
  c(pred[, "lwr"], rev(pred[, "upr"])),
  border = NA,
  col = adjustcolor(col_bef, alpha.f = 0.2)
)
mtext("B", side = 3, line = 1.2, adj = 0, font = 2, cex = 2, col = ax_col)

# ----------  BEF ~ precip (ajustada por explicada-por-rasgos) ----------
BEF.precipitation <- with(datos,
                          BEF - p2["R2.exp"] * R2.exp +
                            mean(p2["R2.exp"] * R2.exp, na.rm = TRUE))
plot(BEF.precipitation ~ datos$lluvia,
     xlab = " ", ylab = "",
     bty = "l", pch = 19, col = pt_bef, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)
mtext("Annual precipitation", side = 1, line = 3.75, cex = 2, col = ax_col)
double_curve(
  fun  = function(x) p2["(Intercept)"] + p2["lluvia"] * x +
    mean(p2["R2.exp"] * datos$R2.exp, na.rm = TRUE),
  from = min(datos$lluvia, na.rm = TRUE),
  to   = max(datos$lluvia, na.rm = TRUE),
  col  = col_bef
)
x.seq <- seq(
  min(datos$lluvia, na.rm = TRUE),
  max(datos$lluvia, na.rm = TRUE),
  length.out = 200
)

newdat <- data.frame(
  lluvia = x.seq,
  R2.exp = mean(datos$R2.exp, na.rm = TRUE)
)

pred <- predict(
  m_bef,
  newdata = newdat,
  interval = "confidence"
)

adj <- mean(p2["R2.exp"] * datos$R2.exp, na.rm = TRUE)

lines(x.seq, pred[, "fit"] - p2["R2.exp"] * mean(datos$R2.exp) + adj,
      col = col_bef, lwd = 3)

lines(x.seq, pred[, "lwr"] - p2["R2.exp"] * mean(datos$R2.exp) + adj,
      col = col_bef, lty = 2)

lines(x.seq, pred[, "upr"] - p2["R2.exp"] * mean(datos$R2.exp) + adj,
      col = col_bef, lty = 2)
polygon(
  c(x.seq, rev(x.seq)),
  c(pred[, "lwr"], rev(pred[, "upr"])),
  border = NA,
  col = adjustcolor(col_bef, alpha.f = 0.2)
)



































######################################################
### ---- Fig. SUPLEMENTARY

# Here only are presented Fig.S2, the rest can be found at the following scripts:

# Fig. S3 "Functional_diversity.R" 

# Fig.S13: "NST.R" -> It contains both NST analyses and the relation between NST and R2.traits.

# Fig.S4-S12: "BEF.R"-> It contains the necessary script for both estimate and visualize BEF slopes 
#                       with different implemented models and the estimate of BEF~R2.traits for each one.

#######################################
### ----- FIG S2. DIVERSITY ------ ####

div<- NULL
for(i in 2005:2025){
  apply(M[which(M$año==i),-c(1:4)],1,sum)-> x
  cbind(mean(x, na.rm=T), round(nn.1$riqueza,0)[which(nn.1$año==i)])-> x.temp
  div<- rbind(div, x.temp)
}
div

div.charco <- NULL
for(i in 2005:2025){
  temp <- br.muestreos_df[br.muestreos_df$año == i,]
  ch <- rowsum(
    temp[, -(1:4)],
    group = temp$Charco
  )
  ch[ch > 0] <- 1
  riqueza.charco <- rowSums(ch)
  div.charco <- rbind(
    div.charco,
    cbind(
      mean(riqueza.charco),
      round(nn.1$riqueza, 0)[nn.1$año == i]
    )  )  }
div.charco
cbind(div[,1], div.charco)-> div


p <- coefficients(lm((div[,1] ~ div[,2])))
summary(lm((div[,1] ~ div[,2])))
par(mfrow=c(1,1), mar=c(5,5,2.5,1))
plot(div[,1] ~ div[,2],
     xlab = "Riqueza de especies de la metacomunidad",
     ylab = "Riqueza de especies promedio UM",
     bty = "l",
     pch = 19, col = pt_B, cex = 2, cex.lab = 1.4,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun = function(x) p[1] + p[2]*x,
  from = min(div[,2], na.rm = TRUE),
  to   = max(div[,2], na.rm = TRUE),
  col = col_B)
text(y = 4, x = 90, "F(1,19): 114.1\np < 0.001\nR²: 0.86",
     cex = 1.3, pos = 1, col = ax_col)



par(mfrow=c(3,1), mar=c(9,7,4,1))
p <- coefficients(lm((div[,1] ~ nn.1$R2.exp)))
df_1 <- data.frame(y = div[,1],R2.exp = nn.1$R2.exp)
m1<- lm(y ~ R2.exp, data=df_1)
summary(lm((div[,1] ~ nn.1$R2.exp)))
plot(div[,1] ~ nn.1$R2.exp,
     xlab = " ",
     ylab = "Mean species richness per sampling unit",
     bty = "l",
     pch = 19, col = pt_B, cex = 2, cex.lab = 2.7, cex.axis=2,
     col.axis = ax_col, col.lab = ax_col)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")"))),
  side = 1, line = 7, cex = 1.65, col = ax_col)
double_curve(
  fun = function(x) p[1] + p[2]*x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col = col_B)
text(y = 3.75, x = 0.25, "F(1,19): 6.85\np: 0.017\nR²: 0.27",
     cex = 2.5, pos = 1, col = ax_col)
mtext("A", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),
             max(nn.1$R2.exp, na.rm = TRUE), length.out = 200)
pred <- predict(m1, newdata = data.frame(R2.exp = x.seq),
                interval = "confidence")
polygon(c(x.seq, rev(x.seq)),
        c(pred[, "lwr"], rev(pred[, "upr"])), border = NA,
        col = adjustcolor(col_B, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_B,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_B,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_B,  lty = 2, lwd = 1.5)

p <- coefficients(lm((div[,2] ~ nn.1$R2.exp)))
df_2 <- data.frame(y = div[,2],R2.exp = nn.1$R2.exp)
m2<- lm(y ~ R2.exp, data=df_2)
summary(lm((div[,2] ~ nn.1$R2.exp)))
plot(div[,2] ~ nn.1$R2.exp,
     xlab = " ",
     ylab = "Mean species richness per pond",
     bty = "l",
     pch = 19, col = pt_B, cex = 2, cex.lab = 2.7, cex.axis=2,
     col.axis = ax_col, col.lab = ax_col)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")"))),
  side = 1, line = 7, cex = 1.65, col = ax_col)
double_curve(
  fun = function(x) p[1] + p[2]*x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col = col_B)
text(y = 8.2, x = 0.25, "F(1,19): 12.1\np: 0.003\nR²: 0.39",
     cex = 2.5, pos = 1, col = ax_col)
mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),
             max(nn.1$R2.exp, na.rm = TRUE), length.out = 200)
pred <- predict(m2, newdata = data.frame(R2.exp = x.seq),
                interval = "confidence")
polygon(c(x.seq, rev(x.seq)),
        c(pred[, "lwr"], rev(pred[, "upr"])), border = NA,
        col = adjustcolor(col_B, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_B,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_B,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_B,  lty = 2, lwd = 1.5)

p <- coefficients(lm((div[,3] ~ nn.1$R2.exp)))
df_3 <- data.frame(y = div[,3],R2.exp = nn.1$R2.exp)
m3<- lm(y ~ R2.exp, data=df_3)
summary(lm((div[,3] ~ nn.1$R2.exp)))
plot(div[,3] ~ nn.1$R2.exp,
     xlab = " ",
     ylab = "Metacommunity species richness",
     bty = "l",
     pch = 19, col = pt_B, cex = 2, cex.lab = 2.7, cex.axis=2,
     col.axis = ax_col, col.lab = ax_col)
mtext(expression(
  atop(
    "Variance in species occurrences",
    paste("explained by traits (", R[plain(Traits)]^2, ")"))),
  side = 1, line = 7, cex = 1.65, col = ax_col)
double_curve(
  fun = function(x) p[1] + p[2]*x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col = col_B)
text(y = 42, x = 0.25, "F(1,19): 5.9\np: 0.03\nR²: 0.24",
     cex = 2.5, pos = 1, col = ax_col)
mtext("C", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)
x.seq <- seq(min(nn.1$R2.exp, na.rm = TRUE),
             max(nn.1$R2.exp, na.rm = TRUE), length.out = 200)
pred <- predict(m3, newdata = data.frame(R2.exp = x.seq),
                interval = "confidence")
polygon(c(x.seq, rev(x.seq)),
        c(pred[, "lwr"], rev(pred[, "upr"])), border = NA,
        col = adjustcolor(col_B, alpha.f = 0.15))
lines(x.seq,pred[, "fit"],col = col_B,lwd = 3)
lines(x.seq,pred[, "lwr"],col = col_B,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_B,  lty = 2, lwd = 1.5)



