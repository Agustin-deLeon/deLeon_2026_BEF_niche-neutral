###############################################################
###############################################################

# de León et al. 
# Biodiversity-function relationship determined
# by the niche-neutral gradient of community assembly

#=========================================================
### ------- NORMALIZED STOCHASTICITY RATIO (NST) ------- ###
#=========================================================
#
# This script quantifies the annual contribution of stochastic
# assembly processes using the Normalized Stochasticity Ratio
# (NST; Ning et al. 2019).
#
# Specifically, it:
#   (1) calculates annual NST values under three null models
#       (PF, PP and FF);
#   (2) evaluates NST using both annual and global regional
#       species pools;
#   (3) summarizes annual NST estimates across all analytical
#       scenarios; 
#   (4) compares annual NST values with the annual position
#       along the niche–neutral gradient estimated by CATS; and
#   (5) evaluates the robustness of the main results by
#       replacing the annual CATS estimates with the
#       Normalized Stochasticity Ratio (NST) as the predictor
#       of annual BEF.
#
# Required objects:
#   - br.muestreos_df: community matrix (sampling units × species)
#     including year, pond and sampling unit identifiers.
#   - nn.1: annual CATS results, including R2.exp.
#
# The object nn.1 must be generated previously by running
# CATS.R.
#=========================================================


################
# --- NST ------

library(NST)

years <- sort(unique(br.muestreos_df$año))

# especies
spp_cols <- colnames(br.muestreos_df)[-(1:4)]

# pool regional global
meta.freq.global <- colSums(
  br.muestreos_df[, spp_cols]
)

names(meta.freq.global) <- spp_cols

# modelos nulos
null_models <- c("PF", "PP", "FF")

# tipos de pool
pools <- c("annual", "global")

NST_out <- NULL

for(pool in pools){
  
  for(nm in null_models){
    
    cat("\n====================\n")
    cat("Pool:", pool, "\n")
    cat("Null model:", nm, "\n")
    cat("====================\n\n")
    
    for(i in years){
      
      cat("Año:", i, "\n")
      
      m <- br.muestreos_df[
        br.muestreos_df$año == i,
      ]
      
      comm <- as.matrix(m[, -(1:4)])
      
      rownames(comm) <- paste0(
        "C", m$Charco,
        "_M", m$Marco
      )
      
      grupo <- data.frame(
        group = rep("year", nrow(comm))
      )
      
      rownames(grupo) <- rownames(comm)
      
      # eliminar muestras vacías
      keep <- rowSums(comm) > 0
      
      comm  <- comm[keep, , drop = FALSE]
      grupo <- grupo[keep, , drop = FALSE]
      
      cat("Muestras vacías eliminadas:",
          sum(!keep), "\n")
      
      if(nrow(comm) < 3){
        cat("Muy pocas muestras, salteando\n")
        next
      }
      
      # pool anual
      if(pool == "annual"){
        
        nst <- tNST(
          comm = comm,
          group = grupo,
          dist.method = "jaccard",
          abundance.weighted = FALSE,
          rand = 200,
          null.model = nm,
          nworker = 4
        )
        
      }
      
      # pool global
      if(pool == "global"){
        
        mf <- matrix(
          meta.freq.global,
          nrow = 1,
          dimnames = list(
            "global",
            names(meta.freq.global)
          )
        )
        
        nst <- tNST(
          comm = comm,
          group = grupo,
          dist.method = "jaccard",
          abundance.weighted = FALSE,
          rand = 200,
          null.model = nm,
          meta.frequency = mf,
          nworker = 4
        )
        
      }
      
      NST_out <- rbind(
        NST_out,
        data.frame(
          año = i,
          pool = pool,
          null_model = nm,
          NST = nst$index.grp$NST.i,
          MST = nst$index.grp$MST.i,
          ST = nst$index.grp$ST.i,
          N_muestras = nrow(comm)
        )
      )
      
    }
  }
}

NST_out

library(ggplot2)

ggplot(
  NST_out,
  aes(año, NST, colour = null_model,
      linetype = pool)) +
  geom_line() +
  geom_point() +
  theme_bw()

aggregate(
  NST ~ pool + null_model,
  data = NST_out,
  range
)


library(tidyr)
library(ggplot2)
tmp <- merge(
  nn.1,
  subset(
    NST_out,
    pool == "global" &
      null_model == "PF"
  ),
  by = "año"
)
tmp <- tmp[, !duplicated(colnames(tmp))]

tmp_long <- rbind(
  data.frame(
    año = tmp$año,
    indice = "CATS",
    valor = tmp$R2.exp
  ),
  data.frame(
    año = tmp$año,
    indice = "NST",
    valor = tmp$NST
  )
)

ggplot(tmp_long,
       aes(x = año,
           y = valor,
           colour = indice)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    x = "Año",
    y = "Índice",
    colour = ""
  )

ggplot(tmp_long,
       aes(x = año,
           y = valor,
           colour = indice)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = "Año",
    y = "Índice",
    colour = ""
  )



library(ggplot2)
m1 <- lm(NST ~ R2.exp, data = tmp)

r2 <- summary(m1)$r.squared
p  <- summary(m1)$coefficients[2,4]
b0 <- coef(m1)[1]
b1 <- coef(m1)[2]

eq <- paste0(
  "NST = ",
  round(b0, 3),
  ifelse(b1 >= 0, " + ", " - "),
  round(abs(b1), 3),
  " × R²\n",
  "R² = ", round(r2, 3), "\n",
  "p = ", signif(p, 3)
)

ggplot(tmp,
       aes(x = R2.exp,
           y = NST)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.1,
    vjust = 1.1,
    label = eq,
    size = 5
  ) +
  theme_bw() +
  labs(
    x = expression(R[traits]^2),
    y = "NST"
  )


##########################################
### ------------ FIG S13. ------------- ##

par(mfrow = c(2,2), mar = c(9,7,5,2))

# ------ Panel A: NST vs CATS ------ #
m_nst <- lm(NST ~ R2.exp, data = tmp)
summary(m_nst)
p <- coefficients(m_nst)
plot(NST ~ R2.exp, data = tmp,
  ylab = "Normalized Stochasticity Ratio (NST)", cex.axis = 2,
  xlab = " ", pch = 19, col = pt_B, cex = 3,
  bty = "l", col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun = function(x) p[1] + p[2]*x,
  from = min(tmp$R2.exp, na.rm = TRUE),
  to = max(tmp$R2.exp, na.rm = TRUE),
  col = col_B)
mtext(xlab_exp, side = 1,line = 7,cex = 2,col = ax_col)
mtext("A", side = 3, line = 1.5, adj = 0, font = 2, cex = 3, col = ax_col)
x.seq <- seq(min(tmp$R2.exp, na.rm = TRUE),max(tmp$R2.exp, na.rm = TRUE),length.out = 200)
pred <- predict(m_nst,
  newdata = data.frame(R2.exp = x.seq),
  interval = "confidence")
polygon(c(x.seq, rev(x.seq)),
  c(pred[, "lwr"], rev(pred[, "upr"])),
  border = NA, col = adjustcolor(col_B, alpha.f = 0.15))
lines(x.seq,pred[, "lwr"],col = col_B,lty = 2,lwd = 1.5)
lines(x.seq,pred[, "upr"],col = col_B,lty = 2,lwd = 1.5)
s <- summary(m_nst)
Fval <- round(s$fstatistic[1], 1)
pval <- pf(s$fstatistic[1],
  s$fstatistic[2], s$fstatistic[3],
  lower.tail = FALSE)
r2 <- round(s$r.squared, 2)
text(x = 0.25, y = 0.35,
  labels = paste0(
    "F(",
    s$fstatistic[2], ",",
    s$fstatistic[3],
    "): ",
    Fval,
    "\n",
    "p: ",
    format.pval(pval, digits = 2),
    "\n",
    "R²: ",
    r2,
    "\n",
    "Cor: -0.65" 
  ),
  pos = 3,
  cex = 2,
  col = ax_col)
plot.new()


###############################################
# ---- Panel B: BEF ~ NST + precipitation --- #

tmp$BEF<- pendientes.lineal[,2] #BEF estimated with GLMM model
col_NST<- "#3A8EC1"
pt_NST<- "#3A8EC199"
m_BEF <- lm(BEF ~ NST + lluvia_anual, data = tmp)
p <- coefficients(m_BEF)

# ---------- BEF ~ NST  ----------
BEF.neutral <- with(tmp,
                    BEF - p["lluvia_anual"] * lluvia_anual +
                      mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))
plot(BEF.neutral ~ tmp$NST,
     ylab = "Richness-biomass relationship (BEF)",
     xlab = "",
     bty = "l", pch = 19, col = pt_NST, cex = 3, cex.axis = 2,
     col.axis = ax_col, col.lab = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * tmp$lluvia_anual, na.rm = TRUE) +
    p["NST"] * x,
  from = min(tmp$NST, na.rm = TRUE),to   = max(tmp$NST, na.rm = TRUE), col  = col_NST)

mtext("Normalized Stochasticity Ratio (NST)", side = 1, line = 3.75, cex = 2, col = ax_col)
mtext("B", side = 3, line = 1.5, adj = 0, font = 2, cex = 2.5, col = ax_col)

text(x = 0.82, y = 1, "F(2,18): 8.2\np: 0.003\nR²: 0.48", cex = 2,
     pos = 1, col = ax_col)
x.seq <- seq(min(tmp$NST, na.rm = TRUE),
             max(tmp$NST, na.rm = TRUE), length.out = 200)
pred <- predict(m_BEF,
                newdata = data.frame(NST = x.seq,
                                     lluvia_anual = mean(tmp$lluvia_anual, na.rm = TRUE)
                ),interval = "confidence")
pred[,c("fit","lwr","upr")] <- pred[,c("fit","lwr","upr")] -
  p["lluvia_anual"]*mean(tmp$lluvia_anual, na.rm=TRUE) +
  mean(p["lluvia_anual"]*tmp$lluvia_anual, na.rm=TRUE)
polygon(c(x.seq, rev(x.seq)),c(pred[,"lwr"], rev(pred[,"upr"])),
        border = NA,col = adjustcolor(col_NST, alpha.f = 0.15))
lines(x.seq, pred[,"fit"], col = col_NST, lwd = 3)
lines(x.seq, pred[,"lwr"], col = col_NST, lty = 2, lwd = 1.5)
lines(x.seq, pred[,"upr"], col = col_NST, lty = 2, lwd = 1.5)

# ---------- BEF ~ precip  ----------
BEF.precipitation <- with(tmp,
                          BEF - p["NST"] * NST +
                            mean(p["NST"] * NST, na.rm = TRUE))

plot(BEF.precipitation ~ tmp$lluvia_anual,
     xlab = " ", ylab = "",
     bty = "l", pch = 19, col = pt_NST, cex = 3, cex.axis = 2,
     col.axis = ax_col, col.lab = ax_col)

mtext("Precipitación anual", side = 1, line = 3.75, cex = 2, col = ax_col)
double_curve(
  fun  = function(x) p["(Intercept)"] + p["lluvia_anual"] * x +
    mean(p["NST"] * tmp$NST, na.rm = TRUE),
  from = min(tmp$lluvia_anual, na.rm = TRUE),
  to   = max(tmp$lluvia_anual, na.rm = TRUE),
  col  = col_NST)
x.seq <- seq(
  min(tmp$lluvia_anual, na.rm = TRUE),
  max(tmp$lluvia_anual, na.rm = TRUE),
  length.out = 200)
pred <- predict(m_BEF,
                newdata = data.frame(
                  NST = mean(tmp$NST, na.rm = TRUE),
                  lluvia_anual = x.seq),interval = "confidence")
pred[,c("fit","lwr","upr")] <-
  pred[,c("fit","lwr","upr")] -
  p["NST"]*mean(tmp$NST, na.rm=TRUE) +
  mean(p["NST"]*tmp$NST, na.rm=TRUE)
polygon(c(x.seq, rev(x.seq)), c(pred[,"lwr"], rev(pred[,"upr"])),
        border = NA, col = adjustcolor(col_NST, alpha.f = 0.15))
lines(x.seq, pred[,"fit"], col = col_NST, lwd = 3)
lines(x.seq, pred[,"lwr"], col = col_NST, lty = 2, lwd = 1.5)
lines(x.seq, pred[,"upr"], col = col_NST, lty = 2, lwd = 1.5)

