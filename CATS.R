###############################################################
###############################################################

# de León et al. 
# Biodiversity-function relationship determined
# by the niche-neutral gradient of community assembly

#=========================================================
### -------- COMMUNITY ASSEMBLY (CATS) ANALYSES ------- ###
#=========================================================
#
# This script estimates the annual position of plant
# metacommunities along the niche–neutral assembly gradient
# using the Community Assembly via Trait Selection (CATS)
# framework (Shipley et al.).
#
# Before running this script, execute:
#
#   Data_preprocessing.R
#
# which converts trait variables
# to the appropriate formats for the CATS analyses.
#
# This script:
#
#   (1) fits annual CATS models for each metacommunity;
#   (2) partitions the variation in species occurrences
#       explained by functional traits while accounting
#       for metacommunity abundances;
#   (3) calculates annual expected (R2.exp) and
#       unexplained (R2.unexp) trait effects, representing
#       each year's position along the niche–neutral
#       gradient; and
#   (4) exports the annual CATS estimates.
#
# Required objects (generated in 01_Prepare_data.R):
#
#   - br
#   - traits.pcoa.ejes
#   - pcoa_axes
#
# Main output:
#
#   - nn.1 (annual CATS estimates; includes R2.exp,
#     R2.unexp and model coefficients)
#
# The object nn.1 is subsequently used by:
#
#   - BEF.R
#   - Functional_diversity.R
#   - NST.R
#   - Main_figures.R
#=========================================================
########################################################
# Set directory
setwd("/Users/agustin/Desktop/Facultad/Maestria")

### Import databases 
library(readxl)
br<-read_excel("br.xlsx")                       # br is the matrix of species occurrences
br<-matrix(unlist(br),byrow = F, nrow=nrow(br), dimnames = list(rownames(br),colnames(br)))       # transform to matrix

bm.nuevo=NULL
for(i in 2005:2025){
  as.matrix(br[which(br[,1]==i),]) -> bm.temp
  if(i==2006){bm.temp[which(bm.temp[,2]==6),]-> bm.temp}
  if(i==2007){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2008){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2009){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  if(i==2010){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  rbind(bm.nuevo, bm.temp)-> bm.nuevo
}
as.matrix(bm.nuevo)-> br

# Import trait databases
# traits: numeric version used for Gower distances and PCoA. Obtained in "Data_preprocessing.R"
traits<- read.csv("traits_numerico_n.csv")
traits[-c(15,35,119),-1]-> traits #when i import it it adds an extra col, i eliminate it

# traits_n: original categorical version used for trait interpretation (e.g. envfit)
traits_n <-read_excel("Traits.xlsx")  

# Edit the format of trait variables (factors)
names.traits<-traits[,1]
traits<-matrix(as.numeric(unlist(traits)),byrow = F, 
               nrow=nrow(traits), dimnames = list(rownames(traits),
                                                  colnames(traits)))
rownames(traits) <- as.vector(unlist(names.traits))
traits <- traits[,-1]
traits<-as.data.frame(traits)
sapply(traits, class) ### ok, numeric
### pass categorical variables to factors
# colnames(traits)
cols <- c(1,3, 4, 5, 7, 11, 12, 13, 14, 15, 16,17) ## factors
cols.ord <- c(2, 6, 8, 18) ## ordered factors :  "Seed size"  "Leaf size"  "Sculthorpe" "Raunkiaer" 
colnames(traits[,cols.ord]); 
colnames(traits[,cols]) ## check
for(i in cols.ord){
  traits[,i] <- ordered(traits[,i])
} 
for(i in cols){
  traits[,i] <- factor(traits[,i])
} 

######################### PCoA ############################
library(ape)
library(cluster)
library(scales)
library(vegan)

### Gower distances
### computing distances among species based on their functional traits
spp_plants<-rownames(traits)
# spp_plants_viejo <- rownames(traits.pcoa)
# rownames(traits)<-spp_plants #### HAY QUE HACER LISTA CORTA DE ESTOS NOMBRES Y TEXTO PARA PONER EN LA FIG
rownames(traits) <- 1:122
D <- daisy(traits, metric = "gower", stand = TRUE) ### species x species distances (100 x 100)
# D1 <- daisy(traits, metric = "gower")
# identical(D, D1); rm(D1) ## ok
PCoA<- pcoa(D, rn = rownames(traits))
traits.pcoa.ejes<-PCoA$vectors
pcoa_axes <- as.data.frame(PCoA$vectors[, 1:5])  # Ejes 1 a 5
colnames(pcoa_axes) <- c("Axis1", "Axis2", "Axis3", "Axis4", "Axis5")
rownames(pcoa_axes)<- names.traits
proporciones <- PCoA$values$Relative_eig
porcentaje_acumulado <- cumsum(proporciones) * 100
print(proporciones)
print(porcentaje_acumulado)


library(glmmTMB)
library(MuMIn)
library(performance)
out.final=NULL
M<- br[,-c(19,39,125)] #Delate spp with 0 occurrences in selected sampling periods
apply(br[,c(19,39,125)],2,sum) #Check 0 anundance: Cirsium_vulgare, Enydra_sessilis, Stellaria_media 
years<-sort(unique(M[,1]))                                       # Ponds observed 
metacomm<-apply(M[,5:ncol(M)],2,sum)
#  id.spp<-match(colnames(M[,5:ncol(M)]),rownames(Traits))
#  Traits<-Traits[id.spp,]
#  metacomm<-metacomm[ii.observed]           
Traits<- traits.pcoa.ejes[,1:5] # reduce pool to species observed on this year
#  Traits<-Traits[ii.observed,]                                     # reduce traits matrix to observed species
it= 200 # lo ideal, lo estoy corriendo con 10 porque demora jeje
for(i in years){ # }, .combine=rbind)%dopar%  {
  charco<-as.numeric(paste(M[which(M[,1]==i),3]))                    # Create a vector indicating survey times
  metacomm<-apply(M[,5:ncol(M)],2,sum)                             # Metacommunity abundances: uso el global, algo a considerar
  metacomm.U<-rep(sum(metacomm)/length(metacomm),length(metacomm)) # Uniform expectation
  print(i)
  out<-NULL
  sigma.obs.nt<-NULL                                               # Null objects for deviance decomposition...
  sigma.obs.ut<-NULL                                               # 
  sigma.null.np<-NULL                                              # 
  sigma.obs.up<-NULL                                               # 
  m<-M[which(M[,1]==i),]                                         # pond matrix
  ss<-charco                                 # sampling times
  m.temp<-NULL                                                   # temporal out
  # if more than a sample unit was observed
  for(j in 1:nrow(m)){                                         # for each sample unit
    cat("row ",c(j, "of year",i),"\n")    
    y<-m[j,5:ncol(m)]
    m.temp.2<-cbind(metacomm,ss[j], Traits, metacomm.U,y)
    rbind(m.temp, m.temp.2)-> m.temp
  }
  #m.temp[1:2214,]-> m.temp
  m.temp <- na.omit(m.temp) 
  colnames(m.temp)[2]<-"charco.id"
  colnames(m.temp)[ncol(m.temp)]<-"y"
  colnames(m.temp)[3:7]<- c("Axis1" ,    "Axis2"    , "Axis3"    , "Axis4"   ,  "Axis5")
  as.data.frame(m.temp)-> m.temp
  # Variables defined in the space bcecause glmmTMB can not use "data=", SI DEJA, solo que tiene que ser un data.frame      
  rm(Axis1, Axis2, Axis3, Axis4, metacomm, metacomm.U)
  attach(m.temp)
  
  print("ok?")
  ######### FIRST MODEL
  # neutral mass effect and observed traits 
  M1.n.t<-glmmTMB(y~Axis1+Axis2+Axis3+Axis4+Axis5+
                    I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                  offset=log(metacomm),
                  data=m.temp, family=binomial)
  print("ok?")
  
  M1.n.t_r2<- r2(M1.n.t)
  r2.n.t<-as.numeric(M1.n.t_r2$R2_marginal)
  # if(adjustr2.to.random==T)if(is.na(M1.n.t_r2$R2_conditional)==F)r2.n.t<-r2.n.t/(1-(M1.n.t_r2$R2_conditional-r2.n.t))
  sel.coeff<-as.vector(unlist(fixef(M1.n.t)$cond))
  print("model1 ok")
  
  ##### SECOND MODEL
  # neutral mass effect and observed traits 
  M1.u.t<-glmmTMB(y~Axis1+Axis2+Axis3+Axis4+Axis5+
                    I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                  offset=log(metacomm.U),
                  data=m.temp, family=binomial)
  sigma.obs.ut<-c(sigma.obs.ut,M1.u.t$sigma) #desviación estandar? para qué la quiere
  #M0.u.t<-glmmTMB(y~1 +(1|year.month),
  #                offset=log(metacomm.U),
  #                data=NULL, family=binomial)
  M1.u.t_r2<-r2(M1.u.t)
  r2.u.t<-as.numeric(M1.u.t_r2$R2_marginal)
  #  if(adjustr2.to.random==T)if(is.na(M1.u.t_r2$R2_conditional)==F)r2.u.t<-r2.u.t/(1-(M1.u.t_r2$R2_conditional-r2.u.t))
  print("model2 ok")
  
  r2.n.p<-NULL #no se pierden si se borran ya acá? ok no, estos son p de permutados (los next)
  r2.u.p<-NULL
  ############
  #####. Randomized traits
  for(rand in 1:it){ #cuantas veces aleatorizo los traits creo que sería
    print(rand)
    y.p<-m.temp[sample(nrow(m.temp)),ncol(m.temp)]
    m.temp.p<-cbind(m.temp[,-ncol(m.temp)],y.p)
    
    #######THIRD MODEL        
    # neutral mass effect and permuted traits         
    M1.n.p<-glmmTMB(y.p~Axis1+Axis2+Axis3+Axis4+Axis5+
                      I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                    offset=log(metacomm),
                    data=NULL, family=binomial)
    M1.n.p_r2<-r2(M1.n.p)
    #r.squaredGLMM(M1.n.p)[1,]-> M1.n.p_r2
    r2.n.p.temp<-as.numeric(M1.n.p_r2$R2_marginal)
    #r2.n.p.temp<-as.numeric(M1.n.p_r2[1])
    # if(adjustr2.to.random==T) if(is.na(M1.n.p_r2$R2_conditional)==F)r2.n.p.temp<-r2.n.p.temp/(1-(M1.n.p_r2$R2_conditional-r2.n.p.temp))
    #   if(adjustr2.to.random==T) if(is.na(M1.n.p_r2[2])==F)r2.n.p.temp<-r2.n.p.temp/(1-(M1.n.p_r2[2]-r2.n.p.temp))
    ######### FOURTH MODEL
    # Uniform metacomm and permuted traits      
    M1.u.p<-glmmTMB(y.p~Axis1+Axis2+Axis3+Axis4+Axis5+
                      I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                    offset=log(metacomm.U),
                    data=NULL, family=binomial)
    M1.u.p_r2<-r2(M1.u.p)
    #r.squaredGLMM(M1.u.p)[1,]-> M1.u.p_r2
    r2.u.p.temp<-as.numeric(M1.u.p_r2$R2_marginal)
    #r2.u.p.temp<-as.numeric(M1.u.p_r2[1])
    #if(adjustr2.to.random==T)if(is.na(M1.u.p_r2$R2_conditional)==F)r2.u.p.temp<-r2.u.p.temp/(1-(M1.u.p_r2$R2_conditional-r2.u.p.temp))
    #   if(adjustr2.to.random==T) if(is.na(M1.u.p_r2[2])==F)r2.u.p.temp<-r2.u.p.temp/(1-(M1.u.p_r2[2]-r2.u.p.temp))
    r2.n.p<-c(r2.n.p,r2.n.p.temp)
    r2.u.p<-c(r2.u.p,r2.u.p.temp)
  }
  r2.n.p<-mean(r2.n.p) #ok, me quedo con la media de todas las veces que lo corrí
  r2.u.p<-mean(r2.u.p)
  
  # Following Shipley 2014 R2 values are adjusted
  r2.u.t<-max(r2.u.t,r2.u.p)  # R2 with observed traits could not be inferior to R2 with randomized traits
  r2.n.p<-max(r2.n.p,r2.u.p)  # R2 with observed metacomm could not be inferior to R2 with Uniform metacomm
  r2.n.t<-max(r2.n.t,r2.u.t)   # metacommunity prior and observed trait adjusted. It should not be less than R2.u.t.
  out<-rbind(out,c(as.numeric(i),r2.n.t,r2.u.t,r2.n.p,r2.u.p, as.vector(sel.coeff)))
  
  
  colnames(out)<-c("year","r2.n.t","r2.u.t","r2.n.p","r2.u.p", "intercept","axis1", "axis2", "axis3", "axis4", "axis5", "axis1_2", "axis2_2", "axis3_2", "axis4_2", "axis5_2") #VER bien que son estos valores
  out-> out.temp
  out.final<- rbind(out.final, out.temp)
  print(date())
}

as.data.frame(out.final)-> nn
nn.1<-cbind(nn[,1:5],(1-(nn[,2]))/(1-nn[,5]),nn[,(6:ncol(nn))])
colnames(nn.1)[6]<-"R2.unexp"        # Estimates DRIFT corrected by expected drift (R2) from randomized trait 
head(nn.1)
nn.1$R2.exp<- 1-(nn.1$R2.unexp)
nn.1-> CATS_salida
write.csv(CATS_salida,"CATS_2025.csv")

####################################
# Armar las matrices a utilizar
####################################

setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
library(openxlsx)
library(glmmTMB)
library(quantreg)  
library(ggplot2) 
library(car)
library(visreg)
library(MuMIn) 
library(ggpubr)
library(gridExtra)
library(caret)
library(broom.mixed)
library(ggeffects)
library(DHARMa)
library(bestglm)
read.xlsx("bm.xlsx")-> bm
read.xlsx("br.xlsx")-> br2

bm[,7:10]<- c(rep(1, nrow(bm)))
colnames(bm)<-  c("año","mes","charco","um","riq","biom","lluvia_muestreo","lluvia_anual","temp_muestreo", "temp_anual")  

read.xlsx("Ambientales_resumen.xlsx")-> Ambiente
read.xlsx("Environmental.xlsx")-> Environment
bm[,11:14]<- c(rep(1, nrow(bm)))
for (i in unique(Environment[,1])) {
  print(i)
  bm_subset <- bm[bm[,3] == i, 11:14]
  env_subset <- Environment[Environment[,1] == i, c(1:18)]
  bm[bm[,3] == i, 11:28]<- env_subset
}
colnames(bm)
colnames(bm)<-  c("año","mes","charco","um","riq","biom","lluvia_muestreo","lluvia_anual","temp_muestreo", "temp_anual",
                  "Pond.id"  ,      "X"           ,    "Y"       ,        "DM"   ,           "ddmm"     ,      
                  "Shape"     ,      "Islands"    ,     "log.Area"     ,   "log.Volumen"    , "Mean.Depth"  ,   
                  "Sd.Depth"    ,    "CV.Depth"   ,     "Hydroperiod"   ,  "Degree"      ,    "log.Betweenness",
                  "Closenness"   ,   "grado.percol" ,   "clos.percol"   )  
#Keep sampling events of interest 
bm.nuevo=NULL
for(i in 2005:2025){
  as.matrix(bm[which(bm[,1]==i),]) -> bm.temp
  if(i==2006){bm.temp[which(bm.temp[,2]==6),]-> bm.temp}
  if(i==2007){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2008){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2009){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  if(i==2010){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  rbind(bm.nuevo, bm.temp)-> bm.nuevo}

as.data.frame(bm.nuevo)-> bm.nuevo
ambientales<- read.xlsx("Ambientales_resumen.xlsx")
ambientales[,c(1:3, 9:16)]-> ambientales 
amb.nuevo=NULL
for(i in 2005:2025){
  as.matrix(ambientales[which(ambientales[,1]==i),]) -> amb.temp
  if(i==2006){amb.temp[which(amb.temp[,2]==6),]-> amb.temp}
  if(i==2007){amb.temp[which(amb.temp[,2]==8),]-> amb.temp}
  if(i==2008){amb.temp[which(amb.temp[,2]==8),]-> amb.temp}
  if(i==2009){amb.temp[which(amb.temp[,2]==7),]-> amb.temp}
  if(i==2010){amb.temp[which(amb.temp[,2]==7),]-> amb.temp}
  rbind(amb.nuevo, amb.temp)-> amb.nuevo
  amb.temp=NULL}
rm(amb.temp, ambientales)
as.data.frame(amb.nuevo)-> amb.nuevo

bm.nuevo[,29:36]<- c(rep(1, nrow(bm.nuevo)))
bm.act=NULL
for (i in unique(amb.nuevo[,1])) {
  print(i)
  bm.temp<- bm.nuevo[which(bm.nuevo[,1]==i),] #divido para cada año
  amb.temp<- amb.nuevo[which(amb.nuevo[,1] == i),]
  for(ii in unique(bm.temp[,3])){ #ii es id de cada charco de ese año
    print(ii)
    bm_subset <- bm.temp[which(bm.temp[,3] == ii), 29:36]
    env_subset <- amb.temp[which(amb.temp[,3] == ii), c(4:11)]
    if(length(which(amb.temp[,3] == ii))>0){
      bm.temp[which(bm.temp[,3] == ii), 29:36]<- env_subset } 
    else{ bm.temp[which(bm.temp[,3] == ii), 29:36]<- NA }} 
  rbind(bm.act, bm.temp)-> bm.act }
colnames(bm.act)
colnames(bm.act)<-  c("año","mes","charco","um","riq","biom","lluvia_muestreo","lluvia_anual","temp_muestreo", "temp_anual",
                      "Pond.id"  ,      "X"           ,    "Y"       ,        "DM"   ,           "ddmm"     ,      
                      "Shape"     ,      "Islands"    ,     "log.Area"     ,   "log.Volumen"    , "Mean.Depth"  ,   
                      "Sd.Depth"    ,    "CV.Depth"   ,     "Hydroperiod"   ,  "Degree"      ,    "log.Betweenness",
                      "Closenness"   ,   "grado.percol" ,   "clos.percol" , 
                      "cortes"   ,  "prof.media", "sd.prof" ,   "CVprof"   ,  "area"  ,
                      "log.area" ,  "vol",   "log.vol"   )  
bm.nuevo.repuesto<- bm.nuevo
bm.nuevo<- bm.act
bm.nuevo.temp<- bm.nuevo 

bm.st.total=NULL
for(i in 2005:2025){
  as.matrix(bm[which(bm[,1]==i),]) -> bm.temp
  if(i==2006){bm.temp[which(bm.temp[,2]==6),]-> bm.temp}
  if(i==2007){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2008){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2009){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  if(i==2010){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  
  for(ii in unique(bm.temp[,3])){
    sum(bm.temp[which(bm.temp[,3]==ii),6])-> biomasa
    max(bm.temp[which(bm.temp[,3]==ii),4])-> um
    (biomasa/um)-> biomasa.st
    c(bm.temp[max(which(bm.temp[,3]==ii)),1:3], um , biomasa, biomasa.st)-> bm.st
    rbind(bm.st.total, bm.st)-> bm.st.total
  }
}

br.muestreos =NULL
for(i in 2005:2025){
  as.matrix(br2[which(br2[,1]==i),]) -> br2.temp
  if(i==2006){br2.temp[which(br2.temp[,2]==6),]-> br2.temp}
  if(i==2007){br2.temp[which(br2.temp[,2]==8),]-> br2.temp}
  if(i==2008){br2.temp[which(br2.temp[,2]==8),]-> br2.temp}
  if(i==2009){br2.temp[which(br2.temp[,2]==7),]-> br2.temp}
  if(i==2010){br2.temp[which(br2.temp[,2]==7),]-> br2.temp}
  rbind(br.muestreos, br2.temp)-> br.muestreos
}

charcos.evaluar <- c("1", "2", "3", "4", "5", "6", "7", "8", 
                     "10", "11", "13", "14", "15", "16", "17", "12", 
                     "21", "24", "25", "26", "27", "29", "30", "32", 
                     "33", "38", "40", "41", "42", "43", "44", "45", 
                     "47", "48", "49", "50", "51", "55", "56", "666", 
                     "91", "10022", "21022", "134", "137")

#I replace the missing data of the environmental variables with the average of that puddle for the total of the samples taken.
bm.nuevo$cortes <- ifelse(is.na(bm.nuevo$cortes), bm.nuevo$Islands, bm.nuevo$cortes)
bm.nuevo$prof.media <- ifelse(is.na(bm.nuevo$prof.media), bm.nuevo$Mean.Depth, bm.nuevo$prof.media)
bm.nuevo$CVprof <- ifelse(is.na(bm.nuevo$CVprof), bm.nuevo$CV.Depth, bm.nuevo$CVprof)
bm.nuevo$log.area <- ifelse(is.na(bm.nuevo$log.area), bm.nuevo$log.Area, bm.nuevo$log.area)
bm.nuevo$sd.prof <- ifelse(is.na(bm.nuevo$sd.prof), bm.nuevo$Sd.Depth, bm.nuevo$sd.prof)


profundidad<- NULL
for(ii in 2005:2025){
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:16,29,34,30:32,23:25)]
  colnames(df.um.ch) <- c("um.año", "ch.id", "um.biom", "um.rich", "lluvia_muestreo", "lluvia_anual", "temp_muestreo", 
                          "temp_anual", "DM", "ddmm", "Shape", "Islands", "log.Area", "Mean.Depth", 
                          "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness")
  df.um.ch <- na.omit(df.um.ch)
  mean(df.um.ch$Mean.Depth)-> media_prof
  sd(df.um.ch$Mean.Depth)-> sd_prof
  (media_prof)/(sd_prof)-> media_cv
  cbind(ii, media_prof, media_cv)-> prof.temp
  rbind(profundidad, prof.temp)-> profundidad
}

