########################################################
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/CATS")

### Import databases 
###.  br is available for review after published will be available in Dyad with embargo time
# import biological database
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

# Import environmental information
Environment<- read_excel("Environmental_variables.xlsx")
Environment<-matrix(unlist(Environment),byrow = F, nrow=nrow(Environment), 
  dimnames = list(rownames(Environment),colnames(Environment)))

# Import species traits data base
traits <- traits %>%
  dplyr::select(
    Species_ID,
    Status,
    `Seed size`,
    `Seed shape`,
    `Dispersal syndrome`,
    `Vegetative sperad`,
    `Leaf size`,
    Habit,
    Sculthorpe,
    `Plant height`,
    `Stem length`,
    Persistence,
    `Reproductive period`,
    `Stem type`,
    `Nitrogen fixation`,
    `Anaerobiosis tolerance`,
    `Drougth tolerance`,
    `Photosynthetic path`,
    Raunkiaer
  )
as.data.frame(traits)-> traits

traits[,1]-> temp.species
unique(traits[,2]) #  "Native"       "Exotic"       "Cosmopolitan"
traits[,2]-> temp.status
ifelse(temp.status=="Native", 1, temp.status)-> temp.status
ifelse(temp.status=="Exotic", 2, temp.status)-> temp.status
ifelse(temp.status=="Cosmopolitan", 3, temp.status)-> temp.status
unique(traits[,3]) #  "1-3mm" "<1mm"  ">5mm"  "3-5mm"
traits[,3]-> temp.seed.size
ifelse(temp.seed.size=="<1mm", 1, temp.seed.size)-> temp.seed.size
ifelse(temp.seed.size=="1-3mm", 2, temp.seed.size)-> temp.seed.size
ifelse(temp.seed.size=="3-5mm", 3, temp.seed.size)-> temp.seed.size
ifelse(temp.seed.size==">5mm", 4, temp.seed.size)-> temp.seed.size
unique(traits[,4]) #  "Winged or pappus" "Compressed"  "Ellipsoid"  "Spherical"  "Elongated"  "Polyhedral"  
traits[,4]-> temp.seed.shape
ifelse(temp.seed.shape=="Winged or pappus", 1, temp.seed.shape)-> temp.seed.shape
ifelse(temp.seed.shape=="Compressed", 2, temp.seed.shape)-> temp.seed.shape
ifelse(temp.seed.shape=="Ellipsoid", 3, temp.seed.shape)-> temp.seed.shape
ifelse(temp.seed.shape=="Spherical", 4, temp.seed.shape)-> temp.seed.shape
ifelse(temp.seed.shape=="Elongated", 5, temp.seed.shape)-> temp.seed.shape
ifelse(temp.seed.shape=="Polyhedral", 6, temp.seed.shape)-> temp.seed.shape
unique(traits[,5]) #  "Anemochory" "Barochory" 
traits[,5]-> temp.disp
ifelse(temp.disp=="Anemochory", 0, temp.disp)-> temp.disp
ifelse(temp.disp=="Barochory", 1, temp.disp)-> temp.disp
unique(traits[,6]) #  "Yes" "No" 
traits[,6]-> temp.veg.sperad
ifelse(temp.veg.sperad=="Yes", 1, temp.veg.sperad)-> temp.veg.sperad
ifelse(temp.veg.sperad=="No", 0, temp.veg.sperad)-> temp.veg.sperad
unique(traits[,7]) #  "2-5cm2"   "5-10cm2"  "<2cm2"    "Leafless" ">10cm2"  
traits[,7]-> temp.leaf.size
ifelse(temp.leaf.size=="Leafless", 0, temp.leaf.size)-> temp.leaf.size
ifelse(temp.leaf.size=="<2cm2", 1, temp.leaf.size)-> temp.leaf.size
ifelse(temp.leaf.size=="2-5cm2", 2, temp.leaf.size)-> temp.leaf.size
ifelse(temp.leaf.size=="5-10cm2", 3, temp.leaf.size)-> temp.leaf.size
ifelse(temp.leaf.size==">10cm2", 4, temp.leaf.size)-> temp.leaf.size
unique(traits[,8]) #  "Decumbent"   "Postrate" "Rhizomatous or Cespitose-stoloniferous" "Supported by water" "Erect" "Erect cespitose" "Rosette"  "Postrate cespitose" 
traits[,8]-> temp.habit
ifelse(temp.habit=="Decumbent", 1, temp.habit)-> temp.habit
ifelse(temp.habit=="Postrate", 2, temp.habit)-> temp.habit
ifelse(temp.habit=="Rhizomatous or Cespitose-stoloniferous", 3, temp.habit)-> temp.habit
ifelse(temp.habit=="Supported by water", 4, temp.habit)-> temp.habit
ifelse(temp.habit=="Erect", 5, temp.habit)-> temp.habit
ifelse(temp.habit=="Erect cespitose", 6, temp.habit)-> temp.habit
ifelse(temp.habit=="Rosette", 7, temp.habit)-> temp.habit
ifelse(temp.habit=="Postrate cespitose", 8, temp.habit)-> temp.habit
unique(traits[,9]) # "No hydrophyte"  "Rooted emergent"  "Amphibious"  "Free floating" "Rooted submerged" "Free submerged" "Rooted floating leaves"
traits[,9]-> temp.sculthorpe
ifelse(temp.sculthorpe=="No hydrophyte", 0, temp.sculthorpe)-> temp.sculthorpe
ifelse(temp.sculthorpe=="Rooted emergent", 1, temp.sculthorpe)-> temp.sculthorpe
ifelse(temp.sculthorpe=="Amphibious", 2, temp.sculthorpe)-> temp.sculthorpe
ifelse(temp.sculthorpe=="Free floating", 3, temp.sculthorpe)-> temp.sculthorpe
ifelse(temp.sculthorpe=="Rooted submerged", 4, temp.sculthorpe)-> temp.sculthorpe
ifelse(temp.sculthorpe=="Free submerged", 5, temp.sculthorpe)-> temp.sculthorpe
ifelse(temp.sculthorpe=="Rooted floating leaves", 6, temp.sculthorpe)-> temp.sculthorpe
traits[,10]-> temp.plant.height
traits[,11]-> temp.stem.length
unique(traits[,12]) # "Perennial" "Annual" "Biennial" 
traits[,12]-> temp.persistence
ifelse(temp.persistence=="Perennial", 1, temp.persistence)-> temp.persistence
ifelse(temp.persistence=="Annual", 2, temp.persistence)-> temp.persistence
ifelse(temp.persistence=="Biennial", 3, temp.persistence)-> temp.persistence
unique(traits[,13]) # "Spring"     "Summer"     "No flowers" "Winter"     "Autumn"  
traits[,13]-> temp.rep.period
ifelse(temp.rep.period=="No flowers", 0, temp.rep.period)-> temp.rep.period
ifelse(temp.rep.period=="Summer", 1, temp.rep.period)-> temp.rep.period
ifelse(temp.rep.period=="Autumn", 2, temp.rep.period)-> temp.rep.period
ifelse(temp.rep.period=="Winter", 3, temp.rep.period)-> temp.rep.period
ifelse(temp.rep.period=="Spring", 4, temp.rep.period)-> temp.rep.period
unique(traits[,14]) # "Herbaceous" "Sufrutice"  "Shrub" 
traits[,14]-> temp.stem.type
ifelse(temp.stem.type=="Herbaceous", 1, temp.stem.type)-> temp.stem.type
ifelse(temp.stem.type=="Sufrutice", 2, temp.stem.type)-> temp.stem.type
ifelse(temp.stem.type=="Shrub", 3, temp.stem.type)-> temp.stem.type
unique(traits[,15]) #  "Yes" "No" 
traits[,15]-> temp.N.fixation
ifelse(temp.N.fixation=="Yes", 1, temp.N.fixation)-> temp.N.fixation
ifelse(temp.N.fixation=="No", 0, temp.N.fixation)-> temp.N.fixation
unique(traits[,16]) # "Low"  "High" "No"
traits[,16]-> temp.anaerobiosis.tol
ifelse(temp.anaerobiosis.tol=="No", 0, temp.anaerobiosis.tol)-> temp.anaerobiosis.tol
ifelse(temp.anaerobiosis.tol=="Low", 1, temp.anaerobiosis.tol)-> temp.anaerobiosis.tol
ifelse(temp.anaerobiosis.tol=="High", 2, temp.anaerobiosis.tol)-> temp.anaerobiosis.tol
unique(traits[,17]) # "Low"  "High" "No"
traits[,17]-> temp.drought.tol
ifelse(temp.drought.tol=="No", 0, temp.drought.tol)-> temp.drought.tol
ifelse(temp.drought.tol=="Low", 1, temp.drought.tol)-> temp.drought.tol
ifelse(temp.drought.tol=="High", 2, temp.drought.tol)-> temp.drought.tol
unique(traits[,18]) # "C3"  "C4"  "CAM"
traits[,18]-> temp.photo
ifelse(temp.photo=="C3", 1, temp.photo)-> temp.photo
ifelse(temp.photo=="C4", 2, temp.photo)-> temp.photo
ifelse(temp.photo=="CAM", 3, temp.photo)-> temp.photo
unique(traits[,19]) #  "Chamaephyte"     "Hemicryptophyte" "Therophyte"      "Cryptophyte"     "Phanerophyte"  
traits[,19]-> temp.raunkiaer
ifelse(temp.raunkiaer=="Chamaephyte", 1, temp.raunkiaer)-> temp.raunkiaer
ifelse(temp.raunkiaer=="Hemicryptophyte", 2, temp.raunkiaer)-> temp.raunkiaer
ifelse(temp.raunkiaer=="Therophyte", 3, temp.raunkiaer)-> temp.raunkiaer
ifelse(temp.raunkiaer=="Cryptophyte", 4, temp.raunkiaer)-> temp.raunkiaer
ifelse(temp.raunkiaer=="Phanerophyte", 5, temp.raunkiaer)-> temp.raunkiaer

cbind(temp.species, temp.status, temp.seed.size, temp.seed.shape,  temp.disp, 
      temp.veg.sperad, temp.leaf.size, temp.habit, temp.sculthorpe, temp.plant.height, temp.stem.length, temp.persistence,
      temp.rep.period, temp.stem.type, temp.N.fixation, temp.anaerobiosis.tol, temp.drought.tol,  temp.photo, temp.raunkiaer)-> traits.nuevos
colnames(traits.nuevos)<- nombres
write.csv(traits.nuevos, "traits_numerico.csv")
traits.nuevos-> traits
traits<- read.csv("traits_numerico.csv")

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
### traits <- traits[order(rownames(traits)), ] ## esto no porque sino queda mal
## respecto a los números que puso Ana

######################### PCoA ############################
# The following lines estimate PCoA axis. as we do for cats analyses
library(ape)
library(cluster)
library(scales)
library(vegan)

### Gower distances
### computing distances among species based on their functional traits
spp_plants<-rownames(traits)
# spp_plants_viejo <- rownames(traits.pcoa)
# rownames(traits)<-spp_plants #### HAY QUE HACER LISTA CORTA DE ESTOS NOMBRES Y TEXTO PARA PONER EN LA FIG
rownames(traits) <- 1:123
D <- dyaisy(traits, metric = "gower", stand = TRUE) ### species x species distances (100 x 100)
# D1 <- daisy(traits, metric = "gower")
# identical(D, D1); rm(D1) ## ok
PCoA<- pcoa(D, rn = rownames(traits))
traits.pcoa.ejes<-PCoA$vectors

###########################################
### axis correlation to traits ####
## trait variables were related to the ordination considering the first 4 axes ###
## factor variables were passed to dummy variables and resacled to range 0-1
## for simplicity, import character data frame
###### traits_character
traits_char<- read.csv("traits_usados.csv")
traits.aux <-read_excel("traits.xlsx", sheet = "Traits_todos")
as.data.frame(cbind(traits.aux[,1],traits_char[,-1]))-> traits_char

traits_char <- as.data.frame(traits_char)
rownames(traits_char) <- traits_char[,1]
# rownames(traits_char) == spp_plants # ok
traits_char <- traits_char[,-1]
library(fastDummies)
traits_char_dummy <- dummy_columns(traits_char, ignore_na = TRUE, 
                                   remove_selected_columns = TRUE)
# a <-dummy.data.frame(traits_char, sep = "_", dummy.class = c("character"))
a.d <- decostand(traits_char_dummy, 2, method = "range", na.rm = TRUE) #Estandariza
# 1. Verificar y filtrar columnas
filtered_a.d <- a.d[, !grepl("_NA", colnames(traits_char_dummy))]
# 2. Comprobar que quedan columnas válidas
if (ncol(filtered_a.d) == 0) {
  stop("No quedan columnas válidas en a.d después del filtrado.")
}
# 3. Ejecutar envfit solo si hay columnas válidas
fit <- envfit(PCoA$vectors, filtered_a.d, na.rm = TRUE, choices = 1:4)
# 4. Revisar los resultados
print(fit)

fit_barplot <- fit$vectors$arrows[which(fit$vectors$pvals <= 0.05), ]
### plot correlations bigger than 0,7 ########
#### axis 1 ###
par(mfrow=c(1,1))
par(mar = c(14, 8, 2, 5))
barplot(fit_barplot[which(abs(fit_barplot[, 1]) >= 0.6),1], las = 2)
#guardar pdf 6 x 4 landscape
#### axis 2 ###
barplot(fit_barplot[which(abs(fit_barplot[, 2]) >= 0.6), 2],
        horiz = FALSE, las = 2, cex.names = 1, cex.axis = 1)
#### axis 3 ### lo mismo pero en distinto orden las barras
par(mar = c(14, 9, 2, 5))
barplot(fit_barplot[which(abs(fit_barplot[, 3]) >= 0.6), 3], las = 2)
## axis 4 ###
par(mar = c(14, 10, 2, 5))
barplot(fit_barplot[which(abs(fit_barplot[, 4]) >= 0.6), 4], las = 2)

pcoa_axes <- as.data.frame(PCoA$vectors[, 1:5])  # Ejes 1 a 4
colnames(pcoa_axes) <- c("Axis1", "Axis2", "Axis3", "Axis4", "Axis5")

