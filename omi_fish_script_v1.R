## R code to analyze the relationship between niche measures and abundance 
## Prepared by Tadeu Siqueira
## Rio Claro, 11/April/2019


library(ade4)
library(vegan)
library(betareg)


# data
# species abundance (columms) data across sites (rows) 
(abund <- read.csv("abundancia_tudo.csv", header = T, row.names = 1, sep = ",", 
                   dec = ","))

# abiotic variables (columms) data across sites (rows)
# these are used to estimate niche metrics
(env <- read.csv("ambientais_tudo.csv", header = T, row.names = 1, sep = ",", 
                 dec = ","))

# species trait data; species names as rows; traits as collums
(traits <- read.csv("dados_modelo.csv", header = T, row.names = 1, sep = ",", 
                    dec = ","))

# check general structure of the data frame
str(env)
apply(env, 1, mean)
apply(env, 2, mean)

# function for resampling 
omi.resamp <- function (bio, env, traits) {

# I want to split the data in two
# One to estimate abundance and occupancy 
# and the other to estimate niche measures (which depend on abund and occur data)

nsites <- 1:nrow(bio)
sites.abund <- sample(nsites, round(nrow(bio)/2), F)
sites.niche <- nsites[-sites.abund]

# Just to make sure we are using sites and variables with information
goodcols.nic <- colSums(bio[sites.niche,])>0
goodrows.nic <- rowSums(bio[sites.niche,])>0

goodcols.abund <- colSums(bio[sites.abund,])>0
goodrows.abund <- rowSums(bio[sites.abund,])>0

# these are use to estimate omi
niche.sample <- bio[sites.niche,][goodrows.nic, goodcols.nic]
env.sample <- env[sites.niche,][goodrows.nic,]

# estimate niche measures 
dudi.lo <- dudi.pca(env.sample, scale = TRUE, scan = FALSE, nf = 3)
# pay attention to the standardization method
# maybe not all variables should be standardized the same way

nic.lo <- niche(dudi.lo, niche.sample, scann = FALSE)

omi.lo <- data.frame("Species"= rownames(niche.param(nic.lo)),
                   "OMI.Local" = niche.param(nic.lo)[,2]) 
tol.lo <- data.frame("Species"= rownames(niche.param(nic.lo)),
                   "Tol.Local" = niche.param(nic.lo)[,3])

# and these are use to estimate abund and occur
occup.sample <- data.frame("Species"= colnames(bio[sites.abund, goodcols.abund]),
                           "Occupancy" = colSums(decostand(bio[sites.abund, 
                                                               goodcols.abund], 
                                                           "pa")) / 
                             nrow(bio[sites.abund, goodcols.abund]))

occup.sample.t <- data.frame("Species"= colnames(bio[sites.abund, 
                                                     goodcols.abund]), 
                             "Occupancy" = colSums(decostand(bio[sites.abund, 
                                                                 goodcols.abund], 
                                                             "pa")))

abund.sample <- data.frame("Species"= colnames(bio[sites.abund, goodcols.abund]),
                           "Abundance" = colSums(bio[sites.abund, 
                                                     goodcols.abund]) / 
                             occup.sample.t[,2])

# I also have a data frame with species traits (e.g., body size)
traits.df <- data.frame("Species"= rownames(traits), traits)
 
# because we can loose some species during the resampling
# we need to make sure that our data frames always have
# the same species in their rows

final.data <- data.frame(merge(
  merge(
    merge(
        merge(omi.lo, tol.lo, by = 1, all.x=F), 
    abund.sample, by = 1, all.x=F),
  occup.sample, by = 1, all.x=F),
  traits.df, by = 1, all.x=F))

# regression models

names(final.data)

# occup model1

occp_mod <- betareg(Occupancy ~ log(Abundance) + OMI.Local + Tol.Local + 
                      log(comprimento_S_medio) + NT_media + 
                      potencial_natatorio, data = final.data)
                      
# occup model2

occp_mod2 <- betareg(Occupancy ~ log(Abundance) + OMI.Local + Tol.Local + 
                      log(peso_medio) + NT_media + 
                      potencial_natatorio, data = final.data)  

# organizing results 
# these are the info I want from the models
resu <- list(
  modComp = c(as.vector(summary(occp_mod)$coefficients$mean[-1,]),
                    summary(occp_mod)$pseudo.r.squared,
                    cor(occp_mod$fitted.values, final.data$Occupancy)),
  modPeso = c(as.vector(summary(occp_mod2)$coefficients$mean[-1,]),
               summary(occp_mod2)$pseudo.r.squared,
               cor(occp_mod2$fitted.values, final.data$Occupancy)))
} 

# end of the function ==========================================================
#===============================================================================

# run it

resam <- list(modComp = matrix(0, nrow=1000, ncol = 26),
              modPeso = matrix(0, nrow=1000, ncol = 26))
              
for(i in 1:2) {
  for (j in 1:1000) {
    resam[[i]][j,] <- omi.resamp (abund, env, traits)[[i]]
  }
}                        

# All the info I extracted from models as columms of the data frame
colnames(resam[[1]]) <- c("Abund_coef", "OMI_coef", "Tol_coef", "Size_coef", 
                     "NT_coef", "Nat_coef", "Abund_se", "OMI_se", "Tol_se", 
                     "Size_se", "NT_se", "Nat_se", "Abund_z", "OMI_z", "Tol_z", 
                     "Size_z", "NT_z", "Nat_z", "Abund_P", "OMI_P", "Tol_P", 
                     "Size_P", "NT_P", "Nat_P", "ocup_pseuR2", "cor_fit_obs")

colnames(resam[[2]]) <- c("Abund_coef", "OMI_coef", "Tol_coef", "Wei_coef", 
                          "NT_coef", "Nat_coef", "Abund_se", "OMI_se", "Tol_se", 
                          "Wei_se", "NT_se", "Nat_se", "Abund_z", "OMI_z", 
                          "Tol_z", "Wei_z", "NT_z", "Nat_z", "Abund_P", "OMI_P", 
                          "Tol_P", "Wei_P", "NT_P", "Nat_P", "ocup_pseuR2", 
                          "cor_fit_obs")


str(resam)










