##############################################
### CODE TO ESTIMATE OP TEMPORAL TRENDS ###
##############################################

# This code was launched on a cluster.

set.seed(1234)

# Loading of packages
library("DHARma")
library("stringr")
library("zoo")
library("lmtest")
library("lme4")
library("Matrix")
library("carData")
library("car")
library("gap")
library("glmmTMB")
library("data.table")
library("estimability")
library("effects")

# Importing the dataset
ukbms <- read.csv2("/1.Temporal_trends_and_interannual_variations/data/ukBMS_abundance_data_may_august_Eng_W_Sud.csv", sep = ";")


lu <- function(x){
  y <- length(unique(x))
  return(y) 
}

overdisp_fun <- function(model) { # creation de la fonction qui teste la surdispersion
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Selection des especes observees suffisamment de fois
nbYEAR <- aggregate(ukBMS$YEAR,by=list(ukBMS$GROUP),lu)
colnames(nbYEAR) <- c("GROUP","nbYEAR")
nbs <- aggregate(ukBMS$TRANSECT_ID,by=list(ukBMS$GROUP),lu) 
colnames(nbs) <- c("GROUP","nbs")
SPECIESobs <- merge(nbYEAR,nbs,by=c("GROUP"),all.x=TRUE, all.y=TRUE)

SPECIESpourtend2 <- subset(SPECIESobs, SPECIESobs$nbYEAR>4 & SPECIESobs$nbs>20)

ukBMS <- subset(ukBMS, ukBMS$GROUP %in% SPECIESpourtend2$GROUP) 

rm(nbYEAR, nbs, SPECIESobs, SPECIESpourtend2)

#e = unique(ukBMS$GROUP)[1] # Pour tester avec un seul groupe

#for (e in unique(ukBMS$GROUP)){
#   print(e)
#   flush.console()
#   espece_e <-  subset(ukBMS, GROUP == e)
#   
#   espece_e$YEAR_sc <- scale(espece_e$YEAR)
#   espece_e$LONGITUDE_sc <- scale(espece_e$LONGITUDE )
#   espece_e$LATITUDE_sc <- scale(espece_e$LATITUDE)
#   espece_e$MONTH_sc <- scale(espece_e$MONTH)
#   
#   mod.TMB.nb <- glmmTMB(COUNT ~ YEAR_sc + LATITUDE_sc + LONGITUDE_sc + MONTH_sc + (1 | TRANSECT_ID), data = espece_e, family="nbinom2")
# 
#   sTMBnb <- summary(mod.TMB.nb)
#   ovd <- sTMBnb$AICtab[4] / sTMBnb$AICtab[5]
#   overdisp_fun(mod.TMB.nb)[4]
#   b <- car::Anova(mod.TMB.nb)
#   
#   ovd.DHARMa <- testDispersion(simulateResiduals(mod.TMB.nb), alternative = "two.sided", plot = F, type = c("DHARMa"))$statistic
#   ovd.DHARMa.pv <- testDispersion(simulateResiduals(mod.TMB.nb), alternative = "two.sided", plot = F, type = c("DHARMa"))$p.value # verifie s'il y a presence de sur ou sous dispersion
#   ks <- testUniformity(simulateResiduals(mod.TMB.nb))$statistic
#   ks.pv <- testUniformity(simulateResiduals(mod.TMB.nb))$p.value # verifie si on s'eloigne de la distribution theorique
#   #outlier <- testOutliers(simulateResiduals(mod.TMB.nb), alternative = c("two.sided"))$statistic # verifie la presence d'exces de residus plus extremes que toutes les simulations
#   outlier.pv <- testOutliers(simulateResiduals(mod.TMB.nb), alternative = c("two.sided"))$p.value
#   zero_inflation <- testZeroInflation(simulateResiduals(mod.TMB.nb))$statistic
#   zero_inflation.pv <- testZeroInflation(simulateResiduals(mod.TMB.nb))$p.value
#   
#   sdYEAR = attributes(mod.TMB.nb$frame$YEAR_sc)$`scaled:scale`
#   estYEAR <- (sTMBnb$coefficients$cond["YEAR_sc", "Estimate"])/sdYEAR
#   
#   sdLATITUDE = attributes(mod.TMB.nb$frame$LATITUDE_sc)$`scaled:scale`
#   estLATITUDE <- (sTMBnb$coefficients$cond["LATITUDE_sc", "Estimate"])/sdLATITUDE
#   
#   sdLONGITUDE = attributes(mod.TMB.nb$frame$LONGITUDE)$`scaled:scale`
#   estLONGITUDE <- (sTMBnb$coefficients$cond["LONGITUDE_sc", "Estimate"])/sdLONGITUDE
#   
#   sdMONTH = attributes(mod.TMB.nb$frame$MONTH)$`scaled:scale`
#   estMONTH <- (sTMBnb$coefficients$cond["MONTH_sc", "Estimate"])/sdMONTH
#   
#   ddTMBnb <- data.frame(Esp = e, Effect = row.names(sTMBnb$coefficients$cond[,0]), Estimate_scaled = sTMBnb$coefficients$cond[,1], 
#                         StdErr_scaled = sTMBnb$coefficients$cond[,2], SdErrVariable = c(NA,sdYEAR, sdLATITUDE, sdLONGITUDE, sdMONTH),
#                         Estimate = c(NA, estYEAR, estLATITUDE, estLONGITUDE, estMONTH), Chisq = c(NA,b[,1]), Df = c(NA,b[,2]), pval = c(sTMBnb$coefficients$cond[1,4],b[,3]), 
#                         sigma = sigma(mod.TMB.nb),Overdisp = sTMBnb$AICtab[4] / sTMBnb$AICtab[5], OverdispBolker = overdisp_fun(mod.TMB.nb)[4],
#                         ovd.DHARMa = ovd.DHARMa,  OverdispDHARMa = ovd.DHARMa.pv, AIC = sTMBnb$AICtab[1],
#                         BIC = sTMBnb$AICtab[2], zero_infl = zero_inflation, zero_infl.pv = zero_inflation.pv, KS = ks, KS.pv = ks.pv ,Outlier.pv = outlier.pv)
#   
#   ddTMBnb$NbPresence <- nrow(subset(espece_e, COUNT != 0))
#   ddTMBnb$Total <- dim(espece_e)[1]
#   
#   # DHARMa residuals
#   p <- paste("Res.TMB.nb/",e,as.character("_DHARMa"),".jpeg",sep="")
#   jpeg(p)
#   DHARMa.res <- plot(simulateResiduals(mod.TMB.nb))
#   dev.off()
#   rm(DHARMa.res)
#   write.table(ddTMBnb, paste0("/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_TT/TMB_nb/resultatTMBnb_",e,".csv", row.names=FALSE), sep = ";")
#   
#   mod.TMB <- glmmTMB(COUNT ~ YEAR_sc + LATITUDE_sc + LONGITUDE_sc + MONTH_sc + (1 | TRANSECT_ID), data = espece_e, family="poisson")
#   sTMB <- summary(mod.TMB)
#   ovd <- sTMB$AICtab[4] / sTMB$AICtab[5] # rapport deviance / df.resid pour evaluer la surdispersion
#   overdisp_fun(mod.TMB)[4] # ou test de Bolker, si pval< 0.05 -> overdisp
#   a <- car::Anova(mod.TMB) # test="F" a utiliser dans l'Anova en cas de surdispersion pour les glm, mais pas applicable ici car glmer
#   
#   ovd.DHARMa <- testDispersion(simulateResiduals(mod.TMB), alternative = "two.sided", plot = F, type = c("DHARMa"))$statistic
#   ovd.DHARMa.pv <- testDispersion(simulateResiduals(mod.TMB), alternative = "two.sided", plot = F, type = c("DHARMa"))$p.value # verifie s'il y a presence de sur ou sous dispersion
#   ks <- testUniformity(simulateResiduals(mod.TMB))$statistic
#   ks.pv <- testUniformity(simulateResiduals(mod.TMB))$p.value # verifie si on s'eloigne de la distribution theorique
#   outlier.pv <- testOutliers(simulateResiduals(mod.TMB), alternative = c("two.sided"))$p.value # verifie la presence d'exces de residus plus extremes que toutes les simulations
#   zero_inflation <- testZeroInflation(simulateResiduals(mod.TMB))$statistic
#   zero_inflation.pv <- testZeroInflation(simulateResiduals(mod.TMB))$p.value
#   
#   sdYEAR = attributes(mod.TMB$frame$YEAR_sc)$`scaled:scale`
#   estYEAR <- (sTMB$coefficients$cond["YEAR_sc", "Estimate"])/sdYEAR
#   
#   sdLATITUDE = attributes(mod.TMB$frame$LATITUDE_sc)$`scaled:scale`
#   estLATITUDE <- (sTMB$coefficients$cond["LATITUDE_sc", "Estimate"])/sdLATITUDE
#   
#   sdLONGITUDE = attributes(mod.TMB$frame$LONGITUDE_sc)$`scaled:scale`
#   estLONGITUDE <- (sTMB$coefficients$cond["LONGITUDE_sc", "Estimate"])/sdLONGITUDE
#   
#   sdMONTH = attributes(mod.TMB$frame$MONTH_sc)$`scaled:scale`
#   estMONTH <- (sTMB$coefficients$cond["MONTH_sc", "Estimate"])/sdMONTH
#   
#   
#   ddTMB <- data.frame(Esp = e, Effect = row.names(sTMB$coefficients$cond[,0]), Estimate_scaled = sTMB$coefficients$cond[,1], StdErr_scaled = sTMB$coefficients$cond[,2], SdErrVariable = c(NA,sdYEAR, sdLATITUDE, sdLONGITUDE, sdMONTH),
#                       Estimate = c(NA, estYEAR, estLATITUDE, estLONGITUDE, estMONTH), Chisq = c(NA,b[,1]), Df = c(NA,b[,2]), pval = c(sTMB$coefficients$cond[1,4],b[,3]), sigma = sigma(mod.TMB),
#                       Overdisp = sTMBnb$AICtab[4] / sTMBnb$AICtab[5], OverdispBolker = overdisp_fun(mod.TMB.nb)[4],ovd.DHARMa = ovd.DHARMa,  OverdispDHARMa = ovd.DHARMa.pv, AIC = sTMBnb$AICtab[1],
#                       BIC = sTMBnb$AICtab[2], zero_infl = zero_inflation, zero_infl.pv = zero_inflation.pv, KS = ks, KS.pv = ks.pv ,Outlier.pv = outlier.pv )
#   
#   
#   ddTMB$NbPresence <- nrow(subset(espece_e, COUNT != 0))
#   ddTMB$Total <- dim(espece_e)[1]
#   
#   
#   write.table(ddTMB, paste0("/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_TT/TMB/resultatTMB_",e,".csv"), sep = ";", row.names=FALSE)
#   
#   
#   if (e == unique(ukBMS$GROUP)[1]) resultFreqTMBnb <- ddTMBnb else resultFreqTMBnb <- rbind(resultFreqTMBnb,ddTMBnb)
#   if (e == unique(ukBMS$GROUP)[1]) resultFreqTMB <- ddTMB else resultFreqTMB <- rbind(resultFreqTMB,ddTMB)
#   
#   rm(a, sTMB, sTMBnb, ddTMBnb, ddTMB, mod.TMB.nb, mod.TMB, ovd)
# 
#   write.table(resultFreqTMB, "/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_TT/TMB/ukBMS_month.TMB.csv", sep = ";", row.names = FALSE)
#   write.table(resultFreqTMBnb, "/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_TT/TMB_nb/ukBMS_month.TMB.NB.csv", sep = ";", row.names = FALSE)
#   
# }

#rm(b, espece_e, resultFreqTMB, resultFreqTMBnb)

# Pour les variations interannuelles
for (e in unique(ukBMS$GROUP)) {
print(e)
espece_e <- subset(ukBMS, ukBMS$GROUP==e)
YEAR2<-as.factor(espece_e$YEAR)

espece_e$LATITUDE_sc <- scale(espece_e$LATITUDE)
espece_e$LONGITUDE_sc <- scale(espece_e$LONGITUDE)
espece_e$MONTH_sc <- scale(espece_e$MONTH)

mod.TMB <- glmmTMB(COUNT ~ YEAR2 +  LATITUDE_sc + LONGITUDE_sc + MONTH_sc + (1 | TRANSECT_ID), data = espece_e, family="poisson")
sTMB <- summary(mod.TMB)

obj<-allEffects(mod.TMB,xlevels=list(YEAR2))

varcwm <- data.frame(obj$YEAR2)
names(varcwm)[1] <- "year"
varcwm$annees <- as.numeric(levels(varcwm$year))[varcwm$year]
varcwm$Espece <- rep(e, nrow(varcwm))
varcwm$Programme <- rep("ukBMS", nrow(varcwm))

# sdLATITUDE = attributes(mod.TMB$frame$LATITUDE_sc)$`scaled:scale`
# estLATITUDE <- (sTMB$coefficients$cond["LATITUDE_sc", "Estimate"])/sdLATITUDE
# 
# sdLONGITUDE = attributes(mod.TMB$frame$LONGITUDE_sc)$`scaled:scale`
# estLONGITUDE <- (sTMB$coefficients$cond["LONGITUDE_sc", "Estimate"])/sdLONGITUDE
# 
# sdMONTH = attributes(mod.TMB$frame$MONTH_sc)$`scaled:scale`
# estMONTH <- (sTMB$coefficients$cond["MONTH_sc", "Estimate"])/sdMONTH

#write.table(resultFreqGLMER, "/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_TT/ukBMS.GLMER.csv", sep = ";", row.names = FALSE)
#write.table(resultFreqGLMERnb, "/home/genouest/mnhn_cesco/sagnoux/ukBMS.GLMER.NB.csv", sep = ";", row.names = FALSE)
write.table(varcwm, paste0("/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_VI/resultFreqVarcwm_",e,".csv"), sep = ";", row.names = FALSE)
write.table(resultFreqVarcwm, "/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_VI/resultFreqVarcwm.csv", sep = ";", row.names = FALSE)

if (e == unique(ukBMS$GROUP)[1]) {resultFreqVarcwm <- varcwm} else {resultFreqVarcwm <- rbind(resultFreqVarcwm,varcwm) }

}
#rm(a, sglmer, sglmernb, ddglmernb, ddglmer, mod.glmer.nb, mod.glmer, ovd)
#rm(a, sglmer, ddglmer, mod.glmer, ovd)

#write.table(resultFreqGLMER, "/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_TT/ukBMS.TMB.csv", sep = ";", row.names = FALSE)
#write.table(resultFreqGLMERnb, "/home/genouest/mnhn_cesco/sagnoux/ukBMS.GLMER.NB.csv", sep = ";", row.names = FALSE)
write.table(varcwm, paste0("/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_VI/resultFreqVarcwm_TMB_",e,".csv"), sep = ";", row.names = FALSE)
write.table(resultFreqVarcwm, "/home/genouest/mnhn_cesco/sagnoux/ukBMS/TMB_VI/resultFreqVarcwm_TMB.csv", sep = ";", row.names = FALSE)

#rm(YEAR2, mod.glmer,sum, ddglmer, obj, a, blower, bupper, varcwm)

