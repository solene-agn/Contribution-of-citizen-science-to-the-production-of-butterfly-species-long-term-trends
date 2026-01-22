library("data.table")
library("speedglm") # fitting linear and generalized linear models to large data sets
library("devtools")
library("glmmTMB")
library("forcats")
library("DHARMa")
library("reshape2")
library("plyr")
library("sf")
library("terra")
library("tmap")
library("DT")
library("mapview")
library("rbms")
library("effects")
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("ggrepel")

set.seed(1234)

# Functions --------------------------------------------------------
lu <- function(x){
  y <- length(unique(x))
  return(y) 
}

overdisp_fun <- function(model) { # to test overdispersion
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Data importation --------------------------------------------------

obs <- data.table::fread("data/data/STERF_groups_corrige_mai_aout_nationale.csv", encoding="UTF-8")

names(obs)[names(obs)=="GROUP"] <- "SPECIES"

obs$LONGITUDE <- as.numeric(as.character(obs$LONGITUDE))
obs$LATITUDE <- as.numeric(as.character(obs$LATITUDE))
obs$DAY <- lubridate::day(obs$DATE)

obs <- subset(obs, obs$REGION_R %in% c("Normandie", "Île-de-France", "Hauts-de-France", "Bretagne", "Grand Est")) # data for north of France

# Step 1 : 
# We consider transects as sites to interpolate missing counts at the transect scale
obs$Site_transect_ID <- paste0(obs$SITE_ID, "_", obs$TRANSECT_ID)
obs <- obs[,c("Site_transect_ID", "DATE", "YEAR", "MONTH", "DAY", "SPECIES", "COUNT", "LONGITUDE", "LATITUDE")]
names(obs)[names(obs)=="Site_transect_ID"] <- "SITE_ID"

# We keep information on sites
info_site <- dplyr::distinct(obs[,c("SITE_ID", "LONGITUDE", "LATITUDE")])

# STEP 2:
# Creation of the time series corresponding to the data
ts_date <- ts_dwmy_table(InitYear = min(obs$YEAR),
                         LastYear = max(obs$YEAR) ,
                         WeekDay1 = 'monday') # Monday is indicated as the first day of the week (default is Sunday)

# Creation of data tables for phenology and temporal analysis
ts_season <- ts_monit_season(ts_date,
                             StartMonth = 5, # 4 for the STERF usually
                             EndMonth = 8, # 9
                             Anchor = TRUE,
                             AnchorLength = 3,
                             AnchorLag = 1,
                             TimeUnit = 'w')

MY_visit_region <- unique(data.frame(SITE_ID = obs$SITE_ID, DATE = obs$DATE))
ts_season_visit <- rbms::ts_monit_site(ts_season, MY_visit_region)

MY_count_region <- obs


# Species selection 
# "Aurores", "Machaons" # not enough data
species <- c("Amaryllis",  "Belle-dame", "Citrons", "Cuivres", "Demi-deuils", "Hesperides orangees",
             "Lycenes bleus", "Megeres", "Myrtil", "Paon du jour", "Petites tortues", "Pierides blanches",
             "Procris", "Robert-le-diable", "Soucis", "Sylvains", "Tabac d'Espagne", "Tircis", "Tristan", "Vulcain")

#s_sp=species[1] # to test the function with the first species

for (s_sp in species) {
  
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, MY_count_region, sp = s_sp)
  
  ts_flight_curve <- flight_curve(ts_season_count, NbrSample = NULL, MinVisit = 3, MinOccur = 2, MinNbrSite = 3,
                                  MaxTrial = 3, GamFamily = 'nb', SpeedGam = FALSE, SelectYear = NULL,
                                  TimeUnit = 'w')
  
  pheno <- ts_flight_curve$pheno

  # Impute & Site index ====
  impt_counts <- rbms::impute_count(ts_season_count = ts_season_count, 
                                    ts_flight_curve = pheno, 
                                    YearLimit = 2, 
                                    TimeUnit = "w")
  
  impt_counts_final <- impt_counts
  impt_counts_final$FINAL_COUNT <- impt_counts_final$COUNT
  
  impt_counts_final <- impt_counts_final[is.na(FINAL_COUNT) & !is.na(IMPUTED_COUNT), FINAL_COUNT := IMPUTED_COUNT]
  
  # We keep one value per month
  permutation <- order(impt_counts_final$FINAL_COUNT, decreasing=TRUE)
  impt_counts_final <- impt_counts_final[permutation,]
  impt_counts_final <- dplyr::distinct(impt_counts_final, SITE_ID, M_YEAR, MONTH, .keep_all=TRUE)
  impt_counts_final <- subset(impt_counts_final, MONTH >= 5 & MONTH <=8)
  
  impt_counts_final <- impt_counts_final[,c("SPECIES", "SITE_ID", "YEAR", "MONTH", "NM", "TOTAL_NM", "COUNT", "IMPUTED_COUNT", "FINAL_COUNT")]
  setnames(impt_counts_final, "SPECIES", "GROUP")
  espece.S <- impt_counts_final
  
  espece.S <- merge(espece.S, info_site, by="SITE_ID", all.x=TRUE, all.y=FALSE)
  espece.S[, c("SITE_ID", "TRANSECT_ID") := tstrsplit(SITE_ID, "_")]
  
  espece.S <- espece.S[complete.cases(espece.S$FINAL_COUNT)]
   nrow(espece.S[espece.S$FINAL_COUNT!=0,])
  
  # espece_e <- espece.S
  # espece_e <- espece_e[!is.na(espece_e$LATITUDE),]
  # 
  # espece_e$YEAR <- scale(espece_e$YEAR)
  # espece_e$LONGITUDE <- scale(espece_e$LONGITUDE)
  # espece_e$LATITUDE <- scale(espece_e$LATITUDE)
  # espece_e$MONTH <- scale(espece_e$MONTH)
  # 

  # To test temporal trend models
  
  # mod.TMB.nb <- glmmTMB(COUNT ~ YEAR*MONTH + LATITUDE+LONGITUDE   + (1 | SITE_ID/TRANSECT_ID), data = espece_e, family="nbinom2")
  # summary(mod.TMB.nb)
  # 
  # mod.TMB.nb <- glmmTMB(COUNT ~ YEAR*MONTH + LATITUDE+LONGITUDE+LATITUDE^2+LONGITUDE^2+LATITUDE:LONGITUDE  + (1 | SITE_ID/TRANSECT_ID), data = espece_e, family="nbinom2")
  # summary(mod.TMB.nb)
  # 
  # mod.TMB.nb <- glmmTMB(COUNT ~ YEAR*MONTH + LATITUDE+LONGITUDE   + (1 | SITE_ID/TRANSECT_ID), data = espece_e, family="poisson")
  # summary(mod.TMB.nb)
  # 
  # mod.TMB.nb <- glmmTMB(COUNT ~ YEAR*MONTH + LATITUDE+LONGITUDE+LATITUDE^2+LONGITUDE^2+LATITUDE:LONGITUDE  + (1 | SITE_ID/TRANSECT_ID), data = espece_e, family="poisson")
  # summary(mod.TMB.nb)
  
  # STEP 3:
  
  flush.console()
  
  # espece.S <-  subset(sterf, GROUP == e)
  
  
#INTERANNUAL VARIATIONS------------------------------------------
  YEAR2.S<-as.factor(espece.S$YEAR)
  espece.S$LONGITUDE <- scale(espece.S$LONGITUDE)
  espece.S$LATITUDE <- scale(espece.S$LATITUDE)
  espece.S$MONTH <- scale(espece.S$MONTH)
  
  #mod.TMB <- glmmTMB(COUNT ~ YEAR2.S + LATITUDE + LONGITUDE + MONTH + (1 | SITE_ID/TRANSECT_ID), data = espece.S, family="poisson")
  # mod.TMB <- glmmTMB(FINAL_COUNT ~ YEAR2.S + LATITUDE + LONGITUDE + (1 | SITE_ID/TRANSECT_ID), weights=TOTAL_NM,  data = espece.S, family="poisson")
  mod.TMB <- glmmTMB(FINAL_COUNT ~ YEAR2.S*MONTH + LATITUDE + LONGITUDE + (1 | SITE_ID/TRANSECT_ID), weights=TOTAL_NM,  data = espece.S, family="poisson")

  summary(mod.TMB)

  obj.S<-allEffects(mod.TMB,xlevels=list(YEAR2.S))

  var.abond.S <- data.frame(obj.S$YEAR2.S)
  names(var.abond.S)[1] <- "year"
  var.abond.S$annees <- as.numeric(levels(var.abond.S$year))[var.abond.S$year]
  var.abond.S$Espece <- rep(s_sp, nrow(var.abond.S))
  var.abond.S$Programme <- rep("STERF", nrow(var.abond.S))
  var.abond.S$AIC <- summary(mod.TMB)$AIC[1]

  # ggplot(var.abond.S, aes(annees,fit))+
  #   #geom_point()+
  #   geom_line(linewidth=1)+
  #   geom_rangeframe() +
  #   #theme(aspect.ratio = 1) +
  #   theme_classic() +
  #   labs(title = paste0("", s_sp), # Tendance temporelle d'abondance des populations d'imagos :
  #        subtitle = ,  x="Year", y = "Index of abundance") + #ajouter text_trend ? subtitle # Indice d'abondance
  #   geom_smooth(formula = y~ x, method=lm, data = var.abond.S ,se=TRUE, fullrange=TRUE, aes(y = fit, x = annees), alpha = 0.1, fill = '#3cd3dc') +
  #  theme(aspect.ratio = 1, plot.title = element_text(size = 10, face = "bold"), axis.text.x= element_text(angle = 30 , size = 2)) +
    # scale_x_continuous(breaks = seq(2006, 2019, 1), lim = c(2006, 2019))

  write.table(var.abond.S, paste("sterf/IV/sterf_result_glmmTMB_", s_sp, ".csv"), sep = ";", row.names = FALSE)
  
  if (s_sp == species[1]) resultFreqvar.abond.S <- var.abond.S else resultFreqvar.abond.S <- rbind(resultFreqvar.abond.S,var.abond.S)

  write.table(resultFreqvar.abond.S, "sterf/IV/sterf_resultFreqvarabond_glmmTMB.csv", sep = ";", row.names = FALSE)

  rm(espece.S, YEAR2.S, mod.glmer.S,sum, ddglmer.S, obj.S, a.S, blower.S, bupper.S, var.abond.S)
# }

}



#   # TEMPORAL TRENDS ------------------------------------------
#   
#   print(s_sp)
#   flush.console()
#   espece_e <-  espece.S
#   
#   espece_e$YEAR_sc <- scale(espece_e$YEAR)
#   espece_e$LONGITUDE_sc <- scale(espece_e$LONGITUDE )
#   espece_e$LATITUDE_sc <- scale(espece_e$LATITUDE)
#   espece_e$MONTH_sc <- scale(espece_e$MONTH)
#   
#   mod.TMB.nb <- glmmTMB(FINAL_COUNT ~ YEAR_sc*MONTH_sc + LATITUDE_sc + LONGITUDE_sc + (1 | SITE_ID/TRANSECT_ID), weights=TOTAL_NM, data = espece_e, family="nbinom2")

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
#   seYEAR <- (sTMBnb$coefficients$cond["YEAR_sc", "Std. Error"])/sdYEAR
#   
#   sdLATITUDE = attributes(mod.TMB.nb$frame$LATITUDE_sc)$`scaled:scale`
#   estLATITUDE <- (sTMBnb$coefficients$cond["LATITUDE_sc", "Estimate"])/sdLATITUDE
#   seLATITUDE<- (sTMBnb$coefficients$cond["LATITUDE_sc", "Std. Error"])/sdLATITUDE
#   
#   sdLONGITUDE = attributes(mod.TMB.nb$frame$LONGITUDE)$`scaled:scale`
#   estLONGITUDE <- (sTMBnb$coefficients$cond["LONGITUDE_sc", "Estimate"])/sdLONGITUDE
#   seLONGITUDE <- (sTMBnb$coefficients$cond["LONGITUDE_sc", "Std. Error"])/sdLONGITUDE
#   
#   sdMONTH = attributes(mod.TMB.nb$frame$MONTH)$`scaled:scale`
#   estMONTH <- (sTMBnb$coefficients$cond["MONTH_sc", "Estimate"])/sdMONTH
#   seMONTH  <- (sTMBnb$coefficients$cond["MONTH_sc", "Std. Error"])/sdLONGITUDE
#   
#   ddTMBnb <- data.frame(Esp = s_sp, Effect = row.names(sTMBnb$coefficients$cond[,0]), Estimate_scaled = sTMBnb$coefficients$cond[,1], 
#                         StdErr_scaled = sTMBnb$coefficients$cond[,2], SdErrVariable = c(NA,seYEAR, seLATITUDE, seLONGITUDE, seMONTH, NA),
#                         Estimate = c(NA, estYEAR, estLATITUDE, estLONGITUDE, estMONTH, NA), Chisq = c(NA,b[,1]), Df = c(NA,b[,2]), pval = c(sTMBnb$coefficients$cond[1,4],b[,3]), 
#                         sigma = sigma(mod.TMB.nb),Overdisp = sTMBnb$AICtab[4] / sTMBnb$AICtab[5], OverdispBolker = overdisp_fun(mod.TMB.nb)[4],
#                         ovd.DHARMa = ovd.DHARMa,  OverdispDHARMa = ovd.DHARMa.pv, AIC = sTMBnb$AICtab[1],
#                         BIC = sTMBnb$AICtab[2], zero_infl = zero_inflation, zero_infl.pv = zero_inflation.pv, KS = ks, KS.pv = ks.pv ,Outlier.pv = outlier.pv)
#   
#   ddTMBnb$NbPresence <- nrow(subset(espece_e, COUNT != 0))
#   ddTMBnb$Total <- dim(espece_e)[1]
#   
#   # # DHARMa residuals
#   # p <- paste0("E:/Papier-Papillons/Scripts_R_GenOuest/Method_2/STERF/TMB_TT/Res.TMB.nb/",s_sp,as.character("_DHARMa"),".jpeg",sep="")
#   # jpeg(p)
#   # DHARMa.res <- plot(simulateResiduals(mod.TMB.nb))
#   # dev.off()
#   # rm(DHARMa.res)
#   #write.table(ddTMBnb, paste0("E:/Papier-Papillons/Scripts_R_GenOuest/Method_2/STERF/TMB_TT/TMB_nb/resultatTMBnb_",s_sp,".csv"), row.names=FALSE, sep = ";")
#   write.table(ddTMBnb, paste0("sterf/CorresultatTMBnb_",s_sp,".csv"), row.names=FALSE, sep = ";")
#   
#   mod.TMB <- glmmTMB(FINAL_COUNT ~ YEAR_sc*MONTH_sc + LATITUDE_sc + LONGITUDE_sc + (1 | SITE_ID/TRANSECT_ID), weights=TOTAL_NM, data = espece_e, family="poisson")
#   #mod.TMB <- glmmTMB(FINAL_COUNT ~ YEAR_sc*MONTH_sc + LATITUDE_sc+LONGITUDE_sc+LATITUDE_sc^2+LONGITUDE_sc^2+LATITUDE_sc:LONGITUDE_sc  + (1 | SITE_ID/TRANSECT_ID), weights=TOTAL_NM, data = espece_e, family="poisson")
#   
#   sTMB <- summary(mod.TMB)
#   ovd <- sTMB$AICtab[4] / sTMB$AICtab[5] # rapport deviance / df.resid pour evaluer la surdispersion
#   overdisp_fun(mod.TMB)[4] # ou test de Bolker, si pval< 0.05 -> overdisp
#   a <- car::Anova(mod.TMB) # test="F" a utiliser dans l'Anova en cas de surdispersion pour les glm, mais pas applicable ici car glmer
#   
#   ovd.DHARMa <- testDispersion(simulateResiduals(mod.TMB), alternative = "two.sided", plot = F, type = c("DHARMa"))$statistic
#   ovd.DHARMa.pv <- testDispersion(simulateResiduals(mod.TMB), alternative = "two.sided", plot = F, type = c("DHARMa"))$p.value
#   ks <- testUniformity(simulateResiduals(mod.TMB))$statistic
#   ks.pv <- testUniformity(simulateResiduals(mod.TMB))$p.value
#   outlier.pv <- testOutliers(simulateResiduals(mod.TMB), alternative = c("two.sided"))$p.value 
#   zero_inflation <- testZeroInflation(simulateResiduals(mod.TMB))$statistic
#   zero_inflation.pv <- testZeroInflation(simulateResiduals(mod.TMB))$p.value
#   
#   sdYEAR = attributes(mod.TMB$frame$YEAR_sc)$`scaled:scale`
#   estYEAR <- (sTMB$coefficients$cond["YEAR_sc", "Estimate"])/sdYEAR
#   seYEAR <- (sTMB$coefficients$cond["YEAR_sc", "Std. Error"])/sdYEAR
#   
#   sdLATITUDE = attributes(mod.TMB$frame$LATITUDE_sc)$`scaled:scale`
#   estLATITUDE <- (sTMB$coefficients$cond["LATITUDE_sc", "Estimate"])/sdLATITUDE
#   seLATITUDE<- (sTMB$coefficients$cond["LATITUDE_sc", "Std. Error"])/sdLATITUDE
#   
#   sdLONGITUDE = attributes(mod.TMB$frame$LONGITUDE_sc)$`scaled:scale`
#   estLONGITUDE <- (sTMB$coefficients$cond["LONGITUDE_sc", "Estimate"])/sdLONGITUDE
#   seLONGITUDE <- (sTMB$coefficients$cond["LONGITUDE_sc", "Std. Error"])/sdLONGITUDE
#   
#   sdMONTH = attributes(mod.TMB$frame$MONTH_sc)$`scaled:scale`
#   estMONTH <- (sTMB$coefficients$cond["MONTH_sc", "Estimate"])/sdMONTH
#   seMONTH  <- (sTMB$coefficients$cond["MONTH_sc", "Std. Error"])/sdLONGITUDE
#   
#   
#   ddTMB <- data.frame(Esp = s_sp, Effect = row.names(sTMB$coefficients$cond[,0]), Estimate_scaled = sTMB$coefficients$cond[,1], StdErr_scaled = sTMB$coefficients$cond[,2], SdErrVariable = c(NA,seYEAR, seLATITUDE, seLONGITUDE, seMONTH, NA),
#                       Estimate = c(NA, estYEAR, estLATITUDE, estLONGITUDE, estMONTH, NA), Chisq = c(NA,b[,1]), Df = c(NA,b[,2]), pval = c(sTMB$coefficients$cond[1,4],b[,3]), sigma = sigma(mod.TMB),
#                       Overdisp = sTMB$AICtab[4] / sTMB$AICtab[5], OverdispBolker = overdisp_fun(mod.TMB)[4],ovd.DHARMa = ovd.DHARMa,  OverdispDHARMa = ovd.DHARMa.pv, AIC = sTMB$AICtab[1],
#                       BIC = sTMB$AICtab[2], zero_infl = zero_inflation, zero_infl.pv = zero_inflation.pv, KS = ks, KS.pv = ks.pv ,Outlier.pv = outlier.pv )
# 
#   ddTMB$NbPresence <- nrow(subset(espece_e, FINAL_COUNT != 0))
#   ddTMB$Total <- dim(espece_e)[1]
#   
#   
#   write.table(ddTMB, paste0("sterf/Corr_resultatTMB_",s_sp,".csv"), sep = ";", row.names=FALSE)
#   
#   
#   if (s_sp == species[1]) resultFreqTMBnb <- ddTMBnb else resultFreqTMBnb <- rbind(resultFreqTMBnb,ddTMBnb)
#   if (s_sp == species[1]) resultFreqTMB <- ddTMB else resultFreqTMB <- rbind(resultFreqTMB,ddTMB)
#   
#   rm(a, sTMB, sTMBnb, ddTMBnb, ddTMB, mod.TMB.nb, mod.TMB, ovd)
#   
#   # write.table(resultFreqTMB, "E:/Papier-Papillons/Scripts_R_GenOuest/Method_2/STERF/TMB_TT/TMB/sterf_month.TMB.csv", sep = ";", row.names = FALSE)
#   # write.table(resultFreqTMBnb, "E:/Papier-Papillons/Scripts_R_GenOuest/Method_2/STERF/TMB_TT/TMB_nb/sterf_month.TMB.NB.csv", sep = ";", row.names = FALSE)
#   # 
#   write.table(resultFreqTMB, "sterf/sterf_month.TMB.csv", sep = ";", row.names = FALSE)
#   write.table(resultFreqTMBnb, "sterf/sterf_month.TMB.NB.csv", sep = ";", row.names = FALSE)
#   
# }



#---------------------------------------------------------------------------------------------------
# To count amount of data (with and without imputation)

for (s_sp in species) {
  
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, MY_count_region, sp = s_sp)
  
  ts_flight_curve <- flight_curve(ts_season_count, NbrSample = NULL, MinVisit = 3, MinOccur = 2, MinNbrSite = 3,
                                  MaxTrial = 3, GamFamily = 'nb', SpeedGam = FALSE, SelectYear = NULL,
                                  TimeUnit = 'w')
  
  pheno <- ts_flight_curve$pheno
  
  # Impute & Site index ====
  impt_counts <- rbms::impute_count(ts_season_count = ts_season_count, 
                                    ts_flight_curve = pheno, 
                                    YearLimit = 2, 
                                    TimeUnit = "w")
  
  impt_counts_final <- impt_counts
  impt_counts_final$FINAL_COUNT <- impt_counts_final$COUNT
  
  impt_counts_final <- impt_counts_final[is.na(FINAL_COUNT) & !is.na(IMPUTED_COUNT), FINAL_COUNT := IMPUTED_COUNT]
  
  # We keep one value per month
  permutation <- order(impt_counts_final$FINAL_COUNT, decreasing=TRUE)
  impt_counts_final <- impt_counts_final[permutation,]
  impt_counts_final <- dplyr::distinct(impt_counts_final, SITE_ID, M_YEAR, MONTH, .keep_all=TRUE)
  impt_counts_final <- subset(impt_counts_final, MONTH >= 5 & MONTH <=8)
  
  impt_counts_final <- impt_counts_final[,c("SPECIES", "SITE_ID", "YEAR", "MONTH", "NM", "TOTAL_NM", "COUNT", "IMPUTED_COUNT", "FINAL_COUNT")]
  setnames(impt_counts_final, "SPECIES", "GROUP")
  espece.S <- impt_counts_final
  
  espece.S <- merge(espece.S, info_site, by="SITE_ID", all.x=TRUE, all.y=FALSE)
  espece.S[, c("SITE_ID", "TRANSECT_ID") := tstrsplit(SITE_ID, "_")]
  
  # With data imputation
  espece.S <- espece.S[complete.cases(espece.S$FINAL_COUNT)]
  total_avec_imp <- nrow(espece.S)
  obs_avec_imp <- nrow(espece.S[espece.S$FINAL_COUNT!=0,])
  
  # Without data imputation
  data_s <- subset(obs, obs$SPECIES==s_sp)
  total_sans_imp <- nrow(data_s)
  obs_sans_imp <- nrow(data_s[data_s$COUNT!=0,])
  
  info_s <- data.frame(Species=s_sp,
                       Total_without_imp=total_sans_imp,
                       Presences_without_imp=obs_sans_imp,
                       Total_with_imp=total_avec_imp,
                       Presences_with_imp=obs_avec_imp)
  
  if(s_sp == species[1]) {Info <- info_s} else {Info <- rbind(Info, info_s)}
  
}

write.table(Info, "data/Count_sterf.csv", sep=";", row.names=FALSE)

