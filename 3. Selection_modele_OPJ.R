library("dplyr")

glmm.TMB <- read.csv2("data/OPJ/TMB_TT/TMB/opj_month.TMB.csv", sep=";")
glmm.TMB[c("Estimate_scaled","StdErr_scaled", "Estimate", "SdErrVariable", "Chisq", "Df", "pval", "Overdisp", "OverdispBolker","ovd.DHARMa", "OverdispDHARMa", "AIC","BIC", "zero_infl", "zero_infl.pv", "KS","KS.pv",  "Outlier.pv", "NbPresence", "Total")] <- sapply(glmm.TMB[c("Estimate_scaled","StdErr_scaled", "Estimate", "SdErrVariable", "Chisq", "Df", "pval", "Overdisp", "OverdispBolker","ovd.DHARMa", "OverdispDHARMa", "AIC","BIC", "zero_infl", "zero_infl.pv", "KS","KS.pv",  "Outlier.pv", "NbPresence", "Total")], as.numeric)
glmm.TMB <- glmm.TMB %>% mutate_if(is.numeric,~round(.,3)) # Arrondir les valeurs des colonnes numeriques
glmm.TMB <- glmm.TMB[(glmm.TMB$Effect=="YEAR"),]
names(glmm.TMB)[names(glmm.TMB)=="Esp"] <- "GROUP"
glmm.TMB[is.na(glmm.TMB)] <- "NA"
glmm.TMB$AIC <- as.numeric(as.character(glmm.TMB$AIC))

glmm.TMB.nb <- read.csv2("data/OPJ/TMB_TT/TMB_nb/opj_month.TMB.NB.csv", sep=";")
glmm.TMB.nb[c("Estimate_scaled","StdErr_scaled", "Estimate", "SdErrVariable", "Chisq", "Df", "pval", "Overdisp", "OverdispBolker","ovd.DHARMa", "OverdispDHARMa", "AIC","BIC", "zero_infl", "zero_infl.pv", "KS","KS.pv",  "Outlier.pv", "NbPresence", "Total")] <- sapply(glmm.TMB.nb[c("Estimate_scaled","StdErr_scaled", "Estimate", "SdErrVariable", "Chisq", "Df", "pval", "Overdisp", "OverdispBolker","ovd.DHARMa", "OverdispDHARMa", "AIC","BIC", "zero_infl", "zero_infl.pv", "KS","KS.pv",  "Outlier.pv", "NbPresence", "Total")], as.numeric)
glmm.TMB.nb <- glmm.TMB.nb %>% mutate_if(is.numeric,~round(.,3)) # Arrondir les valeurs des colonnes numeriques
glmm.TMB.nb <- glmm.TMB.nb[(glmm.TMB.nb$Effect=="YEAR"),]
names(glmm.TMB.nb)[names(glmm.TMB.nb)=="Esp"] <- "GROUP"
glmm.TMB.nb[is.na(glmm.TMB.nb)] <- "NA"

glmm.TMB.nb <- subset(glmm.TMB.nb, glmm.TMB.nb$AIC!="NA")
glmm.TMB.nb$AIC <- as.numeric(as.character(glmm.TMB.nb$AIC))

columns=c("GROUP", "Modele", "Estimate_scaled", "StdErr_scaled", "Estimate", "SdErrVariable", "pval", "OverdispBolker.ratio", "OverdispBolker.pv", "DeltaAIC")
Selection_opj <- data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(Selection_opj) <- columns

for (e in unique(glmm.TMB.nb$GROUP)){
  if (glmm.TMB$AIC[glmm.TMB$GROUP==e] <= glmm.TMB.nb$AIC[glmm.TMB.nb$GROUP==e]) {Selection_opj <- rbind(Selection_opj, data.frame(GROUP=e, 
                                                                                                                                  Modele="Poisson", 
                                                                                                                                  Estimate_scaled=glmm.TMB$Estimate_scaled[glmm.TMB$GROUP==e], 
                                                                                                                                  StdErr_scaled=glmm.TMB$StdErr_scaled[glmm.TMB$GROUP==e], 
                                                                                                                                  Estimate=glmm.TMB$Estimate[glmm.TMB$GROUP==e], 
                                                                                                                                  SdErrVariable=glmm.TMB$SdErrVariable[glmm.TMB$GROUP==e],
                                                                                                                                  pval=glmm.TMB$pval[glmm.TMB$GROUP==e],
                                                                                                                                  # OverdispBolker.ratio=glmm.TMB$[glmm.TMB$GROUP==e],
                                                                                                                                  # OverdispBolker.pv=glmm.TMB$OverdispBolker.pv[glmm.TMB$GROUP==e],
                                                                                                                                  DeltaAIC=glmm.TMB$AIC[glmm.TMB$GROUP==e]-glmm.TMB.nb$AIC[glmm.TMB.nb$GROUP==e],
                                                                                                                                  NbPresence=glmm.TMB$NbPresence[glmm.TMB$GROUP==e],
                                                                                                                                  Total=glmm.TMB$Total[glmm.TMB$GROUP==e])
  )}
  else {Selection_opj <- rbind(Selection_opj, data.frame(GROUP=e, 
                                                         Modele="NB", 
                                                         Estimate_scaled=glmm.TMB.nb$Estimate_scaled[glmm.TMB.nb$GROUP==e], 
                                                         StdErr_scaled=glmm.TMB.nb$StdErr_scaled[glmm.TMB.nb$GROUP==e], 
                                                         Estimate=glmm.TMB.nb$Estimate[glmm.TMB.nb$GROUP==e], 
                                                         SdErrVariable=glmm.TMB.nb$SdErrVariable[glmm.TMB.nb$GROUP==e],
                                                         pval=glmm.TMB.nb$pval[glmm.TMB.nb$GROUP==e],
                                                         # OverdispBolker.ratio=glmm.TMB.nb$OverdispBolker.ratio[glmm.TMB.nb$GROUP==e],
                                                         # OverdispBolker.pv=glmm.TMB.nb$OverdispBolker.pv[glmm.TMB.nb$GROUP==e],
                                                         DeltaAIC=glmm.TMB$AIC[glmm.TMB$GROUP==e]-glmm.TMB.nb$AIC[glmm.TMB.nb$GROUP==e],
                                                         NbPresence=glmm.TMB.nb$NbPresence[glmm.TMB.nb$GROUP==e],
                                                         Total=glmm.TMB.nb$Total[glmm.TMB.nb$GROUP==e])
  )}
}

Selection_opj <- dplyr::distinct(Selection_opj, GROUP, .keep_all=TRUE)

missing <- glmm.TMB[glmm.TMB$GROUP %in% setdiff(glmm.TMB$GROUP, glmm.TMB.nb$GROUP),
                    c("GROUP", "Estimate_scaled", "StdErr_scaled", "Estimate", "SdErrVariable","pval", "NbPresence", "Total")]
missing$Modele <- "Poisson"
missing$DeltaAIC <- "NA"
missing <- missing[,c("GROUP", "Modele", "Estimate_scaled", "StdErr_scaled", "Estimate", "SdErrVariable", "pval", "DeltaAIC", "NbPresence", "Total")]

Selection_opj <- rbind(Selection_opj, missing)

write.table(Selection_opj, "data/OPJ/TMB_TT/Trend_opj.csv", sep=";", row.names=FALSE)
