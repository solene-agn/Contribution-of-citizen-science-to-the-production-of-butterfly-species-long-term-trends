# Libraries
library("lme4")
library("lmerTest") 
library("glmmTMB")
library("ggeffects")
library("GGally")
library("broom.helpers")
library("ggpubr")
library("forcats")

################################################################################

# To update species and species groups names
traits <- read.csv("data/Traits.csv", sep=";")
count_sterf <- read.csv("data/Count_sterf.csv", sep=";")

# Importing temporal trends ----

# OPJ
O.glmm_TMB <- read.csv("data/OPJ/TMB_TT/Trend_opj.csv", sep=";")
O.glmm_TMB <- O.glmm_TMB[,c("GROUP", "Estimate", "SdErrVariable", "pval", "NbPresence", "Total")]
O.glmm_TMB$Scheme <- "OPJ"
O.glmm_TMB <- merge(O.glmm_TMB, traits[,c("French_name", "Scientific_name")], by.x="GROUP", by.y="French_name", all.x=TRUE, all.y=FALSE)

for (e in unique(O.glmm_TMB$GROUP)) {
  est <- O.glmm_TMB$Estimate[O.glmm_TMB$GROUP == e]
  se  <- O.glmm_TMB$SdErrVariable[O.glmm_TMB$GROUP == e]
  
  if ((est - se) <= 0 & (est + se) >= 0) {
    O.glmm_TMB$pvalue[O.glmm_TMB$GROUP == e] <- "NS"
  } else {
    O.glmm_TMB$pvalue[O.glmm_TMB$GROUP == e] <- "S" 
  }
}

O.glmm_TMB$Trend <- with(O.glmm_TMB,
                         ifelse((Estimate - SdErrVariable) > 0, "Increase",
                                ifelse((Estimate + SdErrVariable) < 0, "Decrease", "Stable"))
)

# STERF
S.glmm_TMB <- read.csv("data/STERF/TMB_TT/Trend_sterf_cor.csv", sep=";")
S.glmm_TMB$Scheme <- "STERF"
S.glmm_TMB <- merge(S.glmm_TMB, traits[,c("French_name", "Scientific_name")], by.x="GROUP", by.y="French_name", all.x=TRUE, all.y=FALSE)

for (e in unique(S.glmm_TMB$GROUP)) {
  est <- S.glmm_TMB$Estimate[S.glmm_TMB$GROUP == e]
  se  <- S.glmm_TMB$SdErrVariable[S.glmm_TMB$GROUP == e]
  
  if ((est - se) <= 0 & (est + se) >= 0) {
    S.glmm_TMB$pvalue[S.glmm_TMB$GROUP == e] <- "NS"
  } else {
    S.glmm_TMB$pvalue[S.glmm_TMB$GROUP == e] <- "S" 
  }
}

S.glmm_TMB$Trend <- with(S.glmm_TMB,
                         ifelse((Estimate - SdErrVariable) > 0, "Increase",
                                ifelse((Estimate + SdErrVariable) < 0, "Decrease", "Stable"))
)

# UKBMS
U.glmm_TMB <- read.csv("data/UKBMS/TMB_TT/Trend_ukbms.csv", sep=";")
U.glmm_TMB <- U.glmm_TMB[,c("GROUP", "Estimate", "SdErrVariable", "pval", "NbPresence", "Total")]
U.glmm_TMB$Scheme <- "UKBMS"
U.glmm_TMB <- merge(U.glmm_TMB, traits[,c("French_name", "Scientific_name")], by.x="GROUP", by.y="French_name", all.x=TRUE, all.y=FALSE)

for (e in unique(U.glmm_TMB$GROUP)) {
  est <- U.glmm_TMB$Estimate[U.glmm_TMB$GROUP == e]
  se  <- U.glmm_TMB$SdErrVariable[U.glmm_TMB$GROUP == e]
  
  if ((est - se) <= 0 & (est + se) >= 0) {
    U.glmm_TMB$pvalue[U.glmm_TMB$GROUP == e] <- "NS"
  } else {
    U.glmm_TMB$pvalue[U.glmm_TMB$GROUP == e] <- "S" 
  }
}

U.glmm_TMB$Trend <- with(U.glmm_TMB,
                         ifelse((Estimate - SdErrVariable) > 0, "Increase",
                                ifelse((Estimate + SdErrVariable) < 0, "Decrease", "Stable"))
)

# Selection of shared species between the three schemes
Groupes <- intersect(S.glmm_TMB$GROUP, U.glmm_TMB$GROUP)
Groupes <- data.frame(GROUP=intersect(Groupes, O.glmm_TMB$GROUP))

S.glmm_TMB <- merge(S.glmm_TMB,Groupes, by="GROUP", all.y=TRUE )
U.glmm_TMB <- merge(U.glmm_TMB, Groupes, by="GROUP", all.y=TRUE)
O.glmm_TMB <- merge(O.glmm_TMB, Groupes, by="GROUP", all.y=TRUE)

S.glmm_TMB$Scientific_name[S.glmm_TMB$GROUP=="Sylvains"] <- "white admirals"
O.glmm_TMB$Scientific_name[O.glmm_TMB$GROUP=="Sylvains"] <- "white admirals"
U.glmm_TMB$Scientific_name[U.glmm_TMB$GROUP=="Sylvains"] <- "white admirals"

(100*nrow(O.glmm_TMB[O.glmm_TMB$pvalue=="NS",]))/nrow(O.glmm_TMB)
(100*nrow(U.glmm_TMB[U.glmm_TMB$pvalue=="NS",]))/nrow(U.glmm_TMB)
(100*nrow(S.glmm_TMB[S.glmm_TMB$pvalue=="NS",]))/nrow(S.glmm_TMB)

(100*nrow(O.glmm_TMB[O.glmm_TMB$Trend=="Increase",]))/nrow(O.glmm_TMB)
(100*nrow(U.glmm_TMB[U.glmm_TMB$Trend=="Increase",]))/nrow(U.glmm_TMB)
(100*nrow(S.glmm_TMB[S.glmm_TMB$Trend=="Increase",]))/nrow(S.glmm_TMB)

(100*nrow(O.glmm_TMB[O.glmm_TMB$Trend=="Decrease",]))/nrow(O.glmm_TMB)
(100*nrow(U.glmm_TMB[U.glmm_TMB$Trend=="Decrease",]))/nrow(U.glmm_TMB)
(100*nrow(S.glmm_TMB[S.glmm_TMB$Trend=="Decrease",]))/nrow(S.glmm_TMB)

S.glmm_TMB$Trend_STERF <- S.glmm_TMB$Trend
O.glmm_TMB$Trend_OP <- O.glmm_TMB$Trend
U.glmm_TMB$Trend_UKBMS <- U.glmm_TMB$Trend


Compa_trend <- merge(S.glmm_TMB[,c("GROUP", "Trend_STERF")], O.glmm_TMB[,c("GROUP", "Trend_OP")],
                     by.x="GROUP", by.y="GROUP", all.x=TRUE, all.y=TRUE)

Compa_trend <- merge(Compa_trend, U.glmm_TMB[,c("GROUP", "Trend_UKBMS")],
                     by.x="GROUP", by.y="GROUP", all.x=TRUE, all.y=TRUE)

Compa_trend <- merge(Compa_trend, traits[,c("French_name", "Scientific_name")], 
                     by.x="GROUP", by.y="French_name", all.x=TRUE, all.y=FALSE)

# Appendix S2c
Total <- rbind(O.glmm_TMB, S.glmm_TMB)
Total <- rbind(Total, U.glmm_TMB)
Total <- Total[,c("Scheme", "GROUP", "Scientific_name", "Estimate", "SdErrVariable", "Total", "NbPresence")]
write.table(Total, "1.Temporal_trends_and_interannual_variations/data/Appendix_S2c.csv", sep=";", row.names=FALSE)

# Comparison between schemes
trends <- rbind(O.glmm_TMB, U.glmm_TMB)
trends <- rbind(trends, S.glmm_TMB)

mod_trend <- lmer(Estimate ~ Scheme + (1 | GROUP),
                  data = trends,
                  REML = TRUE)
hist(residuals(mod_trend))
Anova(mod_trend)

pwpm(emmeans(mod_trend,  ~ Scheme),means = FALSE, flip = TRUE,reverse = TRUE)


# Plot
S.glmm_TMB$Scientific_name <- factor(S.glmm_TMB$Scientific_name,
                             levels = c("white Pieridae", "Gonepteryx rhamni", "Colias crocea", # Pieridae
                                        "Vanessa cardui", "Vanessa atalanta", "Pyronia tithonus", "Polygonia c-album", "Pararge aegeria", 
                                        "Melanargia galathea", "Maniola jurtina", "Lasiommata maera", "Coenonympha pamphilus", 
                                        "Argynnis paphia", "Aphantopus hyperantus", "Aglais urticae", "Aglais io", "white admirals", # Nymphalidae
                                        "small blue Lycenidae", "Lycaena phlaeas", # Lycaenidae
                                        "orange skippers" # Hesperiidae
                             ))


jpeg(paste("Engl_Selected_model_Comparaisons_Tendances_temporelles_correct_dalto_noms_corr", "jpeg", sep = "."), width = 10, height =12, units="cm", quality=75, res=300)

p <- ggplot(S.glmm_TMB, aes(x=Scientific_name, y=Estimate)) +
  coord_flip() + theme_bw() + 
  geom_point(data = S.glmm_TMB, aes(x = Scientific_name, y = Estimate, shape=pvalue, col = "#3cd3dc"), size=2, show.legend=FALSE) + #shape=factor(pvalue), size=NbTotal
  geom_point(data = U.glmm_TMB, aes(x = Scientific_name, y = Estimate, shape=pvalue, col = "#a92eb7"), size=2, show.legend=FALSE) + #, shape=factor(pvalue), size=NbTotal
  geom_point(data = O.glmm_TMB, aes(x = Scientific_name, y = Estimate, shape=pvalue, col = "#69c223"), size=2, show.legend=FALSE) + # , shape=factor(pvalue), size=NbTotal
  geom_linerange(data = S.glmm_TMB, aes(ymin = Estimate-SdErrVariable, ymax = Estimate+SdErrVariable,), col = "#3cd3dc") +
  geom_linerange(data = U.glmm_TMB, aes(ymin = Estimate-SdErrVariable, ymax = Estimate+SdErrVariable,), col = "#a92eb7") +
  geom_linerange(data = O.glmm_TMB, aes(ymin = Estimate-SdErrVariable, ymax = Estimate+SdErrVariable,), col = "#69c223") +
  labs(x = "Species/species group",y = "Temporal trends", color = "Legend") +
  scale_color_manual(values = c("#69c223" = "#69c223", "#3cd3dc" = "#3cd3dc", "#a92eb7" = "#a92eb7")) +
  # scale_y_discrete(
  #   labels = c("Pyronia tithonus" = expression(italic("Pyronia tithonus"))) +
  scale_shape_manual(values=c(1,19)) + 
  theme(axis.text.y = element_text(size=8, face="italic"))
p
dev.off()


# Information on the number of data
# Comptage_STERF <- read.csv("E:/Papier-Papillons/Donnees/STERF/Comptage_donnees_STERF.csv", sep = ";")
# S.glmm_TMB <- merge(S.glmm_TMB, Comptage_STERF, by="GROUP")
# 
# Comptage_ukbms <- read.csv("E:/Papier-Papillons/Donnees/ukBMS/Comptage_donnees_ukBMS.csv", sep = ";")
# U.glmm_TMB <- merge(U.glmm_TMB, Comptage_ukbms, by="GROUP")
# 
# Comptage_OPJ <- read.csv("E:/Papier-Papillons/Donnees/OPJ/Comptage_donnees_OPJ.csv", sep = ";")
# O.glmm_TMB$NbPresence <- NULL
# O.glmm_TMB <- merge(O.glmm_TMB, Comptage_OPJ, by="GROUP")
# 
# ens_O <- O.glmm_TMB[,c("GROUP", "Estimate","pval","Scheme","NbTotal","NbPresence")]
# names(ens_O) <- c("GROUP", "Estimate","p_value","Scheme","NbTotal","NbPresence")
# ens_S <- S.glmm_TMB[,c("GROUP", "Estimate","p_value","Scheme","NbTotal","NbPresence")]
# ens_U <- U.glmm_TMB[,c("GROUP", "Estimate","p_value","Scheme","NbTotal","NbPresence")]

################################################################################
#' Calculating the distance to the bisector
#' For each species or group of species, we represent the average estimate of the year effect obtained from a programme
#' as a function of that obtained by another programme
#' The closer the estimates of the two programmes are, the closer the point will be to the line x=y

# Function that calculates the shortest distance between the coordinate point (x0,y0) and the bisector
distance <- function(y0, x0){
  d<-round(abs(y0-x0)/sqrt(2),3)
}

# For each species or group of species
#g = unique(Groupes$Esp)[1]
for(g in unique(Groupes$GROUP)){
  d_OS_g <- distance(O.glmm_TMB[O.glmm_TMB$GROUP==g,"Estimate"],S.glmm_TMB[S.glmm_TMB$GROUP==g,"Estimate"])
  d_OU_g <- distance(O.glmm_TMB[O.glmm_TMB$GROUP==g,"Estimate"],U.glmm_TMB[U.glmm_TMB$GROUP==g,"Estimate"])
  d_US_g <- distance(U.glmm_TMB[U.glmm_TMB$GROUP==g,"Estimate"],S.glmm_TMB[S.glmm_TMB$GROUP==g,"Estimate"])
  distances_g <- data.frame(Groupe=g, Distance=c(d_OS_g, d_OU=d_OU_g, d_US=d_US_g),
                            Schemes=c("OPJ-STERF", "OPJ-UKBMS", "UKBMS-STERF"), 
                            NbTotal=c(sum(O.glmm_TMB$NbTotal[O.glmm_TMB$GROUP==g]+S.glmm_TMB$NbTotal[S.glmm_TMB$GROUP==g]),
                                      sum(O.glmm_TMB$NbTotal[O.glmm_TMB$GROUP==g]+U.glmm_TMB$NbTotal[U.glmm_TMB$GROUP==g]),
                                      sum(U.glmm_TMB$NbTotal[U.glmm_TMB$GROUP==g]+S.glmm_TMB$NbTotal[S.glmm_TMB$GROUP==g])),
                            NbPresence=c(sum(O.glmm_TMB$NbPresence[O.glmm_TMB$GROUP==g]+S.glmm_TMB$NbPresence[S.glmm_TMB$GROUP==g]),
                                         sum(O.glmm_TMB$NbPresence[O.glmm_TMB$GROUP==g]+U.glmm_TMB$NbPresence[U.glmm_TMB$GROUP==g]),
                                         sum(U.glmm_TMB$NbPresence[U.glmm_TMB$GROUP==g]+S.glmm_TMB$NbPresence[S.glmm_TMB$GROUP==g])),
                            NbZero=c(sum(O.glmm_TMB$NbZero[O.glmm_TMB$GROUP==g]+S.glmm_TMB$NbZero[S.glmm_TMB$GROUP==g]),
                                     sum(O.glmm_TMB$NbZero[O.glmm_TMB$GROUP==g]+U.glmm_TMB$NbZero[U.glmm_TMB$GROUP==g]),
                                     sum(U.glmm_TMB$NbZero[U.glmm_TMB$GROUP==g]+S.glmm_TMB$NbZero[S.glmm_TMB$GROUP==g]))
  )
  if(g==unique(Groupes$GROUP)[1]) {Distances <- distances_g} else {Distances <- rbind(Distances,distances_g)}
  rownames(Distances) <- NULL
  
}

rm(Comptage_OPJ, Comptage_STERF, Comptage_ukbms,distances_g, O.glmm_TMB, S.glmm_TMB, U.glmm_TMB, d_OS_g, d_OU_g, d_US_g, g, distance, bd)

# Adding trait information
traits <- read.csv("data/Traits.csv", sep=";")
Distances <- merge(Distances, traits, by.x="Groupe", by.y="French_name")
Distances[is.na(Distances)] <- "NA"

Distances <- Distances[,c("Groupe", "Scientific_name", "Schemes", "Distance")]

write.table(Distances, "data/Distances_bissectrices_traits_final_corr.csv", sep=";", row.names=FALSE)
####
