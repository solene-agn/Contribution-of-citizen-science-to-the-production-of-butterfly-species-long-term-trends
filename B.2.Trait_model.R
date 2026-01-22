library("spaMM")
library("car")
library("ape")
library("DHARMa")
library("glmmTMB")
library("emmeans")
library("PCAmixdata")
library("ggplot2")
library("ggthemes")
library("ggrepel")
library("ggpubr")
library("gridExtra")
library("ggforce")
library("dplyr")
library("tidyr")

# Data importation

traits <- read.csv("data/Traits.csv", sep=";")
traits$Wingspan[traits$Wingspan=="Medium"] <- "Large"

traits$Scientific_name[traits$Scientific_name=="admirals"] <- "white admirals"

count <- read.csv("Appendix_S2c_corr.csv", sep=";")
count$Scientific_name[count$Scientific_name=="admirals"] <- "white admirals"


count <- count %>%
  select(Scheme, Scientific_name, Total)

count <- count %>%
  pivot_wider(names_from = Scheme, values_from = Total)

count <- count %>%
  mutate(
    OP_STERF = OPJ + STERF,
    STERF_UKBMS = STERF + UKBMS,
    UKBMS_OP = UKBMS + OPJ
  ) %>%
  select(Scientific_name, OP_STERF, STERF_UKBMS, UKBMS_OP)

count <- count %>%
  pivot_longer(cols = -Scientific_name,
               names_to = "Schemes",
               values_to = "Total") %>%
  mutate(Schemes = case_when(
    Schemes == "OP_STERF" ~ "OPJ-STERF",
    Schemes == "STERF_UKBMS" ~ "UKBMS-STERF",
    Schemes == "UKBMS_OP" ~ "OPJ-UKBMS",
    TRUE ~ Schemes
  ))

count$Schemes[count$Schemes=="OPJ-STERF"] <- "OP-STERF"
count$Schemes[count$Schemes=="UKBMS-STERF"] <- "STERF-UKBMS"
count$Schemes[count$Schemes=="OPJ-UKBMS"] <- "OP-UKBMS"

# Create a contingency table
tableau_contingence <- table(traits$Identification, traits$Wingspan)
tableau_contingence <- table(traits$Identification, traits$Habitat)
tableau_contingence <- table(traits$Identification, traits$Migratory)
tableau_contingence <- table(traits$Identification, traits$Group_sp)

# Fisher's test
fisher.test(tableau_contingence,simulate.p.value = TRUE)

# Test MCA -----------------------------------------------------------------------------------------------------------

# We keep species for which the information is complete
select_df <- traits[complete.cases(traits),]

select_df$Family <- as.factor(as.character(select_df$Family))
select_df$Identification <- as.factor(as.character(select_df$Identification))
select_df$Wingspan <- as.factor(as.character(select_df$Wingspan))
select_df$Group_sp <- as.factor(as.character(select_df$Group_sp))
select_df$Habitat <- as.factor(as.character(select_df$Habitat))
select_df$Migratory <- as.factor(as.character(select_df$Migratory))

select_df_stock <- select_df[,c("French_name", "English_name", "Scientific_name", "Family")]

select_df$French_name <- NULL
select_df$English_name <- NULL
select_df$Scientific_name <- NULL
select_df$Nb_generations_per_year <- NULL
select_df$Family <- NULL
select_df$Mobility <- NULL
select_df$Group_sp <- NULL
#select_df$Wingspan <- NULL

# Separation of quantitative and qualitative variables, according to the traits selected
X.quanti <- splitmix(select_df)$X.quanti
X.quali <- splitmix(select_df)$X.quali

# PCAmix analysis
pca<-PCAmix(X.quanti=X.quanti,X.quali=X.quali,
            ndim=5,graph=FALSE, rename.level = TRUE)

# Retrieve variable coordinates to plot variables and individuals
# assembly of variable values in a single table
cor_p <- data.frame(rbind(pca$quanti$coord ,pca$levels$coord))
cor_p <- tibble::rownames_to_column(cor_p,"Variable")

# table of coordinates of individuals, their temporal trend and their floral traits
cor_ind <- data.frame(pca$ind$coord)
cor_ind <- cbind(cor_ind,Species_NbTotal)
cor_ind <- cbind(cor_ind, select_df) 
cor_ind <- cbind(cor_ind, select_df_stock)

# Function for displaying qualitative and quantitative variables on the same plot
Corr_circle <- function(coord) {
  ggplot(coord, aes(x=dim.1,y=dim.2, label=Variable)) +
    geom_point() +
    geom_circle(aes(x0=0,y0=0,r=1)) + 
    geom_segment(
      x = 0, y = 0,
      xend = coord$dim.1, yend = coord$dim.2,
      lineend = "round", 
      linejoin = "round",
      size = 0.5, 
      arrow = arrow(length = unit(0.1, "inches")),
      colour = "black"
    ) +
    geom_text(hjust=1, vjust=-1, size=3) +
    geom_hline(yintercept = 0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    labs(x=paste0("Dim 1 (",round(pca$eig["dim 1","Proportion"],2),"%)"),
         y=paste0("Dim 2 (",round(pca$eig["dim 2","Proportion"],2),"%)"),
         title="Correlation circle")
}

# Plot of variables
Corr_circle(cor_p)

# Plot of variables and individuals
coord <- cor_p
rownames(coord) <- coord$Variable

Plot_MCA <- ggplot(coord, aes(x=dim.1,y=dim.2)) +
  geom_point(data = cor_ind, 
             aes(x = dim.1, y = dim.2, fill=Family), size=5, shape = 21, stroke = 0) + 
  
  scale_fill_manual(values=c("Nymphalidae"="#DC8E36",
                             "Pieridae"="#CC6677",
                             "Lycaenidae"="#88CCEE",
                             "Hesperiidae"="#AA4499",
                             "Papilionidae"="#44AA99"),
                    labels = c(
                      expression(italic(Hesperiidae)),
                      expression(italic(Lycaenidae)),
                      expression(italic(Nymphalidae)),
                      expression(italic(Papilionidae)),
                      expression(italic(Pieridae)))) +
  geom_segment(
    x = 0, y = 0,
    xend = coord$dim.1, yend = coord$dim.2,
    
    lineend = "round",
    linejoin = "round",
    size = 0.5, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "black" 
  ) +
  # geom_text(aes(label=Variable),hjust=1, vjust=-1, size=4) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste0("Dim 1 (",round(pca$eig["dim 1","Proportion"],2),"%)"),
       y=paste0("Dim 2 (",round(pca$eig["dim 2","Proportion"],2),"%)"),
       title="") + # too big if we square the standard error
  theme(panel.background=element_rect(fill="white")) +
  geom_text_repel(data = coord, aes(x = dim.1, y = dim.2, label = rownames(coord)))

Plot_MCA

# Addition of information on pollinator dependency
axes <- data.frame( pca[["ind"]][["coord"]])
select_df$PC1 <-axes$dim.1
select_df$PC2 <-axes$dim.2

select_df$French_name <- select_df_stock$French_name
select_df$Scientific_name <- select_df_stock$Scientific_name


#-----------------------------------------------------------------------------------------------------------------------------
# Correlation coefficients (interannual variations)---------------------------------------------------------------------------
similarity_index <- read.csv("data/Tableau_correlation_groupes_traits_final_pivot_cor.csv", sep=";")

similarity_index <- merge(similarity_index, select_df, by.x="Groupe", by.y="Scientific_name", all.x=TRUE, all.y=FALSE)
similarity_index <- merge(similarity_index, traits[,c("Scientific_name", "Family", "Group_sp")], by.x="Groupe", by.y="Scientific_name",
                          all.x=TRUE, all.y=FALSE)
similarity_index$Group_sp <- as.factor(as.character(similarity_index$Group_sp))
similarity_index <- merge(similarity_index, count, by.x=c("Groupe", "Couple"), by.y=c("Scientific_name", "Schemes"),
                          all.x=TRUE, all.y=FALSE)

similarity_index$STERF <- ifelse(similarity_index$Couple %in% c("OP-STERF","UKBMS-STERF"), 1,0)
similarity_index$OPJ <- ifelse(similarity_index$Couple %in% c("OP-STERF","OP-UKBMS"), 1,0)
similarity_index$UKBMS <- ifelse(similarity_index$Couple %in% c("OP-UKBMS","UKBMS-STERF"), 1,0)

similarity_index$Identification <- as.factor(similarity_index$Identification)
similarity_index$Wingspan <- as.factor(similarity_index$Wingspan)
similarity_index$Couple <- as.factor(similarity_index$Couple)

# 1. Ecological model
  model_IV_eco <- glmmTMB(scale(cor)~ Couple * (PC1 + PC2) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  #model <- glmmTMB(scale(cor)~ Couple * (PC1 + PC2) + Family + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  #Anova(model)

  Anova(model_IV_eco)
  summary(model_IV_eco)
  DHARMa.res <- plot(simulateResiduals(model_IV_eco))
  
  # LRT test
  model_IV_eco_b <- glmmTMB(scale(cor)~  PC1 + PC2 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  model_IV_eco_b <- glmmTMB(scale(cor)~ Couple * PC2 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  model_IV_eco_b <- glmmTMB(scale(cor)~ Couple * PC1 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  model_IV_eco_b <- glmmTMB(scale(cor)~  Couple + PC1+ PC2 + Couple:PC2 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  model_IV_eco_b <- glmmTMB(scale(cor)~  Couple + PC2+ PC1 + Couple:PC1 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
  
  anova(model_IV_eco,model_IV_eco_b,test="Chisq")
  
  # Pair of schemes
  pwpm(emmeans(model_IV_eco,  ~ Couple),means = FALSE, flip = TRUE,reverse = TRUE)
  pwpm(emmeans(model_IV_eco,  ~ Couple:PC1),means = FALSE, flip = TRUE,reverse = TRUE)
  
  emtrends(model_IV_eco, ~ Couple:PC1, var = "PC1")
  
# 2. Model based on scheme characteristics
model_IV_sch <- glmmTMB(scale(cor)~ Couple * (Group_sp + scale(Total)) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
#model_IV_sch <- fitme(scale(cor)~Couple *( Group_sp  +  Total) + (1|Groupe) + (1|OPJ) + (1|STERF) + (1|UKBMS), method="ML",data=similarity_index, family=gaussian)

Anova(model_IV_sch, test="Chisq")
summary(model_IV_sch)
DHARMa.res <- plot(simulateResiduals(model_IV_sch))

# LRT test
model_IV_sch_b <- glmmTMB(scale(cor)~  Group_sp + scale(Total) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
model_IV_sch_b <- glmmTMB(scale(cor)~ Couple * scale(Total) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
model_IV_sch_b <- glmmTMB(scale(cor)~ Couple * Group_sp + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
model_IV_sch_b <- glmmTMB(scale(cor)~  Couple + Group_sp + scale(Total)+ Couple:Total + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)
model_IV_sch_b <- glmmTMB(scale(cor)~  Couple + scale(Total)+ Group_sp + Couple:Group_sp + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)


anova(model_IV_sch,model_IV_sch_b,test="Chisq")

# Pair of schemes
pwpm(emmeans(model_IV_sch, ~ Couple),means = FALSE, flip = TRUE,reverse = TRUE)
pwpm(emmeans(model_IV_sch, ~ Group_sp),means = FALSE, flip = TRUE,reverse = TRUE)
pwpm(emmeans(model_IV_sch, ~ Couple:Group_sp),means = FALSE, flip = TRUE,reverse = TRUE)

similarity_index$Total_sc <- scale(similarity_index$Total)
model_IV_sch <- glmmTMB(scale(cor)~ Couple * (Group_sp + Total_sc) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=similarity_index)

emtrends(model_IV_sch, ~ Couple:Total_sc, var = "Total_sc")

# Temporal trends ---------------------------------------------------------------------------
dissimilarity_index <- read.csv("data/Distances_bissectrices_traits_final.csv", sep=";")
dissimilarity_index <- read.csv("data/Distances_bissectrices_traits_final_corr.csv", sep=";")

dissimilarity_index$Scientific_name[dissimilarity_index$Scientific_name=="admirals"] <- "white admirals"

dissimilarity_index <- merge(dissimilarity_index, select_df[,c("PC1", "PC2", "Scientific_name")], by.x="Scientific_name", by.y="Scientific_name", all.x=TRUE, all.y=FALSE)
dissimilarity_index <- merge(dissimilarity_index, traits[,c("Scientific_name", "Family", "Group_sp")], by.x="Scientific_name", by.y="Scientific_name",
                          all.x=TRUE, all.y=FALSE)

dissimilarity_index$Schemes[dissimilarity_index$Schemes=="OPJ-STERF"] <- "OP-STERF"
dissimilarity_index$Schemes[dissimilarity_index$Schemes=="OPJ-UKBMS"] <- "OP-UKBMS"
dissimilarity_index$Schemes[dissimilarity_index$Schemes=="UKBMS-STERF"] <- "STERF-UKBMS"

dissimilarity_index <- merge(dissimilarity_index, count, by.x=c("Scientific_name", "Schemes"), by.y=c("Scientific_name", "Schemes"),
                          all.x=TRUE, all.y=FALSE)

dissimilarity_index$STERF <- ifelse(dissimilarity_index$Schemes %in% c("OP-STERF","UKBMS-STERF"), 1,0)
dissimilarity_index$OPJ <- ifelse(dissimilarity_index$Schemes %in% c("OP-STERF","OPJ-UKBMS"), 1,0)
dissimilarity_index$UKBMS <- ifelse(dissimilarity_index$Schemes %in% c("OP-UKBMS","UKBMS-STERF"), 1,0)

# model
model_TT_eco <- glmmTMB(scale(Distance)~ Schemes * (PC1 + PC2) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
#model_TT_eco <- glmmTMB(scale(Distance)~ Schemes * (PC1 + PC2) + Family + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)

Anova(model_TT_eco)
summary(model_TT_eco)
DHARMa.res <- plot(simulateResiduals(model_TT_eco))

# LRT test
model_TT_eco_b <- glmmTMB(scale(Distance)~  PC1 + PC2 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_eco_b <- glmmTMB(scale(Distance)~ Schemes * PC2 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_eco_b <- glmmTMB(scale(Distance)~ Schemes * PC1 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_eco_b <- glmmTMB(scale(Distance)~  Schemes + PC1+ PC2 + Schemes:PC2 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_eco_b <- glmmTMB(scale(Distance)~  Schemes + PC2+ PC1 + Schemes:PC1 + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)

anova(model_TT_eco,model_TT_eco_b,test="Chisq")

# Pair of schemes
pwpm(emmeans(model_TT_eco, ~ Schemes),means = FALSE, flip = TRUE,reverse = TRUE)
emtrends(model_TT_eco, ~ Schemes:PC1, var = "PC1")

# 2. Model based on scheme characteristics
model_TT_sch <- glmmTMB(scale(Distance)~ Schemes * (Group_sp + scale(Total)) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
#model_IV_sch <- fitme(scale(cor)~Couple *( Group_sp  +  Total) + (1|Groupe) + (1|OPJ) + (1|STERF) + (1|UKBMS), method="ML",data=similarity_index, family=gaussian)

Anova(model_TT_sch, test="Chisq")
summary(model_IV_sch)
DHARMa.res <- plot(simulateResiduals(model_IV_sch))

# LRT test
model_TT_sch_b <- glmmTMB(scale(Distance)~  Group_sp + scale(Total) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_sch_b <- glmmTMB(scale(Distance)~ Schemes * scale(Total) + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_sch_b <- glmmTMB(scale(Distance)~ Schemes * Group_sp + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_sch_b <- glmmTMB(scale(Distance)~  Schemes + Group_sp + scale(Total)+ Schemes:Total + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)
model_TT_sch_b <- glmmTMB(scale(Distance)~  Schemes + scale(Total)+ Group_sp + Schemes:Group_sp + (1|Groupe) + (1|STERF) + (1|OPJ) + (1|UKBMS), data=dissimilarity_index)

anova(model_TT_sch,model_TT_sch_b,test="Chisq")

# Pair of schemes
pwpm(emmeans(model_TT_sch, ~ Schemes),means = FALSE, flip = TRUE,reverse = TRUE)
pwpm(emmeans(model_TT_sch, ~ Group_sp),means = FALSE, flip = TRUE,reverse = TRUE)


# Plot ----

library("glmmTMB")
library("ggplot2")
library("DHARMa")
library("ggpubr")
library("ggpmisc")
library("rstatix")
library("reshape2")
library("ggstatsplot")
library("emmeans")
library("patchwork")

# Figure 4 : Emmeans, glmm --------------------------------------------------

# Pair of schemes

#IV
df_v <- similarity_index[,c("Groupe","cor","Couple")]
df_v <- df_v[complete.cases(df_v),]
v="Couple"

p.val.test<-pwpm(emmeans(model_IV_eco,  ~ Couple),means = FALSE, flip = TRUE,reverse = TRUE)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
stat.test<-subset(melt(p.matx),!is.na(value))
names(stat.test)<-c("group1","group2","p.adj")
stat.test[stat.test$p.adj<=0.001,"p.adj.signif"]<-"***"
stat.test[stat.test$p.adj>0.001 & stat.test$p.adj<=0.01,"p.adj.signif"]<-"**"
stat.test[stat.test$p.adj>0.01 & stat.test$p.adj<=0.05,"p.adj.signif"]<-"*"
stat.test[ stat.test$p.adj>0.05,"p.adj.signif"]<-"ns"
#stat.test<-mc_tribble(stat.test) 

plot_VI_Pair_of_schemes <- ggplot(df_v, aes(x= Couple, y=cor, color=Couple)) + # fill=Variable avec interactions
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  #theme(axis.text=element_text(size=12))+
  #adjust y label positions
  # theme(axis.title.y = element_text(margin = margin (r = 10)),
  #       #change the plot margins
  #       plot.margin = margin(l = 15,r=10)) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) + # ,face="bold"
  labs(x="Pair of schemes", y="Interannual variation similarity index") + # ,title=paste("Correlation vs.",v)
  scale_color_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  stat_pvalue_manual(stat.test,
                     y.position=1.2,step.increase = 0.1,
                     tip.length = 0.01, bracket.shorten = 0.05,
                     #y.position = max(df_v$Distance)+sd(df_v$Distance),
                     color = "black",
                     size=3)
plot_VI_Pair_of_schemes

# Temporal trends
df_v <- dissimilarity_index[,c("Groupe","Distance","Schemes")]
df_v <- df_v[complete.cases(df_v),]
v="Schemes"

p.val.test<-pwpm(emmeans(model_TT_eco,  ~ Schemes),means = FALSE, flip = TRUE,reverse = TRUE)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
stat.test<-subset(melt(p.matx),!is.na(value))
names(stat.test)<-c("group1","group2","p.adj")
stat.test[stat.test$p.adj<=0.001,"p.adj.signif"]<-"***"
stat.test[stat.test$p.adj>0.001 & stat.test$p.adj<=0.01,"p.adj.signif"]<-"**"
stat.test[stat.test$p.adj>0.01 & stat.test$p.adj<=0.05,"p.adj.signif"]<-"*"
stat.test[ stat.test$p.adj>0.05,"p.adj.signif"]<-"ns"
#stat.test<-mc_tribble(stat.test) 
# stat.test <- tribble(
#   ~group1, ~group2, ~p.adj, ~p.adj.signif,
#   "Non conspicuous", "Conspicuous", 0.0022, "**")

plot_TT_Pair_of_schemes <- ggplot(df_v, aes(x= Schemes, y=Distance, color=Schemes)) + # fill=Variable avec interactions
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  #theme(axis.text=element_text(size=12))+
  #adjust y label positions
  # theme(axis.title.y = element_text(margin = margin (r = 10)),
  #       #change the plot margins
  #       plot.margin = margin(l = 15,r=10)) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) + # ,face="bold"
  labs(x="Pair of schemes", y="Temporal trend dissimilarity index") + # ,title=paste("Correlation vs.",v)
  scale_color_manual(name="Schemes",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  stat_pvalue_manual(stat.test,
                     y.position=0.2,step.increase = 0.1,
                     tip.length = 0.01, bracket.shorten = 0.05,
                     #y.position = max(df_v$Distance)+sd(df_v$Distance),
                     color = "black",
                     size=3)
plot_TT_Pair_of_schemes


# PC1 ----

# IV

plot_IV_PC1 <- ggplot(similarity_index, aes(x = PC1, y = cor)) +
  geom_point(size = 2, aes(color=Couple)) +
  geom_smooth(method = "lm", se = TRUE, fill="#dc3c84", color="#dc3c84", alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "PC1",
    y = "Similarity index"
  ) +
  scale_color_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  theme_classic()
plot_IV_PC1

# Temporal trends

plot_TT_PC1 <- ggplot(dissimilarity_index, aes(x = PC1, y = Distance)) +
  geom_point(size = 2, aes(color=Schemes)) +
  geom_smooth(method = "lm", se = TRUE, fill="#dc3c84", color="#dc3c84", linetype = "dotdash", alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "PC1",
    y = "Dissimilarity index"
  ) +
  scale_color_manual(name="Schemes",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  
  theme_classic()
plot_TT_PC1


# PC1:Pair of schemes ----

# IV

plot_IV_PC1_schemes <- ggplot(similarity_index, aes(x = PC1, y = cor, color=Couple, fill=Couple)) +
  geom_point(size = 2, aes(color=Couple, fill=Couple)) +
  geom_smooth(method = "lm", se = TRUE,aes(color=Couple, fill=Couple, linetype = Couple), alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "PC1",
    y = "Similarity index"
  ) +
  scale_color_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  scale_fill_manual(name="Couple",
                    breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                    values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                    guide = "none") +
  scale_linetype_manual(
    name = "Couple",
    breaks = c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
    values = c("OP-UKBMS"="dotdash","STERF-UKBMS"="dotdash","OP-STERF"="longdash"),
    guide = "none") +
  theme_classic()
plot_IV_PC1_schemes

# Temporal trends

plot_TT_PC1_schemes <- ggplot(dissimilarity_index, aes(x = PC1, y = Distance, color=Schemes, fill=Schemes)) +
  geom_point(size = 2, aes(color=Schemes, fill=Schemes)) +
  geom_smooth(method = "lm", se = TRUE,aes(color=Schemes, fill=Schemes, linetype=Schemes), alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "PC1",
    y = "Similarity index"
  ) +
  scale_color_manual(name="Schemes",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  scale_fill_manual(name="Schemes",
                    breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                    values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                    guide = "none") +
  scale_linetype_manual(
    name = "Couple",
    breaks = c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
    values = c("OP-UKBMS"="dotdash","STERF-UKBMS"="dotdash","OP-STERF"="dotdash"),
    guide = "none") +
  theme_classic()
plot_TT_PC1_schemes


# PC2

# IV
plot_IV_PC2 <- ggplot(similarity_index, aes(x = PC2, y = cor)) +
  geom_point(size = 2, aes(color=Couple)) +
  geom_smooth(method = "lm", se = TRUE, fill="#dc3c84", color="#dc3c84", linetype = "dotdash", alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "PC2",
    y = "Similarity index"
  ) +
  scale_color_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  theme_classic()
plot_IV_PC2


# TT

plot_TT_PC2 <- ggplot(dissimilarity_index, aes(x = PC2, y = Distance)) +
  geom_point(size = 2, aes(color=Schemes)) +
  geom_smooth(method = "lm", se = TRUE, fill="#dc3c84", color="#dc3c84", linetype = "dotdash", alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "PC2",
    y = "Dissimilarity index"
  ) +
  scale_color_manual(name="Schemes",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  theme_classic()
plot_TT_PC2


library("patchwork")

# fig_4 <- (plot_IV_PC1|plot_TT_PC1)/(plot_IV_PC2|plot_TT_PC2)/(plot_VI_Pair_of_schemes|plot_TT_Pair_of_schemes)
fig_4 <- (plot_IV_PC1|plot_TT_PC1)/(plot_IV_PC1_schemes|plot_TT_PC1_schemes)/(plot_IV_PC2|plot_TT_PC2)

fig_4

# Save the plot in high resolution
ggsave(
  filename = "figure_4.jpeg",
  plot = fig_4,
  device = "jpeg",
  path = NULL,
  width = 8,
  height = 11,
  dpi = 300,
  bg = "white"
)


# Figure 5 : Emmeans, glmm --------------------------------------------------

# IV

plot_IV_count <- ggplot(similarity_index, aes(x = Total, y = cor, fill = Couple, color = Couple)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE,  alpha=0.2, aes(linetype=Couple), fullrange = TRUE) +
  labs(
    title = "",
    subtitle = "",
    x = "Total amount of data",
    y = "Interannual variation similarity index"
  ) +
  scale_color_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  scale_fill_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  scale_linetype_manual(
    name = "Couple",
    breaks = c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
    values = c("OP-UKBMS"="dotdash","STERF-UKBMS"="dotdash","OP-STERF"="solid"),
    guide = "none") +
  ylim(-0.6,1) +
  theme_classic()
plot_IV_count

# Temporal trends

plot_TT_count <- ggplot(dissimilarity_index, aes(x = Total, y = Distance, fill = Schemes, color = Schemes)) +
  geom_point(size = 2, aes(color=Schemes)) +
  geom_smooth(method = "lm", aes(linetype=Schemes), se = TRUE, alpha=0.2) +
  labs(
    title = "",
    subtitle = "",
    x = "Total amount of data",
    y = "Temporal trend dissimilarity index"
  ) +
  scale_color_manual(name="Schemes",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  scale_fill_manual(name="Couple",
                    breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                    values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                    guide = "none") +
  scale_linetype_manual(
    name = "Schemes",
    breaks = c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
    values = c("OP-UKBMS"="dotdash","STERF-UKBMS"="dotdash","OP-STERF"="dotdash"),
    guide = "none") +
  theme_classic()
plot_TT_count

# Pair of schemes

#IV
df_v <- similarity_index[,c("Groupe","cor","Couple")]
df_v <- df_v[complete.cases(df_v),]
v="Couple"

p.val.test<-pwpm(emmeans(model_IV_sch,  ~ Couple),means = FALSE, flip = TRUE,reverse = TRUE)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
stat.test<-subset(melt(p.matx),!is.na(value))
names(stat.test)<-c("group1","group2","p.adj")
stat.test[stat.test$p.adj<=0.001,"p.adj.signif"]<-"***"
stat.test[stat.test$p.adj>0.001 & stat.test$p.adj<=0.01,"p.adj.signif"]<-"**"
stat.test[stat.test$p.adj>0.01 & stat.test$p.adj<=0.05,"p.adj.signif"]<-"*"
stat.test[ stat.test$p.adj>0.05,"p.adj.signif"]<-"ns"
#stat.test<-mc_tribble(stat.test) 

plot_VI_Pair_of_schemes <- ggplot(df_v, aes(x= Couple, y=cor, color=Couple)) + # fill=Variable avec interactions
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  #theme(axis.text=element_text(size=12))+
  #adjust y label positions
  # theme(axis.title.y = element_text(margin = margin (r = 10)),
  #       #change the plot margins
  #       plot.margin = margin(l = 15,r=10)) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) + # ,face="bold"
  labs(x="Pair of schemes", y="Interannual variation similarity index") + # ,title=paste("Correlation vs.",v)
  scale_color_manual(name="Couple",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  stat_pvalue_manual(stat.test,
                     y.position=1.2,step.increase = 0.1,
                     tip.length = 0.01, bracket.shorten = 0.05,
                     #y.position = max(df_v$Distance)+sd(df_v$Distance),
                     color = "black",
                     size=3)
plot_VI_Pair_of_schemes

# Temporal trends
df_v <- dissimilarity_index[,c("Groupe","Distance","Schemes")]
df_v <- df_v[complete.cases(df_v),]
v="Schemes"

p.val.test<-pwpm(emmeans(model_TT_sch,  ~ Schemes),means = FALSE, flip = TRUE,reverse = TRUE)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
stat.test<-subset(melt(p.matx),!is.na(value))
names(stat.test)<-c("group1","group2","p.adj")
stat.test[stat.test$p.adj<=0.001,"p.adj.signif"]<-"***"
stat.test[stat.test$p.adj>0.001 & stat.test$p.adj<=0.01,"p.adj.signif"]<-"**"
stat.test[stat.test$p.adj>0.01 & stat.test$p.adj<=0.05,"p.adj.signif"]<-"*"
stat.test[ stat.test$p.adj>0.05,"p.adj.signif"]<-"ns"
#stat.test<-mc_tribble(stat.test) 
# stat.test <- tribble(
#   ~group1, ~group2, ~p.adj, ~p.adj.signif,
#   "Non conspicuous", "Conspicuous", 0.0022, "**")

plot_TT_Pair_of_schemes <- ggplot(df_v, aes(x= Schemes, y=Distance, color=Schemes)) + # fill=Variable avec interactions
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  #theme(axis.text=element_text(size=12))+
  #adjust y label positions
  # theme(axis.title.y = element_text(margin = margin (r = 10)),
  #       #change the plot margins
  #       plot.margin = margin(l = 15,r=10)) +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=11)) + # ,face="bold"
  labs(x="Pair of schemes", y="Temporal trend dissimilarity index") + # ,title=paste("Correlation vs.",v)
  scale_color_manual(name="Schemes",
                     breaks=c("OP-UKBMS", "STERF-UKBMS", "OP-STERF"),
                     values=c("OP-UKBMS"="#c6dc3c", "STERF-UKBMS"="#dc3cd3", "OP-STERF"="#1fad86"),
                     guide = "none") +
  stat_pvalue_manual(stat.test,
                     y.position=0.2,step.increase = 0.1,
                     tip.length = 0.01, bracket.shorten = 0.05,
                     #y.position = max(df_v$Distance)+sd(df_v$Distance),
                     color = "black",
                     size=3) 
plot_TT_Pair_of_schemes

fig_5 <- (plot_VI_Pair_of_schemes|plot_TT_Pair_of_schemes)/(plot_IV_count|plot_TT_count)
fig_5

# Save the plot in high resolution
ggsave(
  filename = "figure_5_corr.jpeg",
  plot = fig_5,
  device = "jpeg",
  path = NULL,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

