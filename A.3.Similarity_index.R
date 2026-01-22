library("effects")
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("ggrepel")
library("glmmTMB")
library("forcats")
library("tidyr")

schemes_count <- read.csv("Appendix_S2c_corr.csv", sep=";")
names(schemes_count)[names(schemes_count)=="GROUP"] <- "Group"

nom <- read.csv("data/Traits.csv", sep=";")

#'----------------------------------------------------------------------------
#'---------------------------- STERF -----------------------------------------
#'----------------------------------------------------------------------------

#sterf <- read.csv2("data/STERF/TMB_VI/Interannual_variations_sterf.csv", sep=";")
sterf <- read.csv2("data/STERF/TMB_VI/Interannual_variations_sterf_corr.csv", sep=";")

sterf$mean_fit <- as.numeric(as.character(sterf$mean_fit))
#sterf <- merge(sterf, nom[,c("French_name", "Scientific_name")], by.x="Espece", by.y="French_name")
sterf <- merge(sterf, schemes_count[schemes_count$Scheme=="STERF",c("Group", "Scientific_name", "Total", "NbPresence")], by.x="Espece", by.y="Group")
sterf$Scheme <- "STERF"

#sterf <- merge(sterf, nom[,c("French_name", "Scientific_name")], by.x="Espece", by.y="French_name")

#'----------------------------------------------------------------------------
#'---------------------------- ukBMS -----------------------------------------
#'----------------------------------------------------------------------------

ukbms <- read.csv2("data/UKBMS/TMB_VI/Interannual_variations_ukbms.csv", sep=";")
ukbms$mean_fit <- as.numeric(as.character(ukbms$mean_fit))
#ukbms <- merge(ukbms, nom[,c("French_name", "Scientific_name")], by.x="Espece", by.y="French_name")
ukbms <- merge(ukbms, schemes_count[schemes_count$Scheme=="UKBMS",c("Group", "Scientific_name", "Total", "NbPresence")], by.x="Espece", by.y="Group")
ukbms$Scheme <- "UKBMS"

#'----------------------------------------------------------------------------
#'---------------------------- OPJ -------------------------------------------
#'----------------------------------------------------------------------------

op <- read.csv2("data/OPJ/TMB_VI/Interannual_variations_op.csv", sep=";")
op$mean_fit <- as.numeric(as.character(op$mean_fit))
#op <- merge(op, nom[,c("French_name", "Scientific_name")], by.x="Espece", by.y="French_name")
op <- merge(op, schemes_count[schemes_count$Scheme=="OPJ",c("Group", "Scientific_name", "Total", "NbPresence")], by.x="Espece", by.y="Group")
op$Scheme <- "OP"

#'----------------------------------------------------------------------------
#'-------------------------- COMMON SPECIES -----------------------
#'----------------------------------------------------------------------------

# Selection des groupes communs aux 3 programmes pour tracer les graphiques des comparaisons interannuelles entre les 3 programmes (possible de le faire pour 2 programmes aussi)
op$Scientific_name[op$Scientific_name=="admirals"] <- "white admirals"
sterf$Scientific_name[sterf$Scientific_name=="admirals"] <- "white admirals"
ukbms$Scientific_name[ukbms$Scientific_name=="admirals"] <- "white admirals"

Groupes <- intersect(sterf$Scientific_name, op$Scientific_name)
Groupes <- intersect(Groupes, ukbms$Scientific_name)


#'----------------------------------------------------------------------------
#'---------------------------- GRAPHIQUES ------------------------------------
#'----------------------------------------------------------------------------


for (e  in Groupes) {
  
  # varcwm_e <- rbind(sterf[sterf$Espece==e,],ukbms[ukbms$Espece==e,]) # en français
  # varcwm_e<- rbind(varcwm_e, opj[opj$Espece==e,])
  
  varcwm_e <- rbind(sterf[sterf$Scientific_name==e,],ukbms[ukbms$Scientific_name==e,])
  varcwm_e <- rbind(varcwm_e, op[op$Scientific_name==e,])
  
  # text_trend <- paste(" ukBMS :", Comptage_ukbms$NbTotal[Comptage_ukbms$Scientific_name == e],"data", "(","including", Comptage_ukbms$NbPresence[Comptage_ukbms$Scientific_name == e]," observations)", "\n",
  #                     "OPJ :", Comptage_OPJ$NbTotal[Comptage_OPJ$Scientific_name == e],"data", "(","including", Comptage_OPJ$NbPresence[Comptage_OPJ$Scientific_name == e]," observations)", "\n",
  #                     "STERF :", Comptage_STERF$NbTotal[Comptage_STERF$Scientific_name == e],"data", "(", "including", Comptage_STERF$NbPresence[Comptage_STERF$Scientific_name == e]," observations)","\n")
  text_trend <- paste(" UKBMS :", unique(ukbms$Total[ukbms$Scientific_name == e]),"data","\n",
                      "OP :", unique(op$Total[op$Scientific_name == e]),"data","\n",
                      "STERF :", unique(sterf$Total[sterf$Scientific_name == e]),"data","\n")

  
  jpeg(paste(e, "jpeg", sep = "."), width = 20, height =12, units="cm", quality=75, res=300)
  a<-(ggplot(varcwm_e, aes(year,mean_fit,col = Scheme))+
        #geom_point()+
        geom_line(aes(x = year, y = mean_fit, color = Scheme, group = Scheme))+
        geom_rangeframe() +
        #theme(aspect.ratio = 1) +
        theme_classic() +
        scale_color_manual(values=c('#69c223', '#23b9c2', '#a92eb7')) + #OPJ STERF UKBMS
        labs(title = paste0("", e), # Tendance temporelle d'abondance des populations d'imagos : 
             subtitle = text_trend,  
             x="Year", y = "Index of abundance") + #ajouter text_trend ? subtitle # Indice d'abondance
        geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Scheme=="OP",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#69c223' ) +
        geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Scheme=="STERF",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#3cd3dc') +
        geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Scheme=="UKBMS",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#a92eb7') + ##F7230C
        # geom_smooth(method="glm",method.args = list(family = poisson(link = "log")),  data = varcwm_e[varcwm_e$Programme=="OPJ",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#69c223' ) +
        # geom_smooth(method="glm", method.args = list(family = poisson(link = "log")),  data = varcwm_e[varcwm_e$Programme=="STERF",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#3cd3dc') +
        # geom_smooth(method="glm", method.args = list(family = poisson(link = "log")), data = varcwm_e[varcwm_e$Programme=="ukBMS",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#F7230C') +
        theme(aspect.ratio = 1, plot.title = element_text(size = 10, face = "bold"), axis.text.x= element_text(angle = 30 , size = 2)) +
        # geom_errorbar(data = varcwm.O,aes(ymin=stderrl, ymax=stderru), width=.2, position=position_dodge(0.05))+
        #   geom_errorbar(data = varcwm.P,aes(ymin=stderrl, ymax=stderru), width=.2, position=position_dodge(0.05))+
        #   geom_errorbar(data = varcwm.S,aes(ymin=stderrl, ymax=stderru), width=.2, position=position_dodge(0.05))+
        #theme(aspect.ratio = 1) +
        #geom_ribbon(aes(ymin=stderrl, ymax=stderru),alpha=0.2) +
        #scale_y_continuous(lim = c(round(min(varcwm$cwm),2), round(max(varcwm$cwm),2)), scales::pro) +
        theme(plot.subtitle = element_text(size = 8), aspect.ratio = 1, plot.title = element_text(size = 10, face = "italic"), axis.text.x= element_text(angle = 30, size = 8), axis.text.y= element_text( size = 8)) +
        
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        #theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=18), axis.text.y = element_text(size=18))+
        #labs(title=paste0("Tendances temporelles - ",e),) +
        scale_x_continuous(breaks = seq(2006, 2019, 1), lim = c(2006, 2019)))
  
  print(a)  
  dev.off()
  # rm(p,g)
  rm(varcwm)
}

cvdPlot(a)

sterf$Scientific_name[sterf$Scientific_name=="Small Blue"] <- "Small Blue Lycenidae"
opj$Scientific_name[opj$Scientific_name=="Small Blue"] <- "Small Blue Lycenidae"
ukbms$Scientific_name[ukbms$Scientific_name=="Small Blue"] <- "Small Blue Lycenidae"

sterf$Scientific_name[sterf$Scientific_name=="Large White"] <- "Large White Pieridae"
opj$Scientific_name[opj$Scientific_name=="Large White"] <- "Large White Pieridae"
ukbms$Scientific_name[ukbms$Scientific_name=="Large White"] <- "Large White Pieridae"

Comptage_STERF$Scientific_name[Comptage_STERF$Scientific_name=="Small Blues"] <- "Small Blue Lycenidae"
Comptage_OPJ$Scientific_name[Comptage_OPJ$Scientific_name=="Small Blues"] <- "Small Blue Lycenidae"
Comptage_ukbms$Scientific_name[Comptage_ukbms$Scientific_name=="Small Blues"] <- "Small Blue Lycenidae"

Comptage_STERF$Scientific_name[Comptage_STERF$Scientific_name=="Large White"] <- "Large White Pieridae"
Comptage_OPJ$Scientific_name[Comptage_OPJ$Scientific_name=="Large White"] <- "Large White Pieridae"
Comptage_ukbms$Scientific_name[Comptage_ukbms$Scientific_name=="Large White"] <- "Large White Pieridae"

Comptage_STERF$Scientific_name[Comptage_STERF$Scientific_name=="White Admiral"] <- "Admirals"
Comptage_OPJ$Scientific_name[Comptage_OPJ$Scientific_name=="White Admiral"] <- "Admirals"
Comptage_ukbms$Scientific_name[Comptage_ukbms$Scientific_name=="White Admiral"] <- "Admirals"

graph <- function(e){
  
  varcwm_e <- rbind(sterf[sterf$Scientific_name==e,],ukbms[ukbms$Scientific_name==e,])
  varcwm_e<- rbind(varcwm_e, op[op$Scientific_name==e,])
  
  
  # text_trend <- paste(" ukBMS :", Comptage_ukbms$NbTotal[Comptage_ukbms$Scientific_name == e],"data", "(","including", Comptage_ukbms$NbPresence[Comptage_ukbms$Scientific_name == e]," observations)", "\n",
  #                     "OPJ :", Comptage_OPJ$NbTotal[Comptage_OPJ$Scientific_name == e],"data", "(","including", Comptage_OPJ$NbPresence[Comptage_OPJ$Scientific_name == e]," observations)", "\n",
  #                     "STERF :", Comptage_STERF$NbTotal[Comptage_STERF$Scientific_name == e],"data", "(", "including", Comptage_STERF$NbPresence[Comptage_STERF$Scientific_name == e]," observations)","\n")
  text_trend <- paste(" UKBMS :", unique(ukbms$Total[ukbms$Scientific_name == e]),"data","\n",
                      "OP :", unique(op$Total[op$Scientific_name == e]),"data","\n",
                      "STERF :", unique(sterf$Total[sterf$Scientific_name == e]),"data","\n")
  
  a<-(ggplot(varcwm_e, aes(year,mean_fit,col = Scheme))+
        geom_line(size=0.8)+
        geom_rangeframe() +
        theme_classic() +
        scale_color_manual(values=c('#69c223', '#23b9c2', '#a92eb7')) + #OPJ Propage STERF
        labs(title = paste0("", e), # Tendance temporelle d'abondance des populations d'imagos : 
             subtitle = text_trend,  x="Year", y = "Index of abundance") + #ajouter text_trend ? subtitle # Indice d'abondance
        
        # geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Programme=="OPJ",],se=TRUE, fullrange=TRUE, aes(y = fit, x = year), size=1.2,  alpha = 0.1, fill = '#69c223' ) +
        # geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Programme=="STERF",],se=TRUE, fullrange=TRUE, aes(y = fit, x = year), size=1.2, alpha = 0.1, fill = '#3cd3dc') +
        # geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Programme=="UKBMS",],se=TRUE, fullrange=TRUE, aes(y = fit, x = year), size=1.2, alpha = 0.1, fill = '#a92eb7') + ##F7230C
        # theme(aspect.ratio = 1, plot.title = element_text(size = 14, face = "bold"), axis.text.x= element_text(angle = 30 , size = 14)) +
        # theme(plot.subtitle = element_text(size = 10), aspect.ratio = 1, 
        #       plot.title = element_text(size = 14, face = "italic"), 
        #       axis.text.x= element_text(angle = 30, size = 12, h=1), 
        #       axis.text.y= element_text( size = 12),
        #       axis.title=element_text(size=14,face="bold")) +
        # theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      # scale_x_continuous(breaks = seq(2006, 2019, 1), lim = c(2006, 2019)))
      
      geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Scheme=="OP",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#69c223' ) +
        geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Scheme=="STERF",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#3cd3dc') +
        geom_smooth(formula = y~ x, method=lm, data = varcwm_e[varcwm_e$Scheme=="UKBMS",],se=TRUE, fullrange=TRUE, aes(y = mean_fit, x = year), alpha = 0.1, fill = '#a92eb7') + ##F7230C
        theme(aspect.ratio = 1, plot.title = element_text(size = 14, face = "bold")) +
        theme(plot.subtitle = element_text(size = 10), aspect.ratio = 1, 
              plot.title = element_text(size = 14, face = "italic"), 
              axis.text.x= element_text(angle = 30, size = 10, h=1), 
              axis.text.y= element_text( size = 10),
              axis.title=element_text(size=12,face="bold")) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        scale_x_continuous(breaks = seq(2006, 2019, 1), lim = c(2006, 2019)))
  print(a)
  
}



sort(Groupes)

u <- graph("white admirals")

a <- graph("Aglais io")
b <- graph("Aglais urticae")

#v <- graph("Antocharis cardamines")
c <- graph("Aphantopus hyperantus")

d <- graph("Argynnis paphia")
e <- graph("Coenonympha pamphilus")

f <- graph("Colias crocea")
g <- graph("Gonepteryx rhamni")

h <- graph("white Pieridae")
i <- graph("Lasiommata maera")

j <- graph("Lycaena phlaeas")
k <- graph("Maniola jurtina")

l <- graph("Melanargia galathea")
m <- graph("orange skippers")

#n <- graph("Papilio machaon")
o <- graph("Pararge aegeria")

p <- graph("Polygonia c-album")
q <- graph("Pyronia tithonus")

r <- graph("small blue Lycenidae")
s <- graph("Vanessa atalanta")

t <- graph("Vanessa cardui")

# setwd("E:/Papier-Papillons/Scripts_R/Figures/Variations_interannuelles/Annexes/")

jpeg(paste("Apendix_S0", "jpeg", sep = "."), width = 20, height = 10, units="cm", quality=75, res=300) # pour 4 lignes

plot_Sa <- ggarrange(
  #u,a,
  #b,v,
  #c,d,
  #e,f,
  #g,h,
  #i,j,
  #k,l,
  #m,n,
  #o,p,
  #q,r,
  #s,t,
  # u,v,
  
  t, e,
  ncol=2,
  nrow=1,
  labels=c("a","b"),
  font.label=list(size = 14),
  common.legend = TRUE)
t
plot_Sa
dev.off()

annexe_S4_1 <- (u|a)/(b|c)/(d|e)
annexe_S4_2 <- (f|g)/(h|i)/(j|k)
annexe_S4_3 <- (l|m)/(o|p)/(q|r)
annexe_S4_4 <- (s|t)

# Save the plot in high resolution
ggsave(
  filename = "annexe_S4_1_IV.jpeg",
  plot = annexe_S4_1,
  device = "jpeg",
  path = NULL,
  width = 11,
  height = 15,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Figure_2_Vc_Cp.jpeg",
  plot = plot_Sa,
  device = "jpeg",
  path = NULL,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"
)

jpeg(paste("Vanessa_cardui", "jpeg", sep = "."), width = 20, height = 10, units="cm", quality=75, res=300) # pour 4 lignes
t
dev.off()

jpeg(paste("Aphantopus_hyperantus", "jpeg", sep = "."), width = 20, height = 10, units="cm", quality=75, res=300) # pour 4 lignes
t
dev.off()

jpeg(paste("Test", "jpeg", sep = "."), width = 20, height = 10, units="cm", quality=75, res=300) # pour 4 lignes

plot1 <- ggarrange(t,c,
                   ncol=2,
                   nrow=1,
                   labels=c("a","b"),
                   font.label=list(size = 14),
                   common.legend = TRUE)
plot1
dev.off()


jpeg(paste("Dist_Cor_NbTotal", "jpeg", sep = "."), width = 20, height = 10, units="cm", quality=75, res=300) # pour 4 lignes

plot <- ggarrange(plot_cor, plot_dist,
                  ncol=2, nrow=1,
                  labels=c("a","b"),
                  font.label=list(size = 12),
                  common.legend = TRUE)
plot
dev.off()

cvdPlot(plot)


ggplot(varcwm_e, aes(x = year, y = fit, col = Programme)) +
  geom_point(aes(y = fit)) +
  geom_line(aes(y = fit, col = Programme)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  theme_minimal()

#'----------------------------------------------------------------------------
#'------------------------- TABLEAU SIMILITUDE -------------------------------
#'----------------------------------------------------------------------------

for (e  in Groupes) {
  print(e)
  
  bo <- subset(op, op$Scientific_name==e)
  bo <-  bo[order(bo[,"year"]), ]
  
  bu <- subset(ukbms, ukbms$Scientific_name==e)
  bu <-  bu[order(bu[,"year"]), ]
  
  bs <- subset(sterf, sterf$Scientific_name==e)
  bs <-  bs[order(bs[,"year"]), ]
  
  o <- bo$mean_fit
  u <- bu$mean_fit
  s <- bs$mean_fit
  
  ou <- cor.test(o, u , method = "pearson")
  us <- cor.test(u, s , method = "pearson")
  so <- cor.test(s, o , method = "pearson")
  
  a <- ccf (o, u, lag = 0, correlation = TRUE, pl = FALSE)
  b <- ccf (u, s, lag = 0, correlation = TRUE, pl = FALSE)
  c <- ccf (s, o, lag = 0, correlation = TRUE, pl = FALSE)
  
  cor <- data.frame(Groupe = e, cor.OU = round(ou$estimate,3), cor.US = round(us$estimate,3), cor.SO = round(so$estimate,3) ,pv.OU = round(ou$p.value,3), pv.US = round(us$p.value,3),  pv.SO = round(so$p.value,3))
  cor.ccf <- data.frame(Groupe = e, cor.ccf.OU = a$acf[1], cor.ccf.US = b$acf[1], cor.ccf.SO = c$acf[1]) # Comparer des series annuelles
  cor_e <- merge(cor, cor.ccf, by="Groupe")
  
  
  if (e == Groupes[1]) Total.cor<- cor_e else Total.cor <- rbind(Total.cor,cor_e)
  
}

write.table(Total.cor, "data/Tableau_correlation_groupes_traits_final_cor.csv", sep = ";", row.names = FALSE) # presence de p-value
Total.cor <- read.csv("data/Tableau_correlation_groupes_traits_final_cor.csv", sep=";")

Total.cor$Groupe[Total.cor$Groupe=="admirals"] <- "white admirals"

Total.cor[,c("Sign.OU","Sign.US","Sign.SO")] <- "NA"
for (e in unique(Total.cor$Groupe)){ 
  if (Total.cor$pv.OU[Total.cor$Groupe == e] <= 0.05) {Total.cor$Sign.OU[Total.cor$Groupe == e] = "S"}
  if (Total.cor$pv.OU[Total.cor$Groupe == e] > 0.05) {Total.cor$Sign.OU[Total.cor$Groupe == e] = "NS"}
  
  if (Total.cor$pv.US[Total.cor$Groupe == e] <= 0.05) {Total.cor$Sign.US[Total.cor$Groupe == e] = "S"}
  if (Total.cor$pv.US[Total.cor$Groupe == e] > 0.05) {Total.cor$Sign.US[Total.cor$Groupe == e] = "NS"}
  
  if (Total.cor$pv.SO[Total.cor$Groupe == e] <= 0.05) {Total.cor$Sign.SO[Total.cor$Groupe == e] = "S"}
  if (Total.cor$pv.SO[Total.cor$Groupe == e] > 0.05) {Total.cor$Sign.SO[Total.cor$Groupe == e] = "NS"}
  
}
Total.cor <-  Total.cor[order(Total.cor[,"cor.US"]), ] # on trie les sp par estimates decroissants


Total.cor$Groupe <- factor(Total.cor$Groupe,
                                     levels = c("white Pieridae", "Gonepteryx rhamni", "Colias crocea", # Pieridae
                                                "Vanessa cardui", "Vanessa atalanta", "Pyronia tithonus", "Polygonia c-album", "Pararge aegeria", 
                                                "Melanargia galathea", "Maniola jurtina", "Lasiommata maera", "Coenonympha pamphilus", 
                                                "white admirals", "Argynnis paphia", "Aphantopus hyperantus", "Aglais urticae", "Aglais io",  # Nymphalidae
                                                "small blue Lycenidae", "Lycaena phlaeas", # Lycaenidae
                                                "orange skippers" # Hesperiidae
                                     ))

jpeg(paste("Engl_Selected_model_Comparaisons_Variations_interannuelle_dalto_cor", "jpeg", sep = "."), width = 10, height =12, units="cm", quality=75, res=300)

p <- ggplot(Total.cor, aes(x=Groupe, y=cor.OU)) +
  coord_flip() + theme_bw() + 
  geom_point(data = Total.cor, aes(x = Groupe, y = cor.OU, shape=Sign.OU, col = "#c6dc3c"),size=2,  show.legend=FALSE) + # "#dc603c"
  geom_point(data = Total.cor, aes(x = Groupe, y = cor.US,  shape=Sign.US, col = "#dc3cd3"),size=2, show.legend=FALSE) + 
  geom_point(data = Total.cor, aes(x = Groupe, y = cor.SO,  shape=Sign.SO, col = "#1fad86"),size=2, show.legend=FALSE) +
  geom_hline(aes(yintercept=0.6), color="grey", linetype="dashed")+
  geom_hline(aes(yintercept=-0.3), color="grey", linetype="dashed")+
  ylab("Interannual variation similarity index") + xlab("Species/species group")  +
  scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
  scale_color_manual(name='Schemes',
                     breaks=c('OPJ-UKBMS', 'UKBMS-STERF', 'STERF-OPJ'),
                     values=c('#c6dc3c'='#c6dc3c', '#dc3cd3'='#dc3cd3', '#1fad86'='#1fad86')) +
  scale_shape_manual(values=c(1,19)) + 
  theme(axis.text.y = element_text(size=8, face="italic"))
p
dev.off()

cvdPlot(p)

for (e  in Groupes) {
  print(e)
  
  o <- opj[opj$Espece==e,"fit"]
  s <- sterf[sterf$Espece==e,"fit"]
  
  so <- cor.test(s, o , method = "pearson")
  
  c <- ccf (s, o, lag = 0, correlation = TRUE, pl = FALSE)
  
  # cor <- data.frame(Groupe = e, cor.OU = round(ou$estimate,3), cor.US = round(us$estimate,3), cor.SO = round(so$estimate,3) ,pv.OU = round(ou$p.value,3), pv.US = round(us$p.value,3),  pv.SO = round(so$p.value,3))
  cor.ccf <- data.frame(Groupe = e, cor.SO = c$acf[1]) # Comparer des series annuelles
  # cor_e <- merge(cor, cor.ccf, by="Groupe")
  
  
  if (e == Groupes[1]) Total.cor<- cor.ccf else Total.cor <- rbind(Total.cor,cor.ccf)
  
}

Total.cor <-  Total.cor[order(Total.cor[,"cor.SO"]), ] # on trie les sp par estimates decroissants

ggplot(Total.cor, aes(x=fct_inorder(Groupe), y=cor.SO, col=)) +
  coord_flip() + theme_bw() + 
  geom_point(data = Total.cor, aes(x = fct_inorder(Groupe), y = cor.SO), col = "#1fad86") +
  geom_hline(aes(yintercept=0.6), color="grey", linetype="dashed")+
  geom_hline(aes(yintercept=-0.3), color="grey", linetype="dashed")+
  ylab("Corr?lation") + xlab("Esp?ce ou groupe d'esp?ces")  +
  #scale_x_discrete(limits = as.vector(Bayes_propage$NAME)) + # geom_bar(stat = "identity") +
  theme(axis.text.y = element_text(size=4))


# Table

Total.corr_pivot <- Total.cor %>%
  pivot_longer(
    cols = matches("(cor\\.|pv\\.|cor\\.ccf\\.|Sign\\.)"),  # toutes les colonnes concernées
    names_to = c(".value", "Couple"),  # sépare le nom de la variable et le couple
    names_pattern = "(.*)\\.(OU|US|SO)"  # capture 'cor', 'pv', 'Sign', etc. + le suffixe OU/US/SO
  )

Total.corr_pivot$Couple[Total.corr_pivot$Couple=="OU"] <- "OP-UKBMS"
Total.corr_pivot$Couple[Total.corr_pivot$Couple=="US"] <- "STERF-UKBMS"
Total.corr_pivot$Couple[Total.corr_pivot$Couple=="SO"] <- "OP-STERF"

write.table(Total.corr_pivot, "data/Tableau_correlation_groupes_traits_final_pivot_cor.csv", sep = ";", row.names = FALSE) # presence de p-value
