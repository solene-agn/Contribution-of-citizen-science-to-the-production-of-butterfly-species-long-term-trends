library("data.table")
library("speedglm") # fitting linear and generalized linear models to large data sets
library("devtools")
library("ggplot2")
library("reshape2")
library("plyr")
library("sf")
library("terra")
library("tmap")
library("DT")
library("mapview")
library("rbms")

source("workshop_functions.R")

# Script basé sur le workshop https://butterfly-monitoring.github.io/bms_workshop/index.html , et la documentation annexe https://butterfly-monitoring.github.io/bms_workshop/GAI%20approach.pdf
# Data importation --------------------------------------------------

obs <- data.table::fread("data/data/STERF_groups_corrige_mai_aout_nationale.csv")
names(obs)[names(obs)=="GROUP"] <- "SPECIES"

obs$LONGITUDE <- as.numeric(as.character(obs$LONGITUDE))
obs$LATITUDE <- as.numeric(as.character(obs$LATITUDE))
obs$DAY <- lubridate::day(obs$DATE)

obs <- subset(obs, obs$REGION_R %in% c("Normandie", "Île-de-France", "Hauts-de-France", "Bretagne", "Grand Est")) # data for north of France

# Cas 1 : Si on considère les transects comme des sites :
obs$Site_transect_ID <- paste0(obs$SITE_ID, "_", obs$TRANSECT_ID)
obs <- obs[,c("Site_transect_ID", "DATE", "YEAR", "MONTH", "DAY", "SPECIES", "COUNT", "LONGITUDE", "LATITUDE")]
names(obs)[names(obs)=="Site_transect_ID"] <- "SITE_ID"


# Cas 2 : Si on regroupe les transects à l'échelle du site : 

# Importation des infos sur les sites et les transects
site <- read.csv2("data/data/sites.csv", sep=",")
transect <- read.csv2("data/data/transects.csv", sep = ",")
colnames(transect)[3] <- "id.site.ref"
colnames(transect)[7] <- "geom_transect"
site <- merge(site, transect[,c("id.site.ref", "id_transect")], by="id.site.ref", all.x=FALSE, all.y=TRUE)

#Recuperation des donnees geographiques des transects echantillones dans un site
site_distrib <- transect[,c("id.site.ref", "geom_transect")]

#Spatialisation des donnees
st_geometry(site_distrib) <- st_as_sfc(site_distrib$geom_transect)
st_crs(site_distrib) <- 4326

#Recuperation des coordonnees des centroides de tous les transects d'un site
site_distrib <- cbind(SITE_ID = site_distrib$id.site.ref,st_coordinates(st_centroid(site_distrib)))

#Calcul du centroide unique moyen des transects d'un site pour obtenir une valeur geographique par site
site_distrib <- aggregate(cbind(X,Y) ~ SITE_ID, data = site_distrib, mean)
site_distrib$SITE_ID <- as.factor(as.character(site_distrib$SITE_ID))
obs$SITE_ID <- as.factor(as.character(obs$SITE_ID))
sterf_count <- aggregate(cbind(COUNT = COUNT) ~ SITE_ID + DATE + SPECIES + DAY + MONTH + YEAR, data = obs, sum)
sterf_count <- merge(sterf_count, site_distrib, by.x="SITE_ID", by.y="SITE_ID", all.x=TRUE, all.y=FALSE)

names(sterf_count)[names(sterf_count)=="X"] <- "LONGITUDE"
names(sterf_count)[names(sterf_count)=="Y"] <- "LATITUDE"
obs <- sterf_count


# Computing site and collated indices -----------------------------


# Generalized Abundance Index ====

# Creation de la serie temporelle correspondant aux donnees
ts_date <- ts_dwmy_table(InitYear = min(obs$YEAR),
                         LastYear = max(obs$YEAR) ,
                         WeekDay1 = 'monday') #Lundi est indique comme premier jour de la semaine (par defaut dimanche)

# Creation des tables de donnees pour les phenologies et l'analyse temporelle 
ts_season <- ts_monit_season(ts_date,
                             StartMonth = 5, # 4 pour le sterf habituellement
                             EndMonth = 8, # 9
                             Anchor = TRUE,
                             AnchorLength = 3,
                             AnchorLag = 1,
                             TimeUnit = 'w')

MY_visit_region <- unique(data.frame(SITE_ID = obs$SITE_ID, DATE = obs$DATE))
ts_season_visit <- rbms::ts_monit_site(ts_season, MY_visit_region)

MY_count_region <- obs

# Species list
species <- "Belle-dame" # pour tester pour le moment
s_sp <- "Belle-dame"

# for (s_sp in unique(obs$SPECIES)) {
for (s_sp in species) {
  
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, MY_count_region, sp = s_sp)
  
  ts_flight_curve <- flight_curve(ts_season_count, NbrSample = NULL, MinVisit = 3, MinOccur = 2, MinNbrSite = 3,
                                  MaxTrial = 3, GamFamily = 'nb', SpeedGam = FALSE, SelectYear = NULL,
                                  TimeUnit = 'w')
  
  pheno <- ts_flight_curve$pheno
  
  
  # Impute & Site index ====
  impt_counts <- rbms::impute_count(ts_season_count = ts_season_count, 
                                    ts_flight_curve = pheno, 
                                    YearLimit = NULL, 
                                    TimeUnit = "w")
  
  sindex <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.10)
  
  #saveRDS(sindex, file.path("bms_workshop_data", paste( gsub(" ", "_", s_sp), paste(region_bms, collapse="_"), "sindex.rds", sep="_")))
  
  
  # Collated index ====
  co_index <- collated_index(data = sindex, 
                             s_sp = s_sp,
                             sindex_value = "SINDEX",
                             glm_weights = TRUE,
                             rm_zero = TRUE)
  co_index$col_index
  
  
  # Bootstrap CI ====
  bootsample <- rbms::boot_sample(sindex, boot_n = 500)
  co_index <- list()
  
  pb <- txtProgressBar(min = 0, max = dim(bootsample$boot_ind)[1], initial = 0, char = "*",  style = 3)
  
  for (i in c(0, seq_len(dim(bootsample$boot_ind)[1]))) {
    co_index[[i + 1]] <- rbms::collated_index(data = sindex,
                                              s_sp = s_sp,
                                              sindex_value = "SINDEX",
                                              bootID = i,
                                              boot_ind = bootsample,
                                              glm_weights = TRUE,
                                              rm_zero = TRUE)
    
    ## for progression bar, uncomment the following
    setTxtProgressBar(pb, i)
  }
  
  ## collate and append all the result in a data.table format
  co_index <- rbindlist(lapply(co_index, FUN = "[[", "col_index"))
  co_index$SPECIES <- s_sp
  
  co_index_boot <- co_index
  #    saveRDS(co_index_boot, file.path("bms_workshop_data", paste( gsub(" ", "_", s_sp), bms_id, "co_index_boot.rds", sep="_")))
  
  
  # Figure with CI ==== 
  
  #co_index_boot <- readRDS(file.path("bms_workshop_data", paste( gsub(" ", "_", s_sp), bms_id, "co_index_boot.rds", sep="_")))
  
  boot0 <- co_index_boot[ , mean(COL_INDEX, na.rm=TRUE), by = M_YEAR]
  boot0[, LOGDENSITY0 := log(V1) / log(10)]
  boot0[, meanlog0 := mean(LOGDENSITY0)]
  
  # Add mean of BOOTi == 0 (real data) to all data
  boot_ests <- merge(co_index_boot, boot0[, .(M_YEAR, meanlog0)],
                     by = c("M_YEAR")
  )
  
  boot_ests[, LOGDENSITY := log(COL_INDEX) / log(10)]
  boot_ests[, TRMOBSCI := LOGDENSITY - meanlog0 + 2]
  
  boot_ests[, TRMOBSCILOW := quantile(TRMOBSCI, 0.025), by = M_YEAR]
  boot_ests[, TRMOBSCIUPP := quantile(TRMOBSCI, 0.975), by = M_YEAR]
  
  ylim_ <- c(min(floor(range(boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, TRMOBSCI]))),
             max(ceiling(range(boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, TRMOBSCI]))))
  
  plot_col <- ggplot(data = boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, ], aes(x = M_YEAR, y = TRMOBSCI)) + 
    geom_point(colour = "orange", alpha = 0.2) + ylim(ylim_) +
    geom_line(data = boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, median(TRMOBSCI, na.rm=TRUE), by = M_YEAR], aes(x = M_YEAR, y = V1), linewidth = 1, colour = "dodgerblue4") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "grey50") +
    scale_x_continuous(minor_breaks = waiver(), limits = c(2006, 2019), breaks = seq(2006, 2019, by = 1)) + 
    xlab("Year") + ylab("Abundance Incices (Log/Log(10))") +
    labs(subtitle = paste(s_sp, sep = " "))
  
  plot_col
  
  ggsave(filename = paste0(getwd(),"/Results/Trend_plots/Plot_col/plot_col_",s_sp,".jpg"), plot_col,
         width = 30, height = 20, units = "cm", dpi = 300)
  
  # To center the index on the first year and scale it to 100 instead of 2
  #co_index_boot <- readRDS(file.path("bms_workshop_data", paste( gsub(" ", "_", s_sp), bms_id, "co_index_boot.rds", sep="_")))
  
  boot0 <- co_index_boot[ , mean(COL_INDEX, na.rm=TRUE), by = M_YEAR]
  boot0[, LOGDENSITY0 := log(V1) / log(10)]
  boot0[, first_year_log0 := boot0[order(M_YEAR), LOGDENSITY0][1]]
  
  # Add mean of BOOTi == 0 (real data) to all data
  boot_ests <- merge(co_index_boot, boot0[, .(M_YEAR, first_year_log0)],
                     by = c("M_YEAR")
  )
  
  boot_ests[, LOGDENSITY := log(COL_INDEX) / log(10)]
  boot_ests[, TRMOBSCI := 10^(LOGDENSITY - first_year_log0 + 2)]
  
  boot_ests[, TRMOBSCILOW := quantile(TRMOBSCI, 0.025), by = M_YEAR]
  boot_ests[, TRMOBSCIUPP := quantile(TRMOBSCI, 0.975), by = M_YEAR]
  
  ylim_ <- c(min(floor(range(boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, TRMOBSCI]))),
             max(ceiling(range(boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, TRMOBSCI]))))
  
  plot_col <- ggplot(data = boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, ], aes(x = M_YEAR, y = TRMOBSCI)) + 
    geom_point(colour = "orange", alpha = 0.2) + ylim(ylim_) +
    geom_line(data = boot_ests[TRMOBSCI <= TRMOBSCIUPP & TRMOBSCI >= TRMOBSCILOW, median(TRMOBSCI, na.rm=TRUE), by = M_YEAR], aes(x = M_YEAR, y = V1), linewidth = 1, colour = "dodgerblue4") +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
    scale_x_continuous(minor_breaks = waiver(), limits = c(2006, 2019), breaks = seq(2006, 2019, by = 1)) + 
    xlab("Year") + ylab("Abundance Incices (first year set to 100)") +
    labs(subtitle = paste(s_sp, sep = " "))
  
  plot_col
  ggsave(filename = paste0(getwd(),"/Results/Trend_plots/Plot_col/plot_col_center_1st_year_",s_sp,".jpg"), plot_col,
         width = 30, height = 20, units = "cm", dpi = 300)

  
  # Calculating species trends ----
  spp <- s_sp
  
  # # Read in the bootstrapped collated indices
  # co_index <- readRDS(paste0("./bms_workshop_data/", spp, "_co_index_boot.rds"))
  
  ## Only use COL_INDEX larger then 0 for calculation or logdensity and trmobs
  co_index[COL_INDEX > 0, LOGDENSITY:= log(COL_INDEX)/log(10)][, TRMOBS := LOGDENSITY - mean(LOGDENSITY) + 2, by = .(BOOTi)]
  co_index
  
  sp_trend <- estimate_boot_trends(co_index)
  
  sp_trend
  
  Trend <- data.frame(SPECIES=s_sp,
                      sp_trend)
  
  # Calculate mean log index for original data
  co_index0_mean <- mean(co_index[BOOTi == 0]$LOGDENSITY, na.rm = TRUE)
  # Derive interval from quantiles
  co_index_ci <- merge(co_index[BOOTi == 0, .(M_YEAR, TRMOBS)],
                       co_index[BOOTi != 0, .(
                         LOWER = quantile(LOGDENSITY - co_index0_mean + 2, 0.025, na.rm = TRUE),
                         UPPER = quantile(LOGDENSITY - co_index0_mean + 2, 0.975, na.rm = TRUE)), 
                         by = .(M_YEAR)],
                       by=c("M_YEAR"))
  
  co_index_ci$SPECIES <- s_sp
  
  write.table(co_index_ci, file=paste0(getwd(),"/Results/Interannual_variations_table/IV_",s_sp,".csv"), sep=";", row.names=FALSE)
  
  d <- lm(TRMOBS~M_YEAR, data=co_index_ci)
  Trend$co_index_ci_m_year_est[Trend$SPECIES==s_sp] <- coef(d)[2]
  Trend$co_index_ci_m_year_st_err[Trend$SPECIES==s_sp] <- summary(d)$coefficients["M_YEAR", "Std. Error"]
  
  write.table(Trend, file=paste0(getwd(),"/Results/Temporal_trends_table/Trend_",s_sp,".csv"), sep=";", row.names=FALSE)
  
  Trend_VI_plot <- ggplot(co_index_ci, aes(M_YEAR, TRMOBS))+
    theme(text = element_text(size = 14))+
    geom_line()+
    geom_point()+
    geom_ribbon(aes(ymin = LOWER, ymax = UPPER), alpha = .3)+
    geom_smooth(method="lm", se=FALSE, color="red")+
    xlab("Year")+ylab(expression('log '['(10)']*' Collated Index'))+
    ggtitle(paste("Collated index for", gsub("_", " ", spp)))
  Trend_VI_plot
  
  ggsave(filename = paste0(getwd(),"/Results/Trend_plots/Final/Trend_VI_plot_",s_sp,".jpg"), Trend_VI_plot,
         width = 30, height = 20, units = "cm", dpi = 300)
  
  
  if (s_sp == species[1]) {
    
    Trend_total <- Trend
    IV_total <- co_index_ci
  } 
  
  else {
    Trend_total <- rbind(Trend_total, Trend)
    IV_total <- rbind(IV_total, co_index_ci)
  }
  
  write.table(IV_total, file=paste0(getwd(),"/Results/Interannual_variations_table/Total_IV.csv"), sep=";", row.names=FALSE)
  write.table(Trend_total, file=paste0(getwd(),"/Results/Temporal_trends_table/Total_Trend_.csv"), sep=";", row.names=FALSE)
  
}

