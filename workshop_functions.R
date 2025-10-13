
library(ggplot2)
library(data.table)


# Estimate trends across bootstraps
estimate_boot_trends <- function(boot_data, 
                                 bycols = NULL, # columns to apply trends across
                                 tidy = TRUE,
                                 inparallel = FALSE, ncpus = NULL,
                                 minyear = NULL,
                                 maxyear = NULL){
  
  if(inparallel & is.null(ncpus)){
    # Determine number of cpus to use if not specified
    n_l_cores = detectCores()
    n_p_cores = detectCores(logical = FALSE)
    ncpus = n_l_cores - (n_l_cores/n_p_cores)
  } 
  
  # Estimate a trend for each bootstrap and species
  if(!inparallel){
    boot_trends <-  do.call(rbind, 
                            plyr::dlply(boot_data, c(bycols, "BOOTi"),
                                        trend_func, 
                                        bycols = bycols, 
                                        minyear = minyear,
                                        maxyear = maxyear))
  } else {
    cols <- c(bycols, "BOOTi")
    bycol_boot <- unique(boot_data[, ..cols])
    cl <- makeCluster(ncpus)
    clusterEvalQ(cl, library(data.table))
    clusterExport(cl, varlist = c("bycol_boot", "boot_data", 
                                  "bycols", "minyear", "maxyear",
                                  "trend_func"),
                  envir = environment())
    boot_trends <- do.call(rbind, parLapply(cl, 1:nrow(bycol_boot), boot_trend_func))
    stopCluster(cl)
    
  }  
  
  # Summarise bootstrap trends to get a CI
  if(!is.null(bycols)){
    trends_ci <- plyr::ddply(boot_trends , bycols,
                             trend_ci_func)      
  } else {
    trends_ci <- trend_ci_func(boot_trends)
  }
  
  # Classify long term trends
  trends_ci <- trend_class_func(trends_ci,
                                tcol = "TrendClass_lt")
  # Classify 10 years trends
  if("rate_10y_low" %in% colnames(trends_ci))
    trends_ci <- trend_class_func(trends_ci, 
                                  tcol = "TrendClass_10y", 
                                  low = "rate_10y_low", 
                                  upp = "rate_10y_upp")
  
  # Tidy up for output
  if(tidy) trends_ci <- trends_ci[order(trends_ci$TrendClass_lt),
                                  c(bycols, "data_from", "data_to", "data_nyears", 
                                    "minyear", "maxyear", "n",
                                    "nboot_lt",
                                    "rate_lt", "rate_lt_low", "rate_lt_upp",
                                    "pcn_lt", "TrendClass_lt",
                                    if("nboot_10y" %in% colnames(trends_ci)){c("nboot_10y",
                                                                               "rate_10y", "rate_10y_low", "rate_10y_upp",
                                                                               "pcn_10y", "TrendClass_10y")})]
  
  return(trends_ci)
}


# Minor function for parallel trend fitting
boot_trend_func <- function(iboot){
  boot_datai <- merge(boot_data, bycol_boot[iboot,], by = colnames(bycol_boot))
  trend_func(boot_datai, bycols, minyear, maxyear)
}


# Underlying function for (linear) trend calculation
trend_func <- function(collind_df, bycols = NULL, minyear = NULL, maxyear = NULL){
  setDT(collind_df)
  minyear_ <- collind_df[, min(M_YEAR)]
  maxyear_ <- collind_df[, max(M_YEAR)]
  if (is.null(minyear)) minyear <- collind_df[!is.na(TRMOBS), min(M_YEAR)]
  if (is.null(maxyear)) maxyear <- collind_df[!is.na(TRMOBS), max(M_YEAR)]
  if(!is.null(minyear)) collind_df <- collind_df[M_YEAR >= minyear]
  if(!is.null(maxyear)) collind_df <- collind_df[M_YEAR <= maxyear]
  
  if(nrow(collind_df) > 0){
    # Fit linear model
    lm_obj <- try(lm(TRMOBS ~ M_YEAR, collind_df), silent = TRUE)
    
    trend_df <- data.frame(BOOTi = collind_df$BOOTi[1],
                           data_from = collind_df[!is.na(TRMOBS), min(M_YEAR)],
                           data_to = collind_df[!is.na(TRMOBS), max(M_YEAR)],
                           data_nyears = length(collind_df$M_YEAR),
                           minyear = minyear,
                           maxyear = maxyear,
                           n = length(minyear:maxyear),
                           rate_lt = ifelse(!inherits(lm_obj, "try-error"), 
                                            exp(coef(lm_obj)[2]*2.303), NA),
                           pc1_lt = ifelse(!inherits(lm_obj, "try-error"),
                                           100*(exp(coef(lm_obj)[2]*2.303)-1), NA),
                           pcn_lt = ifelse(!inherits(lm_obj, "try-error"), 
                                           100*(exp(coef(lm_obj)[2]*2.303)^(
                                             maxyear-minyear)-1), NA))
    
    # Fit for last 10 years only if data have more than 10 years
    if(length(minyear_:maxyear_) > 10){
      collind_df10 <- collind_df[M_YEAR >= (maxyear-9)]
      lm_obj10 <- try(lm(TRMOBS ~ M_YEAR, collind_df10), silent = TRUE)
      
      trend_df$rate_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                                  exp(coef(lm_obj10)[2]*2.303), NA)
      trend_df$pc1_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                                 100*(exp(coef(lm_obj10)[2]*2.303)-1), NA)
      trend_df$pcn_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                                 100*(exp(coef(lm_obj10)[2]*2.303)^(
                                   9)-1), NA)
    }
    
    if(!is.null(bycols))
      for(i in 1:length(bycols))
        trend_df[, bycols[i]] <- collind_df[1, get(bycols[i])]
    
    return(trend_df)
  } else {
    return(NULL)
  }
}

# Underlying function to estimate confidence interval for trends from bootstraps
trend_ci_func <- function(boot_trends1){
  trend_ci <- boot_trends1[boot_trends1$BOOTi == 0,]
  if(nrow(trend_ci) > 0){
    boot_trends <- boot_trends1[boot_trends1$BOOTi > 0,]
    trend_ci$nboot_lt <- length(boot_trends[!is.na(boot_trends$rate_lt),]$rate_lt)
    trend_ci$rate_lt_low <- quantile(boot_trends$rate_lt, 0.025, na.rm = TRUE)
    trend_ci$rate_lt_upp <- quantile(boot_trends$rate_lt, 0.975, na.rm = TRUE)
    trend_ci$pc1_lt_low <- quantile(boot_trends$pc1_lt, 0.025, na.rm = TRUE)
    trend_ci$pc1_lt_upp <- quantile(boot_trends$pc1_lt, 0.975, na.rm = TRUE) 
    trend_ci$pcn_lt_low <- quantile(boot_trends$pcn_lt, 0.025, na.rm = TRUE)
    trend_ci$pcn_lt_upp <- quantile(boot_trends$pcn_lt, 0.975, na.rm = TRUE)
    if("rate_10y" %in% colnames(boot_trends)){
      trend_ci$nboot_10y <- length(boot_trends[!is.na(boot_trends$rate_10y),]$rate_10y)
      trend_ci$rate_10y_low <- quantile(boot_trends$rate_10y, 0.025, na.rm = TRUE)
      trend_ci$rate_10y_upp <- quantile(boot_trends$rate_10y, 0.975, na.rm = TRUE)
      trend_ci$pc1_10y_low <- quantile(boot_trends$pc1_10y, 0.025, na.rm = TRUE)
      trend_ci$pc1_10y_upp <- quantile(boot_trends$pc1_10y, 0.975, na.rm = TRUE) 
      trend_ci$pcn_10y_low <- quantile(boot_trends$pcn_10y, 0.025, na.rm = TRUE)
      trend_ci$pcn_10y_upp <- quantile(boot_trends$pcn_10y, 0.975, na.rm = TRUE)
    }
    trend_ci$BOOTi <- NULL
    return(trend_ci)
  } else {
    return(NULL)
  }
}

# Function to classify trends (significance)
trend_class_func <- function(trend_ci, tcol = "TrendClass", 
                             low = "rate_lt_low", upp = "rate_lt_upp"){
  
  trend_ci[,tcol] <- NA
  if(sum(trend_ci[,low] > 1) > 0)
    trend_ci[trend_ci[,low] > 1, tcol] <- "Moderate increase"
  if(sum(trend_ci[,low] > 1.05) > 0)
    trend_ci[trend_ci[,low] > 1.05, tcol] <- "Strong increase"
  if(sum(trend_ci[,upp] < 1) > 0)
    trend_ci[trend_ci[,upp] < 1, tcol] <- "Moderate decline"
  if(sum(trend_ci[,upp]< 0.95) > 0)
    trend_ci[trend_ci[,upp] < 0.95, tcol] <- "Strong decline"
  if(sum(is.na(trend_ci[,tcol]) & 
         (trend_ci[,upp] > 1.05 | trend_ci[,low] < 0.95)) > 0)
    trend_ci[is.na(trend_ci[,tcol]) & 
               (trend_ci[,upp] > 1.05 | trend_ci[,low] < 0.95), tcol] <- "Uncertain"
  if(sum(is.na(trend_ci[,tcol])) > 0)
    trend_ci[is.na(trend_ci[,tcol]), tcol] <- "Stable"
  
  trend_ci[,tcol] <- factor(trend_ci[,tcol], 
                            levels = c("Strong decline", 
                                       "Moderate decline",
                                       "Stable",
                                       "Uncertain",
                                       "Moderate increase",
                                       "Strong increase"))
  return(trend_ci)
}




# Produce indicators
produce_indicator0 <- function(collind_region0, interval = "decreasing"){
  
  if(interval != "decreasing")
    warning("Interval not defined as decreasing - will increase in width over time")
  
  # Rescale collated indices from log10 scale
  collind_region0$TRMOBS100 <- 10^collind_region0$TRMOBS
  # Change year column name
  colnames(collind_region0)[colnames(collind_region0) == "M_YEAR"] <- "year"
  
  # Calculate indicator for real data
  indicator0 <- indicator_func(reshape2::dcast(collind_region0, year ~ SPECIES,
                                               value.var = "TRMOBS100"))[,c('year','indicator')]
  
  # Rescale so last year is 100 (for an interval of decreasing width over time)
  if(interval == "decreasing")
    indicator0[,"indicator"] <- indicator0[,"indicator"]/indicator0[nrow(indicator0),"indicator"]*100
  
  # Fit LOESS to get smoothed indicators
  ind_gam <-  predict(loess(indicator0[,2] ~ indicator0[,1], 
                            span = 0.75, degree = 2, 
                            na.action = na.exclude),
                      se = FALSE)
  
  # Rescale such that the smoothed indicator starts at 100 
  msi <- data.frame(indicator0)
  msi$indicator <- msi$indicator/ind_gam[1]*100
  msi$ind_gam0 <- ind_gam
  msi$SMOOTH <- ind_gam/ind_gam[1]*100
  msi <- merge(msi, 
               collind_region0[, .(NSPECIES = uniqueN(SPECIES)), by = "year"],
               by = "year") 
  
  return(msi)
}               

produce_indicators_boot <- function(collind_region_boot,
                                    interval = "decreasing", 
                                    inparallel = FALSE, ncpus = NULL){
  setDT(collind_region_boot)
  if(inparallel & is.null(ncpus)){
    # Determine number of cpus to use if not specified
    n_l_cores = detectCores()
    n_p_cores = detectCores(logical = FALSE)
    ncpus = n_l_cores - (n_l_cores/n_p_cores)
  } 
  
  if(interval != "decreasing")
    warning("Interval not defined as decreasing - will increase in width over time")
  
  # Rescale collated indices from log10 scale
  collind_region_boot[, TRMOBS100 := 10^TRMOBS]
  # Change year column name
  setnames(collind_region_boot, "M_YEAR", "year")
  
  # Calculate indicator for each bootstrap (if present in input data)
  if(nrow(collind_region_boot) > 0){
    if(!inparallel){
      indicators_boot <- do.call(rbind,
                                 plyr::dlply(collind_region_boot, "BOOTi",
                                             function(x){
                                               indicator_func(reshape2::dcast(x, year ~ SPECIES, 
                                                                              value.var = "TRMOBS100"))[,"indicator"]}))
    } else {
      ind_boot_func <- function(iboot){
        y <- collind_region_boot[BOOTi == iboot]
        indicator_func(reshape2::dcast(y, year ~ SPECIES, 
                                       value.var = "TRMOBS100"))[,"indicator"]}
      
      cl <- parallel::makeCluster(ncpus)
      clusterEvalQ(cl, library(data.table))
      clusterExport(cl, varlist = c("collind_region_boot","indicator_func"),
                    envir = environment())
      indicators_boot <- do.call(rbind, parLapply(cl, unique(collind_region_boot$BOOTi),
                                                  ind_boot_func))
      stopCluster(cl)
    }
    
    # Rescale so last year is 100
    if(interval == "decreasing")
      indicators_boot <- indicators_boot/indicators_boot[,ncol(indicators_boot)]*100
    
    
    indicators_gam <- apply(indicators_boot, 1,
                            function(x){z = data.frame(ind = x,
                                                       y = 1:ncol(indicators_boot));
                            predict(loess(ind ~ y, 
                                          span = 0.75, degree = 2, 
                                          na.action = na.exclude, data = z), 
                                    se = FALSE)
                            })
    
    # Rescale based on smoothed indicator starting at 100
    #msi_eu$LOW <- apply(indicators_boot/ind_gam[1]*100, 2, quantile, 0.025)
    #msi_eu$UPP <- apply(indicators_boot/ind_gam[1]*100, 2, quantile, 0.975)
    return(indicators_gam)
  } else {
    return(NULL)
  }
  
  
}     

# Calculate indicator CI from bootstraps
add_indicator_CI <- function(msi0,
                             msi_boot){
  # Rescale based on smoothed indicator starting at 100
  msi0$LOWsmooth1 <- apply(msi_boot/msi0$ind_gam0[1]*100, 
                           1, quantile, 0.025)
  msi0$UPPsmooth1 <- apply(msi_boot/msi0$ind_gam0[1]*100, 
                           1, quantile, 0.975)
  
  return(msi0)
}



# This is basically the rescale_species function from the BRCindicators package
indicator_func <- function(Data, index = 100, max = 10000, min = 1){ 
  geomean <- function(x) exp(mean(log(x), na.rm = T))
  Data <- data.matrix(Data)
  multipliers <- index/Data[1, 2:ncol(Data)]
  indicator_scaled <- t(t(Data[, 2:ncol(Data)]) * multipliers)
  indicator_scaled[indicator_scaled < min & !is.na(indicator_scaled)] <- min
  indicator_scaled[indicator_scaled > max & !is.na(indicator_scaled)] <- max
  geomean_vals <- apply(X = indicator_scaled, MARGIN = 1, FUN = geomean)
  indicator_scaled <- cbind(indicator_scaled, geomean_vals)
  colnames(indicator_scaled)[ncol(indicator_scaled)] <- "geomean"
  NAtop <- colnames(Data)[is.na(Data[1, ])]
  firstYear <- function(x) min(which(!is.na(x)))
  if (length(NAtop) > 1) 
    NAtop <- names(sort(apply(X = Data[, NAtop], MARGIN = 2, 
                              FUN = firstYear)))
  if (length(NAtop) > 0) {
    for (i in 1:length(NAtop)) {
      temp_col <- data.frame(species = !is.na(Data[, NAtop[i]]), 
                             row = 1:nrow(indicator_scaled))
      first_year <- min(temp_col$row[temp_col$species])
      temp_gm <- indicator_scaled[first_year, "geomean"]
      multi <- temp_gm/Data[first_year, NAtop[i]]
      d <- Data[, NAtop[i]] * multi
      d[d < min & !is.na(d)] <- min
      d[d > max & !is.na(d)] <- max
      indicator_scaled[, NAtop[i]] <- d
      indicator_scaled[, "geomean"] <- apply(X = indicator_scaled[, 
                                                                  !colnames(indicator_scaled) %in% "geomean"], 
                                             MARGIN = 1, FUN = geomean)
    }
  }
  fillTailNAs <- function(x) {
    na_true_false <- is.na(x)
    na_position <- grep(FALSE, na_true_false)
    if (!max(na_position) == length(x)) {
      x[(max(na_position) + 1):length(x)] <- x[max(na_position)]
    }
    return(x)
  }
  temp_indicator_scaled <- try(apply(X = indicator_scaled[, -ncol(indicator_scaled)], 
                                     MARGIN = 2, FUN = fillTailNAs), silent=TRUE)
  if(!inherits(temp_indicator_scaled, "try-error")){
    indicator_scaled <- cbind(temp_indicator_scaled, apply(X = temp_indicator_scaled, 
                                                           MARGIN = 1, FUN = geomean))
    colnames(indicator_scaled)[ncol(indicator_scaled)] <- "indicator"
    indicator_scaled <- cbind(Data[, "year"], indicator_scaled)
    colnames(indicator_scaled)[1] <- "year"} else {
      indicator_scaled = NULL
    }
  return(indicator_scaled)
}


# Indicator trends ----
# Estimate indicator trend including significance from bootstraps
estimate_ind_trends <- function(msi0, msi_boot){
  
  setDT(msi0)
  maxyear <- max(msi0$year)
  minyear <- min(msi0$year)
  
  msi0_10 <- msi0[year %in%  tail(msi0$year,10)]  
  
  # Fit linear model to smoothed indicator
  lm_obj <- try(lm(log(SMOOTH) ~ year, msi0), silent = TRUE)
  
  
  msi_trend <- data.frame(rate_lt = ifelse(!inherits(lm_obj, "try-error"), 
                                           exp(coef(lm_obj)[2]), NA),
                          pc1_lt = ifelse(!inherits(lm_obj, "try-error"),
                                          100*(exp(coef(lm_obj)[2])-1), NA),
                          pcn_lt = ifelse(!inherits(lm_obj, "try-error"), 
                                          100*(exp(coef(lm_obj)[2])^(
                                            maxyear-minyear)-1), NA))
  
  
  if(length(minyear:maxyear) > 10){
    lm_obj10 <- try(lm(log(SMOOTH) ~ year, msi0_10), silent = TRUE)
    
    msi_trend$rate_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                      exp(coef(lm_obj10)[2]), NA)
    msi_trend$pc1_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                     100*(exp(coef(lm_obj10)[2])-1), NA)
    msi_trend$pcn_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                     100*(exp(coef(lm_obj10)[2])^(
                       9)-1), NA)
    
  }
  
  
  msi_boot_trends <-  do.call(rbind, 
                              apply(msi_boot, 2,
                                    indicator_trend_func,
                                    msi0))
  
  msi_boot_trends$BOOTi <- seq_len(nrow(msi_boot_trends))
  
  # Combine real and boot trends
  msi_trend$BOOTi <- 0
  msi_boot_trends <- rbind(msi_trend, msi_boot_trends)
  
  # Summarise bootstrap trends to get a CI
  msi_trends_ci <- trend_ci_func(msi_boot_trends) 
  
  # Classify long term trends
  msi_trends_ci <- trend_class_func(msi_trends_ci,
                                    tcol = "TrendClass_lt")
  # Classify 10 years trends
  if("rate_10y" %in% colnames(msi_trends_ci))
    msi_trends_ci <- trend_class_func(msi_trends_ci, 
                                    tcol = "TrendClass_10y", 
                                    low = "rate_10y_low", 
                                    upp = "rate_10y_upp")
  
  msi_trends_ci$minyear <- minyear
  msi_trends_ci$maxyear <- maxyear
  
  msi_trends_ci <- msi_trends_ci[,
                                 c( "minyear", "maxyear", "nboot_lt",
                                    "rate_lt", "rate_lt_low", "rate_lt_upp",
                                    "pcn_lt", "pcn_lt_low", "pcn_lt_upp",
                                    "TrendClass_lt",
                                    if("nboot_10y" %in% colnames(msi_trends_ci)){
                                      c("nboot_10y", "rate_10y", "rate_10y_low", "rate_10y_upp",
                                    "pcn_10y", "pcn_10y_low", "pcn_10y_upp", "TrendClass_10y")})]
  
  return(msi_trends_ci)
}

# Trend function to apply to each bootstrapped indicator
indicator_trend_func <- function(msi_boot1, msi0){
  
  maxyear <- max(msi0$year)
  minyear <- min(msi0$year)
  
  # Fit linear model
  lm_obj <- try(lm(log(msi_boot1) ~ msi0$year), silent = TRUE)
  
  
  trend_df <- data.frame(rate_lt = ifelse(!inherits(lm_obj, "try-error"), 
                                          exp(coef(lm_obj)[2]), NA),
                         pc1_lt = ifelse(!inherits(lm_obj, "try-error"),
                                         100*(exp(coef(lm_obj)[2])-1), NA),
                         pcn_lt = ifelse(!inherits(lm_obj, "try-error"), 
                                         100*(exp(coef(lm_obj)[2])^(
                                           maxyear-minyear)-1), NA))
  
  if(length(minyear:maxyear) > 10){
    # Create matrix for last ten years too
    msi_boot10 <- tail(msi_boot1,10)
    lm_obj10 <- try(lm(log(msi_boot10) ~ tail(msi0$year,10)), silent = TRUE)
    
    trend_df$rate_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                     exp(coef(lm_obj10)[2]), NA)
    trend_df$pc1_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                    100*(exp(coef(lm_obj10)[2])-1), NA)
    trend_df$ pcn_10y <- ifelse(!inherits(lm_obj10, "try-error"), 
                    100*(exp(coef(lm_obj10)[2])^(
                      9)-1), NA)
  }
  
  return(trend_df)
}


