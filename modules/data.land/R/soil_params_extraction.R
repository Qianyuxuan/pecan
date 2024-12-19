rm(list = ls())
library(dplyr)
library(PEcAn.all)
outdir <- "/projectnb/dietzelab/Cherry/SoilGrids_texture/39NEON"
sand <- readRDS("/projectnb/dietzelab/Cherry/SoilGrids_texture/sand_percent.rds") #sand fraction in percentage
clay <- readRDS("/projectnb/dietzelab/Cherry/SoilGrids_texture/clay_percent.rds") #clay fraction in percentage
silt <- readRDS("/projectnb/dietzelab/Cherry/SoilGrids_texture/silt_percent.rds") #silt fraction in percentage

settings <-
  read.settings("/projectnb/dietzelab/Cherry/xml/pecan_monthly_SDA_soilwater.xml")



##' soil_params_extraction function
##' A function to estimate the soil parameters based on soil texture data and write the parameter paths into settings
##' @title soil_params_extraction
##' @name soil_params_extraction
##' 
##' param settings A multi-site settings
##' param sand A data frame containing sand fraction in percentage from SoilGrids250m v2.0 with columns "Depth", "Quantile", "Siteid", and "Value"
##' param clay A data frame containing clay fraction in percentage from SoilGrids250m v2.0 with columns "Depth", "Quantile", "Siteid", and "Value"
##' param silt A data frame containing silt fraction in percentage from SoilGrids250m v2.0 with columns "Depth", "Quantile", "Siteid", and "Value"
##' param outdir Provide the path to store the parameter file
##' param write_into_settings Whether to write the path of parameter file into the setting. The default is TRUE
##' 
##' @example:
##' \dontrun{
##' 
##' outdir <- "/projectnb/dietzelab/Cherry/SoilGrids_texture/39NEON"
##' sand <- readRDS("/projectnb/dietzelab/Cherry/SoilGrids_texture/sand_percent.rds") #sand fraction in percentage
##' clay <- readRDS("/projectnb/dietzelab/Cherry/SoilGrids_texture/clay_percent.rds") #clay fraction in percentage
##' silt <- readRDS("/projectnb/dietzelab/Cherry/SoilGrids_texture/silt_percent.rds") #silt fraction in percentage
##' settings <-read.settings("/projectnb/dietzelab/Cherry/xml/pecan_monthly_SDA_soilwater.xml")
##' soil_params_extraction(settings,sand,clay,silt,outdir)
##'  }
##'  
##'  @return Ensemble soil parameter files defined in outdir and file paths in xml file
##'  @export 
##'  @author Qianyu Li
##'  @importFrom magrittr %>%
##'  
############## Fit the Dirichlet distribution parameter to approximate the observed quantiles with means as known moments
##### Example for means and quantiles:
### means <- c(0.383,0.213,0.403)
### quantiles <-list(
### q5 = c(0.551,0.055,0.069), # 5th percentile for each category: sand, clay, and silt at one site and one depth
### q50 = c(0.691,0.124,0.156), # 50th percentile (median) for each category: sand, clay, and silt at one site and one depth
### q95 = c(0.786,0.258,0.297))  # 95th percentile for each category: sand, clay, and silt at one site and one depth

# Step 1: Estimate alpha_0
estimate_alpha0 <- function(means, quantiles) {
  # Objective function to minimize the difference between observed and simulated quantiles with means as a known moment
  objective_function <- function(alpha0) {
    if (alpha0 <= 0)
      return(Inf)
    alpha <- means * alpha0
    q5 <- quantiles$q5
    q50 <- quantiles$q50
    q95 <- quantiles$q95
    
    # Compute differences with observed quantiles
    samples <- MCMCpack::rdirichlet(10000, alpha) # Generate samples
    sim_q5 <- apply(samples, 2, quantile, probs = 0.05, na.rm = TRUE)
    sim_q50 <-apply(samples, 2, quantile, probs = 0.50, na.rm = TRUE)
    sim_q95 <-apply(samples, 2, quantile, probs = 0.95, na.rm = TRUE)
    quantile_diff <-sum((sim_q5 - q5) ^ 2 + (sim_q50 - q50) ^ 2 + (sim_q95 - q95) ^ 2)
    return(quantile_diff)
  }
  
  # Optimize alpha_0
  result <- optim(
    par = 1, # Initial guess for alpha_0
    fn = objective_function,
    method = "L-BFGS-B",
    lower = 0.01  # alpha_0 must be positive
  )
  return(result$par)
}

# Step 2: Estimate individual alphas
estimate_dirichlet_parameters <- function(means, quantiles) {
 
  alpha0 <- estimate_alpha0(means, quantiles)
  if (alpha0 <= 0) {
    stop("Estimated alpha0 is non-positive, which is invalid.")
  }
  alphas <- means * alpha0
  
  return(alphas)
}


# Rescale the sum of mean texture fraction to 1 as the original sums are slightly different from 1 for some layers
rescale_sum_to_one <- function(sand, clay, silt) {
  total <- sand + clay + silt
  rescaled_sand <- sand / total
  rescaled_clay <- clay / total
  rescaled_silt <- silt / total
  return(list(
    sand = rescaled_sand,
    clay = rescaled_clay,
    silt = rescaled_silt
  ))
}

# Function to reformat the nested list as inputs to "soil2netcdf" function
reformat_soil_list <- function(samples_all_depth) {
  # Define the fractions we want to extract
  fractions <-
    c(
      "fraction_of_sand_in_soil",
      "fraction_of_clay_in_soil",
      "fraction_of_silt_in_soil"
    )
  
  # Initialize a new list to store reformatted data
  reformatted <-
    setNames(vector("list", length(fractions)), fractions)
  
  # Extract data for each fraction
  for (fraction in fractions) {
    reformatted[[fraction]] <-
      unlist(lapply(samples_all_depth, function(depth_list) {
        depth_list[[fraction]] # Extract the fraction value
      })) %>% purrr::set_names(NULL)
  }
  # Combine depth into a single vector for readability
  reformatted$soil_depth <- as.numeric(names(samples_all_depth))
  return(reformatted)
}

# Write to settings
create_mult_list <- function(list.names, paths) {
  out <- as.list(paths)
  names(out) <- list.names
  out
}


soil_params_extraction <- function(settings,sand,clay,silt,outdir,write_into_settings=TRUE){

# Convert values to proportion from %
if (any(c(sand$Value, clay$Value, silt$Value) > 2)) {
  sand$Value <- if (is.null(sand$Value)) { NULL } else { sand$Value / 100 }
  clay$Value <- if (is.null(clay$Value)) { NULL } else { clay$Value / 100 }
  silt$Value <- if (is.null(silt$Value)) { NULL } else { silt$Value / 100 }
}
ens_n <- settings$ensemble$size
# Merge all soil texture data together
texture_all <-merge(sand, clay, by=c("Depth", "Quantile", "Siteid"))  %>% merge(silt, by=c("Depth", "Quantile", "Siteid")) %>%
  `colnames<-`(c(
      "soil_depth", #soil_depth will be used in "soil2netcdf" function
      "quantile",
      "siteid",
      "fraction_of_sand_in_soil",
      "fraction_of_clay_in_soil",
      "fraction_of_silt_in_soil"
    ))

# Substitute the depth range with the middle depth value in meter
texture_all$soil_depth <-
  gsub("100-200cm", 1.5, gsub("60-100cm", 0.8, gsub(
    "30-60cm", 0.45, gsub("15-30cm", 0.225, gsub(
      "5-15cm", 0.1, gsub("0-5cm", 0.025, texture_all$soil_depth))))))
texture_all$soil_depth <- as.numeric(texture_all$soil_depth)
# Reformat the list based on site id
f1 <- factor(texture_all$siteid, levels = unique(texture_all$siteid))
dat <- split(texture_all, f1)
# Grab Site IDs from settings
settings_id <-lapply(settings, function(x) as.numeric(x$run$site$id))

for (i in seq_along(dat)) {
  samples_ens <- list()
  paths <- c()
  siteid <- as.numeric(unique(dat[[i]]$siteid))
  soil_depth <- unique(dat[[i]]$soil_depth)
  str_ns <- paste0(siteid %/% 1e+09, "-", siteid %% 1e+09)
  temp_outdir <- file.path(outdir, siteid)
  dir.create(temp_outdir)
  
  for (depths in c(0.025, 0.1, 0.225, 0.45, 0.8, 1.5)) {
    # Estimate Dirichlet parameters for each depth at each site
    quantiles <- list(
      q5 = filter(dat[[i]], quantile == "0.05", soil_depth == depths) %>% dplyr::select(
        fraction_of_sand_in_soil,
        fraction_of_clay_in_soil,
        fraction_of_silt_in_soil
      ),
      # 5th percentile for each category
      q50 = filter(dat[[i]], quantile == "0.5", soil_depth == depths) %>% dplyr::select(
        fraction_of_sand_in_soil,
        fraction_of_clay_in_soil,
        fraction_of_silt_in_soil
      ),
      # 50th percentile (median) for each category
      q95 = filter(dat[[i]], quantile == "0.95", soil_depth == depths) %>% dplyr::select(
        fraction_of_sand_in_soil,
        fraction_of_clay_in_soil,
        fraction_of_silt_in_soil
      )  # 95th percentile for each category
    )
    
    # Extract the means
    means <- filter(dat[[i]], quantile == "Mean", soil_depth == depths) %>% dplyr::select(fraction_of_sand_in_soil,fraction_of_clay_in_soil,fraction_of_silt_in_soil)
    soil_rescaled <-rescale_sum_to_one(means$fraction_of_sand_in_soil,means$fraction_of_clay_in_soil,means$fraction_of_silt_in_soil)
   
    # Replace the original means with the rescaled ones
    means$fraction_of_sand_in_soil <- soil_rescaled$sand
    means$fraction_of_clay_in_soil <- soil_rescaled$clay
    means$fraction_of_silt_in_soil <- soil_rescaled$silt
    
    # Estimate Dirichlet parameters
    alpha_est <- estimate_dirichlet_parameters(as.matrix(means), quantiles)

    # Generate the ensemble soil texture data based on the ensemble size (ens_n) defined in the settings
    samples <- MCMCpack::rdirichlet(10000, alpha_est)
    colnames(samples) <-c("fraction_of_sand_in_soil","fraction_of_clay_in_soil","fraction_of_silt_in_soil")
    samples_ens_new <-list(samples[sample(1:10000, ens_n), ]) %>% setNames(depths)
    samples_ens <- append(samples_ens, samples_ens_new)
  }
  
  # Generate soil parameter file for each one in ensemble soil texture data
  for (ens in 1:ens_n) {
    samples_all_depth <- lapply(samples_ens, function(x) x[ens, ])
    # Reformat the nested list as inputs to "soil2netcdf" function
    reformatted_soil_list <- reformat_soil_list(samples_all_depth)
    prefix <- paste0("Soil_params_", str_ns, "_", ens)
    new.file <-  file.path(outdir, siteid, paste0(prefix, ".nc"))
    out.ense <- soil2netcdf(reformatted_soil_list, new.file)
    paths <- c(new.file, paths)
  }
  # Write the parameter paths to settings
    if (write_into_settings) {
      ind <- which(settings_id == siteid)
      settings[[ind]]$run$inputs$soilinitcond$source <- "SoilGrids"
      settings[[ind]]$run$inputs$soilinitcond$output <- "soilinitcond"
      settings[[ind]]$run$inputs$soilinitcond$ensemble <- ens_n
      settings[[ind]]$run$inputs$soilinitcond$path <-create_mult_list(rep("path", ens_n), paths)
      write.settings(settings,outputdir = settings$outdir,outputfile = "pecan.xml")
    }
  }
}



# Generate SoilWHC for whole soil profile in write.config
if (length(settings$run$inputs$soilinitcond$path) > 0) {
  soil_ICs_num <- length(settings$run$inputs$soilinitcond$path)
  soil_ICs_path <- settings$run$inputs$soilinitcond$path[[sample(1:soil_ICs_num, 1)]]
  soil_IC_list <- PEcAn.data.land::pool_ic_netcdf2list(soil_ICs_path)
  #SoilWHC
  if ("volume_fraction_of_water_in_soil_at_saturation" %in% names(soil_IC_list$vals)) {
    # Calculate the soilWHC for each soil layer in cm
    soilWHC_0.025 <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[which(soil_IC_list$dims$depth == 0.025)] * 5
    soilWHC_0.1 <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[which(soil_IC_list$dims$depth == 0.1)] * 10
    soilWHC_0.225 <-unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[which(soil_IC_list$dims$depth == 0.225)] * 15
    soilWHC_0.45 <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[which(soil_IC_list$dims$depth == 0.45)] * 30
    soilWHC_0.8 <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[which(soil_IC_list$dims$depth == 0.8)] * 40
    soilWHC_1.5 <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[which(soil_IC_list$dims$depth == 1.5)] * 100
    # Calculate the soilWHC for the whole soil profile in cm
    soilWHC_total <- soilWHC_0.025 + soilWHC_0.1 + soilWHC_0.225 + soilWHC_0.45 + soilWHC_0.8 + soilWHC_1.5
    param[which(param[, 1] == "soilWHC"), 2] <- soilWHC_total
  }
}
