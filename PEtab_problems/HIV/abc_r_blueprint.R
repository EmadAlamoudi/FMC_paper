# New abc_r_blueprint.R, where we store all functions as the atm we cant work
# with relative paths for automatic job submission.

# Load Libraries
library(MotilityLab)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape))
library(stringi)


# Set parameters
celltypes <- c("target", "infected")
time_End   <- 3900
time_Shift <- 300
use_dist <- "normal" # "ks" or "normal"
dist_to_index <- list("normal" = 1,
                      "ks"     = 2)



setParameters<-function(simulation=TRUE,   
                        samplingStep = 30, # "measure" every 30 seconds 
                        timeEnd = 3720,    # Simulation end
                        timeShift = 120){  # Adaption phase 
  # Returns scale and a vector including all timepoints
  # Timepoints are rounded to 5 digits bc of num instability
  if(simulation)
  {
    time_scale  <- 1/60   # To convert to minutes
    time_vector <- seq(time_scale*timeShift,
                       timeEnd*time_scale,
                       samplingStep*time_scale) %>% round(5)
  } else{
    time_scale  <- 0.5        
    time_vector <- NULL
  }
  return(list("t_scale" = time_scale,
              "t_vector" = time_vector))
}


myReadData<-function(file_name, timeScale = 1, dim = 2){
  # Uses the read.tracks.csv function from motility lab
  # timeScale should be inherited from setParameters
  assertthat::assert_that(dim == 2 | dim == 3, 
                          msg = "Dimension has to be 2 or 3")
  if(dim == 2){
    DF <- read.tracks.csv(file_name,  
                          sep = "\t", 
                          id.column = "cell.id", 
                          time.column = "time", 
                          pos.columns = c(3,4), # c(3,NA) would take all 
                          # colums from 3rd to last
                          header = TRUE, 
                          scale.t = timeScale)}
  if(dim == 3){
    DF <- read.tracks.csv(file_name,  
                          sep = "\t", 
                          id.column = "cell.id", 
                          time.column = "time", 
                          pos.columns = c(3,4,5), # c(3,NA) would take all 
                          # colums from 3rd to last
                          header = TRUE, 
                          scale.t = timeScale)
  }
  
  return(DF)
}


myFilterTracks<-function(trackData.orig, 
                         timeVector = 0, 
                         timeThreshold = 0){
  # Filter timeThreshold data points (is to be given in the correct timescale)
  # or filter for all the given timeVector points
  trackData <- trackData.orig
  
  # 1. take out the first timeThreshold of simulation 
  #(time thought for the sim to reach steady-state)
  if(timeThreshold>0){
    trackData<-as.tracks(lapply(trackData, 
                                function(x) x[x[,1] >= timeThreshold, ]))  
  }
  
  # 2. extract the relevant data points according to the timeVector
  if (length(timeVector) > 1){
    trackData <- as.tracks(lapply(trackData, 
                                  function(x) x[round(x[,1], 5) %in% timeVector, ]))
    # Assert if we lost some timepoints due to numerics
    assertthat::assert_that((trackData[[1]] %>% dim)[1] == (length(timeVector)), 
                            msg = "Filtered data lost some data along the way")
  }
  
  return(trackData)
}


myMotilityAnalysis<-function(trackData, velocityThreshold = 2){
  # returns a list with two data frames: the first one containing the MSD analysis 
  # and the second one containing speed, straightness, turning angles, and arrest coefficient
  

  msd            <- aggregate(trackData, squareDisplacement, FUN = "mean.se")
  inst.velocity  <- as.data.frame(sapply(trackData,speed))
  straight.ness  <- as.data.frame(sapply(trackData, straightness))
  turning.angles <- as.data.frame(sapply(trackData,meanTurningAngle)) %>%
                        mutate_all(coalesce, 0) # nan to 0

  arrest <- sapply(trackData, function(x){tmp  <- sapply(subtracks(x,1),speed)
    if(length(tmp) > 0){  
      arrest <- length(tmp[tmp < velocityThreshold])/length(tmp)
    }  else{
      arrest <- NA
    }
  })
  
  out_else <- data.frame(inst.velocity, 
                         straight.ness, 
                         turning.angles, 
                         arrest)
  
  colnames(out_else) <- c("speed", 
                          "straight", 
                          "angles", 
                          "arrest")
  
  return(list("msd"      = msd,
              "motility" = out_else))
}


check_end_slash <- function(path){
  # Checks that the path ends with a slash. If not a "/" is appended
  #Check that directory ends with "/"
  print(paste0("path_R= ",path))
  if(strsplit(path,"")[[1]][length(strsplit(path,"")[[1]])] == "/"){
    path = path
  } else{
    path = paste0(path, "/")
  }
  return(path)
}


return_logger <- function(celltype){
  # Given the celltype as string it returns a given logger name which is 
  # predefined in this function.
  log_list = list(target = "target_tracker.csv", 
                  infected = "infected_tracker.csv")
  assertthat::assert_that(any(names(log_list) == celltype), 
                          msg = "Celltype has no associated logger file.
                          Please check return_logger function.")
  return(log_list[celltype][[1]])
}


return_names <- function(celltype, sim){
  # Return names for the summaryStatistic List given a celltype(str)
  # If Sim == True, Sim is appended else it's Exp
  if(sim == TRUE){
    mot <- paste0("sumMot", celltype, "Sim")
    msd <- paste0("sumMSD", celltype, "Sim")
  } else {
    mot <- paste0("sumMot", celltype, "Exp")
    msd <- paste0("sumMSD", celltype, "Exp")
  }
  return(c(msd, mot))
}


motility_wrapper <- function(celltype,  # Celltype(str)
                             directory, # directory of the logger file 
                             timeVector,    # passed by setParameters
                             timeThreshold, # adaption phase
                             timeScale,     # passed by setParameters
                             sim = TRUE){
  file_ <- return_logger(celltype = celltype)
  in.file <- paste(directory, file_, sep = "")
  
  tracksData <- myReadData(file_name = in.file,   # Read data
                           timeScale = timeScale) 
  tracksData.filtered <- myFilterTracks(trackData.orig = tracksData, # Filter
                                        timeVector = timeVector,     # adaption
                                        timeThreshold = timeThreshold)  # phase
  tracksData.analyzed <- myMotilityAnalysis(trackData = tracksData.filtered)
  
  # Split output
  msdAnalysis         <- data.frame(tracksData.analyzed[[1]])
  otherMotilityParams <- data.frame(tracksData.analyzed[[2]])
  
  # Fuse to named sumStat
  sumStat        <- list(msdAnalysis, otherMotilityParams)
  names(sumStat) <- return_names(celltype, sim)
  
  return(sumStat)
}


mySummaryStatistics <- function(directory, 
                                celltype = celltypes, 
                                timeEnd=time_End, 
                                timeShift=time_Shift, 
                                samplingStep=30, 
                                simulation=TRUE){
  #' Calculates summary statistics from cell position logger (so far 2D)
  #' 
  #' It is important that the summary statistics have names to store it correctly 
  #' in pyABC's database
  #' @param celltype vector or char that has all the celltypes to be analysed
  #' @param directory Directory where to find logger files
  #' @param timeEnd Until which timepoint to evaluate
  #' @param timeShift Length of adaption phase not evaluated
  #' @param samplingStep Interval length of sampling.
  #' @simulation Data from simulation
  #' @return Named list of summary statistics
  
  Sim <- list()
  
  # timeScale (in hours): if the data is coming from the experiments then each 
  # time point is 0.5; 
  # if the data is from the simulations then 1/60 (assuming each MCS is one sec) 
  # timeVector: vector with the time points measrued in the simulations, 
  # i.e. one measurement every 30 s
  time <- setParameters(samplingStep = samplingStep, 
                        timeEnd      = timeEnd, 
                        timeShift    = timeShift,
                        simulation   = simulation)
  timeScale  <- time[[1]] 
  timeVector <- time[[2]]
  
  adaption_phase <- timeShift*timeScale
  directory <- check_end_slash(directory[['loc']])
  # Create sumStat for all celltypes
  for(t in celltype){
    output <- motility_wrapper(celltype  = t, 
                               directory = directory, 
                               timeVector    = timeVector, 
                               timeThreshold = adaption_phase,
                               timeScale = timeScale)
    
    Sim <- append(Sim, output)
  }
  return(Sim)
}


distance_normal <- function(Sim, Exp){
  # Calculates a distance measure between experimental data and simulation,
  # given the mean + the sigma of experimental/simulated motility distribution.
  # To use this distance measure, we assume the distributions to be normal.
  mean_exp <- Exp %>% mean
  mean_sim <- Sim %>% mean
  sig_exp <- Exp %>% sd
  sig_sim <- Sim %>% sd
  
  dist <- (mean_exp - mean_sim)**2/sig_exp**2 + 
    ((sig_exp/mean_exp)-(sig_sim/mean_exp))**2
  return(dist)
}


distance_kolmogorov <- function(Sim, Exp){
  dist <- ks.test(Sim, Exp)[["statistic"]] %>% unname()
  return(dist)
}

distance_list <- list("normal" = distance_normal,
                      "ks"     = distance_kolmogorov)

dist_fun <- distance_list[[dist_to_index[[use_dist]]]]

distance_speed_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["speed"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["speed"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  print(paste('Spd_tar', dist)) 
 return(dist)
}

distance_speed_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["speed"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["speed"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  print(paste('Spd_inf', dist)) 
 return(dist)
}


distance_str_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["straight"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["straight"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  print(paste('Str_tar', dist)) 
 return(dist)
}

distance_str_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["straight"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["straight"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  print(paste('Str_inf', dist))
 
 return(dist)
}


distance_trn_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]]
  sumMotSim <- Sim[[names_sim[[2]]]]
  
  # Cells that are not moving at all == nan!
  nan_count <- sum(is.nan(sumMotSim[["angles"]]))
  sumMotSim_trn <- sumMotSim %>%
    select(angles) %>%
    filter(!is.nan(angles))
  
  if(dim(sumMotSim_trn)[1] > 0){
    dist <- dist_fn(Sim = sumMotSim_trn[["angles"]],
                    Exp = sumMotExp[["angles"]])
  } else {
    dist <- 0  # Initialize dist
  }
  count_NaN <- sumMotSim %>%
    select(angles) %>%
    filter(is.nan(angles)) %>%
    count() %>% pull()
  fraction_NaN <- count_NaN/dim(sumMotSim)[1]
  dist <- dist + fraction_NaN*pi  # Fraction of NaN add maximal distance
   print(paste('Trn_tar', dist))
  return(dist)
}

distance_trn_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]]
  sumMotSim <- Sim[[names_sim[[2]]]]
  
  # Cells that are not moving at all == nan!
  nan_count <- sum(is.nan(sumMotSim[["angles"]]))
  sumMotSim_trn <- sumMotSim %>%
    select(angles) %>%
    filter(!is.nan(angles))
  
  if(dim(sumMotSim_trn)[1] > 0){
    dist <- dist_fn(Sim = sumMotSim_trn[["angles"]],
                    Exp = sumMotExp[["angles"]])
  } else {
    dist <- 0  # Initialize dist
  }
  count_NaN <- sumMotSim %>%
    select(angles) %>%
    filter(is.nan(angles)) %>%
    count() %>% pull()
  fraction_NaN <- count_NaN/dim(sumMotSim)[1]
  dist <- dist + fraction_NaN*pi  # Fraction of NaN add maximal distance
  print(paste('Trn_inf', dist))
  return(dist)
}



distance_arr_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["arrest"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["arrest"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
   print(paste('Arr_tar', dist))
  return(dist)
}


distance_arr_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["arrest"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["arrest"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  print(paste('Arr_inf', dist))
  return(dist)
}


distance_msd_tar <- function(Sim, Exp, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMSDExp <- Exp[[names_exp[[1]]]]
  sumMSDSim <- Sim[[names_sim[[1]]]]
  
  mn_e <- sumMSDExp[["mean"]]
  mn_s <- sumMSDSim[["mean"]]
  sd_e <- sumMSDExp[["mean"]] - sumMSDExp[["lower"]]
  sd_s <- sumMSDSim[["mean"]] - sumMSDSim[["lower"]]
  
  dist <- ((mn_e-mn_s)/sd_e)**2 %>% mean
   print(paste('MSD_tar', dist))
  return(dist*0.01)
}


distance_msd_inf <- function(Sim, Exp, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMSDExp <- Exp[[names_exp[[1]]]]
  sumMSDSim <- Sim[[names_sim[[1]]]]
  
  mn_e <- sumMSDExp[["mean"]]
  mn_s <- sumMSDSim[["mean"]]
  sd_e <- sumMSDExp[["mean"]] - sumMSDExp[["lower"]]
  sd_s <- sumMSDSim[["mean"]] - sumMSDSim[["lower"]]
  
  dist <- ((mn_e-mn_s)/sd_e)**2 %>% mean
  print(paste('MSD_inf', dist))
  return(dist*0.01)
}


home1 <- paste0(getwd(),"/")
home <- paste0(dirname(dirname(home1)),"/HIV/")
print(home)
sumMSDinfectedExp <- read.table(paste0(home, "loose_MSD_infected.csv"),
                                sep="\t", header=T)
sumMotinfectedExp <- read.table(paste0(home, "loose_Mot_infected.csv"),
                                sep="\t", header=T)
sumMSDtargetExp <- read.table(paste0(home, "loose_MSD_target.csv"),
                                sep="\t", header=T)
sumMottargetExp <- read.table(paste0(home, "loose_Mot_target.csv"),
                                sep="\t", header=T)
Exp <- list(sumMSDtargetExp=sumMSDtargetExp, sumMottargetExp=sumMottargetExp,
            sumMSDinfectedExp=sumMSDinfectedExp, sumMotinfectedExp=sumMotinfectedExp)
