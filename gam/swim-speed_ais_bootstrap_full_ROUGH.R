#################################
### BOOTSTRAP GAM: SPEED, AIS ###
#################################
### 25/04/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

### PACKAGES
packages <- c("tidyverse","RColorBrewer", "scales","survMisc", "Metrics", "lme4", "MASS", "mgcv", "cowplot", "parallel", "foreach", "doSNOW")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


#---------------- FUNCTIONS + THEME --------------------
source("code/functions.R")
source("code/themes.R")


#---------------- DATA --------------------

# AIS
tot <- read.csv("intermediate-products/response-var-dfs/speed_ais.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  mutate(dummy = 1) # important dummy variable to factor out random effect when predicting response

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("meandist", "dist.max")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"), # sqrt
         across(contains(c("diff.time")), .fns = list(log = ~log(.)), .names = "{fn}.{col}")) # sqrt dist max

# bootstrapped speeds
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv") %>%
  mutate(datetime = as.POSIXct(datetime))

# model
gam <- readRDS("intermediate-products/gam/swim-speed/ais/gam_swim-speed_ais_final.rds")



#---------------- BOOTSTRAP --------------------

# formula
var.fac = c("group", "seastate", "surface.feeding", "surface.active")
var.int = c("rib.num.30.1500")
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.10.1500")
var.rand = c("folnum.unique")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# we've looked at different file sizes and will go for 1000 rows for now!

#---------------- 1 to 100 --------------------

# PARALLEL FOR LOOP

# how many cores do we have? 
parallel::detectCores() 
n.cores = 28 # let's use 5 for now

# creating cluster of cores
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

# define parameters first
it = 100 # number of iterations

# create empty data frame to populate

v.line <- names(summary(gam)$p.coeff)[-1]
v.smooth <- gsub("s\\(|\\)","",names(summary(gam)$chi.sq))
speed.ais.boot <- data.frame(matrix(ncol = (length(v.line)+length(v.smooth))*2, nrow = 0)) 
names(speed.ais.boot) = c(paste0(v.line,"_estimate"), paste0(v.line,"_p"),
                          paste0(v.smooth,"_edf"), paste0(v.smooth,"_p")) 


speed.ais.boot <- foreach(i = 1:it, .combine = rbind, .packages = c("tidyverse", "mgcv"), .options.snow = opts) %dopar% {
  
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "GCV.Cp", family = gaussian(link = "log")))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  #speed.ais.boot <- rbind(speed.ais.boot,df)
  return(df)
  
  print(paste0("iteration ",i))
  
}


speed.ais.boot <- as.data.frame(speed.ais.boot)
speed.ais.boot.1to100 <- speed.ais.boot
write.csv(speed.ais.boot.1to100, "intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res_full_1to100.csv", row.names = FALSE)
parallel::stopCluster(cl = my.cluster)


#---------------- 101 to 200 --------------------

# PARALLEL FOR LOOP

# how many cores do we have? 
parallel::detectCores() 
n.cores = 28 # let's use 5 for now

# creating cluster of cores
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

# define parameters first
it = 200 # number of iterations

# create empty data frame to populate

v.line <- names(summary(gam)$p.coeff)[-1]
v.smooth <- gsub("s\\(|\\)","",names(summary(gam)$chi.sq))
speed.ais.boot <- data.frame(matrix(ncol = (length(v.line)+length(v.smooth))*2, nrow = 0)) 
names(speed.ais.boot) = c(paste0(v.line,"_estimate"), paste0(v.line,"_p"),
                          paste0(v.smooth,"_edf"), paste0(v.smooth,"_p")) 


speed.ais.boot <- foreach(i = 101:it, .combine = rbind, .packages = c("tidyverse", "mgcv"), .options.snow = opts) %dopar% {
  
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "GCV.Cp", family = gaussian(link = "log")))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  #speed.ais.boot <- rbind(speed.ais.boot,df)
  return(df)
  
  print(paste0("iteration ",i))
  
}

#speed.ais.boot.1to100 <- speed.ais.boot
speed.ais.boot <- as.data.frame(speed.ais.boot)
speed.ais.boot.101to200 <- speed.ais.boot
write.csv(speed.ais.boot.101to200, "intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res_full_101to200.csv", row.names = FALSE)
parallel::stopCluster(cl = my.cluster)


#---------------- 201 to 300 --------------------

# PARALLEL FOR LOOP

# how many cores do we have? 
parallel::detectCores() 
n.cores = 28 # let's use 5 for now

# creating cluster of cores
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

# define parameters first
it = 300 # number of iterations

# create empty data frame to populate

v.line <- names(summary(gam)$p.coeff)[-1]
v.smooth <- gsub("s\\(|\\)","",names(summary(gam)$chi.sq))
speed.ais.boot <- data.frame(matrix(ncol = (length(v.line)+length(v.smooth))*2, nrow = 0)) 
names(speed.ais.boot) = c(paste0(v.line,"_estimate"), paste0(v.line,"_p"),
                          paste0(v.smooth,"_edf"), paste0(v.smooth,"_p")) 


speed.ais.boot <- foreach(i = 201:it, .combine = rbind, .packages = c("tidyverse", "mgcv"), .options.snow = opts) %dopar% {
  
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "GCV.Cp", family = gaussian(link = "log")))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  #speed.ais.boot <- rbind(speed.ais.boot,df)
  return(df)
  
  print(paste0("iteration ",i))
  
}

speed.ais.boot <- as.data.frame(speed.ais.boot)
speed.ais.boot.201to300 <- speed.ais.boot
write.csv(speed.ais.boot.201to300, "intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res_full_201to300.csv", row.names = FALSE)
parallel::stopCluster(cl = my.cluster)




#---------------- 301 to 400 --------------------

# PARALLEL FOR LOOP

# how many cores do we have? 
parallel::detectCores() 
n.cores = 28 # let's use 5 for now

# creating cluster of cores
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

# define parameters first
it = 400 # number of iterations

# create empty data frame to populate

v.line <- names(summary(gam)$p.coeff)[-1]
v.smooth <- gsub("s\\(|\\)","",names(summary(gam)$chi.sq))
speed.ais.boot <- data.frame(matrix(ncol = (length(v.line)+length(v.smooth))*2, nrow = 0)) 
names(speed.ais.boot) = c(paste0(v.line,"_estimate"), paste0(v.line,"_p"),
                          paste0(v.smooth,"_edf"), paste0(v.smooth,"_p")) 


speed.ais.boot <- foreach(i = 301:it, .combine = rbind, .packages = c("tidyverse", "mgcv"), .options.snow = opts) %dopar% {
  
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "GCV.Cp", family = gaussian(link = "log")))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  #speed.ais.boot <- rbind(speed.ais.boot,df)
  return(df)
  
  print(paste0("iteration ",i))
  
}

speed.ais.boot <- as.data.frame(speed.ais.boot)
speed.ais.boot.301to400 <- speed.ais.boot
write.csv(speed.ais.boot.301to400, "intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res_full_301to400.csv", row.names = FALSE)
parallel::stopCluster(cl = my.cluster)


#---------------- 401 to 500 --------------------

# PARALLEL FOR LOOP

# how many cores do we have? 
parallel::detectCores() 
n.cores = 28 # let's use 5 for now

# creating cluster of cores
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

# define parameters first
it = 500 # number of iterations

# create empty data frame to populate

v.line <- names(summary(gam)$p.coeff)[-1]
v.smooth <- gsub("s\\(|\\)","",names(summary(gam)$chi.sq))
speed.ais.boot <- data.frame(matrix(ncol = (length(v.line)+length(v.smooth))*2, nrow = 0)) 
names(speed.ais.boot) = c(paste0(v.line,"_estimate"), paste0(v.line,"_p"),
                          paste0(v.smooth,"_edf"), paste0(v.smooth,"_p")) 


speed.ais.boot <- foreach(i = 401:it, .combine = rbind, .packages = c("tidyverse", "mgcv"), .options.snow = opts) %dopar% {
  
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "GCV.Cp", family = gaussian(link = "log")))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  #speed.ais.boot <- rbind(speed.ais.boot,df)
  return(df)
  
  print(paste0("iteration ",i))
  
}

speed.ais.boot <- as.data.frame(speed.ais.boot)
speed.ais.boot.401to500 <- speed.ais.boot
write.csv(speed.ais.boot.401to500, "intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res_full_401to500.csv", row.names = FALSE)
parallel::stopCluster(cl = my.cluster)
