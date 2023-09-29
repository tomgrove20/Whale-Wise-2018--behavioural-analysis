#################################
### BOOTSTRAP GAM: SPEED, AIS ###
#################################
### 25/04/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

### PACKAGES
packages <- c("tidyverse","RColorBrewer", "scales","survMisc", "Metrics", "lme4", "MASS", "mgcv", "cowplot", "parallel", "foreach")
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
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.30.1500")
var.rand = c("folnum.unique")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.int,", bs = 'cs', k = 3)", collapse = "+"),"+",
                      paste("s(",var.num,", bs = 'cs', k = 10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

# we've looked at different file sizes and will go for 1000 rows for now!
#n = 1000
#sample <- tot[sample(nrow(tot), n), ]
#system.time(g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "GCV.Cp", family = gaussian(link = "log")))) # 148 seconds. 

# PARALLEL FOR LOOP

# how many cores do we have? 
parallel::detectCores() # 8 cores!
n.cores = 5 # let's use 5 for now

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
n = 1000 # size of random sample
it = 500 # number of iterations

# and run!
speed.ais.boot <- foreach(i = 1:it, .combine = rbind, .errorhandling = 'remove', .packages = c("tidyverse", "mgcv")) %dopar% {
  
  # first take a random sample
  sample <- tot[sample(nrow(tot), n), ]
  
  # then replace speed with iteration
  sample <- sample %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
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
view(speed.ais.boot)
write.csv(speed.ais.boot, "intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res.csv")



#---------------- PLOTTING --------------------

#speed.ais.boot <- read.csv("intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot-res.csv")

# first a df for naming
var.fac = c("grouplone", "seastatechoppy", "surface.feeding1", "surface.active1")
var.int = c("rib.num.30.1500")
var.num = c("day.of.year", "log.diff.time", "sqrt.oak.meandist.10.1500", "sqrt.rib.meandist.30.1500")
name.df <- data.frame(
  var = c(var.fac, var.int, var.num),
  name = c("Group size", "Sea state", "Surface feeding", "Surface active", "# of RIBs (30 min., 1500 m)", "Julian day", "Inter-breath interval", "Oak dist. (10 min., 1500 m)", "RIB dist. (10 min., 1500 m)")) %>%
  mutate(type = ifelse(var %in% var.fac, "fac","num"))

# we need to create a df from the original GAM so we can add lines for full model effect sizes
summary <- summary(gam)
line.coeff <- as.data.frame(summary$p.coeff) %>% add_rownames(var = "var") %>% rename(val = 2)
smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
smooth.edf <- as.data.frame(summary$edf) %>% mutate(var = smooth.names) %>% rename(val = 1)
orig.df <- as.data.frame(bind_rows(line.coeff, smooth.edf)) %>% left_join(name.df)
# note: right now, this doesn't make sense because sample sizes are quite different. We're just seeing how spread out they are

# what percentage of runs had p<0.05?
p.df <- speed.ais.boot %>% gather(key = "var", value = "p", contains("_p")) %>%
  mutate(var = gsub("_p","",var)) %>%
  group_by(var) %>%
  summarise(p.90 = paste0((sum(p<=0.05)/n())*100, "%  p<0.05   ")) %>%
  left_join(name.df)

plot.df <- speed.ais.boot %>% gather(key = "var", value = "val", contains("estimate")) %>%
  mutate(var = gsub("_estimate","",var)) %>% left_join(name.df)

# first plotting factor vars
(pfac <- plot.df %>% drop_na(name) %>% filter(type == "fac") %>% ggplot() +
    geom_histogram(aes(x = val, color = name, fill = name), alpha = 0.4, bins = 60) +
    geom_text(data = p.df %>% filter(type == "fac"), aes(x = Inf, y = Inf, label = p.90),  
              vjust = 1.5,hjust = 1, size = 3) +
    geom_vline(data = orig.df %>% filter(type == "fac"), aes(xintercept = val), linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(length(var.fac), "Set1")[-6]) +
    scale_color_manual(values = RColorBrewer::brewer.pal(length(var.fac), "Set1")[-6]) +  
    facet_wrap(~name, ncol = 1) +
    labs(x = "Estimate", y = "Count", title = "Linear terms") +
    plot.theme + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)))

# then smooth vars
(pnum <- plot.df %>% drop_na(name) %>% filter(type == "num") %>% ggplot() +
    geom_histogram(aes(x = val, color = name, fill = name), alpha = 0.4, bins = 60) +
    geom_text(data = p.df %>% filter(type == "num"), aes(x = Inf, y = Inf, label = p.90),  
              vjust = 1.5,hjust = 1, size = 3) +
    geom_vline(data = orig.df %>% filter(type == "num"), aes(xintercept = val), linetype = "dashed", alpha = 0.5) +
    scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
    facet_wrap(~name, ncol = 1) +
    labs(x = "EDF", y = "Count", title = "Smooth terms") +
    plot.theme + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)))


# now try both together
plot_grid(
  plot_grid(pfac,NULL, ncol = 1, rel_heights = c(14,3)),
  pnum)
ggsave("intermediate-products/gam/swim-speed/ais/swim-speed_ais_boot_res.png", height = 9.5, width = 7)
