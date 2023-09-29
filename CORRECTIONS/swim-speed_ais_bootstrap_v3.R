#################################
### V3: BOOTSTRAP GAM: SPEED, AIS ###
#################################
### 16/09/2023
### Tom Grove
### tomgrove20@yahoo.co.uk

### PACKAGES
packages <- c("tidyverse","RColorBrewer", "scales","survMisc", "Metrics", "lme4", "MASS", "mgcv", "cowplot", "parallel", "foreach", "patchwork")
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
         across(contains(c("diff.time")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"), # sqrt dist max
         log.speed = log(speed))

# now we need to create AR.start (defining the start of each follow to group AR structure) and join this to tot (to enable sampling)
ar.start <- (tot %>% group_by(folnum.unique) %>% mutate(n = row_number()) %>% 
               mutate(ar.start = ifelse(n == min(n), TRUE, FALSE)))$ar.start
tot <- mutate(tot, ar.start = ar.start)

# bootstrapped speeds
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv") %>%
  mutate(datetime = as.POSIXct(datetime))

# model
gam <- readRDS("intermediate-products/gam-v3/swim-speed/ais/gam_swim-speed_ais_final_v3.rds")


#---------------- BOOTSTRAP --------------------

# formula
var.num = c("day.of.year", "log.diff.time")
var.fac = c("group", "year", "surface.feeding")
var.rand = c("folnum.unique")
f <- as.formula(paste("log.speed ~", 
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
n = 3940 # size of random sample
it = 500 # number of iterations

# and run!
speed.ais.boot <- foreach(i = 1:it, .combine = rbind, .errorhandling = 'remove', .packages = c("tidyverse", "mgcv"), .inorder = FALSE) %dopar% {
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i)) %>% 
    mutate(log.speed = log(speed))
  
  # now run the model
  g <- do.call("bam", list(as.formula(f), data=as.name("sample"), method = "REML", 
                             family = gaussian(link = "identity"), rho = -0.05, start = sample$ar.start))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  # progress checker
  write.csv(paste0("intermediate-products/progress-checker/iteration ",i,".csv"))
  
  #speed.ais.boot <- rbind(speed.ais.boot,df)
  return(df)
  
}

speed.ais.boot <- as.data.frame(speed.ais.boot)
view(speed.ais.boot)
write.csv(speed.ais.boot, "intermediate-products/gam-v3/swim-speed/ais/swim-speed_ais_boot-res_v3.csv")



#---------------- PLOTTING --------------------

#speed.ais.boot <- read.csv("intermediate-products/gam-v3/swim-speed/ais/swim-speed_ais_boot-res_v3.csv")

# first a df for naming
var.fac = c("year2019", "year2020", "grouplone", "surface.feeding1")
var.num = c("day.of.year", "log.diff.time")

names <- c("Year (2019)", "Year (2020)", "Group type", "Surface feeding", "Julian day", "Inter-breath interval")

name.df <- data.frame(
  var = c(var.fac, var.num),
  name = factor(names, levels = names)) %>%
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
  summarise(p.90 = paste0((sum(p<=0.004)/n())*100, "%  p<0.004   ")) %>%
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
plot_grid(pfac,
          plot_grid(pnum, NULL, ncol = 1, rel_heights = c(1.27,1)))
ggsave("intermediate-products/gam-v3/swim-speed/ais/swim-speed_ais_boot_res_v3.png", height = 8, width = 7)
