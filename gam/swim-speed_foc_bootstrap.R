#################################
### BOOTSTRAP GAM: SPEED, FOC ###
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

# foc
tot <- read.csv("intermediate-products/response-var-dfs/speed_focal.csv") %>%
  mutate_at(vars(contains("datetime")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.active", "surface.feeding"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  mutate(dummy = 1) %>%
  # important dummy variable to factor out random effect when predicting response
  filter(vess.accel.max.300 < 5) %>% # filtering until properly processed
  filter(vess.accel.60 >-0.2) %>% # filtering until properly processed
  filter(vess.speed.var.60<3) # filtering until properly processed

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("diff.time","vess.speed.var")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(ends_with(c("speed.60", "speed.300", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))

# bootstrapped speeds
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv") %>%
  mutate(datetime = as.POSIXct(datetime))

# model
gam <- readRDS("intermediate-products/gam/swim-speed/foc/gam_swim-speed_foc_final.rds")


#---------------- BOOTSTRAP --------------------

# formula
var.fac = c("seastate", "year", "surface.feeding", "surface.active")
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.60")
var.rand = c("folnum.unique")
f <- as.formula(paste("speed ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

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
it = 500 # number of iterations

# now run!
speed.foc.boot <- foreach(i = 1:it, .combine = rbind, .errorhandling = 'remove', .packages = c("tidyverse", "mgcv")) %dopar% {
  
  # replace speed with iteration
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
  
  return(df)
  
  print(paste0("iteration ",i))
  
}

speed.foc.boot <- as.data.frame(speed.foc.boot)

write.csv(speed.foc.boot, "intermediate-products/gam/swim-speed/foc/swim-speed_foc_boot-res.csv", row.names = FALSE)
# always stop the cluster
parallel::stopCluster(cl = my.cluster)


#---------------- PLOTTING --------------------

#speed.foc.boot <- read.csv("intermediate-products/gam/swim-speed/foc/swim-speed_foc_boot-res.csv")

# first a df for naming
var.fac = c("seastatechoppy", "year2019", "year2020", "surface.feeding1", "surface.active1")
var.num = c("day.of.year", "sqrt.dist", "log.diff.time", "sqrt.vess.speed.60",  "log.vess.speed.var.300", "arcsin.vess.di.60", "vess.accel.60")
name.df <- data.frame(
  var = c(var.fac, var.num),
  name = c("Sea state", "Year 2019", "Year 2020", "Surface feeding", "Surface active", "Julian day", "Distance", "Inter-breath interval", "Vessel speed (60 sec.)", "Vessel speed SD (60 sec.)", "Vessel DI (300 sec.)", "Vessel acceleration (60 sec.)")) %>%
  mutate(type = ifelse(var %in% var.fac, "fac","num"))
nrow = max(length(var.fac), length(var.num))

# we need to create a df from the original GAM so we can add lines for full model effect sizes
summary <- summary(gam)
line.coeff <- as.data.frame(summary$p.coeff) %>% add_rownames(var = "var") %>% rename(val = 2)
smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
smooth.edf <- as.data.frame(summary$edf) %>% mutate(var = smooth.names) %>% rename(val = 1)
orig.df <- as.data.frame(bind_rows(line.coeff, smooth.edf)) %>% left_join(name.df)
# note: right now, this doesn't make sense because sample sizes are quite different. We're just seeing how spread out they are

# what percentage of runs had p<0.05?
p.df <- speed.foc.boot %>% gather(key = "var", value = "p", contains("_p")) %>%
  mutate(var = gsub("_p","",var)) %>%
  group_by(var) %>%
  summarise(p.90 = paste0((sum(p<=0.05)/n())*100, "%  p<0.05   ")) %>%
  left_join(name.df)
  
plot.df <- speed.foc.boot %>% gather(key = "var", value = "val", contains("estimate")) %>%
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
  plot_grid(pfac,NULL, ncol = 1, rel_heights = c(8.45,3)),
  pnum)
ggsave("intermediate-products/gam/swim-speed/foc/swim-speed_foc_boot_res.png", height = 12.5, width = 7)





