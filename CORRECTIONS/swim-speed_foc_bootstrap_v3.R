#################################
### V2: BOOTSTRAP GAM: SPEED, FOC ###
#################################
### 17/09/2023
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
isnt.na = function(x){!is.na(x)}


#---------------- DATA --------------------

# conversion rate from m/s to mph
conv <- 2.23694 # for m/s to mph
conv.sq <- 8052.9692102775 # for m/s2 to mph2

# foc
tot <- read.csv("intermediate-products/response-var-dfs/speed_focal.csv") %>%
  mutate_at(vars(contains("datetime")), as.POSIXct) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.active", "surface.feeding"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
  filter_at(vars(contains("vess")), isnt.na) %>%
  mutate(dummy = 1) %>% # important dummy variable to factor out random effect when predicting response
  # we also need to convert m/s to mph!
  mutate(across(contains("vess.speed"), ~.*conv), # speed and sd speed
         across(contains("accel"), ~.*conv.sq)) %>% # acceleration
  filter(vess.accel.max.60<20000, vess.speed.var.60 < 9)

# now we need to add encounter minute to this (rough!!)
follows.encmin <- read.csv("intermediate-products/follows/follows_start+end_consecutive_rough.csv") %>%
  mutate(folnum.unique = as.factor(folnum.unique))

tot <- tot %>% left_join(follows.encmin) %>%
  group_by(enc.num) %>%
  mutate(encounter.minute = as.numeric(difftime(datetime, start, unit = "mins"))) %>% ungroup()

# creating new transformed columns
tot <- tot %>%
  mutate(across(contains(c("diff.time")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(contains(c("vess.speed", "encounter.minute")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains(c("accel", "dist")), .fns = list(cubrt = ~sign(.)*(abs(.)^(1/3))), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))

# bootstrapped speeds
speed.boot <- read.csv("intermediate-products/bootstrap/bootstrap-speed.csv") %>%
  mutate(datetime = as.POSIXct(datetime))

# model
gam <- readRDS("intermediate-products/gam-v3/swim-speed/foc/gam_swim-speed_foc_final_v3.rds")


#---------------- BOOTSTRAP --------------------

# formula
var.fac = c("surface.feeding") 
var.num = c("log.diff.time", "sqrt.vess.speed.60", "arcsin.vess.di.60", "sqrt.encounter.minute")
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
speed.foc.boot <- foreach(i = 1:it, .combine = rbind, .errorhandling = 'remove', .packages = c("tidyverse", "mgcv"), .inorder = FALSE) %dopar% {
  
  # then replace speed with iteration
  sample <- tot %>% left_join(speed.boot %>% dplyr::select(datetime, move.id.tot, paste0("speed",i))) %>%
    dplyr::select(-speed) %>% rename(speed = paste0("speed",i))
  
  # now run the model
  g <- do.call("gam", list(as.formula(f), data=as.name("sample"), method = "REML", family = gaussian(link = "log")))
  summary <- summary(g)
  
  # now extracting important bits!
  line.coeff <- t(as.data.frame(summary$p.coeff)); rownames(line.coeff) = NULL; colnames(line.coeff) = paste0(colnames(line.coeff),"_estimate")
  line.p <- t(as.data.frame(summary$p.pv)); rownames(line.p) = NULL; colnames(line.p) = paste0(colnames(line.p),"_p")
  
  smooth.names <- gsub("s\\(|\\)","",names(summary$chi.sq))
  smooth.edf <- t(as.data.frame(summary$edf)); rownames(smooth.edf) = NULL; colnames(smooth.edf) = paste0(smooth.names,"_estimate")
  smooth.p <- t(as.data.frame(summary$s.pv)); rownames(smooth.p) = NULL; colnames(smooth.p) = paste0(smooth.names,"_p")
  
  df <- cbind(line.coeff, line.p, smooth.edf, smooth.p)
  
  x <- "hello"
  write.csv(x, paste0("intermediate-products/progress-checker/iteration ",i,".csv"))
  
  return(df)
  
}

speed.foc.boot <- as.data.frame(speed.foc.boot)

write.csv(speed.foc.boot, "intermediate-products/gam-v3/swim-speed/foc/swim-speed_foc_boot-res_v3.csv", row.names = FALSE)
# always stop the cluster
parallel::stopCluster(cl = my.cluster)


#---------------- PLOTTING --------------------

#speed.foc.boot <- read.csv("intermediate-products/gam-v3/swim-speed/foc/swim-speed_foc_boot-res_v3.csv")

# first a df for naming
var.fac = c("surface.feeding1") 
var.num = c("log.diff.time", "sqrt.vess.speed.60", "arcsin.vess.di.60", "sqrt.encounter.minute")
names <- c("Surface feeding", "Inter-breath interval",  "Vessel speed", "Vessel DI", "Encounter minute")

name.df <- data.frame(
  var = c(var.fac, var.num),
  name = factor(names, levels = names)) %>%
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
  summarise(p.90 = paste0((sum(p<=0.004)/n())*100, "%  p<0.004   ")) %>%
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
plot_grid(plot_grid(pfac, NULL, ncol = 1, rel_heights = c(1,1.95)),
          pnum)
ggsave("intermediate-products/gam-v3/swim-speed/foc/swim-speed_foc_boot_res_v3.png", height = 8, width = 7)





