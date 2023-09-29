#################################
### V3: BOOTSTRAP GAM: DI, FOC ###
#################################
### 17/09/2023
### Tom Grove
### tomgrove20@yahoo.co.uk


### PACKAGES
packages <- c("tidyverse", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "corrplot", "lme4", "MASS", "mgcv", "tidymv", "mgcViz", "gridExtra", "gratia", "ggcorrplot", "cowplot", "ggpubr", "parallel", "foreach")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)


### DATA

# conversion rate from m/s to mph
conv <- 2.23694 # for m/s to mph
conv.sq <- 8052.9692102775 # for m/s2 to mph2

# foc
tot <- read.csv("intermediate-products/response-var-dfs/di_focal.csv") %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  mutate_at(.vars = c("group", "seastate", "year", "folnum.unique", "surface.feeding", "surface.active"), as.factor) %>%
  rename_all(function(x) gsub("focal.", "",x)) %>%
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
  mutate(arcsin.DI = asin(DI),
         across(contains(c("difftime")), .fns = list(log = ~log(.)), .names = "{fn}.{col}"),
         across(contains(c("vess.speed", "dist")), .fns = list(sqrt = ~sqrt(.)), .names = "{fn}.{col}"),
         across(contains(c("accel", "encounter.minute")), .fns = list(cubrt = ~sign(.)*(abs(.)^(1/3))), .names = "{fn}.{col}"),
         across(contains("vess.di"), .fns = list(arcsin = ~asin(.)), .names = "{fn}.{col}"))


# bootstrapped DI
di.boot <- read.csv("intermediate-products/bootstrap/bootstrap-di.csv") %>%
  mutate(datetime = as.POSIXct(datetime))

# model
gam <- readRDS("intermediate-products/gam-v3/di/foc/gam_di_foc_final_v3.rds")


#---------------- BOOTSTRAP --------------------

# formula
var.fac = c("group", "seastate", "year", "surface.feeding", "surface.active") 
var.num = c("log.difftime.prev", "sqrt.vess.speed.60",  "arcsin.vess.di.60", "cubrt.encounter.minute")
var.rand = c("folnum.unique")
f <- as.formula(paste("arcsin.DI ~", 
                      paste("s(",var.num,", bs = 'cs', k=10)", collapse= "+"),"+",
                      paste(var.fac, collapse = "+"), "+",
                      paste("s(",var.rand,", bs = 're', by = dummy)", collapse = "+")))

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
it = 500 # number of iterations

# run!
di.foc.boot <- foreach(i = 1:it, .combine = rbind,.errorhandling = 'remove', .packages = c("tidyverse", "mgcv"), .inorder = FALSE) %dopar% {
  
  # replace DI with iteration
  sample <- tot %>% left_join(di.boot %>% dplyr::select(datetime, move.id.tot, paste0("DI",i))) %>%
    dplyr::select(-arcsin.DI) %>% rename(arcsin.DI = paste0("DI",i)) %>%
    mutate(arcsin.DI = asin(arcsin.DI))
  
  # now run the model
  g <- do.call("bam", list(as.formula(f), data=as.name("sample"), method = "REML", 
                           family = gaussian(link = "identity"), rho = 0.1, start = ar.start))
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

di.foc.boot <- as.data.frame(di.foc.boot)
write.csv(di.foc.boot, "intermediate-products/gam-v3/di/foc/di_foc_boot-res_v3.csv", row.names = FALSE)
# always stop the cluster
parallel::stopCluster(cl = my.cluster)


#---------------- PLOTTING --------------------

#di.foc.boot <- read.csv("intermediate-products/gam-v3/di/foc/di_foc_boot-res_v3.csv")

# first a df for naming
var.fac = c("year2019", "year2020", "seastatechoppy", "grouplone", "surface.active1", "surface.feeding1") 
var.num = c("log.difftime.prev", "sqrt.vess.speed.60",  "arcsin.vess.di.60", "cubrt.encounter.minute")

names <- c("Year (2019)", "Year (2020)", "Sea state", "Group type", "Surface active", "Surface feeding",  "IBI before", "Vessel speed", "Vessel DI", "Encounter minute")
  
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
p.df <- di.foc.boot %>% gather(key = "var", value = "p", contains("_p")) %>%
  mutate(var = gsub("_p","",var)) %>%
  group_by(var) %>%
  summarise(p.90 = paste0((sum(p<=0.004)/n())*100, "%  p<0.004   ")) %>%
  left_join(name.df)

plot.df <- di.foc.boot %>% gather(key = "var", value = "val", contains("estimate")) %>%
  mutate(var = gsub("_estimate","",var)) %>% left_join(name.df)

# first plotting factor vars
(pfac <- plot.df %>% drop_na(name) %>% filter(type == "fac") %>% ggplot() +
    geom_histogram(aes(x = val, color = name, fill = name), alpha = 0.4, bins = 60) +
    geom_text(data = p.df %>% filter(type == "fac"), aes(x = Inf, y = Inf, label = p.90),  
              vjust = 1.5,hjust = 1, size = 3) +
    geom_vline(data = orig.df %>% filter(type == "fac"), aes(xintercept = val), linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(length(var.fac), "Set1")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(length(var.fac), "Set1")) +  
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
          plot_grid(pnum,NULL, ncol = 1, rel_heights = c(2.28,1)))
ggsave("intermediate-products/gam-v3/di/foc/di_foc_boot_res_v3.png", height = 11, width = 7)
