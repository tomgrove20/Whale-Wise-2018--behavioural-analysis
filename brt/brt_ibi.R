################
### BRT: IBI ###
################
### 07/02/2022
### Tom Grove
### tomgrove20@yahoo.co.uk

# part A: focal vess model
# part B: AIS vess model

### PACKAGES
packages <- c("tidyverse", "viridis", "ppcor","RColorBrewer", "scales","survMisc", "Metrics", "dismo", "gbm", "corrplot")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

### THEME
source("codes/themes.R")
brewer.pal(n=9,"GnBu")
brewer.pal(n=9, "YlOrRd")


############
### DATA ###
############
# IBI = 

# focal vessel + follows
fol.foc <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_calc.csv")

# ais + follows 
fol.ais <- read.csv("intermediate-products/follows/follows+foc+ais+bath_2018-20_no-other-boats_calc.csv")

##########################
##########################
###### FOCAL VESSEL ######-------------------
##########################
##########################

# first defining which variables we use
vars.foc <- c() # from

############################
### TIME FRAME SELECTION ###
############################
# defining variable combinations
ts = c(10,60,300)
var.s = c() 
tran.s = c()
fun.s = c("log", "log")
var.t = c() 
trans.t = c()
fun.t = c("log", "log") 

# Combining together into a lagged data frame
t.df <- rbind(data.frame(var = rep(var.t, each = length(ts)), t = rep(t,length(var.t))), t.orig = 60) %>%
  mutate(var.orig = paste0(var,".",t.orig), var.target = paste0(var,".",t))
                
(vars <- c(var.s, unique(t.df$var.orig)))

# running the loop
for (i in 1:nrow(t.df)) {

  v <- vars %>% dplyr::recode(!!t.df[i,"var.orig"] := lag.df[i,"var.target"]) # variable list
  
  indices <- vector() ; for(j in v){indices = c(indices, grep(j, colnames(tot), fixed = TRUE))} # column indices
  
  # run the gbm
  set.seed(0)
  t <- gbm.step(data = tot, gbm.x = indices, gbm.y = 7, family = "gaussian",
                tree.complexity = 2, n.trees = 1000, max.trees = 1000,
                learning.rate = 0.01, bag.fraction = 0.75, keep.data = TRUE) # tc = 2: 2D interactions
  
  # now calculate deviance explained (2 metrics) and relative variable influence
  null <- t$self.statistics$mean.null # null total deviance
  res <- t$self.statistics$mean.resid # residual deviance
  cv <- t$cv.statistics$deviance.mean # CV deviance
  t.df[i,"dev.ex.res"] <- 1 -(res/null) # residual deviance explained (proportion)
  t.df[i,"dev.ex.cv"] <- (null - cv)/null # cross-validated deviance explained (proportion)
  rownames(t$contributions) = gsub("`","",rownames(t$contributions))
  t.df[i,"inp"] = t$contributions[t.df[i,"var.target"],"rel.inf"] # relative variable influence
  
  print (t.df[i,"var.target"]) # progress checker
}

t.df
write.csv(t.df, "") # saving data frame

# plotting CV deviance explained and variable importance by var
ggplot(data = t.df, aes(x = ts, y = dev.ex.cv, group = var, color = var)) +
  geom_point(size = 2, alpha = 0.5) + geom_line() + 
  scale_color_brewer(palette = "Paired") +
  labs(x = "Month", color = "Variable", y = "CV deviance explained") + plot.theme
ggsave("")

ggplot(data = lag.df, aes(x = ts, y = inp, group = var, color = var)) +
  geom_point() + geom_line() + 
  scale_color_brewer(palette = "Paired") +
  labs(x = "Month", color = "Variable", y = "Relative importance") + plot.theme
ggsave("")

# figuring out which lag is best!
t.sel <- t.df %>% dplyr::select(var, var.target, dev.ex.cv) %>%
  group_by(var) %>% slice(which.max(dev.ex.cv))
(var.sel <- t.sel$var.target)


###################
### FINAL MODEL ###
###################

vars <- c(var.s, var.sel) # full vars for final model

indices <- vector() ; for(i in vars){indices = c(indices, grep(i, colnames(tot)))}
t <- gbm.step(data = tot, gbm.x = indices, gbm.y = 7, family = "gaussian",
              tree.complexity = 2,
              learning.rate = 0.01, step.size = 100,
              bag.fraction = 0.75,
              keep.data = TRUE)

# removing variables that contribute <3%
filter <- summary(t) %>% filter(rel.inf<3) %>%
  mutate(var = gsub("\\(", "\\\\(",gsub("\\)", "\\\\)", gsub("`","", var))))
vars <- setdiff(vars, filter$var)

# and run again ...
indices <- vector() ; for(i in vars){indices = c(indices, grep(i, colnames(tot)))}
t <- gbm.step(data = tot, gbm.x = indices, gbm.y = 7, family = "gaussian",
              tree.complexity = 2,
              learning.rate = 0.01, step.size = 100,
              bag.fraction = 0.75,
              keep.data = TRUE)

saveRDS(t, "")

# CV deviance explained
(null <- t$self.statistics$mean.null)
(cv <- t$cv.statistics$deviance.mean)
(dev.ex <- (null-cv)/null)

### PLOTS
# variable importance
ggplot(data = summary(t) %>% mutate(var = fct_relevel(var, rev(summary(t)$var))), aes(x = rel.inf, y = var)) +
  geom_bar(stat = "identity", aes(fill = rel.inf)) + barplot.theme + labs(x = "Relative influence (%)", y = "") +
  theme(legend.position = "none") +
  scale_fill_continuous(low = "#7BCCC4", high = "#084081") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
ggsave("")

# partial dependency plots (first need to reformat var names)
for (i in 1:length(t$gbm.call$predictor.names)) {
  if(grepl("\\(",t$gbm.call$predictor.names[i])){
    t$gbm.call$predictor.names[i] <- paste0("`",t$gbm.call$predictor.names[i],"`")
  }
}

gbm.plot(t, smooth = TRUE, rug = FALSE, y.label = "marginal effect", common.scale = FALSE, plot.layout = c(2,3), lwd = 2)
png("", width = 2000, height = 1500, pointsize = 50)
gbm.plot(t, smooth = TRUE, rug = FALSE, y.label = "marginal effect", common.scale = FALSE, plot.layout = c(2,3), lwd = 2)
dev.off()

# strength of interactions
names <- rev(as.vector(rownames(as.data.frame(gbm.interactions(t)$interactions))))
int <- as.data.frame(gbm.interactions(t)$interactions) %>% rownames_to_column(var = "v1") %>% gather(2:ncol(.),key = "v2", value = "int")  %>% mutate_at(vars(v1,v2), ~fct_relevel(., names)) %>% 
  group_by(v2) %>% mutate(sum = sum(int)) %>% filter(sum>0) %>% ungroup() %>%
  group_by(v1) %>% mutate(sum = sum(int)) %>% filter(sum>0) %>% ungroup()
ggplot(data = int, aes(x = v1, y = v2, fill = int)) + geom_tile() + 
  geom_text(aes(label=ifelse(int>0,int,""))) + grid.theme + 
  scale_y_discrete(limits = rev) + scale_x_discrete(limits = rev, position = "top") +
  scale_fill_continuous(low = "white", high = "#FC4E2A", trans = "log", na.value = "lightgrey", breaks = c(0,0.1,1,10)) + 
  labs(x = "", y = "", fill = "Interaction") + coord_cartesian(expand = c(0,0)) +
  theme(axis.text.x = element_text(hjust = 0))
ggsave("")



#################
#################
###### AIS ######-------------------
#################
#################




