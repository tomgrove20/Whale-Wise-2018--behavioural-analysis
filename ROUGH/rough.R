
# ---- NEW COLUMN = DIFFERENCE BETWEEN RF DISTANCE AND CAM DISTANCE

df <- read.csv("intermediate-products/follows_2018-20_cleaned_calc_full.csv") %>%
  mutate(dist.diff = dist.rf-dist.cam)
write.csv(df, "dist-rf-vs-cam.csv", row.names = FALSE)

# --- CHECKING OUT SPEEDS!

df <- st_read("final-products/follows+foc+ais+bath_2018-20_Mn_lines.gpkg") %>%
  as.data.frame() %>% dplyr::select(-geom)
write.csv(df, "test.csv", row.names = FALSE)


df2 <- tot %>% as.data.frame() %>% dplyr::select(-geometry)
write.csv(df2, "test-pos.csv")


### from surface-active_ais
## PARTIALS
# oak
plot(b, allTerms = T, select = 1) + labs(x = "oak boats (10 min, 1500 m)", y = "s(oak boats (10 min, 1500 m))") + 
  partial.theme
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_partial_oak.10.1500.png", height=5, width=7)

# rib
plot(b, allTerms = T, select = 2) + labs(x = "RIB boats (30 min, 1500 m)", y = "s(RIB boats (30 min, 1500 m))") + 
  partial.theme
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_partial_rib.30.1500.png", height=5, width=7)

# folnum.unique
plot(b, allTerms = T, select = 3) + labs(y = "Effects", title = "s(follow number)") + partial.theme
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_partial_folnum.png", height=5, width=7)

# group
plot(b, allTerms = T, select = 4) + partial.theme + labs(y = "Partial effect of group")
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_partial_group.png", height=5, width=7)

## PREDICTED
# first predicting
p <- predict(gam, type = "response")
# need to do extra to get proper 95% CIs. https://stats.stackexchange.com/questions/33327/confidence-interval-for-gam-model
p.se <- predict(gam, type = "link", se.fit = TRUE) 
upr <- gam$family$linkinv(p.se$fit + (2 * p.se$se.fit))
lwr <- gam$family$linkinv(p.se$fit - (2 * p.se$se.fit))

pred <- tot %>%
  mutate(predict = p$fit, # fitted value
         upr = upr, lwr = lwr) # confidence intervals

# now plotting
# oak
pred %>% group_by(oak.num.10.1500) %>% summarise(mean = mean(predict), upr = mean(upr), lwr = mean(lwr)) %>%
  ggplot(aes(x = oak.num.10.1500, y = mean)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
  labs(x = "oak boats (10 min, 1500 m)", y = "Rate of surface activity") + partial.theme
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_pred_oak.num.10.1500.png", height=5, width=7)

# rib
pred %>% group_by(rib.num.30.1500) %>% summarise(mean = mean(predict), upr = mean(upr), lwr = mean(lwr)) %>%
  ggplot(aes(x = rib.num.30.1500, y = mean)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = "dashed", color = "black", fill = "transparent") +
  labs(x = "RIB boats (30 min, 1500 m)", y = "Rate of surface activity") + partial.theme
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_pred_rib.num.30.1500.png", height=5, width=7)

# group
ggplot(pred, aes(x = group, y = pred)) + geom_boxplot() + coord_cartesian(ylim=c(0, 0.02)) + partial.theme + 
  labs(x = "Group type", y = "Rate of surface activity")
ggsave("intermediate-products/gam/surface-active/gam_surface-active_ais_pred_group.png", height=5, width=7)