# 04 Create suitability map
# LJ 2021-11-26

# get mean probability for each tip
p_suit_all$p_suit_mean <- rowMeans(p_suit_all[,2:ncol(p_suit_all)])

# create columns of suitability to use - initially mean
p_suit_all$p_suit_toUse <- p_suit_all$p_suit_mean

# if tip is a known host, set p to 1
p_suit_all$p_suit_toUse[which(p_suit_all$tip %in% fpest_hosts_ourStates_phy$Host)] <- 1

# assign each cell a suitability value equal to the max association probability of any tree species present in the cell

ourStates_suit <- left_join(ourStates_species_data, p_suit_all, by=c("species" = "tip"))

ourStates_suit_summary <- ourStates_suit %>%
    group_by(cell_id) %>%
    summarise(p_suit = max(na.omit(p_suit_toUse)))

tmp3 <- ourStates_suit %>%
    dplyr::select(c(cell_id, x, y))
tmp3 <- unique(tmp3)

ourStates_suit_summary <- left_join(ourStates_suit_summary, tmp3)

ourStates_suit_summary <- ourStates_suit_summary[,c(3,4,2)]

ourStates_suit_summary[which(is.na(ourStates_suit_summary$p_suit)),3] <- 0


d <- rasterize(x=ourStates_suit_summary[,c(1:2)], y=map_ourStates_refRaster,
               field= ourStates_suit_summary$p_suit)

png(paste0("output/ourStates_testSuit_", gsub(" ", "T", Sys.time()),".png"),
    width = 580, height = 480)
d %>%
    tm_shape() +
    tm_raster() +
    tm_layout(legend.outside = TRUE)
dev.off()

