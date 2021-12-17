# 03 Phylogenetic models
# LJ 2021-11-26


# modifications needed 2021-12-16:
# - remove self-pairs from logistic regression
# - figure out repeated coefficient method from gilbert/robles rather than using 'predict'
# - identify 3 good pests to focus on
# - harmonize species names in phy/FIA list


library(ape)
library(castor)
# library(MASS)
# library(pscl)
# library(DAAG)

# read in list of host species for which we have range maps
# list must be a single column of names in format 'Genus_species'
# includes 'var.' and 'ssp.'
trees_withRanges <- read.csv("./output/FIAspp_toUse_2021-12-15T15:33:29.csv")

# read in symbiont-host data
symbHost <- read.csv("./output/pest-host_clean_noSSp_2021-12-15T11:53:40.csv")

symbHost$Symbiont[which(symbHost$Symbiont=="Homaduala_anisocentra")] <- "Homadaula_anisocentra"

# read in plant phylogeny
phy <- read.tree("./output/focalPhy_2021-12-16T11:14:16.nwk")

# identify host plants to use: species must be present in all 3 data sources
ourHosts <- intersect(intersect(trees_withRanges$species,
                                symbHost$HostPlant),
                      phy$tip.label)

# subset symbiont-host data to hosts in ourHosts
ourSymbHosts <- subset(symbHost, HostPlant %in% ourHosts)

# expand to matrix
matrix_ourSymbHosts <- table(ourSymbHosts)

# sort by number of hosts
matrix_ourSymbHosts_sorted <- matrix_ourSymbHosts[order(rowSums(matrix_ourSymbHosts),decreasing=T),]

# hist(rowSums(matrix_ourSymbHosts)) # not many generalists

# choose pest for the analysis
rank <- 1
fpest <- row.names(matrix_ourSymbHosts_sorted)[rank]

# get list of all global hosts of focal pest
fpest_hosts <- subset(symbHost, Symbiont == fpest)
fpest_hosts <- fpest_hosts$HostPlant

# subset to only hosts in phy
fpest_hosts_phy <- intersect(fpest_hosts, phy$tip.label)

# prune phylogeny
phy_pruned <- phy
# only host species in pest db
# phy_pruned <- keep.tip(phy, intersect(unique(symbHost$HostPlant),
#                                       phy$tip.label))

# # minimal clade including all known hosts of focal pest
# phy_pruned <- extract.clade(phy, getMRCA(phy, fpest_hosts_phy))

# # minimal clade including all host species in pest db
# phy_pruned <- extract.clade(phy, getMRCA(phy,
#                                          intersect(unique(symbHost$HostPlant),
#                                                    phy$tip.label)))

# family Sapindaceae
# phy_pruned <- extract.clade(phy, "Sapindaceae")

# Ntip(phy_pruned)

# get all pairwise distances between tips on pruned phylogeny
pd_all <- get_all_pairwise_distances(phy_pruned,
                                     only_clades = 1:Ntip(phy_pruned))

# prep data frame for loop output
p_suit_all <- data.frame(tip = phy_pruned$tip.label)

# set run date/time ID
# now <- gsub(" ", "T", Sys.time())

# check for output directories and create if necessary
# if(!dir.exists("output")){dir.create("output")}
# if(!dir.exists("output/phylo_models")){dir.create("output/phylo_models")}
# if(!dir.exists(paste0("output/phylo_models/", now, "_", Ntip(phy_pruned)))){
#     dir.create(paste0("output/phylo_models/", now, "_", Ntip(phy_pruned)))}

# loop to get predicted probabilities for each tip for all focal hosts
for(i in seq_along(fpest_hosts_phy)){
    
    # choose one host to be 'focal host' for PD calculations
    fhost <- fpest_hosts_phy[i]
    
    # subset to only distances to/from focal host
    pd_fpest <- pd_all[which(phy_pruned$tip.label==fhost),]
    
    # assemble data for logistic regression
    fpest_regdata <- data.frame(tip = phy_pruned$tip.label,
                                pd = pd_fpest,
                                hostStatus = 0)
    
    fpest_regdata[which(fpest_regdata$tip%in%fpest_hosts_phy),]$hostStatus <- 1
    
    fpest_regdata$log_pd <- log10(fpest_regdata$pd + 1)
    
    # fpest_regdata <- fpest_regdata[-which(fpest_regdata$tip==fhost),]
    
    # plot(jitter(fpest_regdata$hostStatus)~fpest_regdata$log_pd)
    
    # logistic regression
    fpest_logfit <- glm(hostStatus ~ log_pd,
                        family = binomial(link = "logit"),
                        data = fpest_regdata)
    
    # summary(fpest_logfit)
    # 
    # DAAG::CVbinary(fpest_logfit)
    # 
    # table(round(predict(fpest_logfit, type="response"), 4))
    # 
    # plot(predict(fpest_logfit, type="response")~fpest_regdata$log_pd)
    
    # tipcols <- rep("black", Ntip(phy_pruned))
    # tipcols[which(phy_pruned$tip.label%in%fpest_hosts_phy)] <- "red"
    # tipcols[which(phy_pruned$tip.label==fhost)] <- "blue"
    # plot(phy_pruned, tip.color = tipcols)
    
    # # NB
    # NB_fit <- glm.nb(hostStatus ~ log_pd,
    #              data = fpest_regdata)
    # 
    # # ZIP
    # ZIP_fit <- zeroinfl(hostStatus ~ log_pd | log_pd,
    #             dist = 'poisson',
    #             data = fpest_regdata)
    # 
    # summary(ZIP_fit)
    #     
    # # ZINB
    # ZINB_fit <- zeroinfl(hostStatus ~ log_pd | log_pd,
    #                dist = 'negbin',
    #                data = fpest_regdata)
    # summary(ZINB_fit)
    
    
    ## residual deviance should be similar to df
    
    # E2 <- resid(ZINB_fit, type = "pearson")
    # N  <- nrow(fpest_regdata)
    # p  <- length(coef(ZINB_fit)) + 1 # '+1' is due to theta
    # sum(E2^2) / (N - p)
    # should be close to 1 (>1:overdispersion, <1:underdispersion)
    
    
    # get predicted probabilities of association for all hosts in phy
    fpest_pred <- fpest_regdata
    fpest_pred$p_suit <- predict(fpest_logfit, type="response")
    
    # hist(fpest_pred$p_suit)
    
    p_suit_all[,ncol(p_suit_all)+1] <- fpest_pred$p_suit
    
    names(p_suit_all)[ncol(p_suit_all)] <- paste0("p_suit_", i) 
    
    # # plot phylogeny with tips coloured by predicted probability
    # fpest_pred$p_suit_round <- round(fpest_pred$p_suit, 1)
    # 
    # phy_plot <- extract.clade(phy_pruned,
    #                           getMRCA(phy_pruned,
    #                                   fpest_pred$tip[which(
    #                                       fpest_pred$p_suit_round!=0.0)]))
    # 
        # plot phylogeny with edges coloured by predicted probability
    # create palette of 10 colours
    # cols <- heat.colors(11)
    # cols[11] <- "black"
    # 
    # create edge colour object, default colour is black
    # edge_colours <- rep("black", Nedge(phy_plot))
    
    # give terminal edges a colour dependent on their predicted suitability value
    # for(j in 1:Nedge(phy_plot)){
    #     edge_colours[which.edge(phy_plot, phy_plot$tip.label[j])] <- cols[11-(fpest_pred[which(fpest_pred$tip==phy_plot$tip.label[j]),6]*10)]
    # }
    
    # colour edge leading to focal host blue
    # edge_colours[which.edge(phy_plot, fhost)] <- "blue"
    
    # plot
#     png(paste0("output/phylo_models/", now, "_", Ntip(phy_pruned), "/fhost-", fpest_hosts_phy[i],
#                "_", Ntip(phy_pruned), ".png"))
#     plot(phy_plot, edge.color=edge_colours)
#     dev.off()
#     
}

# get mean probability for each tip
p_suit_all$p_suit_mean <- rowMeans(p_suit_all[,2:ncol(p_suit_all)])


# write predicted suitabilities to file
suits <- subset(p_suit_all[,c("tip", "p_suit_mean")],
                tip %in% trees_withRanges$species)

write.csv(suits,
          paste0("./output/FIAtrees_predSuitability_",
                 gsub(" ", "T", Sys.time()),
                 ".csv"),
          row.names = FALSE)
