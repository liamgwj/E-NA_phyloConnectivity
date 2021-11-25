# LJ started 2021-11-10
# script for MSc plant symbiont connectivity simulation

# load packages:
library(dplyr)
library(sf)
library(raster)
library(tmap)
library(ape)
library(castor)

# set variable parameters:
    # states to include
    ourStates <- c("Pennsylvania", "New York", "Vermont", "New Hampshire", "Maine", "Connecticut", "Massachusetts")

# prep administrative borders map ---------------------------------------------
# read in borders shapefile
map_borders <- sf::st_read("./indata/admin_shapes/ne_10m_admin_1_states_provinces_lakes.shp")

# filter to USA
map_USA <- map_borders %>%
    filter(admin == "United States of America")

# Get rid of extraneous columns in the object
map_USA <- map_USA[, c(1:60, 84)]
map_USA <- map_USA[, -c(6)]

# filter to chosen state(s)
map_ourStates <- map_USA %>%
    filter(name %in% ourStates)


# get CRS and extent of species raster when cropped to chosen state(s) --------
# read in a species occurrence map
map_species_demo <- sf::st_read("./indata/USDA_data/sp318_hybrid_2100.shp")

# get CRS
crs_speciesMaps <- st_crs(map_species_demo)

# re-project state(s) to the same coordinate system as tree data
map_ourStates_alb <- sf::st_transform(map_ourStates, crs_speciesMaps)

# clip species raster to extent of state(s)
map_species_ourStates_demo <- map_species_demo[map_ourStates_alb, ]

# get extent of clipped raster
extent_ourStates <- raster::extent(map_species_ourStates_demo)


# create reference raster of our state(s) -------------------------------------

# calculate the number of rows and columns in the raster from the resolution and the min/max coordinates
# difference between max x-coord and min x-coord divided by 10000 metres (the resolution of each cell)
extent_numcols <- (extent_ourStates[2] - extent_ourStates[1]) / 10000
extent_numrows <- (extent_ourStates[4] - extent_ourStates[3]) / 10000

# Using this information let's create a new reference raster:
map_ourStates_refRaster <- raster::raster(ncols = extent_numcols, nrows = extent_numrows,
                                        xmn = extent_ourStates[1], ymn = extent_ourStates[3],
                                        ext = extent_ourStates,
                                        resolution = 10000,
                                        crs = crs_speciesMaps)


## Produce raster stack (one layer per species) -------------------------------
# For this procedure we need to cycle through the list of species and: 
# * import the corresponding shapefile
# * clip to the Pennsylvania extent
# * rasterize to a raster layer within a raster stack

# First let's import the CSV file that contains the species information:
data_species <- read.csv("./indata/USDA_data/_Species_List.csv")

# rename variables to remove spaces
names(data_species) <- c("FIA_species_code", "Common_name", "Scientific_name", "Model_reliability")

# we'll loop through each row of the file for which "Model_reliability" is not equal to "FIA Only"
data_species_toUse <- data_species %>%
                        filter(Model_reliability != "FIA Only")

# Let's now create a vector of filenames and paths, and this is what we'll loop through:
filenames_species_toUse <- paste(paste("./indata/USDA_data/sp", 
                                 data_species_toUse$FIA_species_code, sep = ""),
                           "hybrid_2100.shp", sep = "_")

# We'll start with the first species in the list to establish the first layer
# in a raster stack, then we'll loop through remaining species to add to the
# stack.

# **NOTE** we'll create 2 different raster stacks: one for current species
# importance values (RFimp), and one for the averaged prediction for three
# general circulation models for future conditions "GCM45".

# rasterize:
species_imp_current <- raster::rasterize(x = map_species_ourStates_demo, 
                                         y = map_ourStates_refRaster, field = "RFimp")
species_imp_future <- raster::rasterize(x = map_species_ourStates_demo, 
                                        y = map_ourStates_refRaster, field = "GCM45")

# Start big loop
for (i in 2:length(filenames_species_toUse)) {
    
    ## conflicting names here need to be sorted out
    # import file
    map_species <- sf::st_read(filenames_species_toUse[i])
    
    # clip to Pennsylvania
    map_species_ourStates <- map_species[map_ourStates_alb, ]
    
    # rasterize:
    species_imp_current <- raster::stack(species_imp_current, 
                                         raster::rasterize(x = map_species_ourStates,
                                                           y = map_ourStates_refRaster, field = "RFimp"))
    species_imp_future <- raster::stack(species_imp_future,
                                        raster::rasterize(x = map_species_ourStates,
                                                          y = map_ourStates_refRaster, field = "GCM45"))
    print(i)
}

# Name each layer in the stack by its latin name:
names(species_imp_current) <- data_species_toUse$Scientific_name
names(species_imp_future) <- data_species_toUse$Scientific_name

# coords
numvals <- nrow(species_imp_current)*ncol(species_imp_current)
ourStates_xy_coords <- as_tibble(
    data.frame(
        cell_id = 1:numvals,
        raster::xyFromCell(species_imp_current[[1]], 1:numvals)
    )
)

# RFimp
ourStates_RFimp <- as_tibble(
    data.frame(
        raster::extract(species_imp_current, ourStates_xy_coords[, 2:3], cellnumbers = TRUE)
    )
)
# rename first column
names(ourStates_RFimp)[1] <- "cell_id"

# Now we can combine these into long format, keeping only the species present
# in a given cell.

ourStates_species_data <- ourStates_RFimp %>%
    tidyr::pivot_longer(cols = names(ourStates_RFimp)[2:ncol(ourStates_RFimp)],
                        names_to = "species",
                        values_to = "RFimp") %>%
    dplyr::filter(!is.na(RFimp)) %>%  # get rid of NAs
    dplyr::filter(RFimp > 0) %>%  # keep only non-zero importance values
    dplyr::left_join(ourStates_xy_coords, by = "cell_id") %>%  # join with coordinates
    dplyr::select(cell_id, x, y, species, RFimp)  # reorder columns

# Have a look at result:
# ourStates_species_data


# replace spaces in species names with underscores
ourStates_species_data$species <- gsub("\\.", "_",
                                       ourStates_species_data$species)

# read in pest/host db
# North American hosts
pest_NA <- read.csv("./indata/pest_database/Data_Metadata/CSV/Insect x NA Host.csv")

pest_NA$NA_Host <- gsub("\\s", "_", pest_NA$NA_Host)

# non-North American hosts
pest_N <- read.csv("./indata/pest_database/Data_Metadata/CSV/Insect x N Host.csv")

pest_N$N_Host <- gsub("\\s", "_", pest_N$N_Host)

# combine the two hosts lists
tmp1 <- pest_NA[1:3]; names(tmp1)[2] <- "Host"
tmp2 <- pest_N; names(tmp2)[2] <- "Host"
pest <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)
pest$Host_Genus <- gsub("_.*", "", pest$Host)

# identify list of all pests of trees found in our focal states
pest_ourStates <- subset(pest, Host %in% unique(ourStates_species_data$species))

# read in plant megaphylogeny
phy <- read.tree("./indata/PhytoPhylo.tre")

# id host species missing from phylogeny
# setdiff(unique(pest_PA$Host), phy$tip.label)
# [1] "Quercus_coccinea" "Salix_nigra"

# create list of all pests/hosts present in our states and phylogeny
pest_ourStates_phy <- pest_ourStates[-which(!pest_ourStates$Host%in%phy$tip.label),]

# expand to matrix
matrix_host <- table(pest_ourStates_phy[,1:2])

# sort by number of hosts
matrix_host_sorted <- matrix_host[order(rowSums(matrix_host),decreasing=T),]

hist(rowSums(matrix_host)) # not many generalists

# choose top pest for the analysis
fpest <- row.names(matrix_host_sorted)[2]

# get host list
fpest_hosts_ourStates_phy <- subset(pest_ourStates_phy, Insect == fpest)

# determine list of all global hosts of focal pest
fpest_hosts_all <- subset(pest, Insect == fpest)

# check intersection with phy
length(fpest_hosts_all$Host) # 18 global hosts
length(intersect(fpest_hosts_all$Host, phy$tip.label)) # 16 in phy

# prune phylogeny to only host species in pest db
phy_host <- keep.tip(phy, intersect(unique(pest$Host), phy$tip.label))
length(phy_host$tip.label) # 288

# get all pairwise distances between tips on host phylogeny
pd_all <- get_all_pairwise_distances(phy_host,
                                    only_clades = 1:length(phy_host$tip.label))

# prep data frame for loop output
p_suit_all <- data.frame(tip = phy_host$tip.label)

# loop to get predicted probabilities for each tip for all focal hosts
for(i in seq_along(fpest_hosts_ourStates_phy$Host)){

# choose one host to be 'focal host' for PD calculations
fhost <- fpest_hosts_ourStates_phy$Host[i]

# subset to only distances to/from focal host
pd_fpest <- pd_all[which(phy_host$tip.label==fhost),]

# assemble data for logistic regression
fpest_regdata <- data.frame(tip = phy_host$tip.label,
                            pd = pd_fpest,
                            hostStatus = 0)

fpest_regdata[which(fpest_regdata$tip %in% fpest_hosts_all$Host),]$hostStatus <- 1

fpest_regdata$log_pd <- log10(fpest_regdata$pd + 1)

# logistic regression
fpest_logfit <- glm(hostStatus ~ log_pd,
                    family = binomial(link = "logit"),
                    data = fpest_regdata)

# summary(fpest_logfit)

# get predicted probabilities of association for all hosts in phy
fpest_pred <- fpest_regdata
fpest_pred$p_suit <- predict(fpest_logfit, type="response")

# hist(fpest_pred$p_suit)

p_suit_all[,ncol(p_suit_all)+1] <- fpest_pred$p_suit

names(p_suit_all)[ncol(p_suit_all)] <- paste0("p_suit_", i) 

# # plot phylogeny with tips coloured by predicted probability
fpest_pred$p_suit_round <- round(fpest_pred$p_suit, 1)
 
phy_plot <- extract.clade(phy_host, getMRCA(phy_host, fpest_pred$tip[which(fpest_pred$p_suit_round!=0.0)]))
# 
# tip_colours <- rep("black", Ntip(phy_plot))
# 
# cols <- heat.colors(10)
# 
# for(i in 1:Ntip(phy_plot)){
#     tip_colours[i] <- cols[11-fpest_pred$p_suit_round[which(fpest_pred$p_suit_round!=0.0)][i]*10]
# }
# 
# plot(phy_plot, tip.color=tip_colours)

# plot phylogeny with edges coloured by predicted probability
# create palette of 10 colours
cols <- heat.colors(10)

# create edge colour object, default colour is black
edge_colours <- rep("black", Nedge(phy_plot))

# give terminal edges a colour dependent on their predicted suitability value
for(j in 1:Nedge(phy_plot)){
    edge_colours[which.edge(phy_plot, phy_plot$tip.label[j])] <- cols[11-(fpest_pred[which(fpest_pred$tip==phy_plot$tip.label[j]),6]*10)]
}

# colour edge leading to focal host blue
edge_colours[which.edge(phy_plot, fhost)] <- "blue"

# plot
png(paste0("output/PD_fhost_", fpest_hosts_ourStates_phy$Host[i], ".png"))
plot(phy_plot, edge.color=edge_colours)
dev.off()

}

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








pothosts <- fpest_pred$tip[which(fpest_pred$p_suit_round!=0.0)] #29 spp

PA_pothosts <- intersect(pothosts, PA_trees$species)

species_imp_current$Quercus.rubra %>%
    tm_shape() +
    tm_raster() +
    tm_layout(legend.outside = TRUE)

plot(species_imp_current$Quercus.incana)
plot(species_imp_current$Quercus.rubra)
plot(species_imp_current$Quercus.bicolor)
plot(species_imp_current$Quercus.michauxii)
plot(species_imp_current$Quercus.macrocarpa)
plot(species_imp_current$Quercus.alba)
