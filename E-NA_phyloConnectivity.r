# LJ started 2021-11-10
# script for MSc plant symbiont connectivity simulation

# load packages:
library(dplyr)
library(sf)
library(raster)
library(tmap)

# set variable parameters:
    # states to include
    ourStates <- c("Connecticut")

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
map_species_tmp <- sf::st_read("./indata/USDA_data/sp318_hybrid_2100.shp")

# get CRS
crs_speciesMaps <- st_crs(map_species_tmp)

# re-project state(s) to the same coordinate system as tree data
map_ourStates_alb <- sf::st_transform(map_ourStates, crs_speciesMaps)

# clip species raster to extent of state(s)
map_species_ourStates_tmp <- map_species_tmp[map_ourStates_alb, ]

# get extent of clipped raster
extent_ourStates <- raster::extent(map_species_ourStates_tmp)

# remove unneeded maps
rm(c(map_species_tmp, map_species_ourStates_tmp))


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

# We'll start with the first species in the list to establish the first layer in a raster stack, then we'll loop through remaining species to add to the stack.

# **NOTE** we'll create 2 different raster stacks: one for current species importance values (RFimp), and one for the averaged prediction for three general circulation models for future conditions "GCM45".

##### comments consistent to here #####


## replace with reference to earlier demo species map
# # import file
# species_map1 <- sf::st_read(species_filenames[1])
# 
# # clip to Pennsylvania
# ourStates_species1 <- species_map1[map_ourStates_alb, ]

# rasterize:
species_imp_current <- raster::rasterize(x = ourStates_species1, 
                                         y = ourStates_reference_raster, field = "RFimp")
species_imp_future <- raster::rasterize(x = ourStates_species1, 
                                        y = ourStates_reference_raster, field = "GCM45")

# Start big loop chunk:

#options(width = 60)
#for (i in 2:4){
for (i in 2:length(species_filenames)) {
    
    ## conflicting names here need to be sorted out
    # import file
    species_map <- sf::st_read(species_filenames[i])
    
    # clip to Pennsylvania
    ourStates_species <- species_map[map_ourStates_alb, ]
    
    # rasterize:
    species_imp_current <- raster::stack(species_imp_current, 
                                         raster::rasterize(x = ourStates_species,
                                                           y = ourStates_reference_raster, field = "RFimp"))
    species_imp_future <- raster::stack(species_imp_future,
                                        raster::rasterize(x = ourStates_species,
                                                          y = ourStates_reference_raster, field = "GCM45"))
    print(i)
}

# Name each layer in the stack by its latin name:
names(species_imp_current) <- species_info_toUse$Scientific_name
names(species_imp_future) <- species_info_toUse$Scientific_name

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

# Now we can combine these into long format, keeping only the species present in a given cell. 

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



# We these items we can do the following procedure:

# * get species richness per grid cell

ourStates_species_richness <- ourStates_species_data %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(n_species  = n()) %>%
    dplyr::left_join(ourStates_xy_coords, by = "cell_id") %>%
    dplyr::select(cell_id, x, y, n_species)

# Look at result:

ourStates_species_richness


# read in list of PA trees
PA_trees <- ourStates_species_data

# replace spaces in species names with underscores
PA_trees$species <- gsub("\\.", "_", PA_trees$species)

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

# identify list of all pests of PA trees
pest_PA <- subset(pest, Host %in% unique(PA_trees$species))

# read in plant megaphylogeny
library(ape)
phy <- read.tree("./indata/PhytoPhylo.tre")

# id host species missing from phylogeny
# setdiff(unique(pest_PA$Host), phy$tip.label)
# [1] "Quercus_coccinea" "Salix_nigra"

# create list of all pests/hosts present in PA and phylogeny
pest_PA_phy <- pest_PA[-which(!pest_PA$Host%in%phy$tip.label),]

# choose one pest for the analysis
set.seed(2021)
fpest <- sample(pest_PA_phy$Insect, 1)

# get host list
fpest_PA_phy <- subset(pest_PA_phy, Insect == fpest)

# choose one host to be 'focal host' for PD calculations
fhost <- sample(fpest_PA_phy$Host, 1)

# determine list of all global hosts of focal pest
fpest_all <- subset(pest, Insect == fpest)

# check intersection with phy
length(fpest_all$Host) # 5 global hosts
length(intersect(fpest_all$Host, phy$tip.label)) # 5 - all present

# prune phylogeny to only host species in pest db
length(unique(pest$Host)) #416 tree species in db
length(intersect(unique(pest$Host), phy$tip.label)) #288 - many missing from phy

phy_host <- keep.tip(phy, intersect(unique(pest$Host), phy$tip.label))
length(phy_host$tip.label) #288

# get all pairwise distances between tips on host phylogeny
library(castor)

pd_all <- get_all_pairwise_distances(phy_host,
                                    only_clades = 1:length(phy_host$tip.label))

# subset to only distances to/from focal host
pd_fpest <- pd_all[which(phy_host$tip.label==fhost),]

# assemble data for logistic regression
fpest_regdata <- data.frame(tip = phy_host$tip.label,
                            pd = pd_fpest,
                            hostStatus = 0)

fpest_regdata[which(fpest_regdata$tip %in% fpest_all$Host),]$hostStatus <- 1

# logistic regression
fpest_logfit <- glm(hostStatus ~ pd,
                    family = binomial(link = "logit"),
                    data = fpest_regdata)

summary(fpest_logfit)

### runs, but non-significant and data is super zero-inflated... will have to figure this bit out

# get predicted probabilities of association for all hosts in phy
fpest_pred <- fpest_regdata
fpest_pred$p_suit <- predict(fpest_logfit, type="response")

# assign each PA cell a suitability value equal to the max association probability of any tree species present in the cell

PA_trees_suit <- left_join(PA_trees, fpest_pred, by=c("species" = "tip"))

PA_trees_suit_summary <- PA_trees_suit %>%
                            group_by(cell_id) %>%
                            summarise(p_suit = max(p_suit))

tmp3 <- PA_trees_suit %>%
            select(c(cell_id, x, y))
tmp3 <- unique(tmp3)

PA_trees_suit_summary <- left_join(PA_trees_suit_summary, tmp3)

PA_trees_suit_summary <- PA_trees_suit_summary[,c(3,4,2)]

PA_trees_suit_summary[which(is.na(PA_trees_suit_summary$p_suit)),3] <- 0

library(raster)

PA_ref <- ourStates_reference_raster

d <- rasterize(x=PA_trees_suit_summary[,c(1:2)], y=PA_ref,
               field= PA_trees_suit_summary$p_suit)

library(tmap)

png("output/PA_testSuit.png",  width = 580, height = 480)
d %>%
    tm_shape() +
    tm_raster(alpha = 0.7) +
    tm_layout(legend.outside = TRUE)
dev.off()

