# backup script to get species list per cell
# follows Jason's 'mapping.Rmd'

# load packages:
library(dplyr)
library(sf)
library(raster)
library(tmap)


# set variable parameters:
# states to include
ourStates <- c("Alabama", "Georgia", "Mississippi", "South Carolina",
               "North Carolina", "Tennessee", "Kentucky", "Virginia",
               "West Virginia", "Indiana", "Louisiana", "Arkansas", "Missouri",
               "Illinois", "Ohio", "Pennsylvania", "Maryland", "Indiana",
               "Iowa", "Delaware", "New Jersey")

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
map_species_demo <- sf::st_read("./indata/OLD/USDA_data/sp318_hybrid_2100.shp")

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
data_species <- read.csv("./indata/OLD/USDA_data/_Species_List.csv")

# rename variables to remove spaces
names(data_species) <- c("FIA_species_code", "Common_name", "Scientific_name", "Model_reliability")

# we'll loop through each row of the file for which "Model_reliability" is not equal to "FIA Only"
data_species_toUse <- data_species %>%
    filter(Model_reliability != "FIA Only")

# Let's now create a vector of filenames and paths, and this is what we'll loop through:
filenames_species_toUse <- paste(paste("./indata/OLD/USDA_data/sp", 
                                       data_species_toUse$FIA_species_code,
                                       sep = ""),
                                 "hybrid_2100.shp", sep = "_")

# We'll start with the first species in the list to establish the first layer
# in a raster stack, then we'll loop through remaining species to add to the
# stack.

# **NOTE** we'll create 2 different raster stacks: one for current species
# importance values (RFimp), and one for the averaged prediction for three
# general circulation models for future conditions "GCM45".

map_species <- sf::st_read(filenames_species_toUse[1])

map_species_ourStates <- map_species[map_ourStates_alb, ]

# rasterize:
species_imp_current <- raster::rasterize(x = map_species_ourStates, 
                                         y = map_ourStates_refRaster,
                                         field = "RFimp")
# species_imp_future <- raster::rasterize(x = map_species_ourStates_demo, 
#                                         y = map_ourStates_refRaster,
#                                         field = "GCM45")

# Start big loop
for (i in 2:length(filenames_species_toUse)) {
    
    ## conflicting names here need to be sorted out
    # import file
    map_species <- sf::st_read(filenames_species_toUse[i])
    
    # clip to chosen state(s)
    map_species_ourStates <- map_species[map_ourStates_alb, ]
    
    # rasterize:
    species_imp_current <- raster::stack(
                            species_imp_current,
                            raster::rasterize(x = map_species_ourStates,
                                              y = map_ourStates_refRaster,
                                              field = "RFimp"))
    
    # species_imp_future <- raster::stack(
    #                         species_imp_future,
    #                         raster::rasterize(x = map_species_ourStates,
    #                                           y = map_ourStates_refRaster,
    #                                           field = "GCM45"))
    
    print(i)
}

# Name each layer in the stack by its latin name:
names(species_imp_current) <- data_species_toUse$Scientific_name
# names(species_imp_future) <- data_species_toUse$Scientific_name

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
ourStates_species_data

# write stack to file
raster::writeRaster(species_imp_current, 
                    filename = "./output/backup_current_IVs_2021-12-09.tif", 
                    format = "GTiff", overwrite = TRUE)

# write species list to file
write.csv(ourStates_species_data, "./output/backup_species_data_2021-12-09.csv")


