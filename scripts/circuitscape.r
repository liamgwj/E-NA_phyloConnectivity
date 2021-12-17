# Circuitscape connectivity analysis
# LJ 2021-11-26

# input data format
# conductance surface:
# a raster with cell values equal to modelled suitability for the focal pest
# source/ground node file:
# a raster with the same extent as the conductance surface, with focal nodes located at uniform intervals around the perimeter of the map


# # OLD: create source nodes [this will now be a separate routine]
# nodes <- raster::raster(matrix(-9999,
#                                nrow(occurrence[[1]]),
#                                ncol(occurrence[[1]])),
#                         xmn = 0, xmx = nrow(occurrence[[1]]),
#                         ymn = 0, ymx = ncol(occurrence[[1]]))
# 
# nodes@file@nodatavalue <- -9999
# 
# nodes[seq(from = 1, to = nrow(occurrence[[1]]), length.out = 10), 1] <- 1:10
# 
# nodes[seq(from = 1, to = nrow(occurrence[[1]]), length.out = 10), ncol(occurrence[[1]])] <- 11:20
# 
# # write to file
# raster::writeRaster(nodes,
#                     paste0("output/", now, "/input-maps/nodes"),
#                     format = "ascii",
#                     overwrite = TRUE)


## create .ini file for circuitscape run --------------------------------------

# check for/create output directories
## will need to adjust for looping thru pest species
now <- gsub(" ", "T", Sys.time())
if(!dir.exists(paste0("output/circuitscape_", now))){
    dir.create(paste0("output/circuitscape_", now))}

# write .ini file
writeLines(
    c("[Circuitscape Mode]",
      "data_type = raster",
      "scenario = pairwise",
      
      "[Version]",
      "version = 5.0.0",
      
      "[Habitat raster or graph]",
      paste0("habitat_file = /home/liam/Documents/MSc/analysis/MSc_repo/input-maps/phy0_suitability.asc"),
      "habitat_map_is_resistances = resistances",
      
      "[Connection Scheme for raster habitat data]",
      "connect_four_neighbors_only = false",
      "connect_using_avg_resistances = false",
      
      "[Short circuit regions (aka polygons)]",
      "use_polygons = false",
      "polygon_file = False",
      
      "[Options for advanced mode]",
      "ground_file_is_resistances = true",
      "source_file = (Browse for a current source file)",
      "remove_src_or_gnd = keepall",
      "ground_file = (Browse for a ground point file)",
      "use_unit_currents = false",
      "use_direct_grounds = false",
      
      "[Mask file]",
      "use_mask = false",
      "mask_file = None",
      
      "[Options for one-to-all and all-to-one modes]",
      "use_variable_source_strengths = false",
      "variable_source_file = None",
      
      "[Options for pairwise and one-to-all and all-to-one modes]",
      "included_pairs_file = (Browse for a file with pairs to include or exclude)",
      "use_included_pairs = false",
      paste0("point_file = /home/liam/Documents/MSc/analysis/MSc_repo/output/input-maps/nodes.asc"),
      
      "[Calculation options]",
      "solver = cg+amg",
      
      "[Output options]",
      "write_cum_cur_map_only = false",
      "log_transform_maps = false",
      paste0("output_file = /home/liam/Documents/MSc/analysis/E-NA_phyloConnectivity/output/circuitscape_", now, "/out_", now),
      "write_max_cur_maps = false",
      "write_volt_maps = true",
      "set_null_currents_to_nodata = false",
      "set_null_voltages_to_nodata = false",
      "compress_grids = false",
      "write_cur_maps = true"
    ),
    con = paste0("output/circuitscape_", now, "/lastRun.ini"))


## run Circuitscape -----------------------------------------------------------

XRJulia::juliaUsing("Circuitscape")

XRJulia::juliaCommand(paste0("compute(\"/home/liam/Documents/MSc/analysis",
                             "/MSc_repo/output/lastRun.ini\")"))


## check output ---------------------------------------------------------------
# 
# suitPlot <- suitability
# suitPlot[which(suitPlot[]==-9999)] <- NA
# 
# png(paste0("output/", now, "/", now, "_resistance_d",
#            params$dispMax, ".png"),
#     width = 960,
#     height = 960)
# 
# # landscapetools::show_landscape(suitPlot)
# raster::plot(suitPlot, col=viridis::viridis(10))
# 
# dev.off()
# 
# d5_curmap <- raster::raster(paste0("output/", now, "/circuitscape-output/disp_", params$dispMax, "/",
#                                    now, "_cum_curmap.asc"))
