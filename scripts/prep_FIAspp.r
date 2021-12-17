# LJ 2021-12-15 prep FIA tree species list


# define convenience functions
source("/home/liam/Documents/MSc/analysis/misc_blerfLab/checkNames_function.r")


# read in new list from Jason
tr <- read.csv("./indata/studyregion_imputed_FIA_species_by_plot_summary_good_taxonomy.csv")


# get vector of names
sp_tr <- unique(tr$UpdatedFIABinomial)

# remove the cross (since crosses are fundamentally at odds with the notion of a phylogeny as a "graph without cycles"!)
sp_tr <- sp_tr[-checkNames(sp_tr,
                     findWhich = "x.X",
                     type = "index")]

# replace spaces with underscores
sp_tr <- gsub(" ", "_", sp_tr)

# check: ok
checkNames(sp_tr)

# write to file
write.csv(data.frame(species = sp_tr),
          paste0("./output/FIAspp_toUse_",
                 gsub(" ", "T", Sys.time()),
                 ".csv"),
          row.names = FALSE)
