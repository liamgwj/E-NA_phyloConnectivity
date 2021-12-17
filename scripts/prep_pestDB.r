# LJ started 2021-11-28 updated 2021-12-06
# clean up pest/host interaction database and compile clean host species list

# 1st output: data.frame with two columns named "HostPlant" and "Symbiont",
# each containing latin names in format 'Genus_species_var._ssp., where every
# row indicates a compatible host/symbiont pairing and both host and symbiont
# names are potentially repeated

# load checkNames function
source("/home/liam/Documents/MSc/analysis/misc_blerfLab/checkNames_function.r")

# read in pest/host db
# North American hosts
pest_NA <- read.csv("./indata/pest_database/Data_Metadata/CSV/Insect x NA Host.csv")

# non-North American hosts
pest_N <- read.csv("./indata/pest_database/Data_Metadata/CSV/Insect x N Host.csv")

# check
checkNames(pest_NA$NA_Host) # many duplicates, 7 subspecies, 1 hyphen
checkNames(pest_N$N_Host) # error - bad character
pest_N$N_Host[grep("\xd7", pest_N$N_Host)] # 3 unrecognized symbol "\xd7" ('×')
# replace with 'x', adding one space
pest_N$N_Host <- gsub("\xd7europaea", "x europaea", pest_N$N_Host)
pest_N$N_Host <- gsub("\xd7", "x", pest_N$N_Host)
checkNames(pest_N$N_Host) # 2 varietals, 19 subspecies, 5 crosses, 1 hyphen

# combine the two hosts lists
pest_NA <- pest_NA[1:2]; names(pest_NA) <- c("Symbiont", "HostPlant")
pest_N <- pest_N[1:2]; names(pest_N) <- c("Symbiont", "HostPlant")
pest <- rbind(pest_NA, pest_N)

# replace spaces in names with underscores
pest$Symbiont <- gsub("\\s", "_", pest$Symbiont)
pest$HostPlant <- gsub("\\s", "_", pest$HostPlant)

checkNames(pest$HostPlant)
# checkNames(pest$Symbiont)

### for now, we'll remove host subspecies and varietals
# remove subspecies
pest[grep("subsp", pest$HostPlant),2] <- gsub("_subsp\\._.*", "", pest[grep("subsp", pest$HostPlant),2])

# remove var
pest[grep("var\\.", pest$HostPlant),2] <- gsub("_var\\._.*", "", pest[grep("var\\.", pest$HostPlant),2])

# unique cases
pest <- unique(pest)

checkNames(pest$HostPlant)
checkNames(pest$HostPlant, findWhich = "hyphens")


# write to file
write.csv(pest,
          paste0("./output/pest-host_clean_noSSp_",
                 gsub(" ", "T", Sys.time()),
                 ".csv"),
          row.names = FALSE)

# 2nd output - host species list
host <- unique(pest$HostPlant)

checkNames(host)

# remove crosses
host <- host[-checkNames(host,
                     findWhich = "x.X",
                     type = "index")]

# write to file
write.csv(data.frame(host = host),
          paste0("./output/pestDB-hostList_clean_",
                 gsub(" ", "T", Sys.time()),
                 ".csv"),
          row.names = FALSE)
