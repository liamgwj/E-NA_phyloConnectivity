# LJ 2021-12-15 prep FIA tree species list


# define convenience functions
source("/home/liam/Documents/MSc/analysis/misc_blerfLab/checkNames_function.r")


# read in FIA code list
tr_codes <- read.csv("indata/final_SPCD_codes.csv")
# code key
tr_key <- read.csv("indata/FIA_spp_tmp/REF_SPECIES.csv")


# subset to focal species
tr <- subset(tr_key, SPCD %in% tr_codes$SPCD)

# select relevant columns
tr <- tr[,c("SPCD", "GENUS", "SPECIES", "VARIETY", "SUBSPECIES", "SPECIES_SYMBOL")]

# create combined name column
tr$latbi <- paste(tr$GENUS, tr$SPECIES, sep=" ")

for(i in seq_along(tr$latbi)){
    if(!tr$VARIETY[i]==""){
        tr$latbi[i] <- paste(tr$latbi[i], "var.", tr$VARIETY[i], sep=" ")
    }
    if(!tr$SUBSPECIES[i]==""){
        tr$latbi[i] <- paste(tr$latbi[i], "subsp.", tr$SUBSPECIES[i], sep=" ")
    }
}

checkNames(tr$latbi)
# yuck... 8 varietals, 21 "spp", 1 subsp, 1 cross, 1 family (nonFirst cap)

# these might be taken out in Corrina's list?

tr_key_c <- read.csv("./indata/riley_corrina_merged.csv")

tr <- subset(tr_key_c, SPCD %in% tr_codes$SPCD)

# select relevant columns
tr <- tr[,c("SPCD", "GENUS", "SPECIES", "VARIETY", "SUBSPECIES", "SPECIES_SYMBOL")]

# create combined name column
tr$latbi <- paste(tr$GENUS, tr$SPECIES, sep=" ")

for(i in seq_along(tr$latbi)){
    if(!tr$VARIETY[i]==""){
        tr$latbi[i] <- paste(tr$latbi[i], "var.", tr$VARIETY[i], sep=" ")
    }
    if(!tr$SUBSPECIES[i]==""){
        tr$latbi[i] <- paste(tr$latbi[i], "subsp.", tr$SUBSPECIES[i], sep=" ")
    }
}

checkNames(tr$latbi)
# spp and family removed - we'll use this

checkNames(tr$latbi, findWhich = "duplicated")
# unclear why these are duplicated

# there are still 8 varietals and 1 subspecies... for now we'll keep them

# however, remove the cross (since crosses are fundamentally at odds with the notion of a phylogeny as a "graph without cycles"!)
tr <- tr[-checkNames(tr$latbi,
                     findWhich = "x.X",
                     type = "index"),]

# species list in tip format
sp_tr <- gsub(" ", "_", unique(tr$latbi))

checkNames(sp_tr)
# everything seems ok


# write to file
write.csv(data.frame(species = sp_tr),
          paste0("./output/FIAspp_toUse_",
                 gsub(" ", "T", Sys.time()),
                 ".csv"),
          row.names = FALSE)

##############################################################################


# read in new list from Jason
tr <- read.csv("./indata/studyregion_imputed_FIA_species_by_plot_summary_good_taxonomy.csv")

sp_tr <- unique(tr$UpdatedFIABinomial)

sp_tr <- sp_tr[-checkNames(sp_tr,
                     findWhich = "x.X",
                     type = "index")]

sp_tr <- gsub(" ", "_", sp_tr)

checkNames(sp_tr)

write.csv(data.frame(species = sp_tr),
          paste0("./output/FIAspp_toUse_",
                 gsub(" ", "T", Sys.time()),
                 ".csv"),
          row.names = FALSE)
