# LJ 2021-12-09 prepare phylogeny, adding tips where needed


# load packages
library(ape)
library(phytools)
library(V.PhyloMaker)
# install with:
# devtools::install_github("jinyizju/V.PhyloMaker")
library(Taxonstand)


# load species lists
# FIA trees
sp_tr <- read.csv("./output/FIAspp_toUse_2021-12-15T15:33:29.csv")$species

# host trees from pest database
sp_pest  <- read.csv("./output/pestDB-hostList_clean_2021-12-15T15:15:34.csv")$host

# Reassign phylogeny for convenience (loaded with 'V.PhyloMaker' package)
phy <- GBOTB.extended

# tip list
sp_phy <- phy$tip.label


# First step is to run a quick check that the phylogeny is as expected.
# The authors say that the tip list is standardized according to TPL; there are
# 74531 tips, so checking them all individually is inefficient. Instead, we'll
# sample a subset of 500 tips (0.7%) and check them against TPL.
# 
# tip_sample <- sample(GBOTB.extended$tip.label, 500)
# 
# tip_sample_TPL <- Taxonstand::TPL(gsub("_", " ", tip_sample),
#                 corr = TRUE,
#                 diffchar = 2,
#                 max.distance = 1,
#                 version = "1.1",
#                 encoding = "UTF-8",
#                 author = TRUE,
#                 drop.lower.level = FALSE,
#                 file = "",
#                 silent = TRUE,
#                 repeats = 6)
# 
# table(tip_sample_TPL$Taxonomic.status)
# 
# tip_sample_TPL[which(tip_sample_TPL$Taxonomic.status == "Synonym"),]
# 
# # Yup... not 100% 'accepted'. A full correction is out of scope at the
# moment, but we'll keep this problem in mind.


# check intersection of species lists
length(setdiff(sp_pest, sp_phy)) # 90 hosts missing from phy
length(setdiff(sp_tr, sp_phy)) # 46 FIA trees missing from phy


# FIA tree spp are key, so we'll focus on those
sp_tr_mis <- setdiff(sp_tr, sp_phy)

# check against TPL
sp_tr_mis_tpl <- Taxonstand::TPL(gsub("_", " ", sp_tr_mis),
                                 corr = TRUE,
                                 diffchar = 2,
                                 max.distance = 1,
                                 version = "1.1",
                                 encoding = "UTF-8",
                                 author = TRUE,
                                 drop.lower.level = FALSE,
                                 file = "",
                                 silent = TRUE,
                                 repeats = 6)

# get list of updated/corrected names
sp_tr_mis_tpl2 <- paste(sp_tr_mis_tpl$New.Genus,
                        sp_tr_mis_tpl$New.Species,
                        sep = "_")

# check new intersection
length(setdiff(sp_tr_mis_tpl2, sp_phy)) # 34 still missing
# we'll have to add these tips

# list of species to add
sp_tr_add <- setdiff(sp_tr_mis_tpl2, sp_phy)

# ideally, a congener for each of these species will be present in the tree, and we will add the tips sister to their congeners

# check for species whose genera are not present in the tree
noCgen <- sp_tr_add[which(!gsub("_.*", "", sp_tr_add) %in%
                              unique(gsub("_.*", "", sp_phy)))]

# only one species affected - it will have to be added separately, at the base of its family

# remove from list of tips to add
sp_tr_add <- sp_tr_add[-which(sp_tr_add==noCgen)]


# get list of tips to keep in phy
# overlap with pest
t1 <- intersect(sp_pest, sp_phy)

# native overlap with FIA
t2 <- intersect(sp_tr, sp_phy)

# corrected overlap with FIA
t3 <- intersect(sp_tr_mis_tpl2, sp_phy)

# congener from each of missing FIA genera
cGen <- setdiff(unique(gsub("_.*", "", sp_tr_add)),
                 unique(c(gsub("_.*", "", t2), gsub("_.*", "", t3))))

t4 <- vector(length=length(cGen))

for(i in seq_along(cGen)){
    t4[i] <- sp_phy[grep(cGen[i], sp_phy)[1]]
}

# combine
finalTips <- c(t1, t2, t3, t4, noCgen)


# add no-congener tip before pruning
# get family
fam <- subset(sp_tr_mis_tpl, New.Genus == gsub("_.*", "", noCgen))$Family

# ugh, family is "Leguminosae", AKA "Fabaceae" - rename
fam <-  "Fabaceae"

# get mrca node
fnode <- getMRCA(phy, tips.info[which(tips.info$family == fam), "species"])

# add tip
phy <- bind.tip(phy,
                tip.label = noCgen,
                where = fnode,
                position = 0
)


# prune phy
ourPhy <- keep.tip(phy, finalTips)


# add missing tips with congeners
# added as polytomies at base of genus
# if only one congeneric tip is present, new mrca node is created at the 1/2 point along the pendant edge

for(i in seq_along(sp_tr_add)){
    
    node <- getMRCA(ourPhy,
                    ourPhy$tip.label[grep(strsplit(sp_tr_add[i], "_")[[1]][1],
                                          ourPhy$tip.label)]
                    )
    
    if(!is.null(node)){
        
        ourPhy <- bind.tip(ourPhy,
                           tip.label = sp_tr_add[i],
                           where = node,
                           position = 0
                           )
    }else{
        
        ourPhy <- bind.tip(ourPhy,
                           tip.label = sp_tr_add[i],
                           
                           where = grep(strsplit(sp_tr_add[i], "_")[[1]][1],
                                        ourPhy$tip.label),
                           
                           position = ourPhy$edge.length[
                               which(ourPhy$edge[,2] == grep(
                                   strsplit(sp_tr_add[i], "_")[[1]][1],
                                   ourPhy$tip.label))] / 2
                           )
    }
}


# check plot
# tree <- ourPhy
# sister <- noCgen #"Magnolia_montana"
# plot(keep.tip(tree, c(c(which(tree$tip.label==sister)-2):c(which(tree$tip.label==sister)+2))[which(c(c(which(tree$tip.label==sister)-2):c(which(tree$tip.label==sister)+2))>0)]))


# remove t4 tips (congeners used to locate tip addition, no longer needed)
ourPhy <- drop.tip(ourPhy, t4)


# write to file
write.tree(ourPhy,
           paste0("./output/focalPhy_",
                  gsub(" ", "T", Sys.time()),
                  ".nwk"))




###############################################################################
# make key to connect FIA species with tips in cases where we renamed tips

## this needs cleaning up

# library(dplyr)
# 
# fia_tax <- read.csv("./indata/studyregion_imputed_FIA_species_by_plot_summary_good_taxonomy.csv")
# 
# new_tax <- subset(sp_tr_mis_tpl, Taxon %in% gsub("_", " ", setdiff(sp_tr, ourPhy$tip.label)), select = c("Taxon", "New.Genus", "New.Species"))
# 
# all_tax <- left_join(new_tax, fia_tax, by = c("Taxon" = "UpdatedFIABinomial"))
# 
# all_tax <- unique(all_tax[,c("Taxon", "New.Genus", "New.Species", "SPCD")])
# 
# all_tax$tip_name <- paste(all_tax$New.Genus, all_tax$New.Species, sep="_")
# 
# all_tax <- all_tax[,c("Taxon", "tip_name", "SPCD")]
# 
# names(all_tax)[1] <- "UpdatedFIABinomial"
# 
# rest_tax <- subset(fia_tax, !SPCD%in%all_tax$SPCD)
# 
# rest_tax <- unique(rest_tax[,c("UpdatedFIABinomial", "SPCD")])
# 
# rest_tax$tip_name <- gsub(" ", "_", rest_tax$UpdatedFIABinomial)
# 
# rest_tax[grep(" x ", rest_tax$UpdatedFIABinomial), "tip_name"] <- NA
# 
# tax <- rbind(all_tax, rest_tax)
