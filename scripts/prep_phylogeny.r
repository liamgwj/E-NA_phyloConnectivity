# LJ 2021-12-09 prepare phylogeny, adding tips where needed


# load packages
library(ape)
library(phytools)
library(V.PhyloMaker)
# install with:
# devtools::install_github("jinyizju/V.PhyloMaker")
library(Taxonstand)


# define convenience functions
source("/home/liam/Documents/MSc/analysis/misc_blerfLab/checkNames_function.r")


# load species lists
# FIA trees
sp_tr <- read.csv("./output/FIAspp_toUse_2021-12-15T15:33:29.csv")$species

# host trees from pest database
sp_pest  <- read.csv("./output/pestDB-hostList_clean_2021-12-15T15:15:34.csv")$host

# Reassign phylogeny for convenience
phy <- GBOTB.extended

sp_phy <- phy$tip.label


# # First step is to run a quick check that the phylogeny is as expected. The authors say that the tip list is standardized according to TPL; there are  74531 tips, so checking them all individually is inefficient. Instead, we'll sample a subset of 500 tips (0.7%) and check them against TPL.
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
# # Yup... not 100% 'accepted'. A full correction is out of scope at the moment, but we'll keep this problem in mind.


# # check overlap of species lists
# length(unique(c(sp_tr, sp_pest))) #560 tree species
# 
# length(intersect(unique(c(sp_tr, sp_pest)),
#                  sp_phy)) # 434 in phylogeny
# 
# setdiff(unique(c(sp_tr, sp_pest)),
#         sp_phy)
# 
# 
# length(sp_tr) #185
# length(intersect(sp_tr, sp_phy)) #231    46 missing
# 
# length(sp_pest) #405
# length(intersect(sp_pest, sp_phy)) #315    90 missing


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

# updated names
sp_tr_mis_tpl2 <- paste(sp_tr_mis_tpl$New.Genus,
                        sp_tr_mis_tpl$New.Species,
                        sep = "_")

# species to add
sp_tr_add <- setdiff(sp_tr_mis_tpl2, sp_phy)

# check if genera are present
gen_tr_add <- unique(gsub("_.*", "", sp_tr_add))

# d<-vector(length=17)
# for(i in 1:17){
# d[i] <- length(grep(gen_tr_add[i], sp_phy))>=1
# }
# 
# grep(gen_tr_add[14], sp_phy)
# 
# gen_tr_add[14]
# all genera in phylogeny except gymnocladus
sp_tr_add <- sp_tr_add[-which(sp_tr_add=="Gymnocladus_dioica")]



# list of tips to keep in phy
# overlap with pest
t1 <- intersect(sp_pest, sp_phy)

# native overlap with FIA
t2 <- intersect(sp_tr, sp_phy)

# corrected overlap with phy
t3 <- intersect(sp_tr_mis_tpl2, sp_phy)

# congener from each of missing FIA genera
ncgen <- setdiff(gen_tr_add,
                 unique(c(gsub("_.*", "", t2), gsub("_.*", "", t3))))

#remove gymnocladus
ncgen <- ncgen[-4]

t4 <- vector(length=length(ncgen))

for(i in seq_along(ncgen)){
    t4[i] <- sp_phy[grep(ncgen[i], sp_phy)[1]]
}

# combine
finalTips <- c(t1, t2, t3, t4)

# prune phy
ourPhy <- keep.tip(phy, finalTips)


# add missing tips
# added as polytoies at base of genus
# if only one congeneric tip is present, new mrca node is created at the 1/2 point along the pendant edge

for(i in seq_along(sp_tr_add)){
    
    node <- getMRCA(ourPhy,
                    ourPhy$tip.label[grep(strsplit(sp_tr_add[i], "_")[[1]][1],
                                          ourPhy$tip.label)])
    
    if(!is.null(node)){
    
ourPhy<-    bind.tip(ourPhy, tip.label = sp_tr_add[i],
             where = node,
             position = 0
             )
}

if(is.null(node)){
    
    ourPhy<-    bind.tip(ourPhy, tip.label = sp_tr_add[i],
                         where = grep(strsplit(sp_tr_add[i], "_")[[1]][1],
                                      ourPhy$tip.label),
                         position = ourPhy$edge.length[which(ourPhy$edge[,2]==grep(strsplit(sp_tr_add[i], "_")[[1]][1], ourPhy$tip.label))]/2
    )
}
}


# # check plot
# tree <- ourPhy
# sister <- "Magnolia_montana" 
# plot(keep.tip(tree, c(c(which(tree$tip.label==sister)-2):c(which(tree$tip.label==sister)+2))[which(c(c(which(tree$tip.label==sister)-2):c(which(tree$tip.label==sister)+2))>0)]))
# 



# add one remaining tip?


# drop t4 tips
ourPhy <- drop.tip(ourPhy, t4)


## this leads to multi-naming issues! will change names of FIA layers instead
# # change names to match FIA key
# key <- subset(sp_tr_mis_tpl, Taxon %in% gsub("_", " ", setdiff(sp_tr, ourPhy$tip.label)))
# 
# key$newname <- paste(key$New.Genus, key$New.Species, sep = "_")
# 
# key <- key[,c("Taxon", "newname")]
# 
# key$Taxon <- gsub(" ", "_", key$Taxon)
# 
# for(i in seq_along(ourPhy$tip.label)){
#     if(ourPhy$tip.label[i]%in%key$newname){
#         ourPhy$tip.label[i] <- key$Taxon[which(key$newname==ourPhy$tip.label[i])]
#     }
# }


# write to file
write.tree(ourPhy,
           paste0("./output/focalPhy_",
                  gsub(" ", "T", Sys.time()),
                  ".nwk"))
