# LJ 2021-11-08 check intersection of three tree species lists:
# 1 USDA eastern USA tree range maps
# 2 QJ 2016 plant phylogeny
# 3 insect pest database of north american trees

# read in data
usda <- read.csv("./indata/USDA_data/_Species_List.csv")
pest <- read.csv("./indata/pest_database/Data_Metadata/CSV/Host Plant Divergence.csv")
phy <- ape::read.tree("./indata/PhytoPhylo.tre")

# take a look at each in turn, check lengths and modify names where needed
head(usda)
usda$latbi <- gsub("\\s", "_", usda$Scientific.Name)
length(usda$latbi) # 148
length(unique(usda$latbi)) # 148 - no duplicates

head(pest)
pest$NA_Host <- gsub("\\s", "_", pest$NA_Host)
pest$N_Host <- gsub("\\s", "_", pest$N_Host)
length(pest$NA_Host) # 285
length(unique(pest$NA_Host)) # 142
length(pest$N_Host) # 285 - odd that they're the same?
length(unique(pest$N_Host)) # 91
length(c(pest$NA_Host, pest$N_Host)) # 570
length(unique(c(pest$NA_Host, pest$N_Host))) # 233 (142 + 91 - all good)
pest_all <- unique(c(pest$NA_Host, pest$N_Host))

head(phy$tip.label)
length(phy$tip.label) # 31389
length(unique(phy$tip.label)) # 31389

# check intersections
length(intersect(usda$latbi, phy$tip.label)) # 122 - 26 USDA spp missing from phy

length(intersect(pest_all, phy$tip.label)) # 184 - 49 pest spp missing from phy
length(intersect(unique(pest$NA_Host), phy$tip.label)) # 112 - 30missing from phy
length(intersect(unique(pest$N_Host), phy$tip.label)) # 72 - 19 missing


length(intersect(usda$latbi, unique(pest$NA_Host))) # 50 - many missing...
length(intersect(usda$latbi, unique(pest$N_Host))) # 2 - but this is expected
length(intersect(usda$latbi, pest_all)) # 52


# similar numbers of spp missing from phy might reflect naming differences - we can standardize this against TPL (ugh)

# why so many missing from USDA lists / pest DB?

setdiff(usda$latbi, unique(pest$NA_Host))

setdiff(unique(pest$NA_Host), usda$latbi)

# these differences seem legit... so, means that many of the trees with eastern NA occurrence maps are not known hosts of the pests in the db, while many db pests infect trees that don't have maps (reasonable since many are western spp)

# for now, we'll focus on the species that are both known hosts and have maps

length(intersect(intersect(usda$latbi, unique(pest$NA_Host)), phy$tip.label)) # 47

setdiff(intersect(usda$latbi, unique(pest$NA_Host)), phy$tip.label)
# [1] "Quercus_coccinea" "Salix_nigra"      "Ulmus_thomasii"  

grep("Quercus_palustris", phy$tip.label) # synonym in tree
grep("Ulmus_racemosa", phy$tip.label) # not in tree

# so we're looking at 47 - 48 species to start

# get list of spp

omnipresent_spp <- data.frame(species = intersect(intersect(usda$latbi, unique(pest$NA_Host)), phy$tip.label))


write.csv(omnipresent_spp, "output/omnipresent_spp.csv", row.names=FALSE)



###############################################################################
## Below lies a mess, retained for now for useful pieces

# LJ 2021-11-08 
# read in intersecting subset of host list
hosts <- read.csv("multipresent_spp_list.csv")
names(hosts) <- "latbi"

# read in list of PA hosts
PA_hosts <- read.csv("/home/liam/Documents/MSc/analysis/usda_mapping/outdata/penn_species_data.csv")
PA_hosts$species <- gsub("\\.", "_", PA_hosts$species)

length(unique(PA_hosts$species)) # 83

length(intersect(hosts$latbi, unique(PA_hosts$species))) # 36 - 47 PA species aren't known hosts

# read in pest db
pest_NA <- read.csv("/home/liam/Documents/MSc/analysis/Data_Metadata/CSV/Insect x NA Host.csv")
pest_NA$NA_Host <- gsub("\\s", "_", pest_NA$NA_Host)

pest_N <- read.csv("/home/liam/Documents/MSc/analysis/Data_Metadata/CSV/Insect x N Host.csv")
pest_N$N_Host <- gsub("\\s", "_", pest_N$N_Host)

library(dplyr)

pest_NA_summary <- pest_NA %>%
                    group_by(Insect) %>%
                    summarize(n_hosts=length(NA_Host))

png("NA_hist.png", width = 580, height = 480)
hist(pest_NA_summary$n_hosts,
     main="",
     xlab="number of North American hosts")
dev.off()
# most pests are specific to a single host - but up to 12

pest_N_summary <- pest_N %>%
    group_by(Insect) %>%
    summarize(n_hosts=length(N_Host))

hist(pest_N_summary$n_hosts)
# same skew, up to 25

pest_sum <- pest %>%
    group_by(Insect) %>%
    summarize(n_hosts=length(Host),
              n_gen=length(unique(Host_Genus)))


png("all_hist.png", width = 580, height = 480)
hist(pest_sum$n_hosts,
     main="",
     xlab="total number of host species")
dev.off()

png("all_gen_hist.png", width = 580, height = 480)
hist(pest_sum$n_gen,
     main="",
     xlab="total number of host genera")
dev.off()

PA_pest_summary <- pest_NA %>%
    filter(NA_Host%in%PA_hosts$species) %>%
    group_by(Insect) %>%
    summarize(n_hosts=length(NA_Host))

hist(PA_pest_summary$n_hosts) # max of 7 - mostly 1's and 2's

focal_pest <- as.character(PA_pest_summary[which(PA_pest_summary$n_hosts==7),1])


fpest_hosts_NA <- subset(pest_NA, Insect==focal_pest)$NA_Host
length(fpest_hosts_NA) # 12

fpest_hosts_N <- subset(pest_N, Insect==focal_pest)$N_Host
length(fpest_hosts_N) # 6

fpest_hosts <- c(fpest_hosts_NA, fpest_hosts_N)
# these are all known hosts for focal pest - use to determine phylo niche

# read in plant megaphylogeny
library(ape)
phy <- read.tree("/home/liam/Documents/MSc/analysis/PhytoPhylo.tre")

length(intersect(fpest_hosts, phy$tip.label)) # 16

setdiff(fpest_hosts, phy$tip.label)
# Acer_grandidentatum - TPL says that this is a subspecies of Acer saccharum
# Acer_hyrcanum - TPL has a few synonyms, but none fit nicely on the tree - we'll omit for now (Acer_amaliae, Acer_italum, Acer_monspessulanum, Acer_opalus)
# grep("Acer_opalus", phy$tip.label)

fpest_hosts <- fpest_hosts[-which(fpest_hosts%in%c("Acer_grandidentatum", "Acer_hyrcanum"))]

fpest_mrca <- getMRCA(phy, fpest_hosts)

fpest_tree <- extract.clade(phy, fpest_mrca)


fpest_cols <- data.frame(tip=fpest_tree$tip.label)
fpest_cols$col <- "black"
fpest_cols$col[which(fpest_cols$tip%in%fpest_hosts)] <- "red"
fpest_cols$col[which(fpest_cols$tip%in%PA_hosts$species)] <- "blue"
fpest_cols$col[intersect(which(fpest_cols$tip%in%fpest_hosts), which(fpest_cols$tip%in%PA_hosts$species))] <- "purple"

png("fpest_phy", width=480, height=960)
plot(fpest_tree, tip.color = fpest_cols$col)
dev.off()

# combine pest/host lists
tmp1 <- pest_NA[1:3]
names(tmp1)[2] <- "Host"
tmp2 <- pest_N
names(tmp2)[2] <- "Host"
pest <- rbind(tmp1, tmp2)

pest_summary <- pest %>%
    group_by(Insect) %>%
    summarize(n_hosts=length(Host))

hist(pest_summary$n_hosts)

pest$Host_Genus <- gsub("_.*", "", pest$Host)
pest_genus_summary <- pest %>%
    group_by(Insect) %>%
    summarize(n_host_genera=length(unique(Host_Genus)))

hist(pest_genus_summary$n_host_genera)
table(pest_genus_summary$n_host_genera)
# collapsing genera leads to mostly n=1...


# Xyleborus glabratus - 11 host genera, 18 species

xg <- "Xyleborus glabratus"

xg_hosts <- subset(pest, Insect==xg)

xg_mrca <- getMRCA(phy, intersect(xg_hosts$Host, phy$tip.label))
xg_tree <- extract.clade(phy, xg_mrca)

xg_cols <- data.frame(tip=xg_tree$tip.label)
xg_cols$col <- "black"
xg_cols$col[which(xg_cols$tip%in%xg_hosts$Host)] <- "red"

png("xg_phy", width=960, height=960)
plot(xg_tree, type="fan", tip.color = xg_cols$col)
dev.off()


# format data for logistic regression
length(unique(pest$Host)) #416

length(intersect(unique(pest$Host), phy$tip.label)) #288 - many missing

host_phy <- keep.tip(phy, intersect(unique(pest$Host), phy$tip.label))

length(host_phy$tip.label) #288

library(castor)

d <- get_all_pairwise_distances(host_phy,
                                only_clades = 1:length(host_phy$tip.label))

focal_pest <- gsub(" ", "_", focal_pest)

d[which(host_phy$tip.label == focal_pest)]


# carry out logistic regression

# data shape:
# two columns, one 1/0 host status, one numeric PD from focal host

Fit_Logistic <- glm(Tb ~ Length,
                    family = binomial(link = " logit "),
                    data = data)
