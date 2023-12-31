#load in data
setwd('treemetrics/')
library(tidyverse)


annofile <- read.table("annot_Pdavidsonii_1mb.gffread.genes.bed")
annofile <- annofile[,c(1,2,3,6,4)]
colnames(annofile) <- c("chr", "start", "end", "strand", "gene_id")


##setup for RF estimation with 10kb trees:
##################################################
library(ape)
library(phangorn)
library(tidyverse)

#read in the species tree and reroot to correct outgroup
speciestree <- read.tree("astral_CDS_noannotations.tre")

#if tip labels in species tree do not match those in gene trees...
for (i in 1:length(speciestree$tip.label)){
  speciestree$tip.label[i] <-
    str_split(str_split(speciestree$tip.label[i], "CDS_bws_")[[1]][2], ".fixed.fa")[[1]][1]
}

speciestree <- root(speciestree, outgroup = "mon_61-7_S440")
plot(speciestree)

#read in the set of gene trees of interest
intrees <- read.tree("combined_10kbwindowtrees.tre")


#read in the treepath names
treename_matrix <- read.table("numbered_10kbtreepaths.txt")
colnames(treename_matrix) <- c('number','oldnames')


#define function to get information about the scaffold and region each window is from.
#(always give your files descriptive names :) )
addinfo <- function(j, treematrix, addtype) {
  if(addtype == "scaffold"){
    return(paste("scaffold_", 
                 strsplit((strsplit(treename_matrix[j,2], 'scaf_')[[1]][2]), '/')[[1]][1],
                 sep = '')
    )
  }
  if(addtype == "bpstart"){
    return(strsplit((strsplit(treematrix[j,2], 'bp_')[[1]][2]), '-')[[1]][1])
  }
  if(addtype == "bpend"){
    return(strsplit((strsplit(treename_matrix[j,2], '-')[[1]][2]), '.fa.treefile')[[1]][1])
  }
}


#add new columns to the treename matrix with this parsed information
treename_matrix$scaffolds <- sapply((1:nrow(treename_matrix)), addinfo,
                                    treematrix = treename_matrix,
                                    addtype = "scaffold")

treename_matrix$bpstart <- sapply((1:nrow(treename_matrix)), addinfo,
                                  treematrix = treename_matrix,
                                  addtype = "bpstart")

treename_matrix$bpend <- sapply((1:nrow(treename_matrix)), addinfo,
                                treematrix = treename_matrix,
                                addtype = "bpend")


#find RF distance for each tree. First, root trees to the outgroup.
for (i in 1:length(intrees)){
  intrees[[i]] <- root(phy = intrees[[i]], outgroup = "mon_61-7_S440")
}

#Then, estimate normalized RF distance and add to treename_matrix
treename_matrix$RFdist <- RF.dist(tree1 = intrees, tree2 = speciestree, normalize = T, rooted = T)

#write this to a file so we can load in again later
write.csv(treename_matrix, "RFdistance_10kbtrees_astralCDStree.csv", row.names = F)
##################################################



#Plotting RF distance in 10kb stretches along scaffold
#and overlay that onto plot of genic content
##################################################
#read in RF values and annofile info
RFfile <- read.csv("RFdistance_10kbtrees_astral10kbtree.csv")
colnames(RFfile)[3] <- "chromosome"


annofile <- read.table("annot_Pdavidsonii_1mb.gffread.genes.bed")
annofile <- annofile[,c(1,2,3,6,4)]
colnames(annofile) <- c("chromosome", "start", "end", "strand", "gene_id")

#make list of bad scaffolds
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")


#next, add midpoint info to the annotation file and treename file.
#this will be used to determine the number of CDS withiin tree windows.
filter_annofile <- annofile %>%
  filter(!chromosome %in% badscafs) %>%
  mutate(midpoint = ceiling((start+end)/2)) %>%
  rename(bpstart = start, bpend = end)
save(filter_annofile, file = "filtered_annotationfile.obj")

filter_RFfile <- RFfile %>%
  filter(!chromosome %in% badscafs) %>%
  mutate(midpoint = ceiling((bpstart+bpend)/2)) %>%
  mutate(chromosome = gsub("fold", "", chromosome))
save(filter_RFfile, file = "filtered_RFfile_10kb_ASTRAL.obj")




###MAKE PLOT, USING GENE DENSITY
genicfractionfile <- read.delim("genicfraction_10kbwin_10kbslide.bed",
                                sep = ' ', header = F) %>%
  rename(chromosome = V1, bpstart = V2, bpend = V3, genic_fraction = V4) %>%
  mutate(bpstart = bpstart + 1,
         bpend = bpend + 1,
         midpoint = ceiling((bpstart+bpend)/2),
         chromosome = gsub("fold", "", chromosome))
save(genicfractionfile, file = "filtered_genic-fraction_10kb.obj")


#writing the plot
combplot_RF_GF <- inner_join(genicfractionfile, filter_RFfile,
                             by = c("chromosome", "midpoint")) %>%
  pivot_longer(., cols = c(genic_fraction, RFdist))

#this is part of figure 3, so it is c

c <- ggplot(combplot_RF_GF, aes(x = midpoint/1000000, y = value, group = name)) +
  geom_smooth(linewidth = 0.5, aes(col = name), method = "loess", se = F, span = 0.15) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  theme_bw() +
  xlab("Position on scaffold (Mb)") +
  #theme(axis.title.y = element_blank(),legend.title = element_blank()) +
  scale_color_manual(values = c("red", "cyan"),
                     labels = c("genic fraction", "Robinson-Foulds distance"))


