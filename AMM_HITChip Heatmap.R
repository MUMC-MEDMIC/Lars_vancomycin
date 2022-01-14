library("microbiome")
library("gplots")

## Set location of files
setwd("Directory")

## Load in the L2-frpa.tab file
l2file <- "L2-frpa.tab"

## Read file as csv
L2tab <- read.csv(l2file, header = TRUE, sep = '\t', row.names =1, as.is = TRUE)

## Convert first line ("samples") to column name.
## Read the first row, split based on tab
## unlist produces a vector which contains all the atomic components which occurred
colnames(L2tab) <- unlist(strsplit(readLines(l2file,1),"\t"))[-1]

## Convert numbers to numeric need to get row and column out to avoid missing info
rnams <- rownames(L2tab) ## Species
cnams <- colnames(L2tab) ## Samples

## Convert to numeric
L2tab <- apply(L2tab, 2, as.numeric)

## Put species and samples back
rownames(L2tab) <- rnams
colnames(L2tab) <- cnams

## Convert to log10. Notice that no NA in the dataframe any(is.na(L2tab))
L2 <- log10(L2tab)

## Load in metadata file
metaData = read.delim(file="metadata_ABS_mod.txt",header=TRUE)
## Turn sampleID into row names
rownames(metaData) <- metaData$sampleID

## Gets samples appear in metadata and L2 file
allSamples <- as.character(intersect(rownames(metaData),colnames(L2)))
## Gets info of all samples from L2 file
sampdata <- L2[, allSamples]

## Gets metadata from matched samples
metadata <- metaData[allSamples,]

## Sort meta based on time period: pre, post and follow-up
metadata_pre <- subset(metaData[allSamples,], time == "pre")
metadata_post <- subset(metaData[allSamples,], time == "post")

## Collect samples with pre and post time stamp
matching_samples <- intersect(unique(metadata_pre$subject),unique(metadata_post$subject))

## Gets subjects that matched with matching samples for pre and post
matching_samples_pre <- subset(metadata_pre, subject %in% matching_samples)
matching_samples_post <- subset(metadata_post, subject %in% matching_samples)

## Combine objects by rows of all two matching
resulting_samples <- rbind(matching_samples_pre,matching_samples_post)

## Rename rownames as sampleID
rownames(resulting_samples) <- resulting_samples$sampleID

## Gets sampleID that have VanB_presence_baseline0 as intervention at different time point
sampVanB0pre <- as.character(subset(metaData[rownames(resulting_samples),], Intervention == "Vancomycin" & VanB_presence_baseline == "0" & time == "pre")$sampleID)
sampVanB0post <- as.character(subset(metaData[rownames(resulting_samples),], Intervention == "Vancomycin" & VanB_presence_baseline == "0" & time == "post")$sampleID)

## Create data based on samples found previously (Log10)
sampl2log <- L2[, allSamples]

## Calculate Wilcoxon between vanBbase0 pre and post. Gets samples and species
samples <- colnames(sampl2log) ## Samples
levels <- rownames(sampl2log) ## Species

## Creates empty matrix
M <- matrix(data = NA, length(levels), 1)
## write Species as row names
rownames(M) <- levels

for (i in 1:length(levels)) {
  ## species
  lvl <- levels[i]
  ## Gets values based on species in VanB baseline0 (pre)
  l.g1 <- sampl2log[lvl, sampVanB0pre]
  ## Gets values based on species in VanB baseline0 (post)
  l.g2 <- sampl2log[lvl, sampVanB0post]
  
  ## Wilcox test between pre and post  
  p <- wilcox.test(as.numeric(l.g1), as.numeric(l.g2), 
                   paired = TRUE)$p.value
  
  # message(lvl, ' p-value: ', p, '\n')
  ## Write p values to matrix
  M[i, 1] <- p
  
}

vanB_base0_pre_postM <- M

## Adjust P-values for Multiple Comparisons with Benjamini & Hochberg (1995)
## ('BH' or its alias 'fdr')
## Stats packages
p.adjust.method = "BH"

correct_vanBbase0_p <- p.adjust(vanB_base0_pre_postM, method = p.adjust.method)

## add species to corected p values
names(correct_vanBbase0_p) <- rownames(vanB_base0_pre_postM)

## Sort and save to p vals 
correct_vanBbase0_p <- sort(correct_vanBbase0_p)

## Gets sampleID that have VanB_presence_baseline1 as intervention at different time point
sampVanB1pre <- as.character(subset(metaData[rownames(resulting_samples),], Intervention == "Vancomycin" & VanB_presence_baseline == "1" & time == "pre")$sampleID)
sampVanB1post <- as.character(subset(metaData[rownames(resulting_samples),], Intervention == "Vancomycin" & VanB_presence_baseline == "1" & time == "post")$sampleID)

## Creates empty matrix
M2 <- matrix(data = NA, length(levels), 1)
## Species as row names
rownames(M2) <- levels

for (i in 1:length(levels)) {
  ## species
  lvl <- levels[i]
  ## Gets values based on species in VanB baseline1 pre
  l.g3 <- sampl2log[lvl, sampVanB1pre]
  ## Gets values based on species in VanB baseline1 post
  l.g4 <- sampl2log[lvl, sampVanB1post]
  
  ## Wilcox test between pre and post  
  p <- wilcox.test(as.numeric(l.g3), as.numeric(l.g4), 
                   paired = TRUE)$p.value
  
  # message(lvl, ' p-value: ', p, '\n')
  ## Write to matrix
  M2[i, 1] <- p
  
}

vanB_base1_pre_postM <- M2

## Correct p values based on BH method.
correct_vanBbase1_p <- p.adjust(vanB_base1_pre_postM, method = p.adjust.method)

## add species to corrected p values
names(correct_vanBbase1_p) <- rownames(vanB_base0_pre_postM)

## Sort and save to pvals 
correct_vanBbase1_p <- sort(correct_vanBbase1_p)


## Combine vanb0 and vanb1 pre and post
sampVanB01pre <- c(sampVanB0pre,sampVanB0post)
sampVanB01post <- c(sampVanB1pre,sampVanB1post)


## Empty
fcs <- c()
## Loop through species in sample data log (sampdata)
for (tax in rownames(sampdata)) {
  ## Mean of spices from GV4 minus mean of spices in GV5
  fcs[[tax]] <- mean(sampdata[tax, sampVanB01pre]) - mean(sampdata[tax, sampVanB01post])
}

## Species names from vanBbase0 with p values < 0.1
vanB_base0_species <- names(which(correct_vanBbase0_p < 0.1))
## Number of vanBbase0 species
length(vanB_base0_species)

## Species names from vanBbase1 with p values < 0.01
vanB_base1_species <- names(which(correct_vanBbase1_p < 0.1))
## Number of vanBbase1 species
length(vanB_base1_species)

## Combining boh species. Remove duplicates
vanB_base01_species <- unique(c(vanB_base0_species,vanB_base1_species))
## Number of vanB_base01 species
length(vanB_base01_species)

## List with pre and post sample name
vanB_base01_subject_list <- split(as.character(c(sampVanB0pre,sampVanB0post,sampVanB1pre,sampVanB1post,rownames(resulting_samples[resulting_samples$Intervention=="Placebo",]))), as.character(sub("_0.$","", c(sampVanB0pre,sampVanB0post,sampVanB1pre,sampVanB1post,rownames(resulting_samples[resulting_samples$Intervention=="Placebo",])))))

## List of samples for placebo, vanB_base0 and vanB_base1
placeboVanB01 <- as.character(c(sampVanB0pre,sampVanB0post,sampVanB1pre,sampVanB1post,rownames(resulting_samples[resulting_samples$Intervention=="Placebo",])))
## Metadata from placebo, vanB_base0 and vanB_base1
vanBbase01_resulting_samples <- resulting_samples[rownames(resulting_samples) %in% placeboVanB01,]

## The log10 of significant species for placebo, vanB_base0 and vanB_base1 samples
l2logSpecies <- L2[vanB_base01_species,rownames(vanBbase01_resulting_samples)]

## Create matrix of log10 for significant species pre and post
vanB_base01_list <- list()

## Create list of samples with refine species and pre post values
for(i in 1:length(vanB_base01_subject_list)){
  ## Gets subject based on name
  x=vanB_base01_subject_list[[names(vanB_base01_subject_list)[i]]]
  ## Convert to character
  y=as.character(x)
  ## Matrix of species per subject
  temp_matrix=l2logSpecies[,y]
  ## write to vanB_base0_list
  vanB_base01_list[[i]] <- temp_matrix
}

## Name sample for each matrix
names(vanB_base01_list) <- names(vanB_base01_subject_list)

## Calculate log10 differences between pre pro. aka pre columns - post column
vanB_base01_values <- list()

for(i in 1:length(vanB_base01_subject_list)){
  
  ## Name of sample
  x=names(vanB_base01_list)[i]
  ## Bacteria change between pro - pre
  pre_pro=vanB_base01_list[[x]][,2] - vanB_base01_list[[x]][,1]
  ## Combine values
  vanB_base01_values[[i]] <- pre_pro
}

## Name subject for each species value matrix
names(vanB_base01_values) <- names(vanB_base01_subject_list)

## Convert to dataframe
vanB_base01_values <- as.data.frame(vanB_base01_values)

## Convert to matrix
vanB_base01_values_matrix <- as.matrix(vanB_base01_values)
#write.csv(vanB_base0_values.matrix, "matrix_pre_pro_Van0.csv")

## Separate all intervention into amoxicillin, vancomycin and placebo
treatmentOrder <- split(as.character(vanBbase01_resulting_samples$subject),as.character(vanBbase01_resulting_samples$Intervention), drop= TRUE)

## Matrix from pre_pro for each amoxy, vanco and placebo
vanB_base0_samplesID <- resulting_samples[rownames(resulting_samples) %in% sampVanB0pre,]$subject
matrix_pre_pro_vanB_base01_Vanco0 <- vanB_base01_values_matrix[,vanB_base0_samplesID]

vanB_base1_samplesID <- resulting_samples[rownames(resulting_samples) %in% sampVanB1pre,]$subject
matrix_pre_pro_vanB_base01_Vanco1 <- vanB_base01_values_matrix[,vanB_base1_samplesID]

matrix_pre_pro_vanB_base01_Placebo <- vanB_base01_values_matrix[,unique(treatmentOrder$Placebo)]

dev.off()
## Heatmap. This takes order of the significant species
hm1plcbo <- heatmap(matrix_pre_pro_vanB_base01_Placebo)
hm2vanB0 <- heatmap(matrix_pre_pro_vanB_base01_Vanco0)
hm3vanB1 <- heatmap(matrix_pre_pro_vanB_base01_Vanco1)

## Combine and sort the same sample order as the heat map above
vanB_base01_input <- cbind(matrix_pre_pro_vanB_base01_Placebo[,hm1plcbo$colInd],matrix_pre_pro_vanB_base01_Vanco0[,hm2vanB0$colInd],matrix_pre_pro_vanB_base01_Vanco1)

## Get taxonomic names from input
vanB_base01_taxa <- rownames(vanB_base01_input)

## Create dataframe 
vanB_base01_df <- data.frame(list(taxa = vanB_base01_taxa))

## Load in phylogeny information
phyfull <- "phylogeny.full.tab"
phytab <- read.csv(phyfull, header = TRUE, sep = "\t", as.is = TRUE)

## Go from level 2 Order to level 1 Phylum
level.from = "L2"
level.to = "L1"

## Match taxa with L2 to get omap
omap <- phytab[match(as.character(vanB_base01_taxa),phytab[["L2"]]),c("L2", "L1")]

## Convert to factors to preserve order
omap[["L2"]] <- factor(omap[["L2"]])
omap[["L1"]] <- factor(omap[["L1"]])

## Gets L1 from omap
vanB_base01_species_list <- omap[, 2]

## drops repeats and write to L1 of taxa
vanB_base01_df[["L1"]] <- unlist(droplevels(vanB_base01_species_list))

## Rownames from taxa in correct order
rownames(vanB_base01_input) <- vanB_base01_df$taxa

vanB_base01_input_ordered <- vanB_base01_input[order(rownames(vanB_base01_input)),]

## Identify species based on experiments
vanB01common <- intersect(vanB_base0_species,vanB_base1_species)
vanB0specific <- setdiff(vanB_base0_species,vanB01common)
vanB1specific <- setdiff(vanB_base1_species,vanB01common)

## All colour black
cols <- rep('black', nrow(vanB_base01_input_ordered))
## vanB_base0 significant only
cols[row.names(vanB_base01_input_ordered) %in% vanB0specific] <- 'green3'
## vanB_base1 significant only
cols[row.names(vanB_base01_input_ordered) %in% vanB1specific] <- 'red2'

######################################################
dev.off()

heatmap.2(vanB_base01_input_ordered, main="",col=bluered(100), dendrogram="none",
          symbreaks=TRUE,trace="none",margin=c(8,15),Rowv=T, colRow = cols,
          Colv=F,cexRow=0.7,cexCol=1,keysize = 1,density.info="none",key.title="NA",
          colsep= c(ncol(matrix_pre_pro_vanB_base01_Placebo),ncol(matrix_pre_pro_vanB_base01_Placebo)  + ncol(matrix_pre_pro_vanB_base01_Vanco0),ncol(matrix_pre_pro_vanB_base01_Placebo)  + ncol(matrix_pre_pro_vanB_base01_Vanco0) + ncol(matrix_pre_pro_vanB_base01_Vanco1)),
          sepwidth =c(0.1))
