#set working directory
setwd("/Volumes/GoogleDrive/Shared drives/Radish amendment microbes/16S_paired_end/diff_abun_analysis")

#load libraries etc.
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("plyr")
library("vegan")
library("ggpubr")
library("gplots")
library("RColorBrewer")
theme_set(theme_bw())

#import a phyloseq object (biom should have taxonomy attached), start with UNRAREFIED otu table
radish_16S <- import_biom("feature-table-with-taxonomy.biom") #, treefilename="radish_rooted-tree-filtered2-prok.nwk" ,refseqfilename="radish_prok-only.fasta")

#rename columns 
colnames(tax_table(radish_16S)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#save otu table (optional) -- this will save a simple OTU table without taxonomy
#write.csv(as.data.frame(otu_table(radish_16S)),"otu_radish_16S_initial.csv")

#read in metadata
radish_16S_metadata <- read.csv("metadata-radish-16S.csv")

#create another row for "SampleID", create a dataframe for the metadata
row.names(radish_16S_metadata) <- radish_16S_metadata$X.SampleID
radish_16S_metadata$SampleID <- radish_16S_metadata$X.SampleID
sample_data(radish_16S) <- radish_16S_metadata

#View(data.frame(sample_data(radish_16S)))

###De-QIIME-ify the taxa table -- this will separate taxonomic ranks into each separate columns. This _0 table is necessary downstream.
tax_table_radish_16S_0 <- as.data.frame(tax_table(radish_16S))

#OPTIONAL export taxa table
# write.csv(tax_table_radish_16S_0, "tax_table_radish_16S_0.csv")

#Make a copy of the table
tax_table_radish_16S <- tax_table_radish_16S_0

#Renaming the taxonomy if not standard (this may not be necessary depending on the taxonomic database used)
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Acidobacteriota", "Acidobacteria", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Actinobacteriota", "Actinobacteria", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Armatimonadota", "Armatimonadetes", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Bacteroidota", "Bacteroidetes", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Gemmatimonadota", "Gemmatimonadetes", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Halobacterota", "Halobacteria", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Planctomycetota", "Planctomycetes", x)}))
#tax_table_radish_16S <- data.frame(lapply(tax_table_radish_16S, function(x) {gsub("Verrucomicrobiota", "Verrucomicrobia", x)}))

#removing taxonomic notations
tax_table_radish_16S$Kingdom<- gsub("k__", "", tax_table_radish_16S$Kingdom)#sometimes "k" is replaced by "d" so make sure that is is properly removed.
tax_table_radish_16S$Phylum <- gsub("p__", "", tax_table_radish_16S$Phylum)
tax_table_radish_16S$Class <- gsub("c__", "", tax_table_radish_16S$Class)
tax_table_radish_16S$Order <- gsub("o__", "", tax_table_radish_16S$Order)
tax_table_radish_16S$Family <- gsub("f__", "", tax_table_radish_16S$Family)
tax_table_radish_16S$Genus <- gsub("g__", "", tax_table_radish_16S$Genus)
tax_table_radish_16S$Species <- gsub("s__", "", tax_table_radish_16S$Species)
#tax_table_radish_16S$Gen_Fam <- paste(tax_table_radish_16S$Genus, " (", tax_table_radish_16S$Family,")",sep="")

#View(tax_table_radish_16S_0)
#View(tax_table_radish_16S)

row.names(tax_table_radish_16S) <- row.names(tax_table_radish_16S_0)
tax_table(radish_16S) <- as.matrix(tax_table_radish_16S)
# View(data.frame(tax_table(radish_16S)))
# tax_table_radish_16S_1 <- as.data.frame(tax_table(radish_16S))
# write.csv(tax_table_radish_16S_1, "tax_table_radish_16S_1.csv")

#subsetting your datasets (often it will requires a slow narrowing down of each category until you get the samples you want)
#can also use for subsetting taxa: radish_16S_sub0 = subset_taxa(radish_16S, Kingdom=="Bacteria")
radish_16S_treatment = subset_samples(radish_16S, Paper == "n-fert")#includes the following column category
radish_16S_treatment = subset_samples(radish_16S_treatment, Treatment != "Control")#excludes the following column category
#radish_16S = subset_samples(radish_16S, SampleID != "02.2.5.1.16S.a")#exclude the following sample based on ID

radish_16S <- radish_16S_treatment

# View(data.frame(otu_table(radish_16S_treatment)))

# Check rarefaction of the data
# rarecurve(t(otu_table(radish_16S)), step=50, cex=0.5)

#removing samples that didn't work -- low read counts (this would best be done in QIIME)
#radish_16S_0 <- prune_samples(sample_sums(radish_16S) >= 100, radish_16S)#if done here, pass this object onto downstream workflow instead of radish_16S
# #View(data.frame(sample_data(radish_16S_0)))
# #View(data.frame(otu_table(radish_16S_0)))
# rarecurve(t(otu_table(radish_16S_0)), step=100, cex=0.5)
# radish_16S_0_otu <- data.frame(otu_table(radish_16S_0))
# write.csv(radish_16S_0_otu,"radish_16S_0_otu.csv")

##########################################
###DESeq Comparison [Urea (left) vs Compost (right), C1 ]
##########################################
###Note: This approach requires you to make a new phyloseq object for each comparison, which is shown below (using "subset_samples")

###View(data.frame(sample_data(radish_16S_0)))
#Subsetting your dataset to make various comparisons
radish_16S_C1 <- subset_samples(radish_16S, Time.Point == "1")
radish_16S_C6 <- subset_samples(radish_16S, Time.Point == "6")

#remove empty cells due to subsetting
#can use this to filter out lower abundance taxa (e.g. x > 10), although in this dataset it didn't make a difference
radish_16S_C1 <- filter_taxa(radish_16S_C1, function(x) sum(x) > 0, TRUE)
radish_16S_C6 <- filter_taxa(radish_16S_C6, function(x) sum(x) > 0, TRUE)

#An error will occur later on when running DESeq for sample C6 because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the radish_16S_C6 phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
radish_16S_C6 <- phyloseq((otu_table(radish_16S_C6)+1), sample_data(radish_16S_C6), tax_table(radish_16S_C6))

#Convert phyloseq to DESeq Data Set object (dds)
radish_16S_C1dds <- phyloseq_to_deseq2(radish_16S_C1, ~Treatment)
radish_16S_C6dds <- phyloseq_to_deseq2(radish_16S_C6, ~Treatment)

#Determine which level should the dataset be set as the REFERENCE sample class
radish_16S_C1dds$Treatment <- relevel(radish_16S_C1dds$Treatment, "Urea")
radish_16S_C6dds$Treatment <- relevel(radish_16S_C6dds$Treatment, "Urea")

#Perform the contrast using standard and recognized parameters and tests
radish_16S_C1dds = DESeq(radish_16S_C1dds, test="Wald", fitType="parametric")
radish_16S_C6dds = DESeq(radish_16S_C6dds, test="Wald", fitType="local")#parametric fit didn't work well for this dataset

#Performing the final calculations and extracting the points
radish_16S_C1dds_results = results(radish_16S_C1dds, cooksCutoff = TRUE)
radish_16S_C6dds_results = results(radish_16S_C6dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "Waialua Farm", and a negative fold change indicates enrichment in "Urban Garden".
# Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
radish_16S_C1dds_results = lfcShrink(dds=radish_16S_C1dds, contrast = c("Treatment","Compost","Urea"), res=radish_16S_C1dds_results, type="normal")
radish_16S_C6dds_results = lfcShrink(dds=radish_16S_C6dds, contrast = c("Treatment","Compost","Urea"), res=radish_16S_C6dds_results, type="normal")

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#this code basically extracts information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table_C1 = radish_16S_C1dds_results[which(radish_16S_C1dds_results$padj <= alpha), ]
sig_table_C6 = radish_16S_C6dds_results[which(radish_16S_C6dds_results$padj <= alpha), ]#cannot go any further because all samples were non-significant

#Bind taxa names to tables of significant taxa
sig_table_C1 = cbind(as(sig_table_C1, "data.frame"), as(tax_table(radish_16S_C1)[rownames(sig_table_C1), ], "matrix"))

#head(sig_table_C1)
# write.csv(sig_table_C1, "sig_table_C1.csv")

#find those that are not significant, ASVs that are common among the treatments
#notsig_table_C1 = radish_16S_C1dds_results[which(radish_16S_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_C1 = cbind(as(notsig_table_C1, "data.frame"), as(tax_table(radish_16S_C1)[rownames(notsig_table_C1), ], "matrix"))
#head(notsig_table_C1)
# write.csv(notsig_table_C1, "notsig_table_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_C1_Phyla = unique(sig_table_C1$Phylum)
# View(sig_table_C1_Phyla)

#Remove anything that do not have a family or phylum taxonomy
sig_table_C1sub <- subset(sig_table_C1, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="")

#Plots the logfold changes
sig_table_C1subp = ggplot(sig_table_C1sub, aes(x=log2FoldChange, y=reorder(Genus,desc(Genus)), color=Phylum))  + geom_point(size=2, stroke = 1) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("Prokaryote taxa at C1 between Urea (left) vs compost (right)") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

#Plot and beautify by faceting based on phyla
sig_table_C1subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

#Color Brewer Pallette
#RColorBrewer::display.brewer.all()

#Manual color scales
Brewer9 <- c("Acidobacteria" = "#a6cee3", "Actinobacteria" = "#1f78b4", "Bacteroidetes" = "#b2df8a", "Chloroflexi" = "#33a02c", "Crenarchaeota" = "#fb9a99", "Nitrospirae" = "#e31a1c", "Planctomycetes" = "#fdbf6f", "Proteobacteria" = "#ff7f00", "Verrucomicrobia" = "#cab2d6")
Rick <- c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Bdellovibrionota" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Myxococcota" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1")
Dark <- c("Acidobacteria" = "#1c9e77", "Actinobacteria" = "#d95f02", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "#756fb4", "Bdellovibrionota" = "bisque3", "Crenarchaeota" = "#e6ab01", "Chloroflexi" = "#66a621", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Myxococcota" = "darkturquoise", "Nitrospirae" = "#f45c85", "OP3" = "black", "Planctomycetes" = "#a6761d","Proteobacteria" = "#0c7be3", "Thermi" = "hotpink", "Verrucomicrobia" = "#646782", "WS3" = "magenta1")

#The code below will assign a color to most common phyla. You may still need to adjust it for your needs.
sig_table_C1subp2 <- sig_table_C1subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank()) + theme(panel.spacing = unit(0.1, "lines")) + theme(text = element_text(size = 10)) +scale_colour_manual(values = Dark) #+scale_color_brewer(palette="Dark2")
sig_table_C1subp2


#This will make an svg of your file, which I prefer to view and edit in Inkscape.
pdf("sig_table_C1subp2.pdf", width = 8, height = 8, pointsize=12)

svg("sig_table_C1subp2.svg", width = 8, height = 8, pointsize=12)
plot(sig_table_C1subp2)
dev.off()
