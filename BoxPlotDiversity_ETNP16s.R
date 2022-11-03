# Box plot Diversity 
## Madeleine Thompson 
### ETNP 16s data 


### Install Packages 
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
install_phyloseq(branch = "devel")
install.packages("ggplot2", "dplyr", "tidyr", "ape", "grid", "gridExtra", "vegan", "ggpubr", "ecolTest", "BiocManager")

############

## Import data 
zotu16s <- read.csv("data/zotutab.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)

############
## Install Packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(dplyr)
library(tidyr)
library(vegan)
library(ggpubr)
library(ggrepel)
library("viridis")
library(ape)
library(grid)
library(gridExtra)
library(ecolTest)

############

# Functions 

#' Separates taxonomic groups and removes unneeded characters in taxonony assignments 
#'
#' @param data is taxonmy data from qiime
#'
#' @returns a data frame with corrected taxonomy names 
#' @examples
#' taxGroup(tax16s)
taxGroup <- function(data){
  data1 <- tidyr::separate(data = data, col =  Taxon,
                           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                           sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn");
  data1$Domain<-gsub("d__","",as.character(data1$Domain));
  data1$Phylum<-gsub("p__","",as.character(data1$Phylum));
  data1$Class<-gsub("c__","",as.character(data1$Class));
  data1$Order <- gsub("o__", "", as.character(data1$Order));
  data1$Family<-gsub("f__","",as.character(data1$Family));
  data1$Genus<-gsub("g__","",as.character(data1$Genus));
  data1$Species<-gsub("s__","",as.character(data1$Species));
  return(data1)
}


#' Transforms phyloseq data 
#'
#' @param data is merged taxonomy, OTU, and metadata S$ data frame 
#'
#' @returns a phyloseq S4 with the top 500 otu, taxa, and metadata  
#' @examples
#' transformTop500(physeq1)
transformTop500 <- function(data){
  physeq_ns <- phyloseq::filter_taxa(data, function(x) sum(x) > 1, prune = TRUE); 
  phylor = phyloseq::transform_sample_counts(physeq_ns, function(x) x / sum(x) );
  phylor_nz <- phyloseq::filter_taxa(phylor, function(x) mean(x) > 0, prune = TRUE);
  top <- names(sort(phyloseq::taxa_sums(phylor_nz), decreasing=TRUE))[1:500]; 
  phylor.t <- phyloseq::prune_taxa(top, phylor_nz);
  return(phylor.t)
}

############

# Import data 
tax16s <- read.csv("data/taxonomy.csv",header = TRUE, row.names = 1)
metadata.phy = phyloseq::sample_data(data.frame(metadata))

# Separate tax groups and remove all starting characters in taxonomy - custom function
tax16sETNP <- taxGroup(tax16s)

tax16sETNP$Phylum <- sub("Crenarchaeota", "Thaumarchaeota", tax16sETNP$Phylum)
tax16sETNP <- subset(tax16sETNP, Phylum != " Bdellovibrionota")
tax16sETNP$Phylum <- sub(" NB1-j", " NB1-j (Deltaproteobacteria)", tax16sETNP$Phylum)
tax16sETNP$Phylum <- sub(" Marinimicrobia_(SAR406_clade)", " Marinimicrobia(SAR406 clade)", tax16sETNP$Phylum)
tax16sETNP$Phylum <- sub(" SAR324_clade(Marine_group_B)", " Marine group B(SAR324 clade)", tax16sETNP$Phylum)



# Make ZOTU table numeric and a matrix 
# Make Taxonomy a matrix and make the row and column names the same 
mat<- as.matrix(zotu16s)
zotuTab <- matrix(as.numeric(mat),
                  ncol = ncol(mat));
colnames(zotuTab) <- colnames(mat);
rownames(zotuTab) <- rownames(mat);

taxTab <- as.matrix(tax16sETNP)

setdiff(rownames(taxTab), rownames(zotuTab))
all(rownames(zotuTab) == rownames(taxTab))
setdiff(colnames(zotuTab), rownames(taxTab))

############

# Combine Tax and ZOTU and metadata
OTU <- phyloseq::otu_table(zotuTab, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxTab) 

physeq = phyloseq::phyloseq(OTU, TAX) # put metadata in this 
phyloseq::taxa_names(physeq)

physeq1 <-phyloseq::merge_phyloseq(physeq, metadata.phy)

phyloseq::sample_names(physeq)
phyloseq::sample_names(metadata.phy)

#relative abundance and transform data  - custom function 
phylor.top <- transformTop500(physeq1)


oxy1 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$OxygenConc. == "1")
oxy1.stn1 <- phyloseq::subset_samples(oxy1, oxy1@sam_data$Station == "Stn1")
oxy1.stn2 <- phyloseq::subset_samples(oxy1, oxy1@sam_data$Station == "Stn2")
oxy1.stn3 <- phyloseq::subset_samples(oxy1, oxy1@sam_data$Station == "Stn3")

oxy2 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$OxygenConc. == "2")
#oxy2.stn1 <- phyloseq::subset_samples(oxy2, oxy2@sam_data$Station == "Stn1")
oxy2.stn2 <- phyloseq::subset_samples(oxy2, oxy2@sam_data$Station == "Stn2")
oxy2.stn3 <- phyloseq::subset_samples(oxy2, oxy2@sam_data$Station == "Stn3")

oxy3 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$OxygenConc. == "3")
oxy3.stn1 <- phyloseq::subset_samples(oxy3, oxy3@sam_data$Station == "Stn1")
oxy3.stn2 <- phyloseq::subset_samples(oxy3, oxy3@sam_data$Station == "Stn2")
oxy3.stn3 <- phyloseq::subset_samples(oxy3, oxy3@sam_data$Station == "Stn3")

#Shannon Station 1
H1.1 <- vegan::diversity(oxy1.stn1@otu_table)
#H2.1 <- vegan::diversity(oxy2.stn1@otu_table) does not have any oxy 2 in stn 1 
H3.1 <- vegan::diversity(oxy3.stn1@otu_table)

# Richness Station 1
R1.1 <- vegan::specnumber(oxy1.stn1@otu_table)
R3.1 <- vegan::specnumber(oxy3.stn1@otu_table)

#Evenness Station 1
E1.1 <- H1.1/log(R1.1)
E3.1 <- H3.1/log(R3.1)



#Shannon Station 2
H1.2 <- vegan::diversity(oxy1.stn2@otu_table)
H2.2 <- vegan::diversity(oxy2.stn2@otu_table)
H3.2 <- vegan::diversity(oxy3.stn2@otu_table)

# Richness Station 2
R1.2 <- vegan::specnumber(oxy1.stn2@otu_table)
R2.2 <- vegan::specnumber(oxy2.stn2@otu_table)
R3.2 <- vegan::specnumber(oxy3.stn2@otu_table)

#Evenness Station 2

E1.2 <- H1.2/log(R1.2)
E2.2 <- H2.2/log(R2.2)
E3.2 <- H3.2/log(R3.2)


#Shannon Station 3
H1.3 <- vegan::diversity(oxy1.stn3@otu_table)
H2.3 <- vegan::diversity(oxy2.stn3@otu_table)
H3.3 <- vegan::diversity(oxy3.stn3@otu_table)

# Richness Station 3
R1.3 <- vegan::specnumber(oxy1.stn3@otu_table)
R2.3 <- vegan::specnumber(oxy2.stn3@otu_table)
R3.3 <- vegan::specnumber(oxy3.stn3@otu_table)

#Evenness Station 3
E1.3 <- H1.3/log(R1.3)
E2.3 <- H2.3/log(R2.3)
E3.3 <- H3.3/log(R3.3)

#Hutchinson 
hutch1 <- ecolTest::Hutcheson_t_test(
  H1.1,
  H3.1,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch1
#Station 1 oxy 1 and 3 are significantly different 

hutch2.1 <- ecolTest::Hutcheson_t_test(
  H1.2,
  H2.2,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch2.1
#Station 2 oxy 1 and 2 are significantly different 

hutch2.2 <- ecolTest::Hutcheson_t_test(
  H1.2,
  H3.2,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch2.2
#Station 2 oxy 1 and 3 are significantly different 

hutch2.3 <- ecolTest::Hutcheson_t_test(
  H3.2,
  H2.2,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch2.3
#Station 2 oxy 2 and 3 are significantly different 

hutch3.1 <- ecolTest::Hutcheson_t_test(
  H1.3,
  H2.3,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch3.1
#Station 3 oxy 1 and 2 are significantly different 

hutch3.2 <- ecolTest::Hutcheson_t_test(
  H1.3,
  H3.3,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch3.2
#Station 3 oxy 1 and 3 are significantly different 

hutch3.3 <- ecolTest::Hutcheson_t_test(
  H3.3,
  H2.3,
  shannon.base = exp(1),
  alternative = "two.sided",
  difference = 0
)
hutch3.3
#Station 3 oxy 2 and 3 are significantly different 





png("figures/boxplot.png", height = 10, width = 9, units = "in", res = 300)
par(mfrow=c(3,2))
boxplot(H1.1, H3.1,
        ylab = "Alpha Diversity Measurement", 
        names = c("Anoxic", "Oxic"),
        cex.lab = 1.5, cex.axis = 1.5)
title("A      Station 1 Shannon Diversity", adj = 0, cex.main = 2)
boxplot(E1.1, E3.1, 
        ylab = "Alpha Diversity Measurement", 
        names = c("Anoxic", "Oxic"),
        cex.lab = 1.5, cex.axis = 1.5)
title("B      Station 1 Evenness", adj = 0, cex.main = 2)
boxplot(H1.2, H2.2, H3.2,
        ylab = "Alpha Diversity Measurement", 
        names = c("Anoxic", "Hypoxic", "Oxic"),
        cex.lab = 1.5, cex.axis = 1.5)
title("C      Station 2 Shannon Diversity", adj = 0, cex.main = 2)
boxplot(E1.2, E2.2, E3.2, 
        ylab = "Alpha Diversity Measurement", 
        names = c("Anoxic", "Hypoxic", "Oxic"),
        cex.lab = 1.5, cex.axis = 1.5)
title("D      Station 2 Evenness", adj = 0, cex.main = 2)
boxplot(H1.3, H2.3, H3.3,
        ylab = "Alpha Diversity Measurement", 
        names = c("Anoxic", "Hypoxic", "Oxic"),
        cex.lab = 1.5, cex.axis = 1.5)
title("E      Station 3 Shannon Diversity", adj = 0, cex.main = 2)
boxplot(E1.3, E2.3, E3.3, 
        ylab = "Alpha Diversity Measurement", 
        names = c("Anoxic", "Hypoxic", "Oxic"),
        cex.lab = 1.5, cex.axis = 1.5)
title("F      Station 3 Evenness", adj = 0, cex.main = 2)
dev.off()
