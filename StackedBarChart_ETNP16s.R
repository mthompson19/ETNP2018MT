# Stacked Bar Chart - Relative Abundance 
## Madeleine Thompson 
### ETNP 16s data 

### Install Packages 
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
install_phyloseq(branch = "devel")
install.packages("ggplot2", "dplyr", "tidyr", "ape", "grid", "gridExtra", "vegan", "ggpubr", "ecolTest", "BiocManager")

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

# Import data 
zotu16s <- read.csv("data/zotutab.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)
tax16s <- read.csv("data/taxonomy.csv",header = TRUE, row.names = 1)
metadata.phy = phyloseq::sample_data(data.frame(metadata))



#Separate taxa into separate columns and remove beginning string
taxGroup <- tidyr::separate(data = tax16s, col =  Taxon,
                           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                           sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")
taxGroup$Domain<-gsub("d__","",as.character(taxGroup$Domain))
taxGroup$Phylum<-gsub("p__","",as.character(taxGroup$Phylum))
taxGroup$Class<-gsub("c__","",as.character(taxGroup$Class))
taxGroup$Order <- gsub("o__", "", as.character(taxGroup$Order))
taxGroup$Family<-gsub("f__","",as.character(taxGroup$Family))
taxGroup$Genus<-gsub("g__","",as.character(taxGroup$Genus))
taxGroup$Species<-gsub("s__","",as.character(taxGroup$Species))

############




# Make ZOTU table numeric and a matrix 
# Make Taxonomy a matrix and make the row and column names the same 
mat<- as.matrix(zotu16s)
zotuTab <- matrix(as.numeric(mat),
                  ncol = ncol(mat));
colnames(zotuTab) <- colnames(mat);
rownames(zotuTab) <- rownames(mat);

taxTab <- as.matrix(taxGroup)

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

#Phyloseq 
physeq_ns <- phyloseq::filter_taxa(physeq1, function(x) sum(x) > 1, prune = TRUE) 
phylor = phyloseq::transform_sample_counts(physeq_ns, function(x) x / sum(x) )
phylor_nz <- phyloseq::filter_taxa(phylor, function(x) mean(x) > 0, prune = TRUE)
top <- names(sort(phyloseq::taxa_sums(phylor_nz), decreasing=TRUE))[1:500]
phylor.top <- phyloseq::prune_taxa(top, phylor_nz)


############

# Separate by filter size and station
p.02 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$Filter == 0.2)
p1.02 <- phyloseq::subset_samples(p.02, p.02@sam_data$Station == "Stn1")
p2.02 <- phyloseq::subset_samples(p.02, p.02@sam_data$Station == "Stn2")
p3.02 <- phyloseq::subset_samples(p.02, p.02@sam_data$Station == "Stn3")



p.2 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$Filter == 2.0)
p1.2 <- phyloseq::subset_samples(p.2, p.2@sam_data$Station == "Stn1")
p2.2 <- phyloseq::subset_samples(p.2, p.2@sam_data$Station == "Stn2")
p3.2 <- phyloseq::subset_samples(p.2, p.2@sam_data$Station == "Stn3")


p.22 <- phyloseq::subset_samples(phylor.top, phylor.top@sam_data$Filter == 22)
p1.22 <- phyloseq::subset_samples(p.22, p.22@sam_data$Station == "Stn1")
p2.22 <- phyloseq::subset_samples(p.22, p.22@sam_data$Station == "Stn2")
p3.22 <- phyloseq::subset_samples(p.22, p.22@sam_data$Station == "Stn3")


############


mycolors <- c("#F8766D",
              "#BB9D00",
              "#00A5FF",
              "#00C0B8",
              "#00B81F", 
              "#E76BF3", 
              "#FF6C90", 
              
              "#FF61C9",
              "#AC88FF", 
              "#00B8E5",
              "#00BF7D", 
              "#85AD00",
              "#E08B00",
              
              "#7997FF",
              "#00C08D",
              "#DC71FA",
              "#FF689F", 
              "#ED8141")


# Bar plots per filter size per station Phylum 
p1.02fig <- phyloseq::plot_bar(p1.02,  fill="Phylum", title = expression(paste("Station 1 - 0.2 ", mu, "m")))
stat1.02um <- plot(p1.02fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     scale_fill_manual(values = mycolors) + 
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.2fig <- phyloseq::plot_bar(p1.2,  fill="Phylum", title = expression(paste("Station 1 - 2 ", mu, "m")))
stat1.2um <- plot(p1.2fig + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylab("Relative Abundance") +
                    xlab("Depth") +
                    scale_fill_manual(values = mycolors) +
                    scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.22fig <- phyloseq::plot_bar(p1.22,  fill="Phylum", title = expression(paste("Station 1 - 22 ", mu, "m")))
stat1.22um <- plot(p1.22fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_fill_manual(values = mycolors) +
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p2.02fig <- phyloseq::plot_bar(p2.02,  fill="Phylum", title = expression(paste("Station 2 - 0.2 ", mu, "m")))
stat2.02um <- plot(p2.02fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_fill_manual(values = mycolors) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.2fig <- phyloseq::plot_bar(p2.2,  fill="Phylum", title = expression(paste("Station 2 - 2 ", mu, "m")))
stat2.2um <- plot(p2.2fig + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylab("Relative Abundance") +
                    xlab("Depth") +
                    scale_fill_manual(values = mycolors) +
                    scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.22fig <- phyloseq::plot_bar(p2.22,  fill="Phylum", title = expression(paste("Station 2 - 22 ", mu, "m")))
stat2.22um <- plot(p2.22fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +                    
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_fill_manual(values = mycolors) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p3.02fig <- phyloseq::plot_bar(p3.02,  fill="Phylum", title = expression(paste("Station 3 - 0.2 ", mu, "m")))
stat3.02um <- plot(p3.02fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_fill_manual(values = mycolors) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.2fig <- phyloseq::plot_bar(p3.2, fill="Phylum", title = expression(paste("Station 3 - 2 ", mu, "m")))
stat3.2um <- plot(p3.2fig + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylab("Relative Abundance") +
                    xlab("Depth") +
                    scale_fill_manual(values = mycolors) +
                    scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.22fig <- phyloseq::plot_bar(p3.22,  fill="Phylum", title = expression(paste("Station 3 - 22 ", mu, "m")))
stat3.22um <- plot(p3.22fig + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 25, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_fill_manual(values = mycolors) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))

#put together 
figure <- ggpubr::ggarrange(stat1.02um, stat1.2um, stat1.22um, 
                            stat2.02um, stat2.2um, stat2.22um, 
                            stat3.02um, stat3.2um, stat3.22um, 
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                            font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 3,
                            common.legend = TRUE, legend= "bottom")
ggplot2::ggsave("figures/stnPhylum3.png", plot=figure, device="png",
                scale=1, width = 100, height=45, units=c("cm"), dpi=300, limitsize = FALSE)





# Bar plots per filter size per station Proteobacteria

mycolors1 <- c("#00BA42", "#00BFC4")


p1.02proteo <- phyloseq::subset_taxa(p1.02, Phylum==" Proteobacteria")
p1.2proteo <- phyloseq::subset_taxa(p1.2, Phylum==" Proteobacteria")
p1.22proteo <- phyloseq::subset_taxa(p1.22, Phylum==" Proteobacteria")

p2.02proteo <- phyloseq::subset_taxa(p2.02, Phylum==" Proteobacteria")
p2.2proteo <- phyloseq::subset_taxa(p2.2, Phylum==" Proteobacteria")
p2.22proteo <- phyloseq::subset_taxa(p2.22, Phylum==" Proteobacteria")

p3.02proteo <- phyloseq::subset_taxa(p3.02, Phylum==" Proteobacteria")
p3.2proteo <- phyloseq::subset_taxa(p3.2, Phylum==" Proteobacteria")
p3.22proteo <- phyloseq::subset_taxa(p3.22, Phylum==" Proteobacteria")



p1.02fig1 <- phyloseq::plot_bar(p1.02proteo,  fill="Class", title = expression(paste("Station 1 - 0.2 ", mu, "m")))
stat1.02um1 <- plot(p1.02fig1 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                    scale_y_continuous(limits=c(0,0.5)) +
                     ylab("Relative Abundance") +
                      scale_fill_manual(values = mycolors1) +
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.2fig1 <- phyloseq::plot_bar(p1.2proteo,  fill="Class", title = expression(paste("Station 1 - 2 ", mu, "m")))
stat1.2um1 <- plot(p1.2fig1 + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                     scale_y_continuous(limits=c(0,0.5)) +
                    coord_flip() + 
                    ylab("Relative Abundance") +
                    xlab("Depth") +
                     scale_fill_manual(values = mycolors1) +
                    scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.22fig1 <- phyloseq::plot_bar(p1.22proteo,  fill="Class", title = expression(paste("Station 1 - 22 ", mu, "m")))
stat1.22um1 <- plot(p1.22fig1 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_fill_manual(values = mycolors1) +
                     xlab("Depth") +
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p2.02fig1 <- phyloseq::plot_bar(p2.02proteo,  fill="Class", title = expression(paste("Station 2 - 0.2 ", mu, "m")))
stat2.02um1 <- plot(p2.02fig1 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_fill_manual(values = mycolors1) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.2fig1 <- phyloseq::plot_bar(p2.2proteo,  fill="Class", title = expression(paste("Station 2 - 2 ", mu, "m")))
stat2.2um1 <- plot(p2.2fig1 + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylab("Relative Abundance") +
                    xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_fill_manual(values = mycolors1) +
                    scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.22fig1 <- phyloseq::plot_bar(p2.22proteo,  fill="Class", title = expression(paste("Station 2 - 22 ", mu, "m")))
stat2.22um1 <- plot(p2.22fig1 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +                    
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_fill_manual(values = mycolors1) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p3.02fig1 <- phyloseq::plot_bar(p3.02proteo,  fill="Class", title = expression(paste("Station 3 - 0.2 ", mu, "m")))
stat3.02um1 <- plot(p3.02fig1 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_fill_manual(values = mycolors1) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.2fig1 <- phyloseq::plot_bar(p3.2proteo, fill="Class", title = expression(paste("Station 3 - 2 ", mu, "m")))
stat3.2um1 <- plot(p3.2fig1 + 
                    theme_classic()+
                    theme(text=element_text(size=25), 
                          axis.text.x = element_text(size = 25),
                          axis.text.y = element_text(size = 25),
                          legend.title = element_text(size = 40, face = "bold"), 
                          legend.text = element_text(size = 25),
                          legend.key.width = unit(5, "cm")) +
                    coord_flip() + 
                    ylab("Relative Abundance") +
                    xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_fill_manual(values = mycolors1) +
                    scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.22fig1 <- phyloseq::plot_bar(p3.22proteo,  fill="Class", title = expression(paste("Station 3 - 22 ", mu, "m")))
stat3.22um1 <- plot(p3.22fig1 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 25, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                      scale_fill_manual(values = mycolors1) +
                     xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))

#put together 
figure <- ggpubr::ggarrange(stat1.02um1, stat1.2um1, stat1.22um1, 
                            stat2.02um1, stat2.2um1, stat2.22um1, 
                            stat3.02um1, stat3.2um1, stat3.22um1, 
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                            font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 3,
                            common.legend = TRUE, legend= "bottom")
ggplot2::ggsave("figures/stnProteo.png", plot=figure, device="png",
                scale=1, width = 100, height=45, units=c("cm"), dpi=300, limitsize = FALSE)

# Bar plots per filter size per station Alphaproteobacteria


p1.02A <- phyloseq::subset_taxa(p1.02proteo, Class==" Alphaproteobacteria")
p1.2A <- phyloseq::subset_taxa(p1.2proteo, Class==" Alphaproteobacteria")
p1.22A <- phyloseq::subset_taxa(p1.22proteo, Class==" Alphaproteobacteria")

p2.02A <- phyloseq::subset_taxa(p2.02proteo, Class==" Alphaproteobacteria")
p2.2A <- phyloseq::subset_taxa(p2.2proteo, Class==" Alphaproteobacteria")
p2.22A <- phyloseq::subset_taxa(p2.22proteo, Class==" Alphaproteobacteria")

p3.02A <- phyloseq::subset_taxa(p3.02proteo, Class==" Alphaproteobacteria")
p3.2A <- phyloseq::subset_taxa(p3.2proteo, Class==" Alphaproteobacteria")
p3.22A <- phyloseq::subset_taxa(p3.22proteo, Class==" Alphaproteobacteria")




p1.02figA <- phyloseq::plot_bar(p1.02A,  fill="Order", title = expression(paste("Station 1 - 0.2 ", mu, "m")))
stat1.02umA <- plot(p1.02figA + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.2figA <- phyloseq::plot_bar(p1.2A,  fill="Order", title = expression(paste("Station 1 - 2 ", mu, "m")))
stat1.2umA <- plot(p1.2figA + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.22figA <- phyloseq::plot_bar(p1.22A,  fill="Order", title = expression(paste("Station 1 - 22 ", mu, "m")))
stat1.22umA <- plot(p1.22figA + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      xlab("Depth") +
                      scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p2.02figA <- phyloseq::plot_bar(p2.02A,  fill="Order", title = expression(paste("Station 2 - 0.2 ", mu, "m")))
stat2.02umA <- plot(p2.02figA + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.2figA <- phyloseq::plot_bar(p2.2A,  fill="Order", title = expression(paste("Station 2 - 2 ", mu, "m")))
stat2.2umA <- plot(p2.2figA + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.22figA <- phyloseq::plot_bar(p2.22A,  fill="Order", title = expression(paste("Station 2 - 22 ", mu, "m")))
stat2.22umA <- plot(p2.22figA + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +                    
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p3.02figA <- phyloseq::plot_bar(p3.02A,  fill="Order", title = expression(paste("Station 3 - 0.2 ", mu, "m")))
stat3.02umA <- plot(p3.02figA + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.2figA <- phyloseq::plot_bar(p3.2A, fill="Order", title = expression(paste("Station 3 - 2 ", mu, "m")))
stat3.2umA <- plot(p3.2figA + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.22figA <- phyloseq::plot_bar(p3.22A,  fill="Order", title = expression(paste("Station 3 - 22 ", mu, "m")))
stat3.22umA <- plot(p3.22figA + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 25, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))

#put together 
figure <- ggpubr::ggarrange(stat1.02umA, stat1.2umA, stat1.22umA, 
                            stat2.02umA, stat2.2umA, stat2.22umA, 
                            stat3.02umA, stat3.2umA, stat3.22umA, 
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                            font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 3,
                            common.legend = TRUE, legend= "bottom")
ggplot2::ggsave("figures/stnAProteoOrder.png", plot=figure, device="png",
                scale=1, width = 100, height=45, units=c("cm"), dpi=300, limitsize = FALSE)




p1.02G <- phyloseq::subset_taxa(p1.02proteo, Class==" Gammaproteobacteria")
p1.2G <- phyloseq::subset_taxa(p1.2proteo, Class==" Gammaproteobacteria")
p1.22G <- phyloseq::subset_taxa(p1.22proteo, Class==" Gammaproteobacteria")

p2.02G <- phyloseq::subset_taxa(p2.02proteo, Class==" Gammaproteobacteria")
p2.2G <- phyloseq::subset_taxa(p2.2proteo, Class==" Gammaproteobacteria")
p2.22G <- phyloseq::subset_taxa(p2.22proteo, Class==" Gammaproteobacteria")

p3.02G <- phyloseq::subset_taxa(p3.02proteo, Class==" Gammaproteobacteria")
p3.2G <- phyloseq::subset_taxa(p3.2proteo, Class==" Gammaproteobacteria")
p3.22G <- phyloseq::subset_taxa(p3.22proteo, Class==" Gammaproteobacteria")




p1.02figG <- phyloseq::plot_bar(p1.02G,  fill="Order", title = expression(paste("Station 1 - 0.2 ", mu, "m")))
stat1.02umG <- plot(p1.02figG + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.2figG <- phyloseq::plot_bar(p1.2G,  fill="Order", title = expression(paste("Station 1 - 2 ", mu, "m")))
stat1.2umG <- plot(p1.2figG + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.22figG <- phyloseq::plot_bar(p1.22G,  fill="Order", title = expression(paste("Station 1 - 22 ", mu, "m")))
stat1.22umG <- plot(p1.22figG + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p2.02figG <- phyloseq::plot_bar(p2.02G,  fill="Order", title = expression(paste("Station 2 - 0.2 ", mu, "m")))
stat2.02umG <- plot(p2.02figG + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.2figG <- phyloseq::plot_bar(p2.2G,  fill="Order", title = expression(paste("Station 2 - 2 ", mu, "m")))
stat2.2umG <- plot(p2.2figG + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.22figG <- phyloseq::plot_bar(p2.22G,  fill="Order", title = expression(paste("Station 2 - 22 ", mu, "m")))
stat2.22umG <- plot(p2.22figG + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +                    
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p3.02figG <- phyloseq::plot_bar(p3.02G,  fill="Order", title = expression(paste("Station 3 - 0.2 ", mu, "m")))
stat3.02umG <- plot(p3.02figG + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.2figG <- phyloseq::plot_bar(p3.2G, fill="Order", title = expression(paste("Station 3 - 2 ", mu, "m")))
stat3.2umG <- plot(p3.2figG + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.22figG <- phyloseq::plot_bar(p3.22G,  fill="Order", title = expression(paste("Station 3 - 22 ", mu, "m")))
stat3.22umG <- plot(p3.22figG + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 25, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))

#put together 
figure <- ggpubr::ggarrange(stat1.02umG, stat1.2umG, stat1.22umG, 
                            stat2.02umG, stat2.2umG, stat2.22umG, 
                            stat3.02umG, stat3.2umG, stat3.22umG, 
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                            font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 3,
                            common.legend = TRUE, legend= "bottom")
ggplot2::ggsave("figures/stnGProteoOrder.png", plot=figure, device="png",
                scale=1, width = 100, height=45, units=c("cm"), dpi=300, limitsize = FALSE)



# Figure for Cyanobacteria 

p1.02cyano <- phyloseq::subset_taxa(p1.02, Phylum==" Cyanobacteria")
p1.2cyano <- phyloseq::subset_taxa(p1.2, Phylum==" Cyanobacteria")
p1.22cyano <- phyloseq::subset_taxa(p1.22, Phylum==" Cyanobacteria")

p2.02cyano <- phyloseq::subset_taxa(p2.02, Phylum==" Cyanobacteria")
p2.2cyano <- phyloseq::subset_taxa(p2.2, Phylum==" Cyanobacteria")
p2.22cyano <- phyloseq::subset_taxa(p2.22, Phylum==" Cyanobacteria")

p3.02cyano <- phyloseq::subset_taxa(p3.02, Phylum==" Cyanobacteria")
p3.2cyano <- phyloseq::subset_taxa(p3.2, Phylum==" Cyanobacteria")
p3.22cyano <- phyloseq::subset_taxa(p3.22, Phylum==" Cyanobacteria")







p1.02fig2 <- phyloseq::plot_bar(p1.02cyano,  fill="Genus", title = expression(paste("Station 1 - 0.2 ", mu, "m")))
stat1.02um2 <- plot(p1.02fig2 + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      scale_y_continuous(limits=c(0,0.5)) +
                      ylab("Relative Abundance") +
                      scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.2fig2 <- phyloseq::plot_bar(p1.2cyano,  fill="Genus", title = expression(paste("Station 1 - 2 ", mu, "m")))
stat1.2um2 <- plot(p1.2fig2 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     scale_y_continuous(limits=c(0,0.5)) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p1.22fig2 <- phyloseq::plot_bar(p1.22cyano,  fill="Genus", title = expression(paste("Station 1 - 22 ", mu, "m")))
stat1.22um2 <- plot(p1.22fig2 + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      xlab("Depth") +
                      scale_x_discrete(labels = c("120m", "110m", "83m", "70m", "40m")))
p2.02fig2 <- phyloseq::plot_bar(p2.02cyano,  fill="Genus", title = expression(paste("Station 2 - 0.2 ", mu, "m")))
stat2.02um2 <- plot(p2.02fig2 + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.2fig2 <- phyloseq::plot_bar(p2.2cyano,  fill="Genus", title = expression(paste("Station 2 - 2 ", mu, "m")))
stat2.2um2 <- plot(p2.2fig2 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p2.22fig2 <- phyloseq::plot_bar(p2.22cyano,  fill="Genus", title = expression(paste("Station 2 - 22 ", mu, "m")))
stat2.22um2 <- plot(p2.22fig2 + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +                    
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("900m", "120m", "113m", "90m", "30m")))
p3.02fig2 <- phyloseq::plot_bar(p3.02cyano,  fill="Genus", title = expression(paste("Station 3 - 0.2 ", mu, "m")))
stat3.02um2 <- plot(p3.02fig2 + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 40, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.2fig2 <- phyloseq::plot_bar(p3.2cyano, fill="Genus", title = expression(paste("Station 3 - 2 ", mu, "m")))
stat3.2um2 <- plot(p3.2fig2 + 
                     theme_classic()+
                     theme(text=element_text(size=25), 
                           axis.text.x = element_text(size = 25),
                           axis.text.y = element_text(size = 25),
                           legend.title = element_text(size = 40, face = "bold"), 
                           legend.text = element_text(size = 25),
                           legend.key.width = unit(5, "cm")) +
                     coord_flip() + 
                     ylab("Relative Abundance") +
                     xlab("Depth") +
                     scale_y_continuous(limits=c(0,0.5)) +
                     scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))
p3.22fig2 <- phyloseq::plot_bar(p3.22cyano,  fill="Genus", title = expression(paste("Station 3 - 22 ", mu, "m")))
stat3.22um2 <- plot(p3.22fig2 + 
                      theme_classic()+
                      theme(text=element_text(size=25), 
                            axis.text.x = element_text(size = 25),
                            axis.text.y = element_text(size = 25),
                            legend.title = element_text(size = 25, face = "bold"), 
                            legend.text = element_text(size = 25),
                            legend.key.width = unit(5, "cm")) +
                      coord_flip() + 
                      ylab("Relative Abundance") +
                      xlab("Depth") +
                      scale_y_continuous(limits=c(0,0.5)) +
                      scale_x_discrete(labels = c("1000m", "70m", "45m", "40m", "33m", "10m")))

#put together 
figure <- ggpubr::ggarrange(stat1.02um2, stat1.2um2, stat1.22um2, 
                            stat2.02um2, stat2.2um2, stat2.22um2, 
                            stat3.02um2, stat3.2um2, stat3.22um2, 
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                            font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                            ncol = 3, nrow = 3,
                            common.legend = TRUE, legend= "bottom")
ggplot2::ggsave("figures/stnCyanoGen.png", plot=figure, device="png",
                scale=1, width = 100, height=45, units=c("cm"), dpi=300, limitsize = FALSE)