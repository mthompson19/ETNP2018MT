#NMDS Plots 
## Madeleine Thompson 
### ETNP 2018 16s data

## Import data 
zotu16s <- read.csv("data/zotutab.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1)

############
library("ggplot2"); packageVersion("ggplot2")
library(dplyr)
library(tidyr)
library(vegan)
library(ggpubr)
library(ggrepel)
library("viridis")


############
## NMDS 

nmdsmeta <- read.csv("data/nmdsmeta.csv", header = TRUE)

nmdsmeta <- na.omit(nmdsmeta)

rownames(nmdsmeta) <- nmdsmeta[,1]

nmdsmeta <- nmdsmeta[,-1]


zotu <- t(zotu16s)

zotutab.chord <- vegan::decostand(zotu, "normalize")

ord <- vegan::metaMDS(zotutab.chord)
summary(ord)
(fit <- vegan::envfit(ord, nmdsmeta, perm = 999))
vegan::scores(fit, "vectors")
plot(ord)
plot(fit)
plot(fit, p.max = 0.05, col = "red")


### remove the outlying station 
zotu1 <- zotu[-c(44),]

nmdsmeta1 <- nmdsmeta[-c(44),]
zotu_nmds_meta = merge(ord1$points, nmdsmeta1, by = 0)

nmdsmeta2 <- nmdsmeta1[,-c(1:2)]

nmdsmeta3 <- nmdsmeta2[,-c(4:5)]

zotutab1.chord <- vegan::decostand(zotu1, "normalize")

ord1 <- vegan::metaMDS(zotutab1.chord)
summary(ord1)

(fit <- vegan::envfit(ord1, nmdsmeta1, perm = 999))
vegan::scores(fit, "vectors")
plot(ord1)
plot(fit)
plot(fit, p.max = 0.05, col = "red")



en = vegan::envfit(ord1, nmdsmeta3, permutations = 999, na.rm = TRUE)


en_coord_cont = as.data.frame(vegan::scores(en, "vectors")) * vegan::ordiArrowMul(en)
en_coord_cat = as.data.frame(vegan::scores(en, "factors")) * vegan::ordiArrowMul(en)

gg = ggplot2::ggplot(data = zotu_nmds_meta, ggplot2::aes(x = MDS1, y = MDS2)) + 
  ggplot2::geom_point(data = zotu_nmds_meta, ggplot2::aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station, size = Filter2)) +
  ggplot2::geom_segment(data = en_coord_cont, ggplot2::aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size =1, alpha = 0.5, colour = "black") + 
  ggplot2::theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) + 
  labs(colour = "Oxygen") + 
  xlab("NMDS1") + ylab("NMDS2") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 50)) + 
  viridis::scale_color_viridis(option = "C")+
  theme_classic()  

gg

zotu.nmds <- ggplot2::ggplot(data=zotu_nmds_meta, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station, size = Filter2)) +
  xlab("NMDS1") + ylab("NMDS2") + 
  scale_size_area(max_size = 5) + 
  theme(text = element_text(size = 50)) + 
  theme_classic() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "black") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", label = row.names(en_coord_cont)) 
zotu.nmds
zotu.nmds1 <- ggplot2::ggplot(data=zotu_nmds_meta, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station, size = Filter2)) +
  geom_label_repel(aes(label = Depth),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   max.overlaps = 40) + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 50)) + 
  scale_color_viridis(option = "C")+
  theme_classic() 
zotu.nmds1

nmdsfig1 <- ggpubr::ggarrange(zotu.nmds, zotu.nmds1, 
                              labels = c("A", "B"), 
                              font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                              ncol = 2, nrow = 1,
                              common.legend = TRUE, legend = c("right"))

ggplot2::ggsave("figures/zotu.nmds.png", gg, device="png",
                scale=1, width = 18, height=12, units=c("cm"), dpi=300, limitsize = FALSE)
### do the singles oxygen conc. then set filter to continuous, do nmds for each filter size 

nmds.filt.22 <- zotu_nmds_meta %>%
  filter(Filter == 22)
nmds.filt.2 <- zotu_nmds_meta %>%
  filter(Filter == 2)
nmds.filt.02 <- zotu_nmds_meta %>%
  filter(Filter == 0.2)

nmdsmeta1 <- nmdsmeta[-c(44),]
zotu_nmds_meta = merge(ord1$points, nmdsmeta1, by = 0)

nmdsmeta2 <- nmdsmeta1[,-c(1:2)]

nmdsmeta3 <- nmdsmeta2[,-c(4)]

en = envfit(ord1, nmdsmeta3, permutations = 999, na.rm = TRUE)


en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

zotu.nmds22 <- ggplot2::ggplot(data=nmds.filt.22, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station), size = 5) +
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size =1, alpha = 0.5, colour = "black") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "black")) + 
  labs(colour = "Oxygen") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 50)) + 
  scale_color_viridis(option = "C")+
  theme_classic()  
zotu.nmds22


zotu.nmds2 <- ggplot2::ggplot(data=nmds.filt.2, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station), size = 5) +
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size =1, alpha = 0.5, colour = "black") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "black")) + 
  labs(colour = "Oxygen") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 50)) + 
  scale_color_viridis(option = "C")+
  theme_classic()  
zotu.nmds2
zotu.nmds02 <- ggplot2::ggplot(data=nmds.filt.02, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station), size = 5) +
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size =1, alpha = 0.5, colour = "black") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "black")) + 
  labs(colour = "Oxygen") + 
  scale_size_area(max_size = 5) +
  theme(text = element_text(size = 50)) + 
  scale_color_viridis(option = "C")+
  theme_classic()  
zotu.nmds02


zotu.nmds22.1 <- ggplot2::ggplot(data=nmds.filt.22, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station), size = 5) +
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_label_repel(aes(label = Depth),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   max.overlaps = 10) + 
  scale_size_area(max_size = 5) + 
  theme(text=element_text(size=50), 
        strip.text.x = element_text(size=rel(50)),
        axis.text=element_text(size=50), 
        legend.key.width = unit(10, "cm")) +
  scale_color_viridis(option = "C")+
  ggtitle(expression(paste("Filter 22 ", mu, "m"))) +
  theme_classic() 
zotu.nmds22.1
zotu.nmds2.1 <- ggplot2::ggplot(data=nmds.filt.2, aes(x=MDS1, y=MDS2)) +
  ggplot2::geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station), size = 5) +
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_label_repel(aes(label = Depth),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   max.overlaps = 10) + 
  scale_size_area(max_size = 5) + 
  theme(text=element_text(size=50), 
        axis.text.x = element_text(size = 50),
        axis.text.y = element_text(size = 50),
        legend.title = element_text(size = 50, face = "bold"), 
        legend.text = element_text(size = 50),
        legend.key.width = unit(10, "cm")) +
  scale_color_viridis(option = "C")+
  ggtitle(expression(paste("Filter 2 ", mu, "m"))) +
  theme_classic() 
zotu.nmds2.1
zotu.nmds02.1 <- ggplot2::ggplot(data=nmds.filt.02, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color=Oxygen2, shape = Station), size = 5) +
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_label_repel(aes(label = Depth),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   max.overlaps = 10) + 
  scale_size_area(max_size = 5) + 
  theme(text=element_text(size=50), 
        legend.key.width = unit(10, "cm")) +
  scale_color_viridis(option = "C")+
  ggtitle(expression(paste("Filter 0.2 ", mu, "m"))) +
  theme_classic() 
zotu.nmds02.1

nmdsfig <- ggpubr::ggarrange(zotu.nmds02, zotu.nmds2, zotu.nmds22, 
                             zotu.nmds02.1, zotu.nmds2.1, zotu.nmds22.1,
                             labels = c("A", "B", "C", "D", "E", "F"), 
                             font.label = list(size = 15, color = "black", face = "bold", family = NULL),
                             ncol = 3, nrow = 2,
                             common.legend = TRUE, legend = c("right"))

ggplot2::ggsave("figures/filtnmds.png", plot=nmdsfig, device="png",
                scale=1, width = 25, height=15, units=c("cm"), dpi=300, limitsize = FALSE)

############

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



#PERMANOVA 
pseq.rel <- microbiome::transform(phylor.top, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

permanova1 <- vegan::adonis2(t(otu) ~ Oxygen,
                            data = meta, permutations=9999, method = "bray")

permanova2 <- vegan::adonis2(t(otu) ~ Filter,
                            data = meta, permutations=9999, method = "bray")
print(as.data.frame(permanova2$aov.tab))

permanova3 <- vegan::adonis2(t(otu) ~ Station,
                            data = meta, permutations=9999, method = "bray")
print(as.data.frame(permanova3$aov.tab))

permanova4 <- vegan::adonis2(t(otu) ~ Temperture,
                            data = meta, permutations=9999, method = "bray")
print(as.data.frame(permanova4$aov.tab))

permanova5 <- vegan::adonis2(t(otu) ~ Salinity,
                            data = meta, permutations=9999, method = "bray")
print(as.data.frame(permanova5$aov.tab))

permanova6 <- vegan::adonis2(t(otu) ~ Nitrate,
                            data = meta, permutations=9999, method = "bray")
print(as.data.frame(permanova6$aov.tab))


