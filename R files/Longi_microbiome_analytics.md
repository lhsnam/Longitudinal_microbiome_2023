---
title: "Diarrhoea_longitudinal_analytics_fin"
author: "Son-Nam H. Le"
date: "`r Sys.Date()`"
output: html_document
---

# load packages:

```{r}
library(treeio) #read newick
library(phyloseq) #phyloseq object
library(ggplot2) #visualize
library(ggpubr) #add stats
library(ape) #phylogenetic tree
library(microViz) #tax table
library(microbiome)
library(ranacapa)
library(tidyr)
library(philr)
library(zCompositions)
library(mbImpute)
library(Wrench)
library(microbiomeutilities)
library(vegan)
library(cols4all)
library(dendextend)
library(fpc) #clustering
library(cluster) #clustering
library(factoextra)
library(randomForestSRC)
library(ggRandomForests)
library(scales)
library(data.table)
library(ANCOMBC)
library(stringr)
library(DESeq2)
library(zinbwave)
library(tibble)
library(edgeR)
library(ggrepel)
library(NetCoMi) #network
library(knitr)
library(finalfit)
library(kableExtra)
library(Maaslin2)
library(flextable)
library(lme4)
library(lmerTest)
```

# Preparation

```{r}
load(file = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/phyloseq_select_1904-OTUs_fin.RData")
meta.df <- read.table("D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/combined_09AV_04EN_metadata.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
bigtree.unroot <- read.newick(file = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/select_filter_1904-OTUs_pasta_trimal_gt5%gap.fasta.reroot.newick")

#remove some patient that do not include any associated (clinical) information
filter.ps <- prune_samples(!(is.na(sample_data(filter.ps)$sex)), filter.ps)
filter.ps <- prune_samples(!(is.na(sample_data(filter.ps)$wfa_zscore)), filter.ps)
filter.ps <- prune_taxa(taxa_sums(filter.ps) > 0, filter.ps)
filter.ps  #410 samples 1886 taxa

#convert "day" to factor/ "NA" to "Control"
sample_data(filter.ps)$day <- factor(sample_data(filter.ps)$day, levels=c('1', '7', '14', 'Control'))
filter.ps@sam_data[["day"]] <- as.factor(filter.ps@sam_data[["day"]])
sample_data(filter.ps)$day <- sample_data(filter.ps)$day %>% replace_na("Control")

#oral phyloseq
oral.df <- read.csv("D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/Blast/filtered_oral_otus.csv")
otu.oral <- oral.df$Query.id 
otu.df <- as.data.frame(otu_table(filter.ps))
oral.subset <- subset(otu.df, row.names(otu.df) %in% otu.oral)
otus.oral <- otu_table(oral.subset, taxa_are_rows = TRUE)
oral.ps <- merge_phyloseq(otus.oral, tax_table(filter.ps), sample_data(filter.ps))
oral.tre <- read.newick(file = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/Raxml_3rd/oral/oral_tree_152otus.newick")
oral.ps <- merge_phyloseq(oral.ps, oral.tre)

#Change "day" as factor
sample_data(oral.ps)$day <- factor(sample_data(oral.ps)$day, levels=c('1', '7', '14', 'Control'))
oral.ps@sam_data[["day"]] <- as.factor(oral.ps@sam_data[["day"]])
oral.ps@sam_data[["day"]][is.na(oral.ps@sam_data[["day"]])] = "Control"

#gut phyloseq
gut.subset <- subset(otu.df, !(row.names(otu.df) %in% otu.oral))
otus.gut <- otu_table(gut.subset, taxa_are_rows = T)
gut.ps <- merge_phyloseq(otus.gut, tax_table(filter.ps), sample_data(filter.ps))
gut.tre <- read.newick(file = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/Raxml_3rd/gut/gut_1740otus.newick")
gut.ps <- merge_phyloseq(gut.ps, gut.tre)

#Change "day" as factor
sample_data(gut.ps)$day <- factor(sample_data(gut.ps)$day, levels=c('1', '7', '14', 'Control'))
gut.ps@sam_data[["day"]] <- as.factor(gut.ps@sam_data[["day"]])
gut.ps@sam_data[["day"]][is.na(gut.ps@sam_data[["day"]])] = "Control"

#Review gut and oral phyloseq file
gut.ps #1734 taxa
oral.ps #152 taxa
#total 1886 taxa = filter.ps's taxa. OK!
```

Phylogenetic trees (2 branches: oral+gut)

```{r}
#Bind 2 trees (Oral + Gut)
#write.tree(phy_tre, file='D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/R work place/new/bigtree.newick')
#write.tree(gut.ps@phy_tree, file='D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/R work place/new/guttree.newick')
#write.tree(oral.ps@phy_tree, file='D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/R work place/new/oraltree.newick')

#open guttree & oraltree then combine it manually as 2 branches of the new phylogenetic tree -> bigtree_rooted.newick
big_tre <- read.newick(file = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/R work place/new/bigtree_rooted.newick")

keep.tips <- which(big_tre$tip.label %in% rownames(otu_table(filter.ps)) == TRUE)
length(keep.tips) #1886 tips
big_tre <- keep.tip(big_tre, keep.tips)
phy_tree(filter.ps) <- big_tre
tax.clean <- as.matrix(tax_table(filter.ps))
tax.clean[,"Genus"] <- stringr::str_remove_all(tax.clean[,"Genus"],"[\\[\\]]")
tax.clean[,"Species"] <- stringr::str_remove_all(tax.clean[,"Species"],"[\\[\\]]")

for (i in 1:nrow(tax.clean)){
  if(tax.clean[i,"Species"] == ""){
    tax.clean[i,"Species"] <- tax.clean[i,"Species"]
  } else{
  tax.clean[i,"Species"] <- paste(substr(tax.clean[i,"Species"], start = 1, stop = 1),
                                  sep = ". ",
                                  word(tax.clean[i,"Species"],-1))
  }
}

for (i in 1:nrow(tax.clean)){
if (tax.clean[i,2] == ""){
kingdom <- paste("unc_", tax.clean[i,1], sep = "")
tax.clean[i, 2:7] <- kingdom
} else if (tax.clean[i,3] == ""){
phylum <- paste("unc_", tax.clean[i,2], sep = "")
tax.clean[i, 3:7] <- phylum
} else if (tax.clean[i,4] == ""){
class <- paste("unc_", tax.clean[i,3], sep = "")
tax.clean[i, 4:7] <- class
} else if (tax.clean[i,5] == ""){
order <- paste("unc_", tax.clean[i,4], sep = "")
tax.clean[i, 5:7] <- order
} else if (tax.clean[i,6] == ""){
family <- paste("unc_", tax.clean[i,5], sep = "")
tax.clean[i, 6:7] <- family
} else if (tax.clean[i,7] == ""){
genus <- paste("unc_", tax.clean[i,6], sep = "")
tax.clean[i, 7] <- genus
}
}

tax.all <- matrix(tax.clean, nrow = nrow(tax.clean), ncol = 8)
colnames(tax.all) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
tax.all[,"OTU"] <- rownames(tax.clean)
rownames(tax.all) <- rownames(tax.clean)

tax_table(filter.ps) <- tax_table(tax.all)
```

## Remove singleton

```{r}
#sort sample total reads, prune taxa
all.ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(filter.ps) > 0, filter.ps)

#explore singleton
all.tx.pres <- rowSums(abundances(all.ps) > 0)
all.singleton <- as.vector(as.matrix(tax_table(all.ps)[names(which(all.tx.pres == 1)),"Genus"]))
all.singleton #793 singleton taxa

names(all.singleton) <- rownames(tax_table(all.ps)[names(which(all.tx.pres == 1)), "Genus"])
sort(table(all.singleton),decreasing = TRUE) ## most common singleton taxa are Bifidobacterium, Bacteroides, Prevotella, Streptococcus, Escherichia, etc.

#counts of each singleton taxa
sort(taxa_sums(all.ps)[names(all.singleton)], decreasing=TRUE)
summary(taxa_sums(all.ps)[names(all.singleton)]) #median 7191 // 3rd quartile 14639

table(all.singleton[names(which(taxa_sums(all.ps)[names(all.singleton)] < 14639))])
all.taxa.rm <- names(which(taxa_sums(all.ps)[names(all.singleton)] < 14639))

## Remove singleton taxa that have read abundance < 14368 reads
all.fil <- prune_taxa(!(taxa_names(all.ps) %in% all.taxa.rm), all.ps)
all.fil #1292 taxa
summary(sample_sums(all.fil)/sample_sums(all.ps)) #median 98.62%

all.ggrare <- ggrare(all.fil, step=100, label=NULL, color='day')
all.ggrare.plot <- all.ggrare +
  xlim(0,3000) + 
  theme(legend.position = 'bottom') +
  theme_classic()

all.ggrare.plot

sort(sample_sums(all.fil))
```

## Imputation for the highly sparse otu table

```{r}
#Prepare materials
sample_data(all.fil) <- sample_data(all.fil)[ , colSums(is.na(sample_data(all.fil))) == 0]
otutab.imp <- as.matrix(unclass(t(otu_table(all.fil)))) #otu table samples(rows)*compositions(columns)
otutab.imp[1:6, 1:6] #check the data => many zeros

meta_data <- as.data.frame(unclass(sample_data(all.fil)), row.names = rownames(sample_data(all.fil)))
day_condition <- meta_data[,5] #pick the vector of day that will play as the condition
meta_data <- meta_data[,-5] #remove "day" vector
meta_data

##check data before run mbImpute (require long computational time)
rownames(meta_data)
day_condition
rownames(otutab.imp)

D <- cophenetic.phylo(phy_tree(all.fil)) #pairwise phylogenetic distance

##mbImpute (temporarily skip)
mbimp.otutab <- mbImpute(condition = day_condition, metadata = meta_data, otu_tab = otutab.imp, D = D, parallel = TRUE, ncores = 20)
imp.otutab <- t(mbimp.otutab$imp_count_mat_origlibsize)

all.norm <- all.fil
otu_table(all.norm) <- otu_table(as.matrix(unclass(imp.otutab)), taxa_are_rows =  TRUE)
otu_table(all.norm)

sample_data(all.norm) <- sample_data(all.ps)
sample_data(all.fil) <- sample_data(all.ps)

```

## Wrench normalization

```{r}
Wrench.norm <- wrench(otu_table(all.norm), condition = sample_data(all.fil)$day)

norm_factors <- Wrench.norm$ccf
head(norm_factors)

#the sweep function to normalize/scale the raw count matrix and generate normalized counts.
norm_counts <- sweep(as.matrix(data.frame(otu_table(all.norm))), 2, norm_factors, FUN = '/')

#Check the sparsity
length(which(abundances(all.fil) ==0))/(phyloseq::nsamples(all.fil)*phyloseq::ntaxa(all.fil)) #97.71%
length(which(abundances(all.norm) ==0))/(phyloseq::nsamples(all.norm)*phyloseq::ntaxa(all.norm)) #94.97%
```
## phILR

```{r}
#zero replacement before running philr
otutab.all.philr <- cmultRepl(otu_table(all.fil),method = "GBM", output = "p-counts")
otutab.all.philr

#Combine tree and ILR calculation
tree.all.philr <- makeNodeLabel(as.phylo(phy_tree(all.fil)), method='number', prefix='n')
gp.all <- philr(t(otutab.all.philr), tree.all.philr, part.weights = 'enorm.x.gm.counts', ilr.weights= 'blw.sqrt')
gp.all.dist <- dist(gp.all, method = "euclidean")
ord_all <- phyloseq::ordinate(all.fil, "PCoA", distance = gp.all.dist)

#Plot scree plot
phyloseq::plot_scree(ord_all) + 
  geom_bar(stat="identity", fill = "firebrick") +
  scale_x_discrete(limits = c(1:10)) +
  labs(x = "Axis", y = "Proportion of Variance")

(ord_plot <- phyloseq::plot_ordination(all.norm, ord_all, type="samples", axes = c(1,2), color = 'day') +
    geom_point(size = 2) +
    scale_color_manual(values = c("#ef476f", "#ffd166", "#26547c", "grey75")) +
    #stat_ellipse(aes(group = day), linetype = 1) +
    #facet_grid(.~study_ID) +
    labs(title = "a") +
  theme_classic() +
    guides(x.sec = "axis",
         y.sec = "axis") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 2/3,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line = element_line(size = .5, color="grey30"),
          axis.text.x.bottom = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y.left =element_text(size=10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
        axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3),
        axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -1)))

ggplot2::ggsave(filename = "PCoA.png", 
       plot = ord_plot,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Beta", 
       width =7, 
       height = 7, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")

(ord_plot2 <- phyloseq::plot_ordination(all.norm, ord_all, type="samples", axes = c(1,2),  color = "study_ID") +
    geom_point(size = 2) +
    scale_color_manual(values = c("#006e90", "#f18f01")) +
    labs(title = "a", subtitle = "Aitchison Distance + PhILR transformed") +
  theme_pubr()
)

ggplot2::ggsave(filename = "PCoA_study.png", 
       plot = ord_plot2,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Beta", 
       width =7, 
       height = 7, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")

head(ord_all$values$Relative_eig, 5) #34.94% and 18.27%
```

# Diversity

## Alpha diversity

```{r}
ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(all.fil),
                         "observed" = phyloseq::estimate_richness(all.fil, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

#plot alpha diversity
day_cols <- c("#BB6457", "#E19A4f", "#f5ed69", "#9B9B9B")

#Alpha diversity by multiple method
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(all.norm, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(all.norm, measures = "Simpson"),
  "Chao1" = phyloseq::estimate_richness(all.norm, measures = "Chao1"),
  "day" = phyloseq::sample_data(all.norm)$day,
  "study" = sample_data(all.norm)$study_ID)

adiv$Chao1 <- adiv$Chao1.Chao1

comps <- make_pairs(sample_data(all.norm)$day)

(adiv.d1.plot <- adiv %>% 
  filter(day == 1) %>%
    ggplot(aes(x = study, y = Shannon)) +
    geom_violin(aes(fill = study), alpha =0.7, show.legend = F) +
    geom_boxplot(width = 0.1) +
    stat_compare_means(method = "wilcox.test", 
                       paired = F)
)

adiv.plot.df <- adiv %>%
  gather(key = metric, value = value, c("Shannon","Simpson", "Chao1")) %>%
  mutate(metric = factor(metric, levels = c("Shannon","Simpson", "Chao1")))

(adiv.plot <- ggplot(adiv.plot.df, aes(x = day, y = value)) + 
    theme_test() +
  introdataviz::geom_split_violin(data = subset(adiv.plot.df, day == 1),
                  aes(x = day, y = value, fill = study), color = NA, alpha = 0.7,
                  inherit.aes = T) +
  geom_violin(data = subset(adiv.plot.df, day != 1), 
              mapping = aes(x = day, y = value, fill = study), color = NA, alpha = 0.7,
              inherit.aes = T) +
  geom_boxplot(data = subset(adiv.plot.df, day == 1), mapping = aes(x = day, y = value, group = study, fill = study),
               color = "#2A2B2DFF",
               width = 0.4,
               outlier.color = NA,
               inherit.aes = T) +
  geom_boxplot(data = subset(adiv.plot.df, day != 1), 
               mapping = aes(x = day, y = value, fill = study),
               color = "#2A2B2DFF",
               width = 0.2, 
               fill = NA,
               outlier.color = NA, 
               inherit.aes = T) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_x_discrete(limits=c("1","7","14", "Control"), 
                   labels =c("1", "7", "14", "C")) +
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) +
  labs(x = "Day",
       y = "", 
       title = "a") +
  theme(legend.position='top', 
        legend.justification = 'center',
        legend.direction='horizontal',
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.15, "in"),
        axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = 1, hjust = 0.5),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 2),
        plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0), 
        aspect.ratio = 1.5,
        plot.margin = unit(c(0.1,0.2,0.1,0.2),"in")
        ) + 
    guides(fill = guide_legend(title = "Study", title.position = "top")) + 
    scale_fill_manual(values = c("#fc776a", "#ffc667"), labels = c("Longitudinal study", "Cross-sectional study")) +
             stat_compare_means(comparisons = comps, 
                     aes(label = paste0("p = ", ..p.format..)),
                     tip.length = 0.02, size = 2,
                     method = "wilcox.test", 
                     hide.ns = F, 
                     vjust = 0,
                     hjust = 0.5,
                     paired = F)
)
ggplot2::ggsave(filename = "Alpha diversity_updated.png", 
       plot = adiv.plot,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Diversity/Alpha", 
       width =7, 
       height = 7, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")

adiv %>%
  group_by(day) %>%
  dplyr::summarise(mean_shannon = mean(Shannon),
            mean_chao1 = mean(Chao1.Chao1),
            mean_simp = mean(Simpson))
```
Relative Abundance (all)
```{r}
std_mean <- function(x) sd(x)/sqrt(length(x))
se <- function(x) sqrt(var(x) / length(x))

day_abund <- microbiome::aggregate_rare(all.norm, level = "Phylum", detection = 0, prevalence = 0)

day_abund <- microbiomeutilities::phy_to_ldf(day_abund, 
                                         transform.counts = "compositional")

day_abund <- day_abund %>% 
  group_by(Phylum, day) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()

day_abund <- as.data.frame(day_abund)

top.phyla <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes")

day_abund <- day_abund %>% dplyr::mutate(pick_Phylum = ifelse(Phylum %in% top.phyla, yes = Phylum, "Others"))

day_abund.control <- day_abund %>% filter(day == "Control")

day_abund <- day_abund %>% filter(day != "Control")

day_abund.fil <- day_abund %>% 
  group_by(pick_Phylum, day) %>%
  summarise(mean_abundance = sum(mean_abundance)) %>% 
  ungroup()

day_abund.control <- day_abund.control %>% 
  group_by(pick_Phylum, day) %>%
  summarise(mean_abundance = sum(mean_abundance)) %>% 
  ungroup()

day_abund.fil$pick_Phylum <- factor(day_abund.fil$pick_Phylum , levels=c("Bacteroidetes", "Actinobacteria", "Firmicutes", "Proteobacteria", "Others"))

day_abund.control$pick_Phylum <- factor(day_abund.control$pick_Phylum , levels=c("Bacteroidetes", "Actinobacteria", "Firmicutes", "Proteobacteria", "Others"))

day1_abund <- subset(day_abund.fil, day == "1")
day7_abund <- subset(day_abund.fil, day == "7")
day14_abund <- subset(day_abund.fil, day == "14")


phyla.col <- c4a("light", 6)
 
(area.plot <- ggplot(data = day_abund.fil) +
    geom_area(aes(x=day, y=mean_abundance, group = pick_Phylum, fill = pick_Phylum), stat = "identity", position = position_fill(reverse = TRUE), color = "white", linewidth = 0.1) +
    geom_col(data = day_abund.control,
                 aes(x=day, y=mean_abundance, group = pick_Phylum, fill = pick_Phylum), position = position_fill(reverse = TRUE), color = "white", width = 0.5, linewidth = 0.1) + 
    geom_text(data = day1_abund, aes(x = 0.94, y = c(0.19, 0.045, 0.55, 0.99, 0.9), label = formatC(round(mean_abundance,2),2,format="f")), hjust = 1, nudge_x = 0.05,  color = "black", size = 3) +
    geom_text(data = day14_abund, aes(x = 3.01, y = c(0.33, 0.1, 0.7, 1, 0.95), hjust = 0, label = formatC(round(mean_abundance,2),2,format="f")), color = "black", size = 3) +
    geom_text(data = day_abund.control, aes(x = 4.3, y = c(0.33, 0.1, 0.7, 1, 0.95), hjust = 0, label = formatC(round(mean_abundance,2),2,format="f")), color = "black", size = 3) +
    geom_label(data = day1_abund, aes(x= 0.68, y = c(0.19, 0.045, 0.55, 0.99, 0.9), label= pick_Phylum, fill = pick_Phylum), size = 5, hjust = 1, nudge_x = -.2, color = "black", label.size = NA) +
    scale_fill_manual(values = phyla.col) +
    labs(title = "b \n",
       #subtitle = "Top highest abundance phyla",
       x = "Day",
       y = "Mean relative abundance") +
    theme_void() +
  theme(plot.title = element_text(size = 15, face = "bold", vjust = 1, hjust = 0),
        plot.subtitle = element_text(size = 10, color = "grey50", vjust = 2, hjust = 0.5),
        axis.title.x = element_text(hjust = 0.68, vjust = -2, size = 10, color = "grey30", face = "bold"),
        axis.text.x = element_text(vjust = 1, size = 10, color = "grey30"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0,0.2,0.2,0.2),"in"),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill=NA,color= "grey80"),
        legend.position = "none",
        aspect.ratio = 0.7) +
    scale_x_discrete(expand= expand_scale(mult = c(0, 0), add = c(3,1)), limits = c("1", "7", "14", "Control")) +
    scale_y_continuous(expand= expand_scale(mult = c(0, 0), add = c(0.02,0.12)))
)

ggplot2::ggsave(filename = "Relative abundance (area).png", 
       plot = area.plot,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Relative abundance", 
       width = 5, 
       height = 4, 
       units = "in", 
       dpi = 1000, 
       limitsize = TRUE,
       bg = "white")
```
## Tukey for Alpha diversity

```{r}
adiv

#Shannon
shan.div <- adiv[,c("Shannon", "day")]
tukey.shan <- TukeyHSD(aov(lm(Shannon~day, shan.div)))$day
tukey.shan <- as.data.frame(tukey.shan)

tukey.shan[order(tukey.shan$`p adj`),]

#Simpson
simp.div <- adiv[,c("Simpson", "day")]
tukey.simp <- TukeyHSD(aov(lm(Simpson~day, simp.div)))$day
tukey.simp <- as.data.frame(tukey.simp)

tukey.simp[order(tukey.simp$`p adj`),]

#Chao1
chao1.div <- adiv[,c("Chao1", "day")]
tukey.chao1 <- TukeyHSD(aov(lm(Chao1~day, chao1.div)))$day
tukey.chao1 <- as.data.frame(tukey.chao1)

tukey.chao1[order(tukey.chao1$`p adj`),]
```
## Beta-diversity

```{r}
dist.mat <- as(gp.all.dist, "matrix")
days_all <- sample_data(all.norm)$day
sub_dist <- list()

for (group in levels(days_all)) { 
    row_group <- which(days_all == group)
    sample_group <- sample_names(all.norm)[row_group]
    sub_dist[[group]] <- dist.mat[ sample_group, sample_group]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

philrdays<- reshape2::melt(sub_dist)
df.philr <- philrdays[complete.cases(philrdays), ]
df.philr$L1 <- factor(df.philr$L1, levels=names(sub_dist))
head(df.philr)

ggplot(df.philr, aes(x=L1, y=value, fill=L1)) + 
  geom_violin(alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.2) +
  theme(legend.position="top") +
  scale_fill_manual(values = day_cols) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_text(vjust=1,size=12), 
        axis.text.y=element_text(size=12)) +
  theme_classic() +
  labs(title = "Beta diversity: phILR distance",
       subtitle = "between each samples per day",
       x = "Day",
       y = "phILR distance") +
  guides(fill = guide_legend(title = "Day", title.position = "left")) +
  stat_compare_means(comparisons = comps, 
                     #label = "p.signif",
                     tip.length = 0.02, 
                     method = "wilcox.test", 
                     hide.ns = TRUE)
```

Permanova

```{r}
#test
meta.df <- meta.df %>% filter(sample_ID %in% sample_names(all.fil))
meta.df$wfa_zscore <- all.fil@sam_data$wfa_zscore
meta.df$day <- all.fil@sam_data$day
adonis2(gp.all.dist ~ age_month + sex + wfa_zscore + day + breastfeeding, data = meta.df, by = "terms", permutations = 999)
```

# Clustering

```{r}
#Colors
CSTcolors <- scale_fill_manual(values = c("#6388B4", "#FFAE34", "#EF6F6A", "#8CC2CA", "#55AD89"))
CSTcolors.value <- c("#6388B4", "#FFAE34", "#EF6F6A", "#8CC2CA", "#55AD89")

#Clustering
hclust.all <- hclust(gp.all.dist, method = "average")

hc <- as.dendrogram(hclust.all)
cut3 <- cutree(hclust.all, k = 3)
cut4 <- cutree(hclust.all, k = 4)

plot(hc)

#visualize dendrogram
temp_col <- CSTcolors.value[as.numeric(cut4)]
temp_col <- temp_col[order.dendrogram(hc)]
temp_col <- factor(temp_col, unique(temp_col))


hc %>% color_branches(clusters = as.numeric(temp_col), col = levels(temp_col), groupLabels = T) %>%
  plot() + legend("topleft", legend = c(1:4), fill = CSTcolors.value)

sample_data(all.norm)$CST <- as.factor(cut3)

(ord_plot.cst <- phyloseq::plot_ordination(all.norm, ord_all, axes = c(1,2), type="samples", color = 'CST') + 
  geom_point(size = 1) +
  labs(subtitle = "Aitchison Distance // hclust 4 CSTs")  +
  stat_ellipse() + scale_color_manual(values = CSTcolors.value) +
  theme_pubr())

#Other clustering method
asw.y <- pamk(gp.all.dist, krange=2:10, criterion="asw", diss=TRUE,ns=10, critout = TRUE) #3 or 4 clusters
BS.y <- nselectboot(gp.all.dist, B = 50, clustermethod = pamkCBI, classification = "averagedist", krange=2:10) # 4 
ps <- prediction.strength(ord_all$vectors, Gmin = 2, Gmax = 10, M = 50, clustermethod = pamkCBI, cutoff=0.60) #3
```

Extract information from tree balance
```{r}
info_balance <- function(node) {
      bal <- name.balance(tree.all.philr, tax_table(all.fil), node, return.votes=c('up', 'down'))
      return(bal)
}
```
## Silhouette & Random forest

```{r}
sample_data(all.norm)$CST <- as.factor(cut4)
clust <- as.array(setNames(as.numeric(sample_data(all.norm)$CST), rownames(sample_data(all.norm))))
silwidths <- silhouette(clust, gp.all.dist)
rownames(silwidths) <- names(clust)
silwidths
#change cluster manually
clust["AHH20080"] <- 1
clust["AHH20099"] <- 3
clust["AHH20112"] <- 1
clust["AHH20135"] <- 2
clust["AHH20137"] <- 2
clust["AHH20157"] <- 4
clust["AHH20161"] <- 3
clust["AHH20162"] <- 3
clust["AHH20634"] <- 1
clust["AHH20639"] <- 3
clust["AHH20654"] <- 3
clust["AHH20657"] <- 4
clust["AHH20779"] <- 3
clust["AHH20796"] <- 4
clust["AHH20806"] <- 1
clust["AHH20809"] <- 3
clust["PHH2393"] <- 1
clust["PHH2419"] <- 4
clust["PHH2811"] <- 3
clust["PHH2814"] <- 3
clust["PHH2825"] <- 2
clust["PHH2860"] <- 3
clust["PHH2880"] <- 3
clust["PHH2893"] <- 1
clust["PHH2940"] <- 3
clust["PHH2944"] <- 1
clust["PHH3342"] <- 1
clust["PHH3347"] <- 1
clust["PHH3350"] <- 2
clust["PHH3355"] <- 1
clust["PHH3378"] <- 3
clust["PHH3381"] <- 3

silwidths <- silhouette(clust, gp.all.dist)
rownames(silwidths) <- names(clust)
fviz_silhouette(silwidths)

sample_data(all.norm)$CST <- as.factor(clust)
sample_data(all.fil)$CST <- as.factor(clust)

#-----# vimp
gp.all.df <- as.data.frame(gp.all)
identical(rownames(gp.all.df), names(cut4)) # TRUE
gp.all.df$CST <- as.factor(clust)
philr.rf.abund <- rfsrc(CST~., data = gp.all.df, ntree=10000, importance = "permute")
print(philr.rf.abund) #9.5%
plot.rfsrc(philr.rf.abund)
plot(gg_vimp(philr.rf.abund, nvar = 50))

sort(philr.rf.abund$importance[,1], decreasing = TRUE)[1:20] ## All: n695, n903, n1003
sort(philr.rf.abund$importance[,2], decreasing = TRUE)[1:20] ## CST1: n695, n903
sort(philr.rf.abund$importance[,3], decreasing = TRUE)[1:20] ## CST2: n903, n904
sort(philr.rf.abund$importance[,4], decreasing = TRUE)[1:20] ## CST3: n695, n1099
sort(philr.rf.abund$importance[,5], decreasing = TRUE)[1:20] ## CST4: n696

imp.nodes <- convert_to_long(gp.all, get_variable(all.norm, 'CST')) %>% filter(coord %in% c("n1","n695"))
info_balance("n695") #Actinobacteria + Proteobacteria + Firmicutes/ Bacteroidetes
info_balance("n696")
info_balance("n1") #Others / Firmicutes

(balance.plot <- ggplot(imp.nodes, aes(x=labels, y=value, fill = labels)) +
    geom_violin(color = NA) +
    geom_boxplot(width = 0.2,
                 alpha = 0.1) +
    facet_grid(.~coord, scales = "free") +
    ylab('Balance value') + 
    theme_bw() + 
    scale_fill_manual(values = c("#6388B4", "#FFAE34", "#EF6F6A", "#8CC2CA", "#55AD89"))  + geom_hline(yintercept = 0, color = "grey20") + 
      theme(legend.position = "none", 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank()) +
    labs(caption = "\nn1: Others/ Firmicutes \nn695: Actinobacteria + Proteobacteria + Firmicutes/ Bacteroidetes") +
      theme(plot.caption = element_text(hjust = 0))
)

(ord_plot.cst <- phyloseq::plot_ordination(all.fil, ord_all, axes = c(1,2), type="samples", color = 'CST') + 
  geom_point(size = 1) +
  labs(subtitle = "Aitchison Distance // hclust 4 CSTs")  +
  stat_ellipse() + 
  scale_color_manual(values = CSTcolors.value) +
  theme_pubr())

ggarrange(ord_plot.cst, balance.plot)
```

Additional function for heatmap and extract clade
```{r}
## update on heatmap plot functions
make_hcb <- function(data, var, name = NULL, fillScale = NULL, ...) {
      hcb <- ggplot(data=data, aes_string(x="index", y=1, fill=var)) + 
            geom_raster(show.legend = T) +
            scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
            scale_x_continuous(expand=c(0,0)) +
            xlab(NULL) + ylab(NULL) +
            theme(axis.title=element_blank(), axis.ticks=element_blank()) +
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_text(size=12, face="bold")) +
            theme(plot.margin=unit(c(0,0,0,0),"lines")) +  
            #axis.ticks.margin =unit(0,"null"), ...) +
            guides(fill=F)
      if(!is.null(fillScale)) hcb <- hcb + fillScale
      return(hcb)
}

mush <- function(hmap, hcbs) {
      require(gtable)
      require(grid)
  require(gridExtra)
      cbgs <- lapply(hcbs, ggplotGrob)
      hmg <- ggplotGrob(hmap)
      
      # Make sure both plots have the same dimensions in grob objects
      for (i in seq_along(cbgs)) {
            cbgs[[i]] <- gtable_add_cols(cbgs[[i]], widths=unit(1,"null"), pos=8)
            cbgs[[i]] <- gtable_add_cols(cbgs[[i]], widths=unit(1,"null"), pos=8)
      }
      ## Make sure that both plots have the same widths
      cbWidths <- lapply(cbgs, function(x) x$widths)
      maxWidth <- do.call(unit.pmax, cbWidths)
      maxWidth <- unit.pmax(hmg$widths, maxWidth)
      ## replace widths in each grob object with maxWidth
      hmg$widths <- maxWidth
      for (i in seq_along(cbgs)){
            cbgs[[i]]$widths <- maxWidth
      }
      heights <- unit.c(unit(rep(1,length(cbgs)), "lines"), unit(1, "null"))
      rval <- do.call(arrangeGrob, args = c(cbgs, list(hmg), ncol=1, heights=list(heights)))
      return(rval)
}

```

## Heat map

```{r}
#top.genus
genus.ps <- tax_glom(all.fil, "Genus")
genus.abund <- transform_sample_counts(genus.ps, function(OTU) OTU/sum(OTU))
top.genus <- top_taxa(genus.abund, 20)
top.genus.ps <- prune_taxa(top.genus, genus.abund)

taxa.order <- names(sort(taxa_sums(top.genus.ps)))
sample.order <- rownames(sample_data(top.genus.ps)[order(sample_data(top.genus.ps)$antibiotic_trt, sample_data(top.genus.ps)$day)])

#phylum_cols <- c("Actinobacteria" = "#EEDD88", "Bacteroidetes" = "#77AADD", "Firmicutes" = "#EE8866", "Proteobacteria" = "#FFAABB",  "Fusobacteria" = "grey90")

(hm <- plot_heatmap(top.genus.ps, 
                    taxa.label="Species", 
                    sample.order=sample.order,
                    taxa.order=taxa.order) +
    geom_tile(color = "grey10", size = 0.05) +
    scale_fill_gradientn(colors = c("grey10", "#457B9D", "#A8DADC", "#F1FAEE", "#E63946"), 
                         values=c(0, 0.25, 0.5, 0.75, 1),
                         na.value = NA, 
                         guide="colourbar", 
                         name="Relative\nproportion") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", size = 15, color = "grey10"),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"in")
        )
)

hcbdf <- data.frame(sample_data(top.genus.ps))[sample.order,]
hcbdf$index <- seq(1,phyloseq::nsamples(top.genus.ps))
hcb <- make_hcb(hcbdf, "day", name="Day", fillScale = scale_fill_manual(values = c("#ef476f", "#ffd166", "#26547c")))
hcb <- hcb + annotate("text", x=tapply(hcbdf$index, hcbdf[,"day",drop=T], mean), y=1, label=levels(hcbdf[,"day",drop=T]), size=4, color = "white")

#Age component
hcbage <- make_hcb(hcbdf, "age_month", name="age_month", fillScale=scale_fill_gradientn(colours=c("#fff8ab",  "#b6c670", "#6c973c", "#0d690e"), values=rescale(c(0, 16,  25, 60)), na.value="black"))

##Day component
hcbday <- make_hcb(hcbdf, "day", name="Day",fillScale = scale_fill_manual(values=day_cols))


#Antibiotic treatment
hcbabt <- make_hcb(hcbdf, "antibiotic_trt", name="antibiotic_trt",fillScale = scale_fill_manual(values= c("yes" = "#087e8b", "no" =  "#ff5a5f", "NA" = "grey75")))

big.hm <- mush(hm, list(hcbday, hcbabt))
big.heatmap <- plot(ggarrange(big.hm))

ggplot2::ggsave(filename = "Heatmap (unorm).png", 
       plot = big.heatmap,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Heatmap", 
       width = 15, 
       height = 7, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")
```

===

# Filter data of 04EN study for the paired-day 's differential abundance analysis

```{r}
data.04en <- subset_samples(all.norm, study_ID == "04EN") #218 samples
data.09av <- subset_samples(all.norm, study_ID == "09AV") #192 samples

#Patients participated in day 1 and day 7
patient.17.ps <- prune_samples(sample_data(data.04en)$day %in% c("1", "7"), data.04en) #160 samples
patientd7.d <- sample_data(ps_filter(patient.17.ps, day == "7" & diarrhea_stool == "yes", .target = "sample_data", .keep_all_taxa = FALSE))$patient_ID #3 patients still had diarrhoea in day 7
patient.17.ps <- prune_samples(!sample_data(patient.17.ps)$patient_ID %in% patientd7.d, patient.17.ps)
patient.17.ps <- prune_taxa(taxa_sums(patient.17.ps) > 0, patient.17.ps)
sample_data(patient.17.ps) <- sample_data(patient.17.ps)[sample_data(patient.17.ps)$patient_ID %in% sample_data(patient.17.ps)$patient_ID[duplicated(sample_data(patient.17.ps)$patient_ID)],] #remove unique patients that only have day 1 or day 7
patient.17.ps #134 samples = 67 patients // 946 taxa


#Patients participated in day 7 and day 14
patient.714.ps <- prune_samples(sample_data(data.04en)$day %in% c("14", "7"), data.04en) #128 samples
patient.714.ps <- prune_taxa(taxa_sums(patient.714.ps) > 0, patient.714.ps)
sample_data(patient.714.ps) <- sample_data(patient.714.ps)[sample_data(patient.714.ps)$patient_ID %in% sample_data(patient.714.ps)$patient_ID[duplicated(sample_data(patient.714.ps)$patient_ID)],] #remove unique patients that only have day 7 or day 14
patient.714.ps #76 samples = 38 patients // 980 taxa

#Patients participated in day 1 and day 14
patient.114.ps <- prune_samples(sample_data(data.04en)$day %in% c("14", "1"), data.04en) #148 samples
patient.114.ps <- prune_samples(sample_data(data.04en)$day != 7, data.04en) #148 samples

patient.114.ps <- prune_taxa(taxa_sums(patient.114.ps) > 0, patient.114.ps)
sample_data(patient.114.ps) <- sample_data(patient.114.ps)[sample_data(patient.114.ps)$patient_ID %in% sample_data(patient.114.ps)$patient_ID[duplicated(sample_data(patient.114.ps)$patient_ID)],] #remove unique patients that only have day 1 or day 14
patient.114.ps # 116 samples and 973 taxa

#Patients participated full-time or d17
patient.17 <- unique(sample_data(patient.17.ps)$patient_ID)
patient.714 <- unique(sample_data(patient.714.ps)$patient_ID)
patient.3d <- unique(c(patient.17, patient.714))
patient.3d.ps <- subset_samples(data.04en, patient_ID %in% patient.3d)
patient.3d.ps

#Unormalized data
data.04en.un <- subset_samples(all.fil, study_ID == "04EN")
#Patients participated in day 1 and day 7
patient.17.un <- prune_samples(sample_data(data.04en.un)$day %in% c("1", "7"), data.04en.un) #160 samples
patientd7.d.un <- sample_data(ps_filter(patient.17.un, day == "7" & diarrhea_stool == "yes", .target = "sample_data", .keep_all_taxa = FALSE))$patient_ID #3 patients still had diarrhoea in day 7
patient.17.un <- prune_samples(!sample_data(patient.17.un)$patient_ID %in% patientd7.d.un, patient.17.un)
patient.17.un <- prune_taxa(taxa_sums(patient.17.un) > 0, patient.17.un)
sample_data(patient.17.un) <- sample_data(patient.17.un)[sample_data(patient.17.un)$patient_ID %in% sample_data(patient.17.un)$patient_ID[duplicated(sample_data(patient.17.un)$patient_ID)],] #remove unique patients that only have day 1 or day 7
patient.17.un #134 samples = 67 patients // 928 taxa

#Patients participated in day 7 and day 14
patient.714.un <- prune_samples(sample_data(data.04en.un)$day %in% c("14", "7"), data.04en.un) 
patient.714.un <- prune_taxa(taxa_sums(patient.714.un) > 0, patient.714.un)
sample_data(patient.714.un) <- sample_data(patient.714.un)[sample_data(patient.714.un)$patient_ID %in% sample_data(patient.714.un)$patient_ID[duplicated(sample_data(patient.714.un)$patient_ID)],] #remove unique patients that only have day 7 or day 14
patient.714.un #76 samples = 38 patients // 963 taxa

#Patients participated in day 1 and day 14
patient.114.un <- prune_samples(sample_data(data.04en.un)$day %in% c("14", "1"), data.04en.un) #116 samples
patient.114.un <- prune_taxa(taxa_sums(patient.114.un) > 0, patient.114.un)
sample_data(patient.114.un) <- sample_data(patient.114.un)[sample_data(patient.114.un)$patient_ID %in% sample_data(patient.114.un)$patient_ID[duplicated(sample_data(patient.114.un)$patient_ID)],] #remove unique patients that only have day 1 or day 14
patient.114.un # 116 samples and 955 taxa
```

# Genus relative abundance plot

```{r}
dt04_abund <- microbiome::aggregate_rare(data.04en, level = "Genus", detection = 0, prevalence = 0)

dt04_abund <- microbiomeutilities::phy_to_ldf(dt04_abund, 
                                         transform.counts = "compositional")

dt04_abund$antibiotic_trt <- as.factor(dt04_abund$antibiotic_trt)



dt04_abund <- dt04_abund %>% 
  group_by(Genus, day, antibiotic_trt) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()

top.gn.abund <- subset(dt04_abund, Genus %in% c("Bifidobacterium", 
                                                "Streptococcus", 
                                                "Bacteroides", 
                                                "Escherichia", 
                                                "Veillonella",
                                                "Phocaeicola"))

top.gn.abund$Genus <- factor(top.gn.abund$Genus,
                             levels = c("Bifidobacterium", 
                                        "Streptococcus", 
                                        "Bacteroides",
                                        "Escherichia",
                                        "Veillonella",
                                        "Phocaeicola"))

(genus.line <- ggplot(top.gn.abund,
                      aes(x = day,
                          y = mean_abundance,
                          group = antibiotic_trt,
                          color = antibiotic_trt)) + 
    geom_errorbar(aes(ymin=mean_abundance-se, 
                      ymax=mean_abundance+se), 
                  width=0.1, 
                  alpha = 0.75) +
    geom_line(size = 1.6) +
    geom_point(size = 1) +
    facet_wrap(.~Genus, 
               nrow = 1) +
    scale_color_manual(values = c("no" = "grey50", "yes" = "#e07a5f"),
                       guide = guide_legend(title = "Antibiotic \ntreatment")) +
    theme_light() +
    guides(x.sec = "axis",
         y.sec = "axis") +
    labs(x = "Day",
         y = "Mean relative abundance") +
    theme(plot.title = element_text(size = 15, 
                                    face = "bold", 
                                    vjust = 3, 
                                    hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.line = element_line(linewidth = .5),
          axis.text.x.bottom = element_text(vjust = -1, 
                                            size = 10),
          axis.text.y.left =element_text(size = 10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.y = element_text(color = "grey30", 
                                      face = "bold", 
                                      size = 10, 
                                      angle = 90, 
                                      vjust = 3),
          strip.background = element_rect(fill = NA, 
                                          color = "grey30", 
                                          linewidth = 1.5),
          strip.text = element_text(face = "bold.italic", color = "black"))
  )


ggplot2::ggsave(filename = "genus.line.pdf", 
       plot = genus.line,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Heatmap", 
       width = 16, 
       height = 4, 
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")


top.gn.abund <- top.gn.abund %>% 
  group_by(Genus, day, antibiotic_trt) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()

filter(dt04_abund, Genus == "Phocaeicola")
```

# Demographic 04EN:

```{r}
sample.04en <- as.data.frame(unclass(sample_data(data.04en)))
names(sample.04en)

sample.04en <- sample.04en %>% pivot_wider(names_from = day, values_from = diarrhea_stool, values_fill = "NA") %>% as.data.frame()

explanatory <- c("sex", "age_month", "Infection_type", "breastfeeding", "antibiotic_trt", "probiotic_trt",  "hospitalization_days", "timetostopdiarr", "wfa_zscore")

demo.04en <- sample.04en %>% summary_factorlist("patient_ID", 
                                                total_col = TRUE, 
                                                cont = "median",
                                                explanatory, 
                                                p=TRUE, 
                                                digits = c(2,2,2,2,0),
                                                orderbytotal = TRUE,
                                                na_include=TRUE)

demo.04en$levels <- c("Male", "Female",
                      "Median (IQR)",
                      "Virus only", "Virus + Bacteria", "Virus + Parasite", "Unknown", "Bacteria only",
                      "No", "Yes",
                      "No", "Yes",
                      "Yes", "No",
                      "Median (IQR)",
                      "Median (IQR)",
                      "Median (IQR)")

demo.04en <- demo.04en[,c("levels", "Total")]
colnames(demo.04en) <- c(" ", "Total (% or IQR)")

demo_table <-kable(demo.04en, escape = F, 
                   align = "l", 
                   caption = "<b> Table 1: Demographic table of the longitudinal study") %>%
  kable_classic(bootstrap_options = c("hover"), 
                full_width = F, 
                position = "left", html_font = "helvetica") %>%
  group_rows(index = c("Sex" = 2, 
                       "Age months" = 1, 
                       "Infection types" = 5, 
                       "Breast feeding" = 2, 
                       "Antibiotic treatment" = 2, 
                       "Probiotic treatment" = 2, 
                       "Hospitalization duration (days)" = 1, 
                       "Diarrhoea duration (hours)" = 1, 
                       "wfa z-score" = 1),
             color = "black",
             background = "grey75") 

tf <- tempfile(fileext = ".docx")

demo_table <- demo_table %>%
  column_spec(1:2, color = "black", width = "2in")

save_kable(x = demo_table, file = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/demo_table.docx", density = 10000, zoom = 1)
```

## Beta diversity for patient have at least 2 timepoints

```{r}
sample.3d <- rownames(sample_data(patient.3d.ps)) 
dist.mat.3d <- dist.mat[sample.3d, sample.3d]
days_04en <- sample_data(patient.3d.ps)$day
sub_dist.3d <- list()

for (group in levels(days_04en)) { 
    row_group <- which(days_04en == group)
    sample_group <- sample_names(patient.3d.ps)[row_group]
    sub_dist.3d[[group]] <- dist.mat.3d[sample_group, sample_group]
    sub_dist.3d[[group]][!lower.tri(sub_dist.3d[[group]])] <- NA
}

philrdays.3d<- reshape2::melt(sub_dist.3d)
df.philr.3d <- philrdays.3d[complete.cases(philrdays.3d),]
df.philr.3d$L1 <- factor(df.philr.3d$L1, levels=names(sub_dist.3d))

head(df.philr.3d)

colnames(df.philr.3d) <- c("Var1", "Var2", "Bdiff", "day")

df.philr.3d.mean <- df.philr.3d %>%
  group_by(day) %>%
  summarise_at(vars(Bdiff), list(mean_bdiff = mean))

wilcox.test(Bdiff ~ day, data = subset(df.philr.3d, day != 14)) #p-value = 1.605e-08
wilcox.test(Bdiff ~ day, data = subset(df.philr.3d, day != 1)) #p-value = 0.1224
wilcox.test(Bdiff ~ day, data = subset(df.philr.3d, day != 7)) #p-value = 0.01671


(bdiv.btw <- ggplot(df.philr.3d, aes(x=day, y=Bdiff)) + 
  geom_violin(linetype = 1) +
  geom_boxplot(aes(fill = day), width = 0.2, outlier.colour = NA) +
  scale_fill_manual(values = c("#ef476f", "#ffd166", "#26547c")) +
  theme_classic() +
  labs(title = "b",
       x = "Day",
       y = "Aitchison distance (phILR transformed)") +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), 
                     comparisons = comps[c(1,2,4)],
                     tip.length = 0.02,
                     method = "wilcox.test",
                     hide.ns = TRUE,
                     size = 3) +
    ylim(0,50) +
    guides(x.sec = "axis",
         y.sec = "axis") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line = element_line(size = .5, color="grey30"),
          axis.text.x.bottom = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y.left =element_text(size=10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -1),
          axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3)))

ggplot2::ggsave(filename = "bdiv_btw.sample.png", 
       plot = bdiv.btw,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Beta", 
       width = 4, 
       height = 4, 
       units = "in",
       dpi = "retina",
       limitsize = TRUE,
       bg = "white")
```
## Bdiv day: among patiens group by antibiotic_trt

```{r}
#df.philr.d1 <- filter(df.philr.3d, day == 1)
#df.d1 <- filter(meta.df, day == 1)

aby.samp <- meta.df[which(sample_data(patient.3d.ps)$antibiotic_trt == "yes"), "sample_ID"] #65 samples
abn.samp <- meta.df[which(sample_data(patient.3d.ps)$antibiotic_trt == "no"), "sample_ID"] #107 samples

df.philr.3d <- df.philr.3d %>%
  mutate(antibiotic = if_else(df.philr.3d$Var1 %in% aby.samp & df.philr.3d$Var2 %in% aby.samp,
                              "btw_yes", 
                              if_else(df.philr.3d$Var1 %in% abn.samp & df.philr.3d$Var2 %in% abn.samp,
                                      "btw_no",
                                      "yes_and_no")
                              )
         )

#plot boxplot
ab.comp <- make_pairs(df.philr.3d$antibiotic)
(bdiv.ab <- ggplot(data = subset(df.philr.3d, antibiotic != "yes_and_no" & day != "1"),
                   aes(day, Bdiff, fill = antibiotic)) +
    geom_boxplot(data = subset(df.philr.3d, day == "1"),
                 mapping = aes(day, Bdiff), fill ="#e07a5f", width = 0.25) +
    geom_boxplot(outlier.color = NA, width = 0.5) +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), 
                       method = "wilcox.test",
                       label.y = 55, 
                       hide.ns = T,
                       size = 5) +
    scale_fill_manual(values = c("#f4f1de", "#e07a5f"), 
                      labels = c("Between non-usages", "Between usages"),
                      guide = guide_legend(title = "Antibiotic usage pair")) +
    theme_classic() +
    labs(title = "c", 
         x = "Day",
         y = "b-diversity") +
    scale_y_continuous(limits = c(5,59),
                       expand = c(0,0)) +
    scale_x_discrete(limits = c("1", "7", "14")) +
    guides(x.sec = "axis",
         y.sec = "axis") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line = element_line(size = .5, color="grey30"),
          axis.text.x.bottom = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y.left =element_text(size=10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -1),
          axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3),
          legend.position = "none"))

ggplot2::ggsave(filename = "Diff beta & abt group_updated.png", 
       plot = bdiv.ab,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Beta/test", 
       width = 5, 
       height = 4, 
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")

ab.pair.mean.d7 <- subset(df.philr.3d, antibiotic != "yes_and_no" & day == 7) %>%
  group_by(antibiotic) %>%
  summarise_at(vars(Bdiff), list(mean_bdiff = mean))

ab.pair.mean.d14 <- subset(df.philr.3d, antibiotic != "yes_and_no" & day == 14) %>%
  group_by(antibiotic) %>%
  summarise_at(vars(Bdiff), list(mean_bdiff = mean))
```
## Beta diversity per patients through 3 timepoints

```{r}
## within patient pairwise diversity 
## day 7 - 1
gp.all.dist.mat <- as.matrix(gp.all.dist)

patient.17
dist17 <- matrix(rep(NA,length(patient.17)), nrow=length(patient.17), ncol=1)
patient.17.ps
p17.df <- as.data.frame(as.matrix(sample_data(patient.17.ps)))

test1 <-combn(rownames(p17.df[which(p17.df$patient_ID == patient.17[1]),]),2)
gp.all.dist.mat[test1[1], test1[2]]


for(i in 1:length(patient.17)) {
   coord <- combn(rownames(p17.df[which(p17.df$patient_ID == patient.17[i]),]),2)
   dist17[i,1] <- gp.all.dist.mat[coord[1], coord[2]]
   print(paste('Successful for patient ID ', i, collapse=""))
}

dist17
summary(dist17)
colnames(dist17) <- "distance17"
rownames(dist17) <- unique(p17.df$patient_ID)

## within patient pairwise diversity 
## day 14 - 7
patient.714
dist714 <- matrix(rep(NA,length(patient.714)), nrow=length(patient.714), ncol=1)
patient.714.ps
p714.df <- as.data.frame(as.matrix(sample_data(patient.714.ps)))

test2 <-combn(rownames(p714.df[which(p714.df$patient_ID == patient.714[1]),]),2)
gp.all.dist.mat[test2[1], test2[2]]

for(i in 1:length(patient.714)) {
   coord <- combn(rownames(p714.df[which(p17.df$patient_ID == patient.714[i]),]),2)
   dist714[i,1] <- gp.all.dist.mat[coord[1], coord[2]]
   print(paste('Successful for patient ID ', i, collapse=""))
}

dist714
summary(dist714)
colnames(dist714) <- "distance714"
rownames(dist714) <- unique(p714.df$patient_ID)

#combine dist for patients have 3 timepoints (problem here !!!)
dist.3d <- as.matrix(dist17[rownames(dist714),])
colnames(dist.3d) <- "distance17"

#create another subset that include all samples
dist.3d <- dist17
dist.3d <- as.data.frame(dist.3d)
dist.3d$distance714 <- dist714[match(rownames(dist.3d), rownames(dist714))]
dist.3d.p <- dist.3d
dist.3d.p$patient_ID <- rownames(dist.3d.p)
dist.3d.p <- dist.3d.p %>% 
  pivot_longer(cols = c("distance17", "distance714"), names_to = "distance_day", values_to = "distance")
dist.3d.p$antibiotic_trt <- p17.df$antibiotic_trt
dist.3d.p <- as.data.frame(dist.3d.p)
dist.3d.p <- dist.3d.p[which(rowSums(is.na(dist.3d.p))==0),]

dist.3d.plot <- ggplot(dist.3d.p, aes(x = distance_day, y = distance)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", 
                     paired = T, 
                     label.y = 35) +
  theme_classic() +
  labs(title = "d", 
       x = "Day pairs",
       y = "b-diversity") +
  scale_y_continuous(limits = c(7,37)) +
  scale_x_discrete(labels = c("Day 7 - Day 1", "Day 14 - Day 7")) +
  theme(plot.title = element_text(face = "bold", size = 12, margin = ggplot2::margin(0,0,30,0)),
        axis.title = element_text(face = "bold", size = 10, color = "grey40"),
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_line(color = "grey40"),
        aspect.ratio = 1,
        panel.grid.major.y = element_line(color = "grey85"))

(bdiv.patient <- ggplot(dist.3d.p, aes(x = distance_day, y = distance, fill = antibiotic_trt)) +
  geom_boxplot(width = 0.5) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                     method = "wilcox.test",
                     label.y = 45,
                     paired = F, 
                     size =5) +
  scale_fill_manual(values = c("#f4f1de", "#e07a5f"), guide = guide_legend(title = "Antibiotic treatment \n")) +
  theme_classic() +
  labs(title = "d", 
       x = "Day pairs",
       y = "b-diversity") +
  scale_y_continuous(limits = c(7,47)) +
  scale_x_discrete(labels = c("Day 7 - Day 1", "Day 14 - Day 7")) +
  guides(x.sec = "axis",
         y.sec = "axis") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line = element_line(size = .5, color="grey30"),
          axis.text.x.bottom = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y.left =element_text(size=10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -1),
          axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3)))

ggplot2::ggsave(filename = "Bdiv between ab.png", 
       plot = bdiv.patient,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Beta/test", 
       width = 5, 
       height = 4, 
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")

#add timetostopdiar into dist3d
dist.3d.df <- meta.df[which(meta.df$patient_ID %in% rownames(dist.3d)),]
dist.3d$timetostopdiarr <- dist.3d.df[which(dist.3d.df$day == "1"),"timetostopdiarr"]
dist.3d$hos_days <- dist.3d.df[which(dist.3d.df$day == "1"),"hospitalization_days"]
dist.3d$antibiotic_trt <- dist.3d.df[which(dist.3d.df$day == "1"),"antibiotic_trt"]


#summary

bdiv.btw
bdiv.ab
bdiv.patient

legend_1 <- get_legend(ord_plot)
legend_2 <- get_legend(bdiv.ab)
legend_3 <- get_legend(bdiv.patient)
legend.arrange <- ggarrange(legend_1, legend_2, legend_3, nrow = 3)

rm_legend <- function(p){p + theme(legend.position = "none")}

(bdiv.sum <- ggarrange(rm_legend(ord_plot),
                       ggarrange(rm_legend(bdiv.ab), 
                       rm_legend(bdiv.patient), nrow = 1),
                       #align = "h",
                       nrow = 2))

ggplot2::ggsave(filename = "summary_beta.pdf", 
       plot = bdiv.sum,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Beta/", 
       width = 8, 
       height = 10, 
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")
```

# Beta diversity test among VS within, and among-d14 VS among-control

```{r}
#Among vs Within samples

df.philr.3d$Bdiff #among samples
dist.3d.p$distance #within samples

summary(df.philr.3d$Bdiff)
summary(dist.3d.p$distance)

wilcox.test(df.philr.3d$Bdiff, dist.3d.p$distance) #within sample's bdiv is sign. LOWER than among sample's

#pair: day 14 and control
##filter extract control 
control.samples <- rownames(sample_data(ps_filter((data.09av), day == "Control"))) 
d14.samples <- rownames(sample_data(ps_filter((data.04en), day == "14")))
d1.samples <- rownames(sample_data(ps_filter((data.04en), day == "1")))


d14c.dist <- gp.all.dist.mat[control.samples,d14.samples]
d14c.dist <- reshape2::melt(d14c.dist)

d1c.dist <- gp.all.dist.mat[control.samples,d1.samples]
d1c.dist <- reshape2::melt(d1c.dist)


summary(dist17)
summary(dist714)
summary(d14c.dist$value)
summary(d1c.dist$value)

##Wilcoxon test
wilcox.test(dist17, dist714) #no diff
wilcox.test(dist17, d14c.dist$value) #significant different
t.test(dist714, d14c.dist$value) #significant different

wilcox.test(d14c.dist$value, d1c.dist$value) 
```

# Linear Mixed-effects in beta diversity

```{r}
bdiv.04en <- data.frame(sample_data(data.04en))
bdiv.04en <- bdiv.04en[,c("patient_ID", "sex", "wfa_zscore", "age_month", "Infection_type", "antibiotic_trt")]

bdiv.04en <- bdiv.04en[which(bdiv.04en$patient_ID %in% dist.3d.p$patient_ID),]
rownames(bdiv.04en) <- NULL
bdiv.04en <- unique(bdiv.04en)

bdiv.04en <- cbind(dist.3d, bdiv.04en)


bdiv.04en.pivot <- bdiv.04en %>% 
  pivot_longer(cols = c("distance17", "distance714"), names_to = "distance_day", values_to = "distance") 

bdiv.04en.pivot <- bdiv.04en.pivot[which(rowSums(is.na(bdiv.04en.pivot))==0),] %>% as.data.frame()
bdiv.04en.pivot$Infection_type <- as.factor(bdiv.04en.pivot$Infection_type)

levels(bdiv.04en.pivot$Infection_type) <- c("bacteria_only", "mixed", "mixed", "unknown", "virus_only")

#Linear Mixed-Effects Models

  
##Bdiv
lmm_bdiv.04en <- lmerTest::lmer(bdiv.04en.pivot$distance ~ sex + wfa_zscore + age_month + Infection_type + distance_day*antibiotic_trt + (1|patient_ID), bdiv.04en.pivot)

summary(lmm_bdiv.04en)

lmm_bdiv.04en.df <- as.data.frame(summary(lmm_bdiv.04en)$coefficients) %>% rownames_to_column("term")

##Export table
bdiv.table <- as_flextable(lmm_bdiv.04en)

save_as_docx(bdiv.table, path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/bdiv_table.docx")

```

# Tukey Post hoc test

```{r}
df.philr.3d
dist.3d
dist.3d.p

#among Day 1
a.d1 <- df.philr.3d[df.philr.3d$L1 == "1",c("value", "antibiotic")]
a.d1$antibiotic <- "no"

names(a.d1)


#among using antibiotic
a.aby <- df.philr.3d[df.philr.3d$antibiotic == "btw_yes",]
a.aby$antibiotic[a.aby$antibiotic == "btw_yes"] <- "yes"

a.aby.7 <- a.aby[a.aby$L1 == "7",c("value", "antibiotic")]
a.aby.14 <- a.aby[a.aby$L1 == "14",c("value", "antibiotic")]

#among non-using antibiotic
a.abn <- df.philr.3d[df.philr.3d$antibiotic == "btw_no",]
a.abn$antibiotic[a.abn$antibiotic == "btw_no"] <- "no"

a.abn.7 <- a.abn[a.abn$L1 == "7",c("value", "antibiotic")]
a.abn.14 <- a.abn[a.abn$L1 == "14",c("value", "antibiotic")]

#within using atibiotic
w.aby <- dist.3d.p[dist.3d.p$antibiotic_trt == "yes",]
w.aby.1 <- w.aby[w.aby$distance_day == "distance17", c("distance", "antibiotic_trt")]
w.aby.2 <- w.aby[w.aby$distance_day == "distance714", c("distance", "antibiotic_trt")]

names(w.aby.1) <- c("value", "antibiotic")
names(w.aby.2) <- c("value", "antibiotic")

#within non-using atibiotic
w.abn <- dist.3d.p[dist.3d.p$antibiotic_trt == "no",]
w.abn.1 <- w.abn[w.abn$distance_day == "distance17", c("distance", "antibiotic_trt")]
w.abn.2 <- w.abn[w.abn$distance_day == "distance714", c("distance", "antibiotic_trt")]

names(w.abn.1) <- c("value", "antibiotic")
names(w.abn.2) <- c("value", "antibiotic")

#combine files

bdiv.com <- list(a.d1, a.aby.7, a.abn.7, a.aby.14, a.abn.14, w.aby.1, w.abn.1, w.aby.2, w.abn.2)

bdiv.com.df <- melt(bdiv.com)

bdiv.com.df$L1 <- as.factor(bdiv.com.df$L1)
bdiv.com.df$L1 <- factor(bdiv.com.df$L1, labels = c("a.d1", "a.aby.7", "a.abn.7", "a.aby.14", "a.abn.14", "w.aby.1", "w.abn.1", "w.aby.2", "w.abn.2"))

names(bdiv.com.df)[4] <- "type"
bdiv.com.df <- bdiv.com.df[,c(1,3,4)]

head(bdiv.com.df)

tukey.df <- TukeyHSD(aov(lm(value~type, bdiv.com.df)))$type
tukey.df <- as.data.frame(tukey.df)
tukey.df[order(tukey.df$`p adj`),]
```


# CST transition

```{r}
meta.df$CST <- as.factor(clust)

dt04.df <- meta.df %>% filter(study_ID == "04EN")
dt04.17.df <- dt04.df[dt04.df$patient_ID%in% patient.17 & dt04.df$day != "14",]
dt04.714.df <- dt04.df[dt04.df$patient_ID %in% patient.714 & dt04.df$day != "1",]

CST17 <- dt04.17.df[,c(2, 6, 35)]
CST714 <- dt04.714.df[,c(2, 6, 35)]

#Day 1 - 7
matrix_chain_17 <- matrix(rep(0),nrow = 4,ncol = 4)
colnames(matrix_chain_17) <- paste0("CST",1:4)
rownames(matrix_chain_17) <- paste0("CST",1:4)

for(i in unique(CST17$patient_ID))
{
  sub_id_1 <- subset(CST17,patient_ID==i & day == 1)
  sub_id_7 <- subset(CST17,patient_ID==i & day == 7)
  matrix_chain_17[sub_id_7$CST[1],sub_id_1$CST[1]] <- matrix_chain_17[sub_id_7$CST[1],sub_id_1$CST[1]] + 1
}

heat_17 <- reshape2::melt(matrix_chain_17)

(cst17.plot <- ggplot(heat_17, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = value), size =5) +
  scale_fill_gradientn(colors = c("#636363", "#BDBDBD", "#F0F0F0")) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(rownames(matrix_chain_17))) +
  labs(x = "Day 1", y = "Day 7", title = "All") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") 
)


#Day 7 - 14
matrix_chain_714 <- matrix(rep(0),nrow = 4,ncol = 4)
colnames(matrix_chain_714) <- paste0("CST",1:4)
rownames(matrix_chain_714) <- paste0("CST",1:4)

for(i in unique(CST714$patient_ID))
{
  sub_id_7 <- subset(CST714,patient_ID==i & day ==7)
  sub_id_14 <- subset(CST714,patient_ID==i & day == 14)
  matrix_chain_714[sub_id_14$CST[1],sub_id_7$CST[1]] <- matrix_chain_714[sub_id_14$CST[1],sub_id_7$CST[1]] + 1
}


heat_714 <- reshape2::melt(matrix_chain_714)

(cst714.plot <- ggplot(heat_714, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = value), size =5) +
  scale_fill_gradientn(colors = c("#636363", "#BDBDBD", "#F0F0F0")) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(rownames(matrix_chain_714))) +
  labs(x = "Day 7", y = "Day 14", title = "All") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")
)

# Test 2 contigency tables
matrix_chain_17
matrix_chain_714

tab.cst <- array(c(matrix_chain_17, matrix_chain_714), dim = c(4,4,2))
tab.cst <- as.table(tab.cst)

names(dimnames(tab.cst)) = c("after_CST", "before_CST", "days")
dimnames(tab.cst)[[1]] = c("1", "2", "3", "4")
dimnames(tab.cst)[[2]] = c("1", "2", "3", "4")
dimnames(tab.cst)[[3]] = c("17", "714") 
tab.cst

#null model that assumes all cells have the same expected value
summary(tab.cst) #p = 0.003466

#The null is rejected. Next, we can fit a saturated model
m.sat <- loglm(~days*before_CST*after_CST, tab.cst)
m.sat

m1 <- loglm(~days + before_CST*after_CST, tab.cst)
sum(tab.cst[,,1])  # 67
sum(tab.cst[,,2])  # 38
m2 <- loglm(~before_CST*after_CST, tab.cst)
anova(m1, m2)
```

## CST transition and antibiotic_trt

```{r}

#Day 1 - 7
CST17.ab <- dt04.17.df[,c(2, 6, 23, 35)]
CST17.aby <- subset(CST17.ab, antibiotic_trt == "yes")
CST17.abn <- subset(CST17.ab, antibiotic_trt == "no")

##Antibiotic_trt = yes
matrix_chain_17.aby <- matrix(rep(0),nrow = 4,ncol = 4)
colnames(matrix_chain_17.aby) <- paste0("CST",1:4)
rownames(matrix_chain_17.aby) <- paste0("CST",1:4)

for(i in unique(CST17.aby$patient_ID))
{
  sub_id_1 <- subset(CST17.aby,patient_ID==i & day == 1)
  sub_id_7 <- subset(CST17.aby,patient_ID==i & day == 7)
  matrix_chain_17.aby[sub_id_7$CST[1],sub_id_1$CST[1]] <- matrix_chain_17.aby[sub_id_7$CST[1],sub_id_1$CST[1]] + 1
}

heat_17.aby <- reshape2::melt(matrix_chain_17.aby)

(cst17.aby.plot <- ggplot(heat_17.aby, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = value), size =5) +
  scale_fill_gradientn(colors = c("#636363", "#BDBDBD", "#F0F0F0")) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(rownames(matrix_chain_17))) +
  labs(x = "Day 1", y = "Day 7", title = "antibiotic_trt = yes (n = 26)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")
)

##Antibiotic_trt = no
matrix_chain_17.abn <- matrix(rep(0),nrow = 4,ncol= 4)
colnames(matrix_chain_17.abn) <- paste0("CST",1:4)
rownames(matrix_chain_17.abn) <- paste0("CST",1:4)

for(i in unique(CST17.abn$patient_ID))
{
  sub_id_1 <- subset(CST17.abn,patient_ID==i & day == 1)
  sub_id_7 <- subset(CST17.abn,patient_ID==i & day == 7)
  matrix_chain_17.abn[sub_id_7$CST[1],sub_id_1$CST[1]] <- matrix_chain_17.abn[sub_id_7$CST[1],sub_id_1$CST[1]] + 1
}

heat_17.abn <- reshape2::melt(matrix_chain_17.abn)

(cst17.abn.plot <- ggplot(heat_17.abn, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = value), size =5) +
  scale_fill_gradientn(colors = c("#636363", "#BDBDBD", "#F0F0F0")) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(rownames(matrix_chain_17))) +
  labs(x = "Day 1", y = "Day 7", title = "antibiotic_trt = no (n = 41)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")
)

##Test antibiotic_trt on d17

tab.cst17 <- array(c(matrix_chain_17.aby, matrix_chain_17.abn), dim = c(4,4,2))
tab.cst17 <- as.table(tab.cst17)

names(dimnames(tab.cst17)) = c("after_CST", "before_CST", "antibiotic")
dimnames(tab.cst17)[[1]] = c("1", "2", "3", "4")
dimnames(tab.cst17)[[2]] = c("1", "2", "3", "4")
dimnames(tab.cst17)[[3]] = c("yes", "no") 
tab.cst17

#null model that assumes all cells have the same expected value
summary(tab.cst17) #p = 0.2385

#The null is rejected. Next, we can fit a saturated model
m.sat.17 <- loglm(~ antibiotic * before_CST*after_CST, tab.cst17)
m.sat.17

m17.1 <- loglm(~before_CST + after_CST, tab.cst17)
sum(tab.cst17[,,1])  # 26
sum(tab.cst17[,,2])  # 41
m17.2 <- loglm(~before_CST*after_CST, tab.cst17)
anova(m17.1, m17.2)

#Day 714
CST714.ab <- dt04.714.df[,c(2, 6, 23, 35)]
CST714.aby <- subset(CST714.ab, antibiotic_trt == "yes")
CST714.abn <- subset(CST714.ab, antibiotic_trt == "no")

##Antibiotic_trt = yes
matrix_chain_714.aby <- matrix(rep(0),nrow = 4,ncol = 4)
colnames(matrix_chain_714.aby) <- paste0("CST",1:4)
rownames(matrix_chain_714.aby) <- paste0("CST",1:4)

for(i in unique(CST714.aby$patient_ID))
{
  sub_id_7 <- subset(CST714.aby,patient_ID==i & day == 7)
  sub_id_14 <- subset(CST714.aby,patient_ID==i & day == 14)
  matrix_chain_714.aby[sub_id_14$CST[1],sub_id_7$CST[1]] <- matrix_chain_714.aby[sub_id_14$CST[1],sub_id_7$CST[1]] + 1
}

heat_714.aby <- reshape2::melt(matrix_chain_714.aby)

(cst714.aby.plot <- ggplot(heat_714.aby, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = value), size =5) +
  scale_fill_gradientn(colors = c("#636363", "#BDBDBD", "#F0F0F0")) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(rownames(matrix_chain_714.aby))) +
  labs(x = "Day 7", y = "Day 14", title = "antibiotic_trt = yes (n = 13)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")
)

##Antibiotic_trt = no
matrix_chain_714.abn <- matrix(rep(0),nrow = 4,ncol = 4)
colnames(matrix_chain_714.abn) <- paste0("CST",1:4)
rownames(matrix_chain_714.abn) <- paste0("CST",1:4)

for(i in unique(CST714.abn$patient_ID))
{
  sub_id_7 <- subset(CST714.abn,patient_ID==i & day == 7)
  sub_id_14 <- subset(CST714.abn,patient_ID==i & day == 14)
  matrix_chain_714.abn[sub_id_14$CST[1],sub_id_7$CST[1]] <- matrix_chain_714.abn[sub_id_14$CST[1],sub_id_7$CST[1]] + 1
}

heat_714.abn <- reshape2::melt(matrix_chain_714.abn)

(cst714.abn.plot <- ggplot(heat_714.abn, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = value), size =5) +
  scale_fill_gradientn(colors = c("#636363", "#BDBDBD", "#F0F0F0")) +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev(rownames(matrix_chain_714.abn))) +
  labs(x = "Day 7", y = "Day 14", title = "antibiotic_trt = no (n = 25)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")
)

##Test antibiotic_trt on d17

tab.cst714<- array(c(matrix_chain_714.aby, matrix_chain_714.abn), dim = c(4,4,2))
tab.cst714 <- as.table(tab.cst714)

names(dimnames(tab.cst714)) = c("after_CST", "before_CST", "antibiotic")
dimnames(tab.cst17)[[1]] = c("1", "2", "3", "4")
dimnames(tab.cst17)[[2]] = c("1", "2", "3", "4")
dimnames(tab.cst17)[[3]] = c("yes", "no") 
tab.cst714

#null model that assumes all cells have the same expected value
summary(tab.cst714) #p = 0.01139

#The null is rejected. Next, we can fit a saturated model
m.sat.714 <- loglm(~antibiotic*before_CST*antibiotic*after_CST, tab.cst714)
m.sat.714

m714.1 <- loglm(~antibiotic*before_CST + antibiotic*after_CST, tab.cst714)
sum(tab.cst714[,,1])  # 13
sum(tab.cst714[,,2])  # 25
m714.2 <- loglm(~before_CST + after_CST, tab.cst714)
anova(m714.1, m714.2)

#final test for antibiotic_trt and days
anova(m17.1, m17.2)
anova(m714.1, m714.2)
##summary plot
(cst.trans <- ggarrange(cst17.plot, cst17.aby.plot, cst17.abn.plot, cst714.plot, cst714.aby.plot, cst714.abn.plot, nrow = 2, ncol = 3))

ggplot2::ggsave(filename = "CST transition btw days_updated.png", 
       plot = cst.trans,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/CST", 
       width = 12, 
       height = 7, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")
```

## Alpha diversity of 3-timepoint patients

```{r}
#Alpha diversity by multiple method
#----#
patient.cont <- ps_filter(data.04en, patient_ID %in% sample_data(patient.17.ps)$patient_ID)
adiv.3d <- data.frame(
  "patient" = phyloseq::sample_data(patient.cont)$patient_ID,
  "day" = phyloseq::sample_data(patient.cont)$day,
  "age" = phyloseq::sample_data(patient.cont)$age_month,
  "breastfeeding" = sample_data(patient.cont)$breastfeeding,
  "hos_days" = phyloseq::sample_data(patient.cont)$hospitalization_days,
  "timetostopdiarr" = phyloseq::sample_data(patient.cont)$timetostopdiarr,
  "antibiotic_trt" = sample_data(patient.cont)$antibiotic_trt,
  "totaldiarrdur" = sample_data(patient.cont)$totaldiarrdur,
  "Shannon" = phyloseq::estimate_richness(patient.cont, measures = "Shannon"))

comps.3d <- make_pairs(sample_data(patient.cont)$day)

dt04.test <- dplyr::select(as.data.frame(unclass(sample_data(data.04en))), c("patient_ID", "sex", "age_month", "breastfeeding", "wfa_zscore", "antibiotic_trt", "hospitalization_days", "timetostopdiarr", "Infection_type", "fever_bl", "vomit_bl"))

dt04.test <- unique(dt04.test)

levels(dt04.test$Infection_type) <- c("bacteria_only", "mixed", "mixed", "unknown", "virus_only")

dt04.test$Infection_type <- factor(dt04.test$Infection_type, levels = c("virus_only", "bacteria_only", "mixed", "unknown"))

dt04.test[which(dt04.test$fever_bl == "unknown"),"fever_bl"] <- NA

md.pattern(dt04.test)

dt04.test <- data.frame(dt04.test,
                        impute_fever = complete(mice(dt04.test, method = "pmm"))$fever_bl)

dt04.test[,c("fever_bl", "impute_fever")]

hosdays_model <- glm(hospitalization_days ~ age_month + breastfeeding + antibiotic_trt + Infection_type + wfa_zscore + vomit_bl + impute_fever, data = dt04.test)

summary(hosdays_model)

ttstop_model <- glm(timetostopdiarr ~ age_month + breastfeeding + antibiotic_trt + Infection_type + wfa_zscore + vomit_bl + impute_fever, data = dt04.test)

summary(ttstop_model)

tbl_uvregression(dt04.test[,!names(dt04.test) %in% c("patient_ID", "timetostopdiarr", "fever_bl")],
                 method = lm,
                 y = hospitalization_days,
                 exponentiate = T)

tbl_uvregression(dt04.test[,!names(dt04.test) %in% c("patient_ID", "hospitalization_days", "fever_bl")],
                 method = glm,
                 y = timetostopdiarr,
                 method.args = list(family = "poisson"),
                 exponentiate = T)

#hospitalization days +1.5 day in antibiotic_yes. add this idea into METHODS

#antibiotic_trt:
##YES: 35 cases = 20 vir (57%), 6 vir_bac (17%), 5 vir_par (14%), 2 bac (6%), 2 unk
##NO: 55 cases = 27 vir (49%), 15 vir_bac (27%), 4 vir_par (7%), 3 bac (5%), 6 unk

comps.3d <- make_pairs(sample_data(patient.3d.ps)$day)

#Pair Wilcoxon test
##Shannon
adiv.3d
adiv.shan <- adiv.3d %>%
  pivot_wider(names_from = "day",
              values_from = "Shannon", names_prefix = "day") %>%
  as.data.frame()
adiv.shan

wilcox.test(adiv.shan$day1,adiv.shan$day7, paired = TRUE, data = adiv.shan) #p = 0.01172
wilcox.test(na.omit(adiv.shan)$day7, na.omit(adiv.shan)$day14, paired = TRUE, data = na.omit(adiv.shan)) #p = 0.04609
wilcox.test(na.omit(adiv.shan)$day1, na.omit(adiv.shan)$day14, paired = TRUE, data = na.omit(adiv.shan)) #p = 0.0004415

adiv.3d$day <- as.numeric(as.character(adiv.3d$day))
adiv.3d$antibiotic_trt <- as.factor(adiv.3d$antibiotic_trt)
levels(adiv.3d$antibiotic_trt) <- c("Antibiotic: No", "Antibiotic: Yes")

#Line plot for the change of alpha diversity in each patient
(line_plot.3d <- ggplot(adiv.3d, aes(x = day, y = Shannon)) +
  geom_line(aes(group = patient, color = hos_days), alpha = 0.7, size = .7) +
  geom_smooth(fill = "grey50", alpha = 0.5, color = "black", method = "loess", linetype = "dashed") +
  scale_color_gradientn(colors = c( "#a8dadc", "#1d3557", "#E63946"), 
                        values = c(0, 0.75, 1), 
                        limits = c(1,14),
                        breaks = c(1, 4, 7, 11, 14),
                        guide = "colorbar", 
                        name = "Hospitalization \nday\n")  +
  labs(x = "Day", y = "Shannon Alpha-diversity", title = "d") +
  theme_test() +
    guides(alpha = F,
           color = guide_colourbar(barwidth = 1, barheight = 5)) +
  facet_grid(.~antibiotic_trt) +
    scale_x_continuous(breaks = c(1, 7, 14), expand = expansion(mult = c(0.01,0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) +
    theme(legend.position='right', 
          legend.justification='right',
          legend.direction='vertical',
          legend.title = element_text(size = 10, vjust = 0, hjust = 0),
          legend.text = element_text(size = 10, color = "grey30", vjust = 0.2, hjust = 0),
          plot.title = element_text(size = 15, face = "bold", vjust = 5, hjust = 0),
          strip.background = element_rect(fill = "white", color = "black", size = 1),
          strip.text = element_text(face = "bold", color = "grey30", hjust = 0.5),
          aspect.ratio = 0.95,
          plot.margin = ggplot2::margin(10,0.5,10,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.ticks.y = element_line(color="grey30"),
          axis.line.x = element_line(size = .1, color="grey30"),
          axis.text.x = element_text(color = "grey30", vjust = -1, size = 8),
          axis.text.y= element_text(color = "grey30", size = 8),
          axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -2),
        axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3))
)

ggplot2::ggsave(filename = "Alpha diversity per patients.png", 
       plot = line_plot.3d,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Alpha", 
       width = 6, 
       height = 3, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")


chisq.test(dt04.test$antibiotic_trt, dt04.test$Infection_type, correct = FALSE) #p = 0.573

hos.plot <- ggplot(dt04.test, 
                   aes(x = antibiotic_trt, 
                       y = hospitalization_days, 
                       fill = antibiotic_trt)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("#e07a5f", "#f4f1de"), guide = guide_legend(title = "Antibiotic treatment \n")) +
  stat_compare_means(method = "wilcox.test") + ylim(0, 17) +
  theme_classic()

dtime.plot <- ggplot(dt04.test, 
                     aes(x = antibiotic_trt, 
                         y = timetostopdiarr, 
                         fill = antibiotic_trt)) +
  geom_boxplot() + scale_fill_manual(values = c("#e07a5f", "#f4f1de"), guide = guide_legend(title = "Antibiotic treatment \n")) + ylim(0, 400) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic()

clin.plot <- ggarrange(hos.plot, dtime.plot, align = "h", common.legend = T)

ggplot2::ggsave(filename = "atb_outcome.pdf", 
       plot = clin.plot,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")

adiv.3d.test <- pivot_wider(data = adiv.3d, names_from = day, values_from = Shannon, names_glue = "Shannon_{day}")
adiv.3d.test <- as.data.frame(adiv.3d.test)
adiv.3d.test <- adiv.3d.test %>% mutate(Shannon17 = Shannon_7 - Shannon_1)
adiv.3d.test <- adiv.3d.test %>% mutate(Shannon714 = Shannon_14 - Shannon_7)
adiv.3d.test <- adiv.3d.test %>% mutate(Shannon114 = Shannon_14 - Shannon_1)
adiv.3d.test$antibiotic_trt <- as.factor(adiv.3d.test$antibiotic_trt)


with(adiv.3d.test, cor(timetostopdiarr, Shannon17, method = "kendall"))
with(na.omit(adiv.3d.test), cor(timetostopdiarr, Shannon714, method = "kendall"))

wilcox.test(Shannon17 ~ antibiotic_trt, data = adiv.3d.test) #p-value = 0.003713
wilcox.test(Shannon714 ~ antibiotic_trt, adiv.3d.test) #p-value = 0.04445

wilcox.test(timetostopdiarr ~ antibiotic_trt, data = adiv.3d.test) #p-value = 0.1388

#not significant, but using antibiotic could increase the timetostopdiarr
ggplot(adiv.3d.test, aes(x = antibiotic_trt, y = hos_days, fill = antibiotic_trt)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

ggplot(adiv.3d.test, aes(x = antibiotic_trt, y = timetostopdiarr, fill = antibiotic_trt)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

#
(adiv17.ab <- ggplot(adiv.3d.test, aes(x = antibiotic_trt, y = Shannon17)) +
    geom_boxplot(outlier.size = 1, width = 0.5) +
    stat_compare_means(method = "wilcox.test", label.y = 2.2, size = 3) +
    labs(x = "Antibiotic treatment",
         y = "Shannon D7 - Shannon D1") +
    theme_classic() +
    ylim(-2,2.5) +
    labs(title = "c") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line.x = element_line(size = .1, color="grey30"),
          axis.text.x = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y=element_text(size=10),
          axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -3),
        axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3))
)

(adiv714.ab <- ggplot(na.omit(adiv.3d.test), aes(x = antibiotic_trt, y = Shannon714)) +
    geom_boxplot(outlier.size = 1, width = 0.5) +
    stat_compare_means(method = "wilcox.test", label.y = 2.2, size  = 3) +
    labs(x = "Antibiotic treatment",
         y = "Shannon D14 - Shannon D7") +
    theme_classic() +
    ylim(-2,2.5) +
    labs(title = "d") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line.x = element_line(size = .1, color="grey30"),
          axis.text.x = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y=element_text(size=10),
          axis.title.x = element_text(color = "grey30", face = "bold", size = 10, vjust = -3),
        axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3))
  )

adiv.3d.test.1 <- adiv.3d.test %>%
  pivot_longer(cols = c("Shannon17", "Shannon714", "Shannon114"),
               names_to = "Week",
               values_to = "Shannon_diff",
               values_drop_na = F) %>% as.data.frame()

adiv.ab <- ggplot(adiv.3d.test.1[which(adiv.3d.test.1$Week != "Shannon114"),], 
       aes(x = Week, 
           y = Shannon_diff, 
           fill = antibiotic_trt)) +
    geom_boxplot(outlier.size = 1, width = 0.5) +
  scale_fill_manual(values = c("#e07a5f", "#f4f1de"), guide = guide_legend(title = "Antibiotic \ntreatment \n")) +
  scale_x_discrete(labels = c("Week 1", "Week 2")) +
    stat_compare_means(method = "wilcox.test", 
                       label.y = 2.2, 
                       size = 3, 
                       label = "p.format") +
    labs(x = "",
         y = "Difference of Shannon index") +
    theme_classic() +
    ylim(-2,2.5) +
    labs(title = "c") +
  guides(x.sec = "axis",
         y.sec = "axis") +
    theme(plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.ticks.x = element_line(color="grey30"),
          axis.line = element_line(size = .5, color="grey30"),
          axis.text.x.bottom = element_text(color = "grey30", vjust = -1, size = 10),
          axis.text.y.left =element_text(size=10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
        axis.title.y = element_text(color = "grey30", face = "bold", size = 10, angle = 90, vjust = 3))


#(adiv.ab <- ggarrange(adiv17.ab, adiv714.ab))

#(adiv.ab <- annotate_figure(adiv.ab, top = text_grob("Figure 3: Alpha diversity changes and Antibiotic treatment", 
               #color = "black", face = "bold", size = 15, hjust = 0.5))
#)

ggplot2::ggsave(filename = "adivdiff.ab.png", 
       plot = adiv.ab,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Alpha/test/", 
       width = 5, 
       height = 4, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")

meta.df$breastfeeding <- as.factor(meta.df$breastfeeding)
meta.df$antibiotic_trt <- as.factor(meta.df$antibiotic_trt)

#filter distance matrix for the data 04EN
dt04.dist <- gp.all.dist.mat[which(rownames(gp.all.dist.mat) %in% sample_names(data.04en)), which(colnames(gp.all.dist.mat) %in% sample_names(data.04en))]
dt04.dist <- as.dist(dt04.dist)

#adonis test for breastfeeding*antibiotic*timetostopdiarr
adonis2(dt04.dist ~ breastfeeding*antibiotic_trt*timetostopdiarr, subset(meta.df, study_ID == "04EN"), by = "term", permutations = 999)
adonis2(dt04.dist ~ breastfeeding + antibiotic_trt + timetostopdiarr, subset(meta.df, study_ID == "04EN"), by = "term", permutations = 999)

##Arrange all a-div plot
adiv.plot
adiv.ab
line_plot.3d

(adiv.summary.plot <- ggarrange(ggarrange(adiv.plot, area.plot, widths = c(1.3,1), align = "hv"), ggarrange(adiv.ab, line_plot.3d, ncol = 2, align = "h"), nrow = 2, heights = c(5,3), align = "hv"))

ggplot2::ggsave(filename = "summary_alpha.pdf", 
       plot = adiv.summary.plot,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Diversity/Alpha", 
       width = 10, 
       height =8, 
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")

```

# Linear Mixed Effects

```{r}
data.04en

#Alpha
df.04en <- data.frame(sample_data(data.04en))
df.04en$Infection_type <- as.factor(df.04en$Infection_type)

levels(df.04en$Infection_type) <- c("bacteria_only", "mixed", "mixed", "unknown", "virus_only")

adiv.04en <- adiv[rownames(df.04en),]

#Linear Mixed-Effects Models

##Shannon
shan.04en <- lmerTest::lmer(adiv.04en$Shannon ~ sex + wfa_zscore + age_month + Infection_type + day*antibiotic_trt + (1|patient_ID), df.04en)

summary(shan.04en)

shan.04en.df <- as.data.frame(summary(shan.04en)$coefficients) %>% rownames_to_column("term")

##Richness
chao1.04en <- lme4::glmer(adiv.04en$Chao1 ~ sex + age_month + wfa_zscore + Infection_type + day*antibiotic_trt + (1|patient_ID), df.04en, family = poisson())

summary(chao1.04en)

##Simpson
simp.04en <- lmerTest::lmer(adiv.04en$Simpson ~ sex + wfa_zscore + age_month + Infection_type + day*antibiotic_trt + (1|patient_ID), df.04en)

summary(simp.04en)

simp.04en.df <- as.data.frame(summary(simp.04en)$coefficients) %>% rownames_to_column("term")

#Export table
shan.table <- as_flextable(shan.04en)

save_as_docx(shan.table, path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/shannon_table.docx")

chao1.table <- as_flextable(chao1.04en)

save_as_docx(chao1.table, path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/chao1_table.docx")

simp.table <- as_flextable(simp.04en)

save_as_docx(simp.table, path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/simp_table.docx")
```


---
>Note of ANCOMBC when qval (or pval) = 0 (from github): These taxa were considered to have structural zeros. For example, when comparing group A vs. B for taxon X, if all samples in group A have 0 abundance, while the abundances of taxon X in group B are not all 0s, we will consider taxon X has structural zeros in group A, and it will be declared to be significant (one group has something while another has nothing) no matter what the test statistic is. P-value and q-value will also be manually set to be 0.

>In this markdown, for the better visualization, all zero qval are changed to the new qval = min(q)^1.2 

# Differential Abundance Analysis:

*Function output ancombc and output deseq2*

```{r}
## create function to output ancombc results
output_ancombc <- function(physeq, ancom.out, group){
      ancom.res <- ancom.out$res
      if (ncol(ancom.res$diff_abn) == 1) {
            difftx <- rownames(ancom.out$feature_table[which(ancom.res$diff_abn[,group] ==TRUE),])
      } else {
            difftx <- ancom.res$diff_abn[which(ancom.res$diff_abn[,group] ==TRUE),"taxon"]
      }
      rs <- data.frame(difftx, ancom.res$W[ancom.res$W[,"taxon"] %in% difftx, group])
      colnames(rs) <- c("taxa", "test_stat")
      rs$qvalue <- as.numeric(ancom.res$q_val[ancom.res$W[,"taxon"] %in% difftx, group])
      rs$pvalue <- as.numeric(ancom.res$p_val[ancom.res$W[,"taxon"] %in% difftx, group])
      rs$lfc <- ancom.res$lfc[ancom.res$W[,"taxon"] %in% difftx, group]
      rs$Species <- unclass(phyloseq::tax_table(all.norm)[difftx,'Species'])
      rownames(rs) <- rs$taxa
      rs <- rs[,c(-1)]
      rs <- rs[order(rs$test_stat, decreasing=TRUE),]
      return(rs)
}

## create function to summarize DESeq2 output
output_deseq2 <- function(physeq, ds, contrast, alpha, coef, testtype) {
      rs <- results(ds, name=contrast, alpha=alpha, pAdjustMethod ="BH")
      rs_shrink <- lfcShrink(dds=ds, res=rs, coef=coef, type=testtype) ## correcting log2foldchange
      #rs_shrink <- rs_shrink[order(rs_shrink$padj, na.last=T),]
      rs_shrink$prank <- seq(1,nrow(rs_shrink))
      rs_shrink$OTUs <- rownames(rs_shrink)
      rs_shrink <- as.data.frame(rs_shrink)
      rs_shrink <- rs_shrink[!is.na(rs_shrink$padj),]
      rs_shrink <- rs_shrink[rs_shrink$padj < alpha,]
      rs_shrink$Species <- unclass(phyloseq::tax_table(physeq)[rownames(rs_shrink),'Species'])
      rs_shrink <- rs_shrink[order(rs_shrink$log2FoldChange, decreasing=TRUE),]
      return(rs_shrink)
}
```

## ANCOMBC

### Day 1 and Day 7
```{r}
patient.17.ps
patient.17.keep <- prune_taxa(rowSums(abundances(patient.17.ps)>0)>=5, patient.17.ps) #at least 5 samples = 3.7% of the pool
patient.17.keep #134 samples and 222 taxa
summary(sample_sums(patient.17.keep)/sample_sums(patient.17.ps)) #92.45% data is covered

patient.17.un
patient.17.un.keep <- prune_taxa(rowSums(abundances(patient.17.un)>0)>=5, patient.17.un)
patient.17.un.keep
summary(sample_sums(patient.17.un.keep)/sample_sums(patient.17.un))

#filter aby and abn
aby.17.un.keep <- ps_filter(patient.17.un.keep, antibiotic_trt == "yes")
abn.17.un.keep <- ps_filter(patient.17.un.keep, antibiotic_trt == "no")

#----------#
#ANCOMBC
ancom_17.aby <- ancombc(data = aby.17.un.keep, 
                     formula = "patient_ID + day",
                     tax_level = NULL,
                     p_adj_method = "fdr", 
                     lib_cut = 1000, 
                     prv_cut=0.001,
                     group = "day", 
                     struc_zero = TRUE, 
                     neg_lb = TRUE,
                     tol = 1e-05,
                     alpha = 0.05, 
                     global = FALSE) ## consider conserve=TRUE -> "sample size > 30" = large

dim(ancom_17.aby$feature_table) 
ancom_res.17.aby <- ancom_17.aby$res

length(which(ancom_res.17.aby$q_val[,'day7'] < 0.05))
length(which(ancom_res.17.aby$diff_abn[,'day7'] == TRUE))

d17.ancom.out.aby <- output_ancombc(aby.17.un.keep, ancom.out = ancom_17.aby, group="day7")

d17.ancom.out.aby

ancom_17.abn <- ancombc(data = abn.17.un.keep, 
                     formula = "patient_ID + day",
                     tax_level = NULL,
                     p_adj_method = "fdr", 
                     lib_cut = 1000, 
                     prv_cut=0.001,
                     group = "day", 
                     struc_zero = TRUE, 
                     neg_lb = TRUE,
                     tol = 1e-05,
                     alpha = 0.05, 
                     global = FALSE) ## consider conserve=TRUE -> "sample size > 30" = large

dim(ancom_17.abn$feature_table) 
ancom_res.17.abn <- ancom_17.abn$res

length(which(ancom_res.17.abn$q_val[,'day7'] < 0.05))
length(which(ancom_res.17.abn$diff_abn[,'day7'] == TRUE))

d17.ancom.out.abn <- output_ancombc(abn.17.un.keep, ancom.out = ancom_17.abn, group="day7")

d17.ancom.out.abn

#transform log e to log2fc
d17.ancom.out.aby$lfc <- log2(exp(d17.ancom.out.aby$lfc))

d17.ancom.out.abn$lfc <- log2(exp(d17.ancom.out.abn$lfc))
```

### Day 7 and Day 14
```{r}
patient.714.ps
patient.714.keep <- prune_taxa(rowSums(abundances(patient.714.ps)>0)>=5, patient.714.ps) #at least 5 samples = 3.7% of the pool
patient.714.keep #76 samples and 189 taxa
summary(sample_sums(patient.714.keep)/sample_sums(patient.714.ps)) #90.49% data is covered

patient.714.un
patient.714.un.keep <- prune_taxa(rowSums(abundances(patient.714.un)>0)>=5, patient.714.un)
patient.714.un.keep
summary(sample_sums(patient.714.un.keep)/sample_sums(patient.714.un))

#filter aby and abn
aby.714.un.keep <- ps_filter(patient.714.un.keep, antibiotic_trt == "yes")
abn.714.un.keep <- ps_filter(patient.714.un.keep, antibiotic_trt == "no")

#----------#
#ANCOMBC
ancom_714.aby <- ancombc(data = aby.714.un.keep, 
                     formula = "patient_ID + day",
                     tax_level = NULL,
                     p_adj_method = "fdr", 
                     lib_cut = 1000, 
                     prv_cut=0.001,
                     group = "day", 
                     struc_zero = TRUE, 
                     neg_lb = TRUE,
                     tol = 1e-05,
                     alpha = 0.05, 
                     global = FALSE) ## consider conserve=TRUE -> "sample size > 30" = large

dim(ancom_714.aby$feature_table) 
ancom_res.714.aby <- ancom_714.aby$res

length(which(ancom_res.714.aby$q_val[,'day14'] < 0.05))
length(which(ancom_res.714.aby$diff_abn[,'day14'] == TRUE))

d714.ancom.out.aby <- output_ancombc(aby.714.un.keep, ancom.out = ancom_714.aby, group="day14")

d714.ancom.out.aby

ancom_714.abn <- ancombc(data = abn.714.un.keep, 
                     formula = "patient_ID + day",
                     tax_level = NULL,
                     p_adj_method = "fdr", 
                     lib_cut = 1000, 
                     prv_cut=0.001,
                     group = "day", 
                     struc_zero = TRUE, 
                     neg_lb = TRUE,
                     tol = 1e-05,
                     alpha = 0.05, 
                     global = FALSE) ## consider conserve=TRUE -> "sample size > 30" = large

dim(ancom_714.abn$feature_table) 
ancom_res.714.abn <- ancom_714.abn$res

length(which(ancom_res.714.abn$q_val[,'day14'] < 0.05))
length(which(ancom_res.714.abn$diff_abn[,'day14'] == TRUE))

d714.ancom.out.abn <- output_ancombc(abn.714.un.keep, ancom.out = ancom_714.abn, group="day14")

d714.ancom.out.abn

#transform log e to log2fc
d714.ancom.out.aby$lfc <- log2(exp(d714.ancom.out.aby$lfc))

d714.ancom.out.abn$lfc <- log2(exp(d714.ancom.out.abn$lfc))
```

---

## DESeq2


### Day 1 and Day 7
```{r}
#----------#
#DeSeq2
length(which(abundances(patient.17.keep) ==0))/(phyloseq::nsamples(patient.17.keep)*phyloseq::ntaxa(patient.17.keep)) ## 74.74% are zeros => highly sparse

#filter aby and abn
aby.17.keep <- ps_filter(patient.17.keep, antibiotic_trt == "yes")
abn.17.keep <- ps_filter(patient.17.keep, antibiotic_trt == "no")

#Antibiotic_yes:
## Consider applying the zero-inflated NB model 
sample_data(aby.17.keep)$day ## Day 1 is reference level
ds.17.aby <- phyloseq_to_deseq2(aby.17.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.17.aby <- zinbwave(ds.17.aby, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 
assay(ds.17.aby, "weights")[1:5, 1:5]
ds.17.aby.dt <- DESeqDataSet(ds.17.aby, design=~0 + patient_ID + day)
DESeq2::sizeFactors(ds.17.aby.dt) <- norm_factors[rownames(sample_data(aby.17.keep))]

## Use size factor normalization by wrench
ds.17.aby.dt <- DESeq(ds.17.aby.dt, test="LRT", reduced=~patient_ID, minmu=0.5,
                     minReplicatesForReplace = 7, fitType = "local", betaPrior = FALSE, useT=TRUE)

## Observe dispersion plot
DESeq2::plotMA(ds.17.aby.dt)
plotDispEsts(ds.17.aby.dt)

## Extract results from DESeq2 analysis
resultsNames(ds.17.aby.dt)
results(ds.17.aby.dt)

rs.17.aby <- results(ds.17.aby.dt, pAdjustMethod = "BH")
rs.17.aby.df <- as.data.frame(lfcShrink(ds.17.aby.dt, res = rs.17.aby, contrast = "day7", coef = 26, lfcThreshold = 0, type = "ashr"))
rs.17.aby.fil <- filter(rs.17.aby.df, baseMean >= 20)
rs.17.aby.fil$Species <- unclass(phyloseq::tax_table(patient.17.keep)[rownames(rs.17.aby.fil),'Species'])

rs.17.aby.fil

#Antibiotic_no
## Consider applying the zero-inflated NB model 
sample_data(abn.17.keep)$day ## Day 1 is reference level
ds.17.abn <- phyloseq_to_deseq2(abn.17.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.17.abn <- zinbwave(ds.17.abn, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 
assay(ds.17.abn, "weights")[1:5, 1:5]
ds.17.abn.dt <- DESeqDataSet(ds.17.abn, design=~0 + patient_ID + day)
DESeq2::sizeFactors(ds.17.abn.dt) <- norm_factors[rownames(sample_data(abn.17.keep))]

## Use size factor normalization by wrench
ds.17.abn.dt <- DESeq(ds.17.abn.dt, test="LRT", reduced=~patient_ID, minmu=0.5,
                     minReplicatesForReplace = 7, fitType = "local", betaPrior = FALSE, useT=TRUE)

## Observe dispersion plot
DESeq2::plotMA(ds.17.abn.dt)
plotDispEsts(ds.17.abn.dt)

## Extract results from DESeq2 analysis
resultsNames(ds.17.abn.dt)
results(ds.17.abn.dt)

rs.17.abn <- results(ds.17.abn.dt, pAdjustMethod = "BH")
rs.17.abn.df <- as.data.frame(lfcShrink(ds.17.abn.dt, res = rs.17.abn, contrast = "day7", coef = 26, lfcThreshold = 0, type = "ashr"))
rs.17.abn.fil <- filter(rs.17.abn.df, !baseMean <= 20)
rs.17.abn.fil$Species <- unclass(phyloseq::tax_table(patient.17.keep)[rownames(rs.17.abn.fil),'Species'])

rs.17.abn.fil
```

### Day 7 and Day 14
```{r}
#----------#
#DeSeq2
length(which(abundances(patient.714.keep) ==0))/(phyloseq::nsamples(patient.714.keep)*phyloseq::ntaxa(patient.714.keep)) ## 69.30% are zeros => highly sparse

#filter aby and abn
aby.714.keep <- ps_filter(patient.714.keep, antibiotic_trt == "yes")
abn.714.keep <- ps_filter(patient.714.keep, antibiotic_trt == "no")

#Antibiotic_yes:
## Consider applying the zero-inflated NB model 
sample_data(aby.714.keep)$day ## Day 7 is reference level
ds.714.aby <- phyloseq_to_deseq2(aby.714.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.714.aby <- zinbwave(ds.714.aby, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 
assay(ds.714.aby, "weights")[1:5, 1:5]
ds.714.aby.dt <- DESeqDataSet(ds.17.aby, design=~0 + patient_ID + day)
DESeq2::sizeFactors(ds.714.aby.dt) <- norm_factors[rownames(sample_data(aby.714.keep))]

## Use size factor normalization by wrench
ds.714.aby.dt <- DESeq(ds.714.aby.dt, test="LRT", reduced=~patient_ID, minmu=0.5,
                     minReplicatesForReplace = 7, fitType = "local", betaPrior = FALSE, useT=TRUE)

## Observe dispersion plot
DESeq2::plotMA(ds.714.aby.dt)
plotDispEsts(ds.714.aby.dt)

## Extract results from DESeq2 analysis
resultsNames(ds.714.aby.dt)
results(ds.714.aby.dt)

rs.714.aby <- results(ds.714.aby.dt, pAdjustMethod = "BH")
rs.714.aby.df <- as.data.frame(lfcShrink(ds.714.aby.dt, res = rs.714.aby, contrast = "day14", coef = 26, lfcThreshold = 0, type = "ashr"))
rs.714.aby.fil <- filter(rs.714.aby.df, !baseMean <= 20)
rs.714.aby.fil$Species <- unclass(phyloseq::tax_table(patient.714.ps)[rownames(rs.714.aby.fil),'Species'])

rs.714.aby.fil

#Antibiotic_no
## Consider applying the zero-inflated NB model 
sample_data(abn.714.keep)$day ## Day 7 is reference level
ds.714.abn <- phyloseq_to_deseq2(abn.714.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.714.abn <- zinbwave(ds.714.abn, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 
assay(ds.714.abn, "weights")[1:5, 1:5]
ds.714.abn.dt <- DESeqDataSet(ds.714.abn, design=~0 + patient_ID + day)
DESeq2::sizeFactors(ds.17.abn.dt) <- norm_factors[rownames(sample_data(abn.17.keep))]

## Use size factor normalization by wrench
ds.714.abn.dt <- DESeq(ds.714.abn.dt, test="LRT", reduced=~patient_ID, minmu=0.5,
                     minReplicatesForReplace = 7, fitType = "local", betaPrior = FALSE, useT=TRUE)

## Observe dispersion plot
DESeq2::plotMA(ds.714.abn.dt)
plotDispEsts(ds.714.abn.dt)

## Extract results from DESeq2 analysis
resultsNames(ds.714.abn.dt)
results(ds.714.abn.dt)

rs.714.abn <- results(ds.714.abn.dt, pAdjustMethod = "BH")
rs.714.abn.df <- as.data.frame(lfcShrink(ds.714.abn.dt, res = rs.714.abn, contrast = "day14", coef = 26, lfcThreshold = 0, type = "ashr"))
rs.714.abn.fil <- filter(rs.714.abn.df, !baseMean <= 20)
rs.714.abn.fil$Species <- unclass(phyloseq::tax_table(patient.714.keep)[rownames(rs.714.abn.fil),'Species'])

rs.714.abn.fil
```

---

## Limma Voom

### Day 1 and Day 7
```{r}
#Antibiotic_yes
#Limma voom with TMM norm and ZINB
ds.17.lm.aby <- phyloseq_to_deseq2(aby.17.un.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.17.lm.aby <- zinbwave(ds.17.lm.aby, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 

ds.17.tmm.aby <- calcNormFactors(ds.17.lm.aby, method = "TMMwsp")           #store TMM norm factor
ds.17.tmm.aby$samples$norm.factors

mm.17.aby <- model.matrix(~ day, ds.17.tmm.aby$samples)              #construct model matrix          
head(mm.17.aby)

y.17.aby <- voom(ds.17.tmm.aby, mm.17.aby, plot = T)                          #obtain Voom weights
fit.17.aby <- lmFit(y.17.aby, mm.17.aby)                                   #fit lm with limma
fit.17.aby <- eBayes(fit.17.aby)
head(coef(fit.17.aby))

limma_res_df.17.aby <- data.frame(topTable(fit.17.aby, coef = "day7", number = Inf))    #extract results
limma_res_df.17.aby$Change <- 10^limma_res_df.17.aby$logFC
limma_res_df.17.aby$log2FoldChange <- log(limma_res_df.17.aby$Change, 2)
limma_res_df.17.aby$Species <-  unclass(phyloseq::tax_table(all.fil)[rownames(limma_res_df.17.aby),'Species'])

fdr_limma.17.aby <- limma_res_df.17.aby %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    rownames_to_column(var = "OTUs")


fdr_limma.17.aby

#Antibiotic_no
#Limma voom with TMM norm and ZINB
ds.17.lm.abn <- phyloseq_to_deseq2(abn.17.un.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.17.lm.abn <- zinbwave(ds.17.lm.abn, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 

ds.17.tmm.abn <- calcNormFactors(ds.17.lm.abn, method = "TMMwsp")           #store TMM norm factor
ds.17.tmm.abn$samples$norm.factors

mm.17.abn <- model.matrix(~ day, ds.17.tmm.abn$samples)              #construct model matrix          
head(mm.17.abn)

y.17.abn <- voom(ds.17.tmm.abn, mm.17.abn, plot = T)                          #obtain Voom weights
fit.17.abn <- lmFit(y.17.abn, mm.17.abn)                                   #fit lm with limma
fit.17.abn <- eBayes(fit.17.abn)
head(coef(fit.17.abn))

limma_res_df.17.abn <- data.frame(topTable(fit.17.abn, coef = "day7", number = Inf))    #extract results
limma_res_df.17.abn$Change <- 10^limma_res_df.17.abn$logFC
limma_res_df.17.abn$log2FoldChange <- log(limma_res_df.17.abn$Change, 2)
limma_res_df.17.abn$Species <-  unclass(phyloseq::tax_table(all.fil)[rownames(limma_res_df.17.abn),'Species'])

fdr_limma.17.abn <- limma_res_df.17.abn %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    rownames_to_column(var = "OTUs")


fdr_limma.17.abn
```

### Day 7 and Day 14
```{r}
#Antibiotic_yes
#Limma voom with TMM norm and ZINB
ds.714.lm.aby <- phyloseq_to_deseq2(aby.714.un.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.714.lm.aby <- zinbwave(ds.714.lm.aby, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 

ds.714.tmm.aby <- calcNormFactors(ds.714.lm.aby, method = "TMMwsp")           #store TMM norm factor
ds.714.tmm.aby$samples$norm.factors

mm.714.aby <- model.matrix(~ day, ds.714.tmm.aby$samples)              #construct model matrix          
head(mm.714.aby)

y.714.aby <- voom(ds.714.tmm.aby, mm.714.aby, plot = T)                          #obtain Voom weights
fit.714.aby <- lmFit(y.714.aby, mm.714.aby)                                   #fit lm with limma
fit.714.aby <- eBayes(fit.714.aby)
head(coef(fit.714.aby))

limma_res_df.714.aby <- data.frame(topTable(fit.714.aby, coef = "day14", number = Inf))    #extract results
limma_res_df.714.aby$Change <- 10^limma_res_df.714.aby$logFC
limma_res_df.714.aby$log2FoldChange <- log(limma_res_df.714.aby$Change, 2)
limma_res_df.714.aby$Species <-  unclass(phyloseq::tax_table(all.fil)[rownames(limma_res_df.714.aby),'Species'])

fdr_limma.714.aby <- limma_res_df.714.aby %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    rownames_to_column(var = "OTUs")


fdr_limma.714.aby

#Antibiotic_no
#Limma voom with TMM norm and ZINB
ds.714.lm.abn <- phyloseq_to_deseq2(abn.714.un.keep, design = ~0 + patient_ID + day) ## design without intercept
ds.714.lm.abn <- zinbwave(ds.714.lm.abn, epsilon=1e10, verbose=TRUE, K=0, observationalWeights=TRUE, BPPARAM=BiocParallel::SerialParam(), X="~1") 

ds.714.tmm.abn <- calcNormFactors(ds.714.lm.abn, method = "TMMwsp")           #store TMM norm factor
ds.714.tmm.abn$samples$norm.factors

mm.714.abn <- model.matrix(~ day, ds.714.tmm.abn$samples)              #construct model matrix          
head(mm.714.abn)

y.714.abn <- voom(ds.714.tmm.abn, mm.714.abn, plot = T)                          #obtain Voom weights
fit.714.abn <- lmFit(y.714.abn, mm.714.abn)                                   #fit lm with limma
fit.714.abn <- eBayes(fit.714.abn)
head(coef(fit.714.abn))

limma_res_df.714.abn <- data.frame(topTable(fit.714.abn, coef = "day14", number = Inf))    #extract results
limma_res_df.714.abn$Change <- 10^limma_res_df.714.abn$logFC
limma_res_df.714.abn$log2FoldChange <- log(limma_res_df.714.abn$Change, 2)
limma_res_df.714.abn$Species <-  unclass(phyloseq::tax_table(all.fil)[rownames(limma_res_df.714.abn),'Species'])

fdr_limma.714.abn <- limma_res_df.714.abn %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    rownames_to_column(var = "OTUs")


fdr_limma.714.abn
```

---

## MaAslin2

>MaAsLin2 with Wrench normalization

### Day 1 and Day 7
```{r}
#antibiotic_yes
mas_17.aby <- Maaslin2(
  input_data = data.frame(otu_table(aby.17.keep)),
  input_metadata = data.frame(sample_data(aby.17.keep)),
  output = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/aby.17",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "NONE", 
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = "day", 
  random_effects = "all",
  reference = "patient_ID", 
  plot_heatmap = FALSE, 
  plot_scatter = FALSE, 
  standardize = FALSE,
  cores = 20)

dfmas.res.17.aby <- mas_17.aby$results

fdr_mas.17.aby <- dfmas.res.17.aby %>%
    dplyr::filter(pval < 0.05)

fdr_mas.17.aby$Species <-  phyloseq::tax_table(aby.17.keep)[fdr_mas.17.aby$feature,'Species']

dim(fdr_mas.17.aby) #7 OTUs
fdr_mas.17.aby[order(fdr_mas.17.aby$coef, decreasing=TRUE),]

fc17.aby<- exp(fdr_mas.17.aby$coef) ## Antilog coef
fdr_mas.17.aby$log2FoldChange <- log2(fc17.aby)
rownames(fdr_mas.17.aby) <- fdr_mas.17.aby$feature

#antibiotic_no
mas_17.abn <- Maaslin2(
  input_data = data.frame(otu_table(abn.17.keep)),
  input_metadata = data.frame(sample_data(abn.17.keep)),
  output = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/abn.17",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "NONE", 
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = "day", 
  random_effects = "all", 
  plot_heatmap = FALSE, 
  plot_scatter = FALSE, 
  standardize = FALSE,
  cores = 20)

dfmas.res.17.abn <- mas_17.abn$results

fdr_mas.17.abn <- dfmas.res.17.abn %>%
    dplyr::filter(pval < 0.05)

fdr_mas.17.abn$Species <-  phyloseq::tax_table(abn.17.keep)[fdr_mas.17.abn$feature,'Species']

dim(fdr_mas.17.abn) #28 OTUs
fdr_mas.17.abn[order(fdr_mas.17.abn$coef, decreasing=TRUE),]

fc17.abn<- exp(fdr_mas.17.abn$coef) ## Antilog coef
fdr_mas.17.abn$log2FoldChange <- log2(fc17.abn)
rownames(fdr_mas.17.abn) <- fdr_mas.17.abn$feature

```

### Day 7 and Day 14

```{r}
#antibiotic_yes
mas_714.aby <- Maaslin2(
  input_data = data.frame(otu_table(aby.714.keep)),
  input_metadata = data.frame(sample_data(aby.714.keep)),
  output = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/aby.714",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "NONE", 
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = "day", 
  random_effects = "all",
  reference = "patient_ID", 
  plot_heatmap = FALSE, 
  plot_scatter = FALSE, 
  standardize = FALSE,
  cores = 20)

dfmas.res.714.aby <- mas_714.aby$results

fdr_mas.714.aby <- dfmas.res.714.aby %>%
    dplyr::filter(pval < 0.05)

fdr_mas.714.aby$Species <-  phyloseq::tax_table(aby.714.keep)[fdr_mas.714.aby$feature,'Species']

dim(fdr_mas.714.aby) #4 OTUs
fdr_mas.714.aby[order(fdr_mas.714.aby$coef, decreasing=TRUE),]

fc714.aby<- exp(fdr_mas.714.aby$coef) ## Antilog coef
fdr_mas.714.aby$log2FoldChange <- log2(fc714.aby)
rownames(fdr_mas.714.aby) <- fdr_mas.714.aby$feature

#antibiotic_no
mas_714.abn <- Maaslin2(
  input_data = data.frame(otu_table(abn.714.keep)),
  input_metadata = data.frame(sample_data(abn.714.keep)),
  output = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/abn.714",
  min_abundance = 0,
  min_prevalence = 0,
  normalization = "NONE", 
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = "day", 
  random_effects = "all", 
  reference = "patient_ID",
  plot_heatmap = FALSE, 
  plot_scatter = FALSE, 
  standardize = FALSE,
  cores = 20)

dfmas.res.714.abn <- mas_714.abn$results

fdr_mas.714.abn <- dfmas.res.714.abn %>%
    dplyr::filter(pval < 0.05)

fdr_mas.714.abn$Species <-  phyloseq::tax_table(abn.714.keep)[fdr_mas.714.abn$feature,'Species']

dim(fdr_mas.714.abn) #28 OTUs
fdr_mas.714.abn[order(fdr_mas.714.abn$coef, decreasing=TRUE),]

fc714.abn<- exp(fdr_mas.714.abn$coef) ## Antilog coef
fdr_mas.714.abn$log2FoldChange <- log2(fc714.abn)
rownames(fdr_mas.714.abn) <- fdr_mas.714.abn$feature

```

---

#Data visualization (for the differential abundance analysis)

## Visualize. Day 1 7

### antibiotic_yes:
```{r}
#upset
color.upsetR <- c("#2C69B0", "#F02720", "#6BA3D6", "#EA6B73")

otu.select17.aby <- list(DeSeq2 = rownames(rs.17.aby.fil), 
                         ANCOMBC = rownames(d17.ancom.out.aby),
                         LimmaVoom = fdr_limma.17.aby$OTUs,
                         Maaslin2 = fdr_mas.17.aby$feature)

otu.upset17.aby <- fromList(otu.select17.aby)
color.upsetR.17.aby <- c("#EA6B73", "#F02720", "#6BA3D6", "#2C69B0")

(upset.plot17.aby <- UpSetR::upset(as.data.frame(otu.upset17.aby), main.bar.color = "grey75", text.scale = 2, point.size = 4,
      sets.bar.color = color.upsetR.17.aby, matrix.color = "grey75", 
      order.by = 'freq', empty.intersections = 'on', 
      queries = list(list(query = intersects, params = list('LimmaVoom'), 
                          color = color.upsetR[1], active = T), 
                     list(query = intersects, params = list('DeSeq2'), 
                          color = color.upsetR[2], active = T),
                     list(query = intersects, params = list('Maaslin2'), 
                          color = color.upsetR[3], active = T),
                     list(query = intersects, params = list('ANCOMBC'), 
                          color = color.upsetR[4], active = T))))


#Otus are significant in at least 2 methods
count.17.aby <- table(unlist(lapply(otu.select17.aby, unique)))
count.17.aby <- as.data.frame(count.17.aby)
count.17.aby <- filter(count.17.aby, Freq >= 2)
count.17.aby <- count.17.aby[order(count.17.aby$Freq, decreasing = T),]
count.17.aby

phyloseq::tax_table(all.fil)[as.character(count.17.aby$Var1), "Species"] #37 taxa 


#Filter only the OTUs that are significant in at least 2 methods in each DA results
d17.ancom.fil.aby <- d17.ancom.out.aby[which(rownames(d17.ancom.out.aby) %in% count.17.aby$Var1),]
rs.17.aby.fil <- rs.17.aby.fil[which(rownames(rs.17.aby.fil) %in% count.17.aby$Var1),]
d17.lm.fil.aby <- limma_res_df.17.aby[which(rownames(limma_res_df.17.aby) %in% count.17.aby$Var1),]
fdr_mas.17.fil.aby <- filter(fdr_mas.17.aby, fdr_mas.17.aby$feature %in% count.17.aby$Var1)


#select only the rownames/log2fc/qval/species:
d17.ancom.fil.aby <- d17.ancom.fil.aby[c(4,2,5)]
rs.17.aby.fil <- rs.17.aby.fil[c(2,5,6)]
d17.lm.fil.aby <- d17.lm.fil.aby[c(8,5,9)]
fdr_mas.17.fil.aby <- fdr_mas.17.aby[c(12, 8, 11)]

#Change column names
d17.ancom.fil.aby <- setNames(d17.ancom.fil.aby, c("log2FC", "p.adj", "Species"))
rs.17.aby.fil <- setNames(rs.17.aby.fil, c("log2FC", "p.adj", "Species"))
d17.lm.fil.aby <- setNames(d17.lm.fil.aby, c("log2FC", "p.adj", "Species"))
fdr_mas.17.fil.aby <- setNames(fdr_mas.17.fil.aby, c("log2FC", "p.adj", "Species"))

#add column to classify the method before bind into 1 dataframe
d17.ancom.fil.aby$method <- rep("ANCOMBC", nrow(d17.ancom.fil.aby))
rs.17.aby.fil$method <- rep("DeSeq2", nrow(rs.17.aby.fil))
d17.lm.fil.aby$method <- rep("Limma-Voom", nrow(d17.lm.fil.aby))
fdr_mas.17.fil.aby$method <- rep("MaAsLin2", nrow(fdr_mas.17.fil.aby))

#add Otu column to avoid the duplicate rownames in combination
d17.ancom.fil.aby$otu <- rownames(d17.ancom.fil.aby)
rs.17.aby.fil$otu <- rownames(rs.17.aby.fil)
d17.lm.fil.aby$otu <- rownames(d17.lm.fil.aby)
fdr_mas.17.fil.aby$otu <- rownames(fdr_mas.17.fil.aby)

rownames(d17.ancom.fil.aby) <- NULL
rownames(rs.17.aby.fil) <- NULL
rownames(d17.lm.fil.aby) <- NULL
rownames(fdr_mas.17.fil.aby) <- NULL
#combine data frame
da.d17.aby <- rbind(d17.ancom.fil.aby, rs.17.aby.fil, d17.lm.fil.aby, fdr_mas.17.fil.aby)
da.d17.aby <- na.omit(da.d17.aby)

#remove duplicate otus

methods_otu <- unique(da.d17.aby$method)
new_da.d17.aby <- data.frame()
for(i in unique(da.d17.aby$otu))
{
  each_otu <- subset(da.d17.aby,otu == i)
  new_da.d17.aby <- rbind(new_da.d17.aby,each_otu[1,])
}
new_da.d17.aby$p.adj <- replace(new_da.d17.aby$p.adj, new_da.d17.aby$p.adj == 0, min(new_da.d17.aby$p.adj[new_da.d17.aby$p.adj > 0])^1.2)

new_da.d17.aby$Phylum <- as.factor(unclass(tax_table(all.fil)[new_da.d17.aby$otu,"Phylum"]))
new_da.d17.aby$Family <- as.factor(unclass(tax_table(all.fil)[new_da.d17.aby$otu,"Family"]))
new_da.d17.aby$Genus <- as.factor(unclass(tax_table(all.fil)[new_da.d17.aby$otu,"Genus"]))


#Volcano Plot
#classify for the significant p.adj
new_da.d17.aby$sig <- "NS"
# if pvalue < 0.05 and abs FC >1, set as "S" 
new_da.d17.aby$sig[new_da.d17.aby$p.adj <= 0.05 & abs(new_da.d17.aby$log2FC) > 0.6] <- "S"

linetype <- c("dotted", "solid")

new_da.d17.aby$Freq <- count.17.aby[order(match(as.character(count.17.aby$Var1), new_da.d17.aby$otu)),]$Freq

(volcano.d17.aby <- ggplot(new_da.d17.aby, aes(x = log2FC, y = -log10(p.adj), shape = method)) +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = max(-log10(p.adj)), ymax = Inf), fill = "grey95") + 
  geom_hline(aes(yintercept = max(-log10(p.adj))), 
             linetype = "dashed",
             color = "grey30") +
  geom_point(data = subset(new_da.d17.aby, sig == "NS"), 
             aes(fill = sig, size = Freq),
             stroke = NA) +
  geom_hline(yintercept=-log10(0.05), col="grey30", 
             linetype = 1) +
  geom_vline(xintercept = c(-0.6, 0.6), col="grey30", 
             linetype = 1) +
  geom_point(data = subset(new_da.d17.aby, sig == "S"), aes(group = Phylum, fill = Phylum, size = Freq)) +
    scale_size_continuous(range = c(4,7), breaks = c(2:4)) +
  scale_fill_manual(values = c("Actinobacteria" = "#EEDD88", "Bacteroidetes" = "#77AADD", "Firmicutes" = "#EE8866", "Proteobacteria" = "#FFAABB",  "NS" = "grey90")) +
  scale_shape_manual(values = c(21:24)) +
  geom_text_repel(data = subset(new_da.d17.aby, sig == "S" & Freq >= 3),
                  aes(label = Species), 
                  size = 5, 
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  box.padding = 1,
                  fontface = "italic") +
  geom_text(aes(x = 8, y = max(-log10(p.adj) + 0.07), label = "p.adj = 0"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  geom_text(aes(x = 8, y = -log10(0.05) + 0.07, label = "p.adj = 0.05"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(expand = c(0,0), limits = c(-8,8))  +
  theme_classic() +
    labs(title = "c",
         subtitle = "Day 7 VS Day 1. Antibiotic_yes") +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
 theme(plot.title = element_text(size = 10, face = "bold", vjust = 0, hjust =0),
        plot.subtitle = element_text(size = 10, color = "grey50", vjust = 0, hjust = 0),
        #plot.margin = unit(c(0.5,0.5,0.5,1.5),"in"),
        axis.title.x = element_text(vjust = 0, hjust = 0.5, size = 10),
        axis.text.x = element_text( size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
       axis.line.y.left = element_blank()) +
  geom_segment(aes(x = -8, y = -Inf, xend = -8, yend = max(-log10(p.adj))), size = 1, linetype = "solid") +
  geom_segment(aes(x = -8, y = max(-log10(p.adj)), xend = -8, yend = Inf), size = 1, linetype = "dotted")
)
```

### antibiotic_no:
```{r}
#upset
otu.select17.abn <- list(DeSeq2 = rownames(rs.17.abn.fil), 
                         ANCOMBC = rownames(d17.ancom.out.abn),
                         LimmaVoom = fdr_limma.17.abn$OTUs,
                         Maaslin2 = fdr_mas.17.abn$feature)

otu.upset17.abn <- fromList(otu.select17.abn)
color.upsetR.17.abn <- c("#F02720", "#EA6B73", "#6BA3D6", "#2C69B0")

(upset.plot17.abn <- UpSetR::upset(as.data.frame(otu.upset17.abn), main.bar.color = "grey75", text.scale = 2, point.size = 4,
      sets.bar.color = color.upsetR.17.abn, matrix.color = "grey75", 
      order.by = 'freq', empty.intersections = 'on', 
      queries = list(list(query = intersects, params = list('LimmaVoom'), 
                          color = color.upsetR[1], active = T), 
                     list(query = intersects, params = list('DeSeq2'), 
                          color = color.upsetR[2], active = T),
                     list(query = intersects, params = list('Maaslin2'), 
                          color = color.upsetR[3], active = T),
                     list(query = intersects, params = list('ANCOMBC'), 
                          color = color.upsetR[4], active = T))))
#Otus are significant in at least 2 methods
count.17.abn <- table(unlist(lapply(otu.select17.abn, unique)))
count.17.abn <- as.data.frame(count.17.abn)
count.17.abn <- filter(count.17.abn, Freq >= 2)
count.17.abn <- count.17.abn[order(count.17.abn$Freq, decreasing = T),]
count.17.abn

phyloseq::tax_table(all.fil)[as.character(count.17.abn$Var1), "Species"] #62 taxa 


#Filter only the OTUs that are significant in at least 2 methods in each DA results
d17.ancom.fil.abn <- d17.ancom.out.abn[which(rownames(d17.ancom.out.abn) %in% count.17.abn$Var1),]
rs.17.abn.fil <- rs.17.abn.fil[which(rownames(rs.17.abn.fil) %in% count.17.abn$Var1),]
d17.lm.fil.abn <- limma_res_df.17.abn[which(rownames(limma_res_df.17.abn) %in% count.17.abn$Var1),]
fdr_mas.17.fil.abn <- filter(fdr_mas.17.abn, fdr_mas.17.abn$feature %in% count.17.abn$Var1)

#select only the rownames/log2fc/qval/species:
d17.ancom.fil.abn <- d17.ancom.fil.abn[c(4,2,5)]
rs.17.abn.fil <- rs.17.abn.fil[c(2,5,6)]
d17.lm.fil.abn <- d17.lm.fil.abn[c(8,5,9)]
fdr_mas.17.fil.abn <- fdr_mas.17.fil.abn[c(12, 8, 11)]

#Change column names
d17.ancom.fil.abn <- setNames(d17.ancom.fil.abn, c("log2FC", "p.adj", "Species"))
rs.17.abn.fil <- setNames(rs.17.abn.fil, c("log2FC", "p.adj", "Species"))
d17.lm.fil.abn <- setNames(d17.lm.fil.abn, c("log2FC", "p.adj", "Species"))
fdr_mas.17.fil.abn <- setNames(fdr_mas.17.fil.abn, c("log2FC", "p.adj", "Species"))

#add column to classify the method before bind into 1 dataframe
d17.ancom.fil.abn$method <- rep("ANCOMBC", nrow(d17.ancom.fil.abn))
rs.17.abn.fil$method <- rep("DeSeq2", nrow(rs.17.abn.fil))
d17.lm.fil.abn$method <- rep("Limma-Voom", nrow(d17.lm.fil.abn))
fdr_mas.17.fil.abn$method <- rep("MaAsLin2", nrow(fdr_mas.17.fil.abn))

#add Otu column to avoid the duplicate rownames in combination
d17.ancom.fil.abn$otu <- rownames(d17.ancom.fil.abn)
rs.17.abn.fil$otu <- rownames(rs.17.abn.fil)
d17.lm.fil.abn$otu <- rownames(d17.lm.fil.abn)
fdr_mas.17.fil.abn$otu <- rownames(fdr_mas.17.fil.abn)

rownames(d17.ancom.fil.abn) <- NULL
rownames(rs.17.abn.fil) <- NULL
rownames(d17.lm.fil.abn) <- NULL
rownames(fdr_mas.17.fil.abn) <- NULL

#combine data frame
da.d17.abn <- rbind(d17.ancom.fil.abn, rs.17.abn.fil, d17.lm.fil.abn, fdr_mas.17.fil.abn)
da.d17.abn <- na.omit(da.d17.abn)

#remove duplicate otus

methods_otu <- unique(da.d17.abn$method)
new_da.d17.abn <- data.frame()
for(i in unique(da.d17.abn$otu))
{
  each_otu <- subset(da.d17.abn,otu == i)
  new_da.d17.abn <- rbind(new_da.d17.abn,each_otu[1,])
}
new_da.d17.abn$p.adj <- replace(new_da.d17.abn$p.adj, new_da.d17.abn$p.adj == 0, min(new_da.d17.abn$p.adj[new_da.d17.abn$p.adj > 0])^1.2)

new_da.d17.abn$Phylum <- as.factor(unclass(tax_table(all.fil)[new_da.d17.abn$otu,"Phylum"]))
new_da.d17.abn$Family <- as.factor(unclass(tax_table(all.fil)[new_da.d17.abn$otu,"Family"]))
new_da.d17.abn$Genus <- as.factor(unclass(tax_table(all.fil)[new_da.d17.abn$otu,"Genus"]))

#Volcano Plot
#classify for the significant p.adj
new_da.d17.abn$sig <- "NS"
# if pvalue < 0.05 and abs FC >1, set as "S" 
new_da.d17.abn$sig[new_da.d17.abn$p.adj <= 0.05 & abs(new_da.d17.abn$log2FC) > 0.6] <- "S"

linetype <- c("dotted", "solid")

new_da.d17.abn$Freq <- count.17.abn[order(match(as.character(count.17.abn$Var1), new_da.d17.abn$otu)),]$Freq
new_da.d17.abn <- filter(new_da.d17.abn, p.adj <= 0.05)


(volcano.d17.abn <- ggplot(new_da.d17.abn, aes(x = log2FC, y = -log10(p.adj), shape = method)) +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = max(-log10(p.adj)), ymax = Inf), fill = "grey95") + 
  geom_hline(aes(yintercept = max(-log10(p.adj))), 
             linetype = "dashed",
             color = "grey30") +
  geom_point(data = subset(new_da.d17.abn, sig == "NS"), 
             aes(fill = sig, size = Freq),
             stroke = NA) +
  geom_hline(yintercept=-log10(0.05), col="grey30", 
             linetype = 1) +
  geom_vline(xintercept = c(-0.6, 0.6), col="grey30", 
             linetype = 1) +
  geom_point(data = subset(new_da.d17.abn, sig == "S"), aes(group = Phylum, fill = Phylum, size = Freq)) +
    scale_size_continuous(range = c(4,7), breaks = c(2:4)) +
  scale_fill_manual(values = c("Actinobacteria" = "#EEDD88", "Bacteroidetes" = "#77AADD", "Firmicutes" = "#EE8866", "Proteobacteria" = "#FFAABB",  "NS" = "grey90")) +
  scale_shape_manual(values = c(21:24)) +
  geom_text_repel(data = subset(new_da.d17.abn, sig == "S" & Freq >= 3),
                  aes(label = Species), 
                  size = 5, 
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  box.padding = 1.3,
                  fontface = "italic") +
  geom_text(aes(x = 8, y = max(-log10(p.adj) + 0.07), label = "p.adj = 0"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  geom_text(aes(x = 8, y = -log10(0.05) + 0.07, label = "p.adj = 0.05"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  scale_y_continuous(expand = c(0,0), limits = c(1,8)) +
  scale_x_continuous(expand = c(0,0), limits = c(-8,8))  +
  theme_classic() +
    labs(title = "a",
         subtitle = "Day 7 VS Day 1. Antibiotic_no") +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
 theme(plot.title = element_text(size = 10, face = "bold", vjust = 0, hjust =0),
        plot.subtitle = element_text(size = 10, color = "grey50", vjust = 0, hjust = 0),
        #plot.margin = unit(c(0.5,0.5,0.5,1.5),"in"),
        axis.title.x = element_text(vjust = 0, hjust = 0.5, size = 10),
        axis.text.x = element_text( size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
       axis.line.y.left = element_blank()) +
  geom_segment(aes(x = -8, y = -Inf, xend = -8, yend = max(-log10(p.adj))), size = 1, linetype = "solid") +
  geom_segment(aes(x = -8, y = max(-log10(p.adj)), xend = -8, yend = Inf), size = 1, linetype = "dotted")
)

volcano.d17 <- ggarrange(volcano.d17.aby, volcano.d17.abn, align = "v", common.legend = T)

ggplot2::ggsave(filename = "Volcano_d17.png", 
       plot = volcano.d17,
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/Volcano/", 
       width = 20, 
       height = 10, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")
```

---

## Visualize. Day 7 and Day 14

### antibiotic_yes:
```{r}
#upset
otu.select714.aby <- list(DeSeq2 = rownames(rs.714.aby.fil), 
                         ANCOMBC = rownames(d714.ancom.out.aby),
                         LimmaVoom = fdr_limma.714.aby$OTUs,
                         Maaslin2 = fdr_mas.714.aby$feature)
otu.upset714.aby <- fromList(otu.select714.aby)
color.upsetR.714.aby <- c("#F02720", "#EA6B73", "#6BA3D6", "#2C69B0")

(upset.plot714.aby <- UpSetR::upset(as.data.frame(otu.upset714.aby), main.bar.color = "grey75", text.scale = 2, point.size = 4,
      sets.bar.color = color.upsetR.714.aby, matrix.color = "grey75", 
      order.by = 'freq', empty.intersections = 'on', 
      queries = list(list(query = intersects, params = list('LimmaVoom'), 
                          color = color.upsetR[1], active = T), 
                     list(query = intersects, params = list('DeSeq2'), 
                          color = color.upsetR[2], active = T),
                     list(query = intersects, params = list('Maaslin2'), 
                          color = color.upsetR[3], active = T),
                     list(query = intersects, params = list('ANCOMBC'), 
                          color = color.upsetR[4], active = T))))
#Otus are significant in at least 2 methods
count.714.aby <- table(unlist(lapply(otu.select714.aby, unique)))
count.714.aby <- as.data.frame(count.714.aby)
count.714.aby <- filter(count.714.aby, Freq >= 2)
count.714.aby <- count.714.aby[order(count.714.aby$Freq, decreasing = T),]
count.714.aby

phyloseq::tax_table(all.fil)[as.character(count.714.aby$Var1), "Species"] #41 taxa 


#Filter only the OTUs that are significant in at least 2 methods in each DA results
d714.ancom.fil.aby <- d714.ancom.out.aby[which(rownames(d714.ancom.out.aby) %in% count.714.aby$Var1),]
rs.714.aby.fil <- rs.714.aby.fil[which(rownames(rs.714.aby.fil) %in% count.714.aby$Var1),]
d714.lm.fil.aby <- limma_res_df.714.aby[which(rownames(limma_res_df.714.aby) %in% count.714.aby$Var1),]
fdr_mas.714.fil.aby <- filter(fdr_mas.714.aby, fdr_mas.714.aby$feature %in% count.714.aby$Var1)


#select only the rownames/log2fc/qval/species:
d714.ancom.fil.aby <- d714.ancom.fil.aby[c(4,2,5)]
rs.714.aby.fil <- rs.714.aby.fil[c(2,5,6)]
d714.lm.fil.aby <- d714.lm.fil.aby[c(8,5,9)]
fdr_mas.714.fil.aby <- fdr_mas.714.fil.aby[c(12, 8, 11)]

#Change column names
d714.ancom.fil.aby <- setNames(d714.ancom.fil.aby, c("log2FC", "p.adj", "Species"))
rs.714.aby.fil <- setNames(rs.714.aby.fil, c("log2FC", "p.adj", "Species"))
d714.lm.fil.aby <- setNames(d714.lm.fil.aby, c("log2FC", "p.adj", "Species"))
fdr_mas.714.fil.aby <- setNames(fdr_mas.714.fil.aby, c("log2FC", "p.adj", "Species"))

#add column to classify the method before bind into 1 dataframe
d714.ancom.fil.aby$method <- rep("ANCOMBC", nrow(d714.ancom.fil.aby))
rs.714.aby.fil$method <- rep("DeSeq2", nrow(rs.714.aby.fil))
d714.lm.fil.aby$method <- rep("Limma-Voom", nrow(d714.lm.fil.aby))
fdr_mas.714.fil.aby$method <- rep("MaAsLin2", nrow(fdr_mas.714.fil.aby))

#add Otu column to avoid the duplicate rownames in combination
d714.ancom.fil.aby$otu <- rownames(d714.ancom.fil.aby)
rs.714.aby.fil$otu <- rownames(rs.714.aby.fil)
d714.lm.fil.aby$otu <- rownames(d714.lm.fil.aby)
fdr_mas.714.fil.aby$otu <- rownames(fdr_mas.714.fil.aby)

rownames(d714.ancom.fil.aby) <- NULL
rownames(rs.714.aby.fil) <- NULL
rownames(d714.lm.fil.aby) <- NULL
rownames(fdr_mas.714.fil.aby) <- NULL

#combine data frame
da.d714.aby <- rbind(d714.ancom.fil.aby, rs.714.aby.fil, d714.lm.fil.aby, fdr_mas.714.fil.aby)
da.d714.aby <- na.omit(da.d714.aby)

#remove duplicate otus

methods_otu <- unique(da.d714.aby$method)
new_da.d714.aby <- data.frame()
for(i in unique(da.d714.aby$otu))
{
  each_otu <- subset(da.d714.aby,otu == i)
  new_da.d714.aby <- rbind(new_da.d714.aby,each_otu[1,])
}
new_da.d714.aby$p.adj <- replace(new_da.d714.aby$p.adj, new_da.d714.aby$p.adj == 0, min(new_da.d714.aby$p.adj[new_da.d714.aby$p.adj > 0])^1.2)

new_da.d714.aby$Phylum <- as.factor(unclass(tax_table(all.fil)[new_da.d714.aby$otu,"Phylum"]))
new_da.d714.aby$Family <- as.factor(unclass(tax_table(all.fil)[new_da.d714.aby$otu,"Family"]))
new_da.d714.aby$Genus <- as.factor(unclass(tax_table(all.fil)[new_da.d714.aby$otu,"Genus"]))


#Volcano Plot
#classify for the significant p.adj
new_da.d714.aby$sig <- "NS"
# if pvalue < 0.05 and abs FC >1, set as "S" 
new_da.d714.aby$sig[new_da.d714.aby$p.adj <= 0.05 & abs(new_da.d714.aby$log2FC) > 0.6] <- "S"
linetype <- c("dotted", "solid")

new_da.d714.aby$Freq <- count.714.aby[order(match(as.character(count.714.aby$Var1),new_da.d714.aby$otu)),]$Freq

(volcano.d714.aby <- ggplot(new_da.d714.aby, aes(x = log2FC, y = -log10(p.adj), shape = method)) +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = max(-log10(p.adj)), ymax = Inf), fill = "grey95") + 
  geom_hline(aes(yintercept = max(-log10(p.adj))), 
             linetype = "dashed",
             color = "grey30") +
  geom_point(data = subset(new_da.d714.aby, sig == "NS"), 
             aes(fill = sig, size = Freq),
             stroke = NA) +
  geom_hline(yintercept=-log10(0.05), col="grey30", 
             linetype = 1) +
  geom_vline(xintercept = c(-0.6, 0.6), col="grey30", 
             linetype = 1) +
  geom_point(data = subset(new_da.d714.aby, sig == "S"), aes(group = Phylum, fill = Phylum, size = Freq)) +
    scale_size_continuous(range = c(4,7), breaks = c(2:4)) +
  scale_fill_manual(values = c("Actinobacteria" = "#EEDD88", "Bacteroidetes" = "#77AADD", "Firmicutes" = "#EE8866", "Proteobacteria" = "#FFAABB",  "NS" = "grey90")) +
  scale_shape_manual(values = c(21:24)) +
  geom_text_repel(data = subset(new_da.d714.aby, sig == "S" & Freq >= 3),
                  aes(label = Species), 
                  size = 5, 
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  box.padding = 1,
                  fontface = "italic") +
  geom_text(aes(x = 8, y = max(-log10(p.adj) + 0.07), label = "p.adj = 0"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  geom_text(aes(x = 8, y = -log10(0.05) + 0.07, label = "p.adj = 0.05"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_x_continuous(expand = c(0,0), limits = c(-8,8))  +
  theme_classic() +
    labs(title = "d",
         subtitle = "Day 14 VS Day 7. Antibiotic_yes") +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
 theme(plot.title = element_text(size = 10, face = "bold", vjust = 0, hjust =0),
        plot.subtitle = element_text(size = 10, color = "grey50", vjust = 0, hjust = 0),
        #plot.margin = unit(c(0.5,0.5,0.5,1.5),"in"),
        axis.title.x = element_text(vjust = 0, hjust = 0.5, size = 10),
        axis.text.x = element_text( size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
       axis.line.y.left = element_blank()) +
  geom_segment(aes(x = -8, y = -Inf, xend = -8, yend = max(-log10(p.adj))), size = 1, linetype = "solid") +
  geom_segment(aes(x = -8, y = max(-log10(p.adj)), xend = -8, yend = Inf), size = 1, linetype = "dotted")
)
```

### antibiotic_no:

```{r}
#upset
otu.select714.abn <- list(DeSeq2 = rownames(rs.714.abn.fil), 
                         ANCOMBC = rownames(d714.ancom.out.abn),
                         LimmaVoom = fdr_limma.714.abn$OTUs,
                         Maaslin2 = fdr_mas.714.abn$feature)
otu.upset714.abn <- fromList(otu.select714.abn)
color.upsetR.714.abn <- c("#F02720", "#EA6B73", "#6BA3D6", "#2C69B0")

(upset.plot714.abn <- UpSetR::upset(as.data.frame(otu.upset714.abn), main.bar.color = "grey75", text.scale = 2, point.size = 4,
      sets.bar.color = color.upsetR.714.aby, matrix.color = "grey75", 
      order.by = 'freq', empty.intersections = 'on', 
      queries = list(list(query = intersects, params = list('LimmaVoom'), 
                          color = color.upsetR[1], active = T), 
                     list(query = intersects, params = list('DeSeq2'), 
                          color = color.upsetR[2], active = T),
                     list(query = intersects, params = list('Maaslin2'), 
                          color = color.upsetR[3], active = T),
                     list(query = intersects, params = list('ANCOMBC'), 
                          color = color.upsetR[4], active = T))))
#Otus are significant in at least 2 methods
count.714.abn <- table(unlist(lapply(otu.select714.abn, unique)))
count.714.abn <- as.data.frame(count.714.abn)
count.714.abn <- filter(count.714.abn, Freq >= 2)
count.714.abn <- count.714.abn[order(count.714.abn$Freq, decreasing = T),]
count.714.abn

phyloseq::tax_table(all.fil)[as.character(count.714.abn$Var1), "Species"] #62 taxa 


#Filter only the OTUs that are significant in at least 2 methods in each DA results
d714.ancom.fil.abn <- d714.ancom.out.abn[which(rownames(d714.ancom.out.abn) %in% count.714.abn$Var1),]
rs.714.abn.fil <- rs.714.abn.fil[which(rownames(rs.714.abn.fil) %in% count.714.abn$Var1),]
d714.lm.fil.abn <- limma_res_df.714.abn[which(rownames(limma_res_df.714.abn) %in% count.714.abn$Var1),]
fdr_mas.714.fil.abn <- filter(fdr_mas.714.abn, fdr_mas.714.abn$feature %in% count.714.abn$Var1)

#select only the rownames/log2fc/qval/species:
d714.ancom.fil.abn <- d714.ancom.fil.abn[c(4,2,5)]
rs.714.abn.fil <- rs.714.abn.fil[c(2,5,6)]
d714.lm.fil.abn <- d714.lm.fil.abn[c(8,5,9)]
fdr_mas.714.fil.abn <- fdr_mas.714.fil.abn[c(12,8,11)]

#Change column names
d714.ancom.fil.abn <- setNames(d714.ancom.fil.abn, c("log2FC", "p.adj", "Species"))
rs.714.abn.fil <- setNames(rs.714.abn.fil, c("log2FC", "p.adj", "Species"))
d714.lm.fil.abn <- setNames(d714.lm.fil.abn, c("log2FC", "p.adj", "Species"))
fdr_mas.714.fil.abn <- setNames(fdr_mas.714.fil.abn, c("log2FC", "p.adj", "Species"))

#add column to classify the method before bind into 1 dataframe
d714.ancom.fil.abn$method <- rep("ANCOMBC", nrow(d714.ancom.fil.abn))
rs.714.abn.fil$method <- rep("DeSeq2", nrow(rs.714.abn.fil))
d714.lm.fil.abn$method <- rep("Limma-Voom", nrow(d714.lm.fil.abn))
fdr_mas.714.fil.abn$method <- rep("MaAsLin2", nrow(fdr_mas.714.fil.abn))

#add Otu column to avoid the duplicate rownames in combination
d714.ancom.fil.abn$otu <- rownames(d714.ancom.fil.abn)
rs.714.abn.fil$otu <- rownames(rs.714.abn.fil)
d714.lm.fil.abn$otu <- rownames(d714.lm.fil.abn)
fdr_mas.714.fil.abn$otu <- rownames(fdr_mas.714.fil.abn)
rownames(d714.ancom.fil.abn) <- NULL
rownames(rs.714.abn.fil) <- NULL
rownames(d714.lm.fil.abn) <- NULL
rownames(fdr_mas.714.fil.abn) <- NULL

#combine data frame
da.d714.abn <- rbind(d714.ancom.fil.abn, rs.714.abn.fil, d714.lm.fil.abn, fdr_mas.714.fil.abn)
da.d714.abn <- na.omit(da.d714.abn)

#remove duplicate otus

methods_otu <- unique(da.d714.abn$method)
new_da.d714.abn <- data.frame()
for(i in unique(da.d714.abn$otu))
{
  each_otu <- subset(da.d714.abn,otu == i)
  new_da.d714.abn <- rbind(new_da.d714.abn,each_otu[1,])
}
new_da.d714.abn$p.adj <- replace(new_da.d714.abn$p.adj, new_da.d714.abn$p.adj == 0, min(new_da.d714.abn$p.adj[new_da.d714.abn$p.adj > 0])^1.2)

new_da.d714.abn$Phylum <- as.factor(unclass(tax_table(all.fil)[new_da.d714.abn$otu,"Phylum"]))
new_da.d714.abn$Family <- as.factor(unclass(tax_table(all.fil)[new_da.d714.abn$otu,"Family"]))
new_da.d714.abn$Genus <- as.factor(unclass(tax_table(all.fil)[new_da.d714.abn$otu,"Genus"]))

#Volcano Plot
#classify for the significant p.adj
new_da.d714.abn$sig <- "NS"
# if pvalue < 0.05 and abs FC >1, set as "S" 
new_da.d714.abn$sig[new_da.d714.abn$p.adj <= 0.05 & abs(new_da.d714.abn$log2FC) > 0.6] <- "S"

new_da.d714.abn$Freq <- count.714.abn[order(match(as.character(count.714.abn$Var1),new_da.d714.abn$otu)),]$Freq

(volcano.d714.abn <- ggplot(new_da.d714.abn, aes(x = log2FC, y = -log10(p.adj), shape = method)) +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = max(-log10(p.adj)), ymax = Inf), fill = "grey95") + 
  geom_hline(aes(yintercept = max(-log10(p.adj))), 
             linetype = "dashed",
             color = "grey30") +
  geom_point(data = subset(new_da.d714.abn, sig == "NS"), 
             aes(fill = sig, size = Freq),
             stroke = NA) +
  geom_hline(yintercept=-log10(0.05), col="grey30", 
             linetype = 1) +
  geom_vline(xintercept = c(-0.6, 0.6), col="grey30", 
             linetype = 1) +
  geom_point(data = subset(new_da.d714.abn, sig == "S"), aes(group = Phylum, fill = Phylum, size = Freq)) +
    scale_size_continuous(range = c(4,7), breaks = c(2:4)) +
  scale_fill_manual(values = c("Actinobacteria" = "#EEDD88", "Bacteroidetes" = "#77AADD", "Firmicutes" = "#EE8866", "Proteobacteria" = "#FFAABB",  "NS" = "grey90")) +
  scale_shape_manual(values = c(21:24)) +
  geom_text_repel(data = subset(new_da.d714.abn, sig == "S" & Freq >= 3),
                  aes(label = Species), 
                  size = 5, 
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  box.padding = 1,
                  fontface = "italic") +
  geom_text(aes(x = 8, y = max(-log10(p.adj) + 0.07), label = "p.adj = 0"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  geom_text(aes(x = 8, y = -log10(0.05) + 0.07, label = "p.adj = 0.05"), 
             hjust = 1, vjust = 0, 
             color = "grey30",
             size = 5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
  scale_x_continuous(expand = c(0,0), limits = c(-8,8))  +
  theme_classic() +
    labs(title = "b",
         subtitle = "Day 14 VS Day 7. Antibiotic_no") +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
 theme(plot.title = element_text(size = 10, face = "bold", vjust = 0, hjust =0),
        plot.subtitle = element_text(size = 10, color = "grey50", vjust = 0, hjust = 0),
        #plot.margin = unit(c(0.5,0.5,0.5,1.5),"in"),
        axis.title.x = element_text(vjust = 0, hjust = 0.5, size = 10),
        axis.text.x = element_text( size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
       axis.line.y.left = element_blank()) +
  geom_segment(aes(x = -8, y = -Inf, xend = -8, yend = max(-log10(p.adj))), size = 1, linetype = "solid") +
  geom_segment(aes(x = -8, y = max(-log10(p.adj)), xend = -8, yend = Inf), size = 1, linetype = "dotted")
)

volcano.d714 <- ggarrange(volcano.d714.aby, volcano.d714.abn, align = "v", common.legend = T)

(da.plot <- ggarrange(volcano.d17.abn, 
                      volcano.d714.abn, 
                      volcano.d17.aby, 
                      volcano.d714.aby, 
                      nrow = 2, 
                      ncol = 2, 
                      align = "hv", 
                      common.legend = T))

ggplot2::ggsave(filename = "Differential Abundance.pdf", 
       plot = da.plot,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/Volcano/", 
       width = 15, 
       height = 15, 
       units = "in", 
       dpi = "retina", 
       limitsize = TRUE,
       bg = "white")

#UpsetR
upset.plot17.abn
upset.plot17.aby
upset.plot714.abn
upset.plot714.aby

#Data
new_da.d17.aby
new_da.d17.abn
new_da.d714.aby
new_da.d714.abn

write_xlsx(new_da.d17.aby, 
           "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/Volcano/da.d17.aby.xlsx")

write_xlsx(new_da.d17.abn, 
           "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/Volcano/da.d17.abn.xlsx")

write_xlsx(new_da.d714.aby, 
           "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/Volcano/da.d714.aby.xlsx")

write_xlsx(new_da.d714.abn, 
           "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/Volcano/da.d714.abn.xlsx")

```

# Different log2fc

```{r}
new_da.d17.abn
new_da.d17.aby

otu17.dup <- intersect(new_da.d17.abn$otu, new_da.d17.aby$otu)

fil.17.aby <- subset(new_da.d17.aby, otu %in% otu17.dup)
fil.17.abn <- subset(new_da.d17.abn, otu %in% otu17.dup)

fil.17.aby$atb <- rep("Yes", nrow(fil.17.aby))
fil.17.abn$atb <- rep("No", nrow(fil.17.abn))

dup.17 <- rbind(fil.17.abn, fil.17.aby)


dup.17$otu <- factor(dup.17$otu, levels = unique(dup.17$otu[order(dup.17$log2FC, decreasing = T)]))

dup.17.species <- dup.17$Species[c(1:nrow(dup.17))]

(dup.17.plot <- ggplot(dup.17, aes(x = otu, y = log2FC)) + 
    geom_hline(yintercept = 0, color = "grey75") +
    geom_line(size = 0.7, color = "black") +
    geom_point(aes(fill = atb), 
               shape = 21, 
               color = "black",
               size = 3,
               stroke = 0.7) +
    scale_x_discrete(labels = setNames(dup.17$Species, dup.17$otu)) +
    #scale_y_continuous(limits = c(-5,7), breaks = seq(-5, 7 ,2)) +
    scale_fill_manual(values = c("#f4f1de", "#e07a5f"), guide = guide_legend(title = "Antibiotic treatment")) +
    theme_bw() +
    labs(subtitle = "Week 1",
         x = "Species") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     face = "italic",
                                     color = "black"),
          plot.subtitle = element_text(color = "black",
                                       face = "bold"),
          axis.line = element_line(linewidth = 0.5),
          axis.title.x = element_blank()
          )
)

#---#

new_da.d714.abn
new_da.d714.aby

otu714.dup <- intersect(new_da.d714.abn$otu, new_da.d714.aby$otu)

fil.714.aby <- subset(new_da.d714.aby, otu %in% otu714.dup)
fil.714.abn <- subset(new_da.d714.abn, otu %in% otu714.dup)

fil.714.aby$atb <- rep("Yes", nrow(fil.714.aby))
fil.714.abn$atb <- rep("No", nrow(fil.714.abn))

dup.714 <- rbind(fil.714.abn, fil.714.aby)

dup.714$otu <- factor(dup.714$otu, levels = unique(dup.714$otu[order(dup.714$log2FC, decreasing = T)]))

dup.714.species <- dup.714$Species[c(1:nrow(dup.714))]

(dup.714.plot <- ggplot(dup.714, aes(x = otu, 
                                     y = log2FC)) + 
    geom_hline(yintercept = 0, color = "grey75") +
    geom_line(color = "black" ,size = 0.7) +
    geom_point(aes(fill = atb), 
               shape = 21, 
               color = "black",
               size = 3,
               stroke = 0.7) +
    scale_x_discrete(labels = setNames(dup.714$Species, dup.714$otu)) +
    #scale_y_continuous(limits = c(-5,7), breaks = seq(-5, 7 ,2)) +
    scale_fill_manual(values = c("#f4f1de", "#e07a5f"), guide = guide_legend(title = "Antibiotic treatment")) +
    theme_bw() +
    labs(subtitle = "Week 2",
         x = "Species") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     face = "italic",
                                     color = "black"),
          plot.subtitle = element_text(color = "black",
                                       face = "bold"),
          axis.line = element_line(linewidth = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
          )
)

intersect(dup.17$otu, dup.714$otu)

(dup.da.plot <- ggarrange(dup.17.plot, 
                          dup.714.plot, 
                          nrow = 1, 
                          widths = c(22,44), 
                          align = "hv",
                          common.legend = TRUE,
                          legend = "bottom")
)

ggplot2::ggsave(filename = "Dup_DA.pdf", 
       plot = dup.da.plot,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Differential Abundance/", 
       width = 8, 
       height = 4, 
       units = "in", 
       dpi = 500, 
       limitsize = TRUE,
       bg = "white")
```

# Heat map for DA
```{r}
#combine otu from all 4 DA subsets
chosen_otu<- unique(c(new_da.d17.abn$otu,
                      new_da.d17.aby$otu,
                      new_da.d714.abn$otu,
                      new_da.d714.aby$otu))

#top.genus
genus.da.ps <- tax_glom(data.04en.un, "Genus")
genus.da.abund <- transform_sample_counts(genus.da.ps, function(OTU) OTU/sum(OTU))

#da.fil.ps <- prune_taxa(chosen_otu, genus.da.abund)

top.da.genus <- top_taxa(genus.da.abund, 20)
da.fil.ps <- prune_taxa(top.da.genus, genus.da.abund)

taxa.da.order <- names(sort(taxa_sums(da.fil.ps)))
sample.da.order <- rownames(sample_data(da.fil.ps)[order(sample_data(da.fil.ps)$day, sample_data(da.fil.ps)$antibiotic_trt)])

#phylum_cols <- c("Actinobacteria" = "#EEDD88", "Bacteroidetes" = "#77AADD", "Firmicutes" = "#EE8866", "Proteobacteria" = "#FFAABB",  "Fusobacteria" = "grey90")

(hm <- plot_heatmap(da.fil.ps, 
                    taxa.label="Genus", 
                    sample.order=sample.da.order,
                    taxa.order=taxa.da.order) +
    geom_tile(color = "grey10", size = 0.05) +
    scale_fill_gradientn(colors = c("grey10", "#457B9D", "#A8DADC", "#E63946"), 
                         values=c(0, 0.1, 0.5, 1),
                         na.value = NA, 
                         guide="colourbar", 
                         name="Relative\nproportion") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", size = 15, color = "grey10"),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"in")
        )
)

hcbdf <- data.frame(sample_data(da.fil.ps))[sample.da.order,]
hcbdf$index <- seq(1,phyloseq::nsamples(da.fil.ps))
hcb <- make_hcb(hcbdf, "day", name="Day", fillScale = scale_fill_manual(values=c("#ef476f", "#ffd166", "#26547c")))
hcb <- hcb + annotate("text", x=tapply(hcbdf$index, hcbdf[,"day",drop=T], mean), y=1, label=levels(hcbdf[,"day",drop=T]), size=4)

#Age component
hcbage <- make_hcb(hcbdf, "age_month", name="age_month", fillScale=scale_fill_gradientn(colours=c("#fff8ab",  "#b6c670", "#6c973c", "#0d690e"), values=rescale(c(0, 16,  25, 60)), na.value="black"))

#Day component
hcbday <- make_hcb(hcbdf, "day", name="Day",fillScale = scale_fill_manual(values=c("#ef476f", "#ffd166", "#26547c")))


#Antibiotic treatment
hcbdf$antibiotic_trt <- as.factor(hcbdf$antibiotic_trt)
hcbabt <- make_hcb(hcbdf, "antibiotic_trt", name="antibiotic_trt",fillScale = scale_fill_manual(values= c("no" = "#f4f1de", "yes" =  "#e07a5f", "NA" = "grey75")))

big.hm <- mush(hm, list(hcb, hcbabt))
big.heatmap <- plot(ggarrange(big.hm))

ggplot2::ggsave(filename = "Heatmap (unorm).pdf", 
       plot = big.heatmap,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Heatmap", 
       width = 16, 
       height = 8, 
       units = "in", 
       dpi = 1000, 
       limitsize = TRUE,
       bg = "white")

fil.ps <- prune_taxa(chosen_otu, data.04en)

mean_abund <- microbiomeutilities::phy_to_ldf(data.04en, 
                                              transform.counts = "compositional")

mean_abund <- mean_abund %>% 
  group_by(antibiotic_trt, OTUID, day) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()


#Bacteroides genus
dt04.genus <- aggregate_taxa(data.04en, "Genus")

mean_abund.genus <- microbiomeutilities::phy_to_ldf(dt04.genus, 
                                              transform.counts = "compositional")

mean_abund.genus <- mean_abund.genus %>% 
  group_by(antibiotic_trt, Family, Genus, day) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()

bac.fam <- mean_abund.genus[which(mean_abund.genus$Family == "Bacteroidaceae" & mean_abund.genus$Genus != "Bacteroidaceae_unclassified"),] %>% as.data.frame()

ggplot(bac.fam, 
       aes(x = day,
           y = mean_abundance,
           group = antibiotic_trt)) +
  geom_line() +
  facet_wrap(.~Genus)

###
mean.17 <- filter(mean_abund, OTUID %in% dup.17$otu & day == 1)
mean.17 <- as.data.frame(mean.17)

dup.17 <- dup.17[order(dup.17$otu),]
mean.17 <- mean.17[order(mean.17$OTUID),]

dup.17$mean_abund.7 <- mean.17$mean_abundance

dup.17$otu <- factor(dup.17$otu, levels = unique(dup.17$otu[order(dup.17$log2FC, decreasing = T)]))

dup.17.species <- dup.17$Species[c(1:nrow(dup.17))]

(dup.17.plot <- ggplot(dup.17, aes(x = otu, y = log2FC)) + 
    geom_hline(yintercept = 0, color = "grey75") +
    geom_line(linewidth = 0.7, color = "black") +
    geom_point(aes(fill = atb, size = log2(mean_abund.7)), 
               shape = 21, 
               color = "black",
               stroke = 0.7) +
    scale_x_discrete(labels = setNames(dup.17$Species, dup.17$otu)) +
    #scale_y_continuous(limits = c(-5,7), breaks = seq(-5, 7 ,2)) +
    scale_fill_manual(values = c("#e07a5f", "#f4f1de"), guide = guide_legend(title = "Antibiotic treatment")) +
    theme_bw() +
    labs(subtitle = "Week 1",
         x = "Species") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     face = "italic",
                                     color = "black"),
          plot.subtitle = element_text(color = "black",
                                       face = "bold"),
          axis.line = element_line(linewidth = 0.5),
          axis.title.x = element_blank()
          )
)
```

# line plot for specfic otus

```{r}
fil.otu <- unique(dup.17$otu,
                  dup.714$otu) %>% as.character()
fil.otu <- append(fil.otu[fil.otu != "Otu00225"], c("Otu00050", "Otu00006", "Otu00004", "Otu00002", "Otu00012", "Otu00065"))
fil.abund <- mean_abund[which(mean_abund$OTUID %in% fil.otu),] %>% 
  as.data.frame()

fil.abund$Species <- unclass(phyloseq::tax_table(all.fil)[fil.abund$OTUID,'Species']) %>% as.factor()


(line.abund <- ggplot(fil.abund, aes(x = day,
                                     y = mean_abundance,
                                     group = antibiotic_trt,
                                     color = antibiotic_trt)) + 
    geom_errorbar(aes(ymin=mean_abundance-se, 
                      ymax=mean_abundance+se), 
                  width=0.1, 
                  alpha = 0.75) +
    geom_line(size = 1.6) +
    geom_point(size = 1) +
    facet_wrap(.~Species, 
               scales = "free_y") +
    scale_color_manual(values = c("no" = "grey50", "yes" = "#e07a5f"),
                       guide = guide_legend(title = "Antibiotic \ntreatment")) +
    theme_light() +
    guides(x.sec = "axis",
         y.sec = "axis") +
    labs(x = "Day",
         y = "Mean relative abundance") +
    theme(plot.title = element_text(size = 15, 
                                    face = "bold", 
                                    vjust = 3, 
                                    hjust = 0),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(10,0.5,0.5,10),
          axis.line = element_line(linewidth = .5),
          axis.text.x.bottom = element_text(vjust = -1, 
                                            size = 10),
          axis.text.y.left =element_text(size = 10),
          axis.text.x.top = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.y = element_text(color = "grey30", 
                                      face = "bold", 
                                      size = 10, 
                                      angle = 90, 
                                      vjust = 3),
          strip.background = element_rect(fill = NA, 
                                          color = "grey30", 
                                          linewidth = 1.5),
          strip.text = element_text(face = "bold.italic", color = "black"))
  )

ggplot2::ggsave(filename = "rel_abund.pdf", 
       plot = line.abund,
       device = "pdf", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Relative abundance", 
       width = 10, 
       height = 10, 
       units = "in", 
       dpi = 800, 
       limitsize = TRUE,
       bg = "white")


#C. difficile

c.diff <- mean_abund[which(mean_abund$OTUID == "Otu00015"),]

```

---

# Correlation network:

> CClasso function

```{r}
##############################################################################################
#-----------------------------------------------------------------------------------------------------
## Taken from CCLASSO methodology
# The CCLasso function from github
################################################################################
# File: cclasso.R
# Aim : Correlation inference for compositional data through lasso
#-------------------------------------------------------------------------------
# Author : Fang Huaying (Peking University)
# Email  : hyfang@pku.edu.cn
# Date   : 2016-01-08
# Version: 2.0
#-------------------------------------------------------------------------------
# Main function: cclasso(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
#                        lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) 
#
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#                    n samples & p compositional variables
#      counts ------ Is the compositional data matrix a count matrix? 
#                    Default: FALSE
#      pseudo ------ pseudo count if counts = TRUE
#                    Default: 0.5
#        k_cv ------ folds of cross validation
#                    Default: 3     
#     lam_int ------ tuning parameter interval
#                    Default: [1e-4, 1]
#       k_max ------ maximum iterations for golden section method
#                    Default: 20
#      n_boot ------ Bootstrap times
#                    Default: 20
#  Output: 
#      A list structure contains:
#       var_w ------ variance estimation
#       cor_w ------ correlation estimation
#      p_vals ------ p-values for elements of cor_w equal 0 or not
#      lambda ------ final tuning parameter
#     info_cv ------ information for cross validation
#-------------------------------------------------------------------------------
#CClasso function
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {
  n <- nrow(x);
  p <- ncol(x);

  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  x <- log(x);
  vx2 <- stats::var(x);

  # Diagonal weight for loss function
  rmean_vx2 <- rowMeans(vx2);
  wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
  wd2 <- sqrt(wd);
  
  # Some global parameters for optimization with single lambda
  rho <- 1;
  u_f <- eigen(diag(p) - 1/p)$vectors;
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
  wd_u_eig <- eigen(wd_u);
  d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) / 
    (2 * rho) + 1 );
  u0_wd <- wd_u_eig$vectors;

  # Golden section method for the selection of lambda (log10 scale)
  sigma <- vx2;
  lam_int2 <- log10(range(lam_int));
  a1 <- lam_int2[1]; 
  b1 <- lam_int2[2];
  # Store lambda and corresponding cross validation's loss
  lams <- NULL; 
  fvals <- NULL;
  # Two trial points in first 
  a2 <- a1 + 0.382 * (b1 - a1); 
  b2 <- a1 + 0.618 * (b1 - a1);
  fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
    sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
  lams <- c(lams, b2); 
  fvals <- c(fvals, fb2$cv_loss);
  fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
    sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
  lams <- c(lams, a2); 
  fvals <- c(fvals, fa2$cv_loss);
  # Error tolerance for convergence
  err_lam2 <- 1e-1 * max(1, lam_int2);
  err_fval <- 1e-4;
    
  err <- b1 - a1;
  k <- 0;
  while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);

    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;      
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
        sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);

      lams <- c(lams, b2); 
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
        sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);

      lams <- c(lams, a2); 
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);      

    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
  }
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2, 
    lam_int = 10^c(a1, b1)); 
  if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
      10^b1, "])\n", sep = "");
  }

  lambda <- 10^((a2 + b2)/2);
  # Bootstrap for cclasso
  lambda2 <- lambda / rho;
  info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2, 
    n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);

  return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w, 
    p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}

## Bootstrap for cclasso
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20, 
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x);
  p <- ncol(x);

  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);

  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
    ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]), 
      lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
        
    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);

  #----------------------------------------  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  #----------------------------------------
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art), 
    lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * 
    (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);
  # Modification of estimation for cclasso
  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );    
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) + 
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;
  #----------------------------------------

  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}

# cross validation's loss of cclasso for single lambda
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma, 
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x);
  p <- ncol(x);

  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);

    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2, 
      wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }

  return(list(cv_loss = cv_loss, sigma = sigma));
}
#-------------------------------------------------------------------------------
# cclasso for single lambda
cclasso_sub <- function(sigma, vx, lambda2, 
                        wd, u_f, u0_wd, d0_wd, 
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;

  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) * 
      d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);
    # Update sigma2
    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) * 
      (A + lambda2);
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);

    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1), 
      abs(sigma2_new - sigma2)/(abs(sigma2) + 1) ); 
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }

  if(k >= k_max) {
    cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max, 
      "&& Relative error:", err, "!\n");
  }
  
  return(sigma);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
plot_cclasso<- function(cor.map, adj.matrix, phyloseq){
      # Loading library
      require(igraph)
      require(plyr)
      require(RColorBrewer)
      require(cols4all)
      # create a graph object from the adjacency matrix
      g <- graph.adjacency(adj.matrix) 
      
      # get the edges as well as the original 
      edgelist <- get.edgelist(g)
      
      edge.cor <- apply(edgelist, MARGIN = 1, FUN = function(x) {cor.map[x[1],x[2]]})
      edgelist <- data.frame(edgelist, edge.cor)
      colnames(edgelist) <- c("from","to", "cor")
      
      # filter edges that link the same node
      filt <- edgelist$from == edgelist$to 
      print("There are some edges that link the same node... now deleting...")
      edgelist <- edgelist[!filt,]
      
      edgelist[,1:2] <- t(apply(edgelist, MARGIN = 1, FUN = function(x){sortEdges(x[1:2])}))
      dupl <- duplicated(paste(edgelist$from, edgelist$to)); sum(dupl)
      edgelist <- edgelist[!dupl,]
      # Take only one direction

      # get the annotation for the nodes that are implicated in the network.
      all.nodes <- unique(c(as.character(edgelist$from), as.character(edgelist$to)))
      # Get the taxonomic classification and abundance of OTUs from original data.
      feat.annot <- taxa_names(phyloseq)[as.numeric(all.nodes)] #OTU names
      feat.annot <- as.data.frame(feat.annot)
      rownames(feat.annot) <- feat.annot$feat.annot
      feat.annot$name <- all.nodes
      feat.annot <- feat.annot[,c(2,1)]
      feat.annot$size <- as.vector(taxa_sums(phyloseq)[as.numeric(all.nodes)]) # abundance
      feat.annot$Family <- as.vector(tax_table(phyloseq)[as.numeric(all.nodes),"Family"]) #Assign taxonomy of the Family level
      print("Finished first part")

      # create the igraph object
      gD <- graph.data.frame(d = edgelist, directed = FALSE, vertices = feat.annot)
      degAll <- igraph::degree(gD, v = V(gD), mode = "all")
      V(gD)$name.origin <- tax_table(phyloseq)[as.numeric(all.nodes), "Species"]
      V(gD)$name <- tax_table(phyloseq)[as.numeric(all.nodes), "Species"]
      V(gD)$otu <- rownames(tax_table(phyloseq)[as.numeric(all.nodes),])
      gD <- set.vertex.attribute(gD, "degree", index = V(gD)$name, value = degAll)
      #gD <- set.vertex.attribute(gD, "degree", index = V(gD)$otu, value = degAll)
      
      # Calculate betweenness for all nodes
      betAll <- igraph::betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
      betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll)); rm(betAll)
      
      # Add new node/edge attributes based on the calculated node properties/similarities
      gD <- set.vertex.attribute(gD, "betweenness", index = V(gD)$name, value = betAll.norm)

      ##Calculate edge properties and add to the network
      #Calculate Dice similarities between all pairs of nodes
      dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
      # The following function will transform a square matrix to an edge driven one and add 
      # values to each edge
      F1 <- function(x) {data.frame(dice = dsAll[which(V(gD)$name == as.character(x$from)), 
                                           which(V(gD)$name == as.character(x$to))])}
      edgelist.ext <- ddply(edgelist, .variables=c("from", "to"), 
                      function(x) data.frame(F1(x))); dim(edgelist.ext)
      gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
      E(gD)[as.character(edgelist.ext$from) %--% as.character(edgelist.ext$to)]$similarity <- as.numeric(edgelist.ext$dice)
      
      # Edge color based on type of correlations
      E(gD)$color <- ifelse(E(gD)$cor < 0, "#EE6677", "#4477AA") # Red for Negative, Blue for Positive
      # Vertex color based on the Order of the OTU
      V(gD)$color <- c(c4a("classic20", 19), c4a("okabe", 7), c4a("rainbow", 23))[as.factor(V(gD)$Family)]

      # resume the network
      print(paste("There are",vcount(gD),"nodes and",ecount(gD),"edges"))
      summary(gD)

      # compute the layout
      l <- layout.auto(gD)
      par(oma=c(2,0.2,0.2,0.2))

      pl <- plot(gD,
            vertex.label = V(gD)$name, 
            vertex.label.cex=1.2,
            #vertex.size=(log10(V(gD)$size))^1.5/2,
            vertex.size= 5,
            vertex.label.dist=0.75, vertex.label.color='black',
            vertex.label.font=3,
            edge.arrow.size=2, asp=TRUE, 
            #edge.width = abs(E(gD)$cor)*7,
            edge.width = 4,
            rescale=TRUE, layout=l*1.2
           )
      par(mfrow=c(1,1), fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
      ## legend 
      legend(1.5, 1, legend = levels(factor(tax_table(phyloseq)[as.numeric(all.nodes),"Family"])), 
       fill = c(c4a("classic20", 19), c4a("okabe", 7), c4a("rainbow", 23)), xpd=TRUE,
      bty='n', cex= 1, xjust = 0, ncol=1, text.width = 0.5 ,x.intersp = 0.1, title = "Family")
      return(pl)
}

## To beautify the cclasso plot 
sortEdges <- function(var.name = c(node1, node2)){
      if(length(var.name)!=2) stop("Two variable names are needed")
      if(var.name[2]<var.name[1])
      {
            b <- var.name[1]
            var.name[1] <- var.name[2]
            var.name[2] <- b
      }
      return(var.name)
}
#----------------------------------------------------------#
```

## D17:

### antibiotic_yes
```{r}
aby.17.keep

#NetCoMi
net_17aby <- netConstruct(data = aby.17.keep,
                        measure = "cclasso",
                        measurePar = list(counts = TRUE, K = 3, kmax = 5000),
                        normMethod = "none",
                        zeroMethod = "multRepl",
                        sparsMethod = "none",
                        dissFunc = "signed",
                        verbose = 3)

props.net_17aby <- netAnalyze(net_17aby,
                             centrLCC = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             normDeg = FALSE, 
                             weightDeg = FALSE)

#summary
summary(props.net_17aby)

tax_table(all.fil)[names(sort(props.net_17aby$centralities$eigenv1, decreasing = TRUE)[1:10]),]
tax_table(all.fil)[props.net_17aby$hubs$hubs1,]

#plot network
net_17.aby.cor.map <- net_17aby$assoEst1
colnames(net_17.aby.cor.map) <- NULL
rownames(net_17.aby.cor.map) <- NULL
net_17.aby.cor.map <- as.matrix(net_17.aby.cor.map)

diag(net_17.aby.cor.map) <- 0 ## all self-interactions are set to 0

#Exploration
plot(density(as.numeric(net_17.aby.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(net_17.aby.cor.map[abs(net_17.aby.cor.map)>0.2])), main="filtered")


net_17aby.adj.mat <- (net_17.aby.cor.map >=0.2| net_17.aby.cor.map <= -0.2) + 0 

length(which(net_17aby.adj.mat == 1))

plot_cclasso(net_17.aby.cor.map, net_17aby.adj.mat, aby.17.keep)

ggplot2::ggsave(filename = "17aby_net.png", 
       plot = plot_cclasso(net_17.aby.cor.map, net_17aby.adj.mat, aby.17.keep),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")
```

### antibiotic_no
```{r}
abn.17.keep

#NetCoMi
net_17abn <- netConstruct(data = abn.17.keep,
                          measure = "cclasso",
                          measurePar = list(counts = TRUE, K = 3, kmax = 5000),
                          normMethod = "none",
                          zeroMethod = "multRepl",
                          sparsMethod = "none",
                          dissFunc = "signed",
                          verbose = 3,
                          cores = 20L)

props.net_17abn <- netAnalyze(net_17abn,
                             centrLCC = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             normDeg = FALSE, 
                             weightDeg = FALSE)

#summary
summary(props.net_17abn)

tax_table(all.fil)[names(sort(props.net_17abn$centralities$eigenv1, decreasing = TRUE)[1:10]),]
tax_table(all.fil)[props.net_17abn$hubs$hubs1,]

#plot network
net_17.abn.cor.map <- net_17abn$assoEst1
colnames(net_17.abn.cor.map) <- NULL
rownames(net_17.abn.cor.map) <- NULL
net_17.abn.cor.map <- as.matrix(net_17.abn.cor.map)

diag(net_17.abn.cor.map) <- 0 ## all self-interactions are set to 0

#Exploration
plot(density(as.numeric(net_17.abn.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(net_17.abn.cor.map[abs(net_17.abn.cor.map)>0.07])), main="filtered")


net_17abn.adj.mat <- (net_17.abn.cor.map >=0.15| net_17.abn.cor.map <= -0.15) + 0 

length(which(net_17abn.adj.mat == 1))

plot_cclasso(net_17.abn.cor.map, net_17abn.adj.mat, abn.17.keep)

ggplot2::ggsave(filename = "17abn_net.png", 
       plot = plot_cclasso(net_17.abn.cor.map, net_17abn.adj.mat, abn.17.keep),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")
```

### Differential network

```{r}
net_17_cclasso <- netConstruct(data = abn.17.keep, 
                               data2 = aby.17.keep, 
                               measurePar = list(counts = TRUE, 
                                                 K = 3, 
                                                 kmax = 5000),
                               measure = "cclasso", 
                               zeroMethod = "pseudo",
                               normMethod = "none",
                               sparsMethod = "none", 
                               dissFunc = "signed",
                               verbose = 3,
                               cores = 20L)

# Differential network construction
diff_17 <- diffnet(net_17_cclasso,
                   diffMethod = "permute", 
                   adjust = "lfdr", 
                   cores = 20L,
                   )

# Differential network plot
plot(diff_17, 
     cexNodes = 0.8, 
     cexLegend = 3,
     cexTitle = 4,
     mar = c(2,2,8,5),
     legendGroupnames = c("group 'no'", "group 'yes'"),
     legendPos = c(0.7,1.6))

props_17_cclasso <- netAnalyze(net_17_cclasso, 
                               clustMethod = "cluster_fast_greedy",
                               weightDeg = TRUE,
                               normDeg = FALSE,
                               gcmHeat = FALSE)
```

## D714:

### antibiotic_yes:

```{r}
aby.714.keep

#NetCoMi
net_714aby <- netConstruct(data = aby.714.keep,
                        measure = "cclasso",
                        measurePar = list(counts = TRUE, K = 3, kmax = 5000),
                        normMethod = "none",
                        zeroMethod = "multRepl",
                        sparsMethod = "none",
                        dissFunc = "signed",
                        verbose = 3)

props.net_714aby <- netAnalyze(net_714aby,
                             centrLCC = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             normDeg = FALSE, 
                             weightDeg = FALSE)

#summary
summary(props.net_714aby)

tax_table(all.fil)[names(sort(props.net_714aby$centralities$eigenv1, decreasing = TRUE)[1:10]),]
tax_table(all.fil)[props.net_714aby$hubs$hubs1,]

#plot network
net_714.aby.cor.map <- net_714aby$assoEst1
colnames(net_714.aby.cor.map) <- NULL
rownames(net_714.aby.cor.map) <- NULL
net_714.aby.cor.map <- as.matrix(net_714.aby.cor.map)

diag(net_714.aby.cor.map) <- 0 ## all self-interactions are set to 0

#Exploration
plot(density(as.numeric(net_714.aby.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(net_714.aby.cor.map[abs(net_714.aby.cor.map)>0.3])), main="filtered")


net_714aby.adj.mat <- (net_714.aby.cor.map >=0.38| net_714.aby.cor.map <= -0.25) + 0 

length(which(net_714aby.adj.mat == 1))

plot_cclasso(net_714.aby.cor.map, net_714aby.adj.mat, aby.714.keep)

ggplot2::ggsave(filename = "714aby_net.png", 
       plot = plot_cclasso(net_714.aby.cor.map, net_714aby.adj.mat, aby.714.keep),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")
```

### antibiotic_no:

```{r}
abn.714.keep

#NetCoMi
net_714abn <- netConstruct(data = abn.714.keep,
                        measure = "cclasso",
                        measurePar = list(counts = TRUE, K = 3, kmax = 5000),
                        normMethod = "none",
                        zeroMethod = "multRepl",
                        sparsMethod = "none",
                        dissFunc = "signed",
                        verbose = 3)

props.net_714abn <- netAnalyze(net_714abn,
                             centrLCC = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             normDeg = FALSE, 
                             weightDeg = FALSE)

#summary
summary(props.net_714abn)

tax_table(all.fil)[names(sort(props.net_714abn$centralities$eigenv1, decreasing = TRUE)[1:10]),]
tax_table(all.fil)[props.net_714abn$hubs$hubs1,]

#plot network
net_714.abn.cor.map <- net_714abn$assoEst1
colnames(net_714.abn.cor.map) <- NULL
rownames(net_714.abn.cor.map) <- NULL
net_714.abn.cor.map <- as.matrix(net_714.abn.cor.map)

diag(net_714.abn.cor.map) <- 0 ## all self-interactions are set to 0

#Exploration
plot(density(as.numeric(net_714.abn.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(net_714.abn.cor.map[abs(net_714.abn.cor.map)>0.1])), main="filtered")


net_714abn.adj.mat <- (net_714.abn.cor.map >=0.08| net_714.abn.cor.map <= -0.08) + 0 

length(which(net_714abn.adj.mat == 1))

plot_cclasso(net_714.abn.cor.map, net_714abn.adj.mat, abn.714.keep)

ggplot2::ggsave(filename = "714abn_net.png", 
       plot = plot_cclasso(net_714.abn.cor.map, net_714abn.adj.mat, abn.714.keep),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")
```

## Week 1:

```{r}
patient.17.un

#NetCoMi
net_17 <- netConstruct(data = patient.17.un,
                       measure = "cclasso",
                       measurePar = list(counts = TRUE, K = 3, kmax = 5000),
                       normMethod = "none",
                       zeroMethod = "multRepl",
                       filtTax = c("numbSamp", "totalReads"),
                       filtTaxPar = list(numbSamp = 5, totalReads = 20),
                       sparsMethod = "none",
                       dissFunc = "signed",
                       verbose = 3, 
                       cores = 20L)

props.net_17 <- netAnalyze(net_17,
                           centrLCC = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           normDeg = FALSE, 
                           weightDeg = FALSE)

#summary
summary(props.net_17)

tax_table(all.fil)[names(sort(props.net_17$centralities$eigenv1, decreasing = TRUE)[1:10]),]
tax_table(all.fil)[props.net_17$hubs$hubs1,]

#plot network
net_17.cor.map <- net_17$assoEst1
colnames(net_17.cor.map) <- NULL
rownames(net_17.cor.map) <- NULL
net_17.cor.map <- as.matrix(net_17.cor.map)

diag(net_17.cor.map) <- 0 ## all self-interactions are set to 0

#Exploration
plot(density(as.numeric(net_17.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(net_17.cor.map[abs(net_17.cor.map)>0.1])), main="filtered")


net_17.adj.mat <- (net_17.cor.map >=0.05| net_17.cor.map <= -0.3) + 0 

length(which(net_17.adj.mat == 1))

plot_cclasso(net_17.cor.map, net_17.adj.mat, patient.17.un)

ggplot2::ggsave(filename = "net_17.png", 
       plot = plot_cclasso(net_17.cor.map, net_17.adj.mat, patient.17.un),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")



### Apply original cclasso function for week1
w1.cclasso <- cclasso(t(otu_table(patient.17.un.keep)), counts=TRUE, k_cv = 3, n_boot = 200)

w1.cor.map <- w1.cclasso$cor_w
diag(w1.cor.map) <- 0 ## all self-interactions are set to 0
w1.cor.pvalue <- w1.cclasso$p_vals 
w1.cor.map[which(w1.cor.pvalue <= 0.01)] 
#Exploration
plot(density(as.numeric(w1.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(w1.cor.map[abs(w1.cor.map)>0.3])), main="filtered")

w1.adj.mat <- ((w1.cor.map >=0.27 | w1.cor.map <= -0.27) & (w1.cor.pvalue <=0.01)) + 0 
length(which(w1.adj.mat == 1))

plot_cclasso(w1.cor.map, w1.adj.mat, patient.17.un.keep)

ggplot2::ggsave(filename = "w1_cclasso.png", 
       plot = plot_cclasso(w1.cor.map, w1.adj.mat, patient.17.un.keep),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 25, 
       height = 14, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")
```

## Week 2:

```{r}
patient.714.un

#NetCoMi
net_714 <- netConstruct(data = patient.714.un,
                       measure = "cclasso",
                       measurePar = list(counts = TRUE, K = 3, kmax = 5000),
                       normMethod = "none",
                       zeroMethod = "multRepl",
                       filtTax = c("numbSamp", "totalReads"),
                       filtTaxPar = list(numbSamp = 5, totalReads = 20),
                       sparsMethod = "none",
                       dissFunc = "signed",
                       verbose = 3, 
                       cores = 20L)

props.net_714 <- netAnalyze(net_714,
                           centrLCC = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           normDeg = FALSE, 
                           weightDeg = FALSE)

#summary
summary(props.net_714)

tax_table(all.fil)[names(sort(props.net_714$centralities$eigenv1, decreasing = TRUE)[1:10]),]
tax_table(all.fil)[props.net_714$hubs$hubs1,]

#plot network
net_714.cor.map <- net_714$assoEst1
colnames(net_714.cor.map) <- NULL
rownames(net_714.cor.map) <- NULL
net_714.cor.map <- as.matrix(net_714.cor.map)

diag(net_714.cor.map) <- 0 ## all self-interactions are set to 0

#Exploration
plot(density(as.numeric(net_714.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(net_714.cor.map[abs(net_714.cor.map)>0.1])), main="filtered")


net_714.adj.mat <- (net_714.cor.map >= 0.1| net_714.cor.map <= -0.1) + 0 

length(which(net_714.adj.mat == 1))

plot_cclasso(net_714.cor.map, net_714.adj.mat, patient.714.un)

ggplot2::ggsave(filename = "net_714.png", 
       plot = plot_cclasso(net_714.cor.map, net_714.adj.mat, patient.714.un),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")


### Apply original cclasso function for week1
w2.cclasso <- cclasso(t(otu_table(patient.714.un.keep)), counts=TRUE, k_cv = 3, n_boot = 200)

w2.cor.map <- w2.cclasso$cor_w
diag(w2.cor.map) <- 0 ## all self-interactions are set to 0
w2.cor.pvalue <- w2.cclasso$p_vals 
w2.cor.map[which(w2.cor.pvalue <= 0.01)] 
#Exploration
plot(density(as.numeric(w2.cor.map)), main="All interactions") ## most interaction has strength centering zero
plot(density(as.numeric(w2.cor.map[abs(w2.cor.map)>0.2])), main="filtered")

w2.adj.mat <- ((w2.cor.map >=0.1| w2.cor.map <= 1) & (w2.cor.pvalue <=0.01)) + 0 
length(which(w2.adj.mat == 1))

plot_cclasso(w2.cor.map, w2.adj.mat, patient.714.un.keep)

ggplot2::ggsave(filename = "w2_cclasso.png", 
       plot = plot_cclasso(w2.cor.map, w2.adj.mat, patient.714.un.keep),
       device = "png", 
       path = "D:/Study/Oucru/SonNam_project/diarrhea_dataset/diarrhea_dataset/updated_images/Final/Network", 
       width = 20, 
       height = 12, 
       units = "in", 
       dpi = 400, 
       limitsize = TRUE,
       bg = "white")
```

# Combine plot:

```{r}
plot_cclasso(net_17.aby.cor.map, net_17aby.adj.mat, aby.17.keep)
plot_cclasso(net_17.abn.cor.map, net_17abn.adj.mat, abn.17.keep)
plot_cclasso(net_714.aby.cor.map, net_714aby.adj.mat, aby.714.keep)
plot_cclasso(net_714.abn.cor.map, net_714abn.adj.mat, abn.714.keep)
```
