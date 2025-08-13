library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(patchwork)
library(vegan)
library(dplyr)


path <- "~/Documents/seqs/carpobrotus_16s/Sequences_noadapt"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(basename(fnFs), function(x) sub("\\..*$", "", x))

sample.names

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[11:12])


filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, 
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) #ho tolto out perchè si è saltata la fase "Filter and Trim"
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/Users/danielacosta/Documents/seqs/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)


species_fasta <- '/Users/danielacosta/Documents/seqs/tax/silva_species_assignment_v138.1.fa.gz'

chunk.size <- 4000
chunks <- split(c(1:nrow(taxa)),
                sort(c(1:nrow(taxa))%%ceiling(nrow(taxa)/chunk.size)))

chunks.species <- lapply(chunks,
                         function(x){
                           return(addSpecies(taxa[x,],
                                             refFasta = species_fasta, verbose = TRUE))
                         })
taxa.species <- do.call(rbind, chunks.species)



save.image("~/Documents/seqs/carpobrotus_16s/Sequences_noadapt/dada2_16s_Giglio.RData")



taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
dim(taxa.print)


#Creo i raggruppamenti "Control", "Invaded", "Eradicated"
sample_df <- data.frame(
  Sample = c("C1","C11","C12","C13","C15","C16","C18","C19","C21",
             "C22","C24","C25","C27","C28","C3","C30","C31","C32",
             "C33","C34","C35","C37","C39","C4","C40","C41","C43",
             "C44","C46","C47","C49","C50","C52","C6","C7","C9"),
  status = c("Native","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated",
             "Native","Native","Native","Native","Native","Native","Native","Native","Native",
             "Eradicated","Eradicated","Invaded","Invaded","Invaded","Native","Invaded","Invaded","Invaded",
             "Invaded","Invaded","Invaded","Invaded","Invaded","Invaded","Native","Eradicated","Eradicated"),
  location = c("G", "G", "G","G","G","G","G","G","G","G","G","G",
               "G","G","G","G","G","G","G","G","G","G","G","G",
               "G","G","G","G","G","G","G","G","G","G","G","G") ,
  code = c("GN", "GE", "GE", "GE", "GE", "GE", "GE", "GE", "GE",
           "GN", "GN", "GN", "GN", "GN", "GN", "GN", "GN", "GN",
           "GE", "GE", "GC", "GC", "GC", "GN", "GC","GC","GC",
           "GC","GC","GC","GC","GC","GC","GN","GE","GE"),
  name =  c("T130_new","T129_new", "T129_new", "T131M", "T131M", "T136_new", "T136_new", "T135", "T135",
            "T147_new", "T147_new", "T134", "T134", "T149", "T130_new", "T149", "T150", "T150",
            "T132_new", "T132_new", "82", "82", "80", "T133_new", "80","81_new","81_new",
            "INV_1_new","INV_1_new","INV_2_new","INV_2_new","84_new","84_new","T133_new","T148_new","T148_new"),
  year  = c("2023", "2023", "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023"))

rownames(sample_df) <- c("C1","C11","C12","C13","C15","C16","C18","C19","C21","C22","C24","C25","C27","C28","C3","C30","C31","C32","C33","C34","C35","C37","C39","C4","C40","C41","C43","C44","C46","C47","C49","C50","C52","C6","C7","C9")

rownames(seqtab.nochim) <- c("C1","C11","C12","C13","C15","C16","C18","C19","C21","C22","C24","C25","C27","C28","C3","C30","C31","C32","C33","C34","C35","C37","C39","C4","C40","C41","C43","C44","C46","C47","C49","C50","C52","C6","C7","C9")


#BUilding a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa.species))


#rename asvs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))


###Build tree
library(phangorn)
library(DECIPHER)

#Get the DNA sequences
seqs <- refseq(ps)  # or DNAStringSet object with your sequences

# Align sequences using DECIPHER
alignment <- DECIPHER::AlignSeqs(seqs, anchor = NA)

# Build distance matrix
?phangorn::dist.ml
dists <- phangorn::dist.ml(as.phyDat(alignment))

# Build neighbor-joining tree
treeNJ <- phangorn::NJ(dists)
tree_rooted <- midpoint(treeNJ)

ps <- merge_phyloseq(ps, tree_rooted)



##Rarefaction curve
otu_mat <- as.matrix(otu_table(ps))
otu_mat = t(otu_mat)

#  Convert matrix columns to abundance vectors
otu_list <- apply(otu_mat, 2, function(x) {
  x <- x[x > 0]                     # remove 0s
  if (length(x) >= 2 && sum(x) > 0) return(x) else return(NULL)  # only keep valid samples
})

#  Remove NULLs
otu_list <- otu_list[!sapply(otu_list, is.null)]

# Run iNEXT
inext_out <- iNEXT(otu_list, q = 0, datatype = "abundance")
ggiNEXT(inext_out, type = 1) + 
  theme_minimal() +
  labs(title = "Rarefaction Curves", x = "Sample Size", y = "Species Richness")



###Tree like figure

ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus_comp <- ps_genus %>% microbiome::transform("compositional")

ps_melted <- psmelt(ps_genus_comp)

# Calculate total abundance per genus
top50_genus <- ps_melted %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 50) %>%
  pull(Genus)

genus_abund <- ps_melted %>%
  filter(Genus %in% top50_genus) %>%
  group_by(OTU, Genus, status) %>%  
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Get OTUs associated with top 50 genera
keep_otus <- unique(genus_abund$OTU)
ps_top50 <- prune_taxa(keep_otus, ps_genus)

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(tidyverse)

#Reshape tree

tree_circ <- phy_tree(ps_top50)
taxa_top50 <- as.data.frame(tax_table(ps_top50)) %>%
  rownames_to_column("label")

tree_circ <- ggtree(tree_circ, layout = "circular")

genus_native <- genus_abund %>% filter(status == "Native")
genus_invaded <- genus_abund %>% filter(status == "Invaded")
genus_eradicated <- genus_abund %>% filter(status == "Eradicated")

genus_native$status     <- as.factor(genus_native$status)
genus_eradicated$status <- as.factor(genus_eradicated$status)
genus_invaded$status    <- as.factor(genus_invaded$status)

Phylum_tile_df <- tree_circ$data %>%
  dplyr::select(label) %>%
  left_join(taxa_top50, by = "label") %>%
  dplyr::mutate(
    y = label,
    Phylum_strip = factor(Phylum)
  ) %>%
  dplyr::select(y, Phylum_strip)

phylum_colors <- c( "#0072B2",  # blue
                           "#D55E00",  # vermilion
                           "#009E73",  # green
                           "#CC79A7",  # magenta
                           "#F0E442",  # yellow
                           "#56B4E9",  # sky blue
                           "#8C564B",  # brown
                           "#E69F00",  # orange
                           "#000000",  # black
                           "#999999"   # grey
)


taxa_top50$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", taxa_top50$Genus)


tree_circ %<+% taxa_top50 +
  geom_fruit(
    data = Phylum_tile_df,
    geom = geom_tile,
    mapping = aes(y = y, fill = Phylum_strip),
    width = 0.45,
    offset = 0.05,
    color = NA,
    alpha = 0.3
  )+
  geom_tiplab(aes(label = Genus), size = 2,
              align = TRUE) + 
   geom_fruit(
    data = genus_eradicated,
    geom = geom_bar,
    mapping = aes(x = Abundance, y = OTU),
    stat = "identity",
    orientation = "y",
    fill = "#4DAF4A",  
    width = 0.4,
    offset = 0.85
  ) +
    geom_fruit(
    data = genus_invaded,
    geom = geom_bar,
    mapping = aes(x = Abundance, y = OTU),
    stat = "identity",
    orientation = "y",
    fill = "#E41A1C", 
    width = 0.4,
    offset = 0.01
  ) +
  geom_fruit(
    data = genus_native,
    geom = geom_bar,
    mapping = aes(x = Abundance, y = OTU),
    stat = "identity",
    orientation = "y",
    fill = "#377EB8",  
    width = 0.4,
    offset = 0.01
  ) +
  scale_color_manual(values= c("Native" = "#377EB8", "Eradicated" = "#4DAF4A", "Invaded"= "#E41A1C")) +
  scale_fill_manual(values = phylum_colors) +
  theme(legend.position = "right")



##top 20

ps_phy <- tax_glom(ps, taxrank = "Phylum")
ps_phy_comp <- ps_phy %>% microbiome ::transform ("compositional") 

top15 <- names(sort(taxa_sums(ps_phy_comp), decreasing=TRUE))[1:15]
ps.top15 <- prune_taxa(top15, ps_phy_comp)

library(speedyseq)
sample_data(ps_phy)
ps_phy_plot <- merge_samples2(ps_phy, group = "status")
ps_pl_comp <- ps_phy_plot %>% microbiome ::transform ("compositional") 

ps_phy_st <- merge_samples2(ps_phy, group = "status")
ps_st_comp <- ps_phy_st %>% microbiome ::transform ("compositional") 


top15 <- names(sort(taxa_sums(ps_st_comp), decreasing=TRUE))[1:55]
ps.top15 <- prune_taxa(top15, ps_st_comp)

keep9 <- c(
  "Acidobacteriota","Actinobacteriota","Chloroflexi" , "Campylobacterota" , "Bacteroidota","Chloroflexi",
  "Firmicutes","Myxococcota","Patescibacteria","Proteobacteria","Verrucomicrobiota"
)

df_top <- psmelt(ps.top15)
df_top$Phylum<- ifelse(
  df_top$Phylum %in% keep9,
  df_top$Phylum,
  "Other")

phylum_colors_pastel <- c(
  "Acidobacteriota" = "#CCE3F0",  
  "Actinobacteriota" ="#F7DFCC",  
  "Bacteroidota" = "#CCECE3",
  "Chloroflexi" = "lavender",
  "Campylobacterota" = "#F5E4ED",  
  "Firmicutes" = "#FCFAD9",  
  "Myxococcota" ="#DDF0FB",  
  "Patescibacteria" = "#E8DDDB",  
  "Proteobacteria" ="#FAECCC",  
  "Verrucomicrobiota" = "#CCCCCC", 
  "Other" = "mintcream" 
)

# Plot with ggplot
ggplot(df_top, aes(x = status, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors_pastel) +
  theme_minimal() +
  theme(
    legend.position = "right"
  ) +
  labs(x = "Status", y = "Abundance", fill = "Phylum")

x##Diversity indexes 

ps_plot<- merge_samples(ps, "name", fun=mean)
index <- estimate_richness(ps, measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
index$name<- rownames(index)
index <- merge(index, sample_df, by = "name")
rownames(index) <- index$Sample
table(index$status)
table(sample_df$name)

obs <- aov(Observed ~ status, data = index)
summary(obs)
shapiro.test(resid(obs))
leveneTest(resid(obs) ~ index$status)

obs_plot <-  ggplot(index, aes(x=status, y=Observed, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c( "#4DAF4A", "#E41A1C", "#377EB8"))+
  ylab("Observed Richness")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.text.x = element_blank()) 

obs_plot

chao <- aov(Chao1 ~ status, data = index)

summary(chao)


chao_plot <-  ggplot(index, aes(x=status, y=Chao1, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c( "#4DAF4A", "#E41A1C", "#377EB8"))+
  ylab("Chao1 Index")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none" ,
          axis.text.x = element_blank()) 
chao_plot



shannon <- aov(Shannon ~ status, data = index)
summary(shannon)
TukeyHSD(shannon)
shapiro.test(resid(shannon))
leveneTest(resid(shannon) ~ index$status)


sh_plot <-  ggplot(index, aes(x=status, y=Shannon, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c("#4DAF4A", "#E41A1C", "#377EB8"))+
  ylab("Shannon Index")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.text.x = element_blank()
  ) +
  stat_compare_means(
    aes(label = ..p.signif..),     
    method = "t",            # or "t.test" if you prefer
    comparisons = list(c("Native", "Invaded")),  # specify the pair
    label.y = max(index$Shannon) * 0.96      # position above the max
  )
sh_plot

simpson <- aov(InvSimpson ~ status, data = index)
summary(simpson)
TukeyHSD(simpson)
shapiro.test(resid(simpson))
leveneTest(resid(simpson) ~ index$status)

library("ggpubr")

simp_plot <-  ggplot(index, aes(x=status, y=InvSimpson, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c( "#4DAF4A", "#E41A1C", "#377EB8"))+
  ylab("InvSimpson Index")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none") +
  stat_compare_means(
    aes(label = ..p.signif..),     
    method = "t",            # or "t.test" if you prefer
    comparisons = list(c("Native", "Invaded")),  # specify the pair
    label.y = max(index$InvSimpson) * 0.8       # position above the max
  )
simp_plot

(obs_plot / sh_plot  / simp_plot)

library(vegan)

##NMDS beta div
ps_comp  <- ps %>% microbiome::transform("compositional")
nmds_its <-  ordinate(ps_comp, method="NMDS", distance="bray")
sample_data(ps_comp)
stress_val <- round(nmds_its$stress, 3)

nmds_df <- as.data.frame(scores(nmds_its, display = "sites"))
nmds_df$Sample <- rownames(nmds_df)

meta<- as.data.frame(sample_data(ps_comp))
nmds_df <- left_join(nmds_df, meta, by = "Sample")

bray_dist <- phyloseq::distance(ps_comp, method = "bray")
bray_mat <- as.matrix(bray_dist)


hulls <- nmds_df %>%
  group_by(status) %>%
  slice(chull(NMDS1, NMDS2)) %>%
  ungroup()


plot_nmds <- plot_ordination(ps_comp, nmds_its, color = "status") +
  geom_point(size = 3) +  
  geom_polygon(data = hulls,
               aes(x = NMDS1, y = NMDS2, fill = status, group = status),
               alpha = 0.2, color = NA) +
  scale_color_manual(values= c("Native" = "#377EB8", "Eradicated" = "#4DAF4A", "Invaded"= "#E41A1C")) + 
  scale_fill_manual(values = c(Native     = "#377EB8",  Eradicated = "#4DAF4A", Invaded    = "#E41A1C"
  )) +
  annotate("text", x = Inf, y = Inf, label = paste("Stress =", stress_val),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic")+
  theme_minimal()  
plot_nmds


perma_bray <- adonis2(bray_mat ~ meta$status, permutations = 9999, method = bray)
print(perma_bray)


save.image("~/Documents/seqs/carpobrotus_16s/Sequences_noadapt/dada2_16s_Giglio.RData")

saveRDS(ps, "ps_Gi_16s.RDS")


###

ps_phylum <- tax_glom(ps_comp, taxrank = "Phylum")
taxa_names(ps_phylum) <- as.character(tax_table(ps_phylum)[, "Phylum"])


# Extract abundance data
abund_data <- psmelt(ps_phylum)

# Group and summarize relative abundance per sample and phylum
phylum_abund <- abund_data %>%
  group_by(Sample, Phylum, status) %>%  
  summarise(RelAbundance = sum(Abundance), .groups = "drop")

# Loop through each phylum and run ANOVA
anova_results <- phylum_abund %>%
  group_by(Phylum) %>%
  summarise(
    F = {
      m <- aov(RelAbundance ~ status, data = cur_data())
      summary(m)[[1]]$`F value`[1]
    },
    p.value = {
      m <- aov(RelAbundance ~ status, data = cur_data())
      summary(m)[[1]]$`Pr(>F)`[1]
    }
  ) %>%
  arrange(p.value)

print(anova_results)

significant_phyla <- anova_results %>%
  filter(p.value < 0.05)
anova_results$Phylum <- gsub("p__", "", anova_results$Phylum)

tukey_results <- list()

# Loop over significant phyla and run TukeyHSD
for (phy in significant_phyla$Phylum) {
  
  # Subset data for that phylum
  df_phy <- phylum_abund %>%
    filter(Phylum == phy)
  
  # Fit ANOVA model
  aov_model <- aov(RelAbundance ~ status, data = df_phy)
  
  # Run Tukey HSD
  tukey_out <- TukeyHSD(aov_model)
  
  # Store results
  tukey_results[[phy]] <- tukey_out
}

# View all results
tukey_results



# Calculate mean relative abundance per phylum
top_phyla <- abund_data %>%
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

# View top 10 most abundant phyla
head(top_phyla, 20)


