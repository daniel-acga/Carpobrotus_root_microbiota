library(dada2)
library(ShortRead)
library(Biostrings)
library(microbiome)
library(dplyr)
library(tidyverse)


path <- "~/Documents/seqs/carpobrotus_its/Sequences_renamed"
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(basename(fnFs), function(x) sub("_1\\..*$", "", x))
sample.names



#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])




##Identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


#cutadapt

cutadapt <- "/opt/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") 

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = 
        sapply(FWD.orients,
               primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                            fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))




# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fq.gz", full.names = TRUE))

#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", paste0(sample.names, "_1_filt.fq.gz"))
filtRs <- file.path(path.cut, "filtered", paste0(sample.names, "_2_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
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
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

save.image("~/Documents/seqs/carpobrotus_its/Sequences_renamed/dada2_its_Giglio.RData")

#Assign taxonomy 
taxa <- assignTaxonomy(seqtab.nochim, '/Users/danielacosta/Documents/seqs/tax/sh_general_release_dynamic_04.04.2024_dev.fasta', 
                       multithread = TRUE, tryRC = TRUE, verbose=TRUE)


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
            "T147_new", "T147_new", "T134", "T134", "T134", "T130_new", "T134", "T150", "T150",
            "T132_new", "T132_new", "82", "82", "80", "T133_new", "80","81_new","81_new",
            "INV_1_new","INV_1_new","INV_2_new","INV_2_new","84_new","84_new","T133_new","T148_new","T148_new"),
  year  = c("2023", "2023", "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023"))

rownames(sample_df) <- c("C1","C11","C12","C13","C15","C16","C18","C19","C21","C22","C24","C25","C27","C28","C3","C30","C31","C32","C33","C34","C35","C37","C39","C4","C40","C41","C43","C44","C46","C47","C49","C50","C52","C6","C7","C9")


library(phyloseq)
rownames(seqtab.nochim) <- c("C1","C11","C12","C13","C15","C16","C18","C19","C21","C22","C24","C25","C27","C28","C3","C30","C31","C32","C33","C34","C35","C37","C39","C4","C40","C41","C43","C44","C46","C47","C49","C50","C52","C6","C7","C9")

#BUilding a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa))


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
save.image("~/Documents/seqs/carpobrotus_its/Sequences_renamed/dada2_its_Giglio.RData")

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
  group_by(OTU, Genus, status) %>%  # Keep OTU (or ASV) info to match tree tip labels
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Get OTUs associated with top 50 genera
keep_otus <- unique(genus_abund$OTU)
ps_top50 <- prune_taxa(keep_otus, ps_genus)

library(ggtree)
library(ggtreeExtra)
library(ggplot2)

#Reshape tree

tree_circ <- phy_tree(ps_top50)
taxa_top50 <- as.data.frame(tax_table(ps_top50)) %>%
  rownames_to_column("label")

tree_circ <- ggtree(tree_circ, layout = "circular")


taxa_top50$Phylum <- sub("p__", "", taxa_top50$Phylum)
taxa_top50$Genus <- sub("g__", "", taxa_top50$Genus)


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

phylum_colors <- c(
  "#332288",  # dark indigo
  "#88CCEE",  # light cyan
  "#117733",  # forest green
  "#DDCC77",  # goldenrod
  "#CC6677",  # rose red
  "#882255",  # maroon
  "#AA4499"   # mauve
)

taxa_top50$Genus <- gsub("gen_Incertae_sedis", "Incertae_sedis", taxa_top50$Genus)


tree_circ %<+% taxa_top50 +
  geom_fruit(
    data = Phylum_tile_df,
    geom = geom_tile,
    mapping = aes(y = y, fill = Phylum_strip),
    width = 0.35,
    offset = 0.65,
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
    offset = 0.89
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


top15 <- names(sort(taxa_sums(ps_st_comp), decreasing=TRUE))[1:22]
ps.top15 <- prune_taxa(top15, ps_st_comp)

keep6 <- c(
  "Ascomycota","Basidiomycota","Chytridiomycota" , "Glomeromycota" , "Mortierellomycota","Rozellomycota"
)

df_top <- psmelt(ps.top15)
df_top$Phylum <- gsub("p__", "", df_top$Phylum)
df_top$Phylum<- ifelse(
  df_top$Phylum %in% keep6,
  df_top$Phylum,
  "Other")


phylum_colors <- c(
  "Ascomycota"       = "#D6D3E7",
  "Basidiomycota"    = "#E7F5FC",
  "Chytridiomycota"  = "#CFE4D6",
  "Glomeromycota"    = "#F8F5E4",
  "Mortierellomycota" = "#F5E0E4",
  "Rozellomycota"    = "#E7D3DD",
  "Other" = "mintcream" 
 )


# Plot with ggplot
ggplot(df_top, aes(x = status, y = Abundance, fill = Phylum)) +
  scale_fill_manual(values = phylum_colors)+
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    legend.position = "right"
  ) +
  labs(x = "Status", y = "Abundance", fill = "Phylum")

df_top <- psmelt(ps.top15)

# Plot with ggplot
ggplot(df_top, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~status, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Hide long sample names
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Abundance", fill = "Phylum")

##Diversity indexes 

index <- estimate_richness(ps, measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
index$Sample <- rownames(index)
index <- merge(index, sample_df, by = "Sample")
rownames(index) <- index$Sample

index$Obs_st <- decostand(index$Observed, method="log")
index$InvSimpson_st <- decostand(index$InvSimpson, method="log")

obs <- aov (Obs_st ~ status, data = index)
summary(obs)
shapiro.test(resid(obs))
leveneTest(resid(obs) ~ index$status)

obs_plot <-  ggplot(index, aes(x=status, y=Observed, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c("#4DAF4A", "#E41A1C", "#377EB8"))+
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
  theme(legend.position = "none",
        axis.text.x = element_blank()) 
chao_plot


shannon <- aov(Shannon ~ status, data = index)
summary(shannon)
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
        axis.text.x = element_blank() ) 
sh_plot

simpson <- aov(InvSimpson_st ~ status, data = index)
summary(simpson)
TukeyHSD(simpson)
shapiro.test(resid(simpson))
leveneTest(resid(simpson) ~ index$status)


simp_plot <-  ggplot(index, aes(x=status, y=InvSimpson_st, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c( "#4DAF4A", "#E41A1C", "#377EB8"))+
  ylab("InvSimpson Index")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none") 
simp_plot

( obs_plot / sh_plot / simp_plot)

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

save.image("~/Documents/seqs/carpobrotus_its/Sequences_renamed/dada2_its_Giglio.RData")

saveRDS(ps, "ps_Gi_its.RDS")


perma_bray <- adonis2(bray_mat ~ meta$status, permutations = 9999, method = bray)
print(perma_bray)


#########
# Extract abundance data
abund_data <- psmelt(ps_phy_comp)

# Group and summarize relative abundance per sample and phylum
phylum_abund <- abund_data %>%
  group_by(Sample, Phylum, status) %>%  # Assuming 'status' is in sample_data
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

significant_phyla

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

head(top_phyla, 20)



phylum_by_status <- abund_data %>%
  group_by(status, Phylum) %>% 
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(status, desc(mean_abundance))


phylum_by_status_wide <- phylum_by_status %>%
  pivot_wider(
    names_from   = status,
    values_from  = mean_abundance,
    values_fill  = 0
  ) %>%
  arrange(desc(Invaded))  # or whichever status you want to sort by

# View wide-form
View(phylum_by_status_wide)


write.csv(phylum_by_status_wide, "rel_ab_phyla_its.csv")


