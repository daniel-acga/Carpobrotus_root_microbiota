################################################################################
# ASV + TAX to Faprotax input
################################################################################

library(dplyr)
library(stringr)
library(tidyverse)

# Preparing data ----
asv <- read.table("../data/networking/its/keystone_otu_nat.csv", sep=",", header=TRUE, row.names = 'id')
tax <- read.table("../data/networking/its/keystone_tax_nat.csv", sep=",", header=TRUE)

# Exact number of rows
asv <- as.data.frame(t(asv))
tax <- tax[,-1]
tax$Species[ is.na(tax$Species) | tax$Species == "" ] <- paste(tax$Genus, "sp.")

small_tax <- tax[tax$ids==rownames(asv),]

tax$taxonomy <- apply(tax[, 2:7], 1, function(x) paste(x[!is.na(x) & x!=""], collapse=";"))

##Nat 8 Citrobacter sp.  27 Rhodoplanes sp 71. UBA6140
#Inv 3 Rhodopila 14 Aminobacter 78 Sporocytophaga. 93 Bosea 114 Paracoccus 


# Creating the data for FAPROTAX
tax_4_funguild <- data.frame(asv, taxonomy = tax$taxonomy)

rownames(tax)
rownames(tax_4_funguild) <- paste("OTU",rownames(tax))
write.table(tax_4_funguild, "../data/funguild/tax_4_fun_nat.tsv", sep="\t", quote = FALSE)

################################################################################
# This part is in bash!!
################################################################################
#Check the paths of the collapse_table.py and the FAPROTAX.txt

collapse_table.py -i tax_4_faprotax_inv.tsv -o func_table.tsv -g FAPROTAX.txt -d "taxonomy" -c "#" -r report.txt -v

  python Guilds_v1.0.py -otu tax_4_fun_inv.tsv -db fungi -m -u


################################################################################
# Functional Analysis using Faprotax output
################################################################################
# This analysis takes into account ONLY the functional groups, but not the specific taxa.

library(phyloseq)
library(ggplot2)
library(paletteer)
library(tidyr)

setwd("path")

# Faprotax output
func_table_inv <- read.csv("../data/faprotax/func_table_inv.tsv",header = T, row.names = 1, sep = '\t')
func_table_inv <- func_table_inv[,-1]

func_table_nat <- read.csv("../data/faprotax/func_table_nat.tsv",header = T, row.names = 1, sep = '\t')
func_table_nat <- func_table_nat[,-1]

func_table_inv$Group <- "Invaded"
func_table_nat$Group<- "Native"

native_long <- func_table_nat %>%
  rownames_to_column("Function") %>%
  pivot_longer(cols=where(is.numeric), names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Group = "Native")

invaded_long <- func_table_inv %>%
  rownames_to_column("Function") %>%
  pivot_longer(cols=where(is.numeric), names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Group = "Invaded")

plot_df <- bind_rows(native_long, invaded_long) %>%
  filter(Abundance > 0) %>%  # Remove zero-abundance values
  group_by(Group, Function) %>%
  summarise(RelativeAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")


ggplot(plot_df, aes(x = Function, y = Group, size = RelativeAbundance, fill = Group)) +
  geom_point(shape = 21, alpha = 0.8) +
  scale_size(range = c(1, 10)) +
  scale_fill_manual(values = c("Native" = "#377EB8", "Invaded" = "#E41A1C")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(size = "Abundance", fill = "Group", x = "Function", y = "")

# Matching the sample names
metadata <- read.delim("metadata.tsv", header=TRUE, sep = "\t",  row.names = 1, check.names = FALSE)
colnames(func_table) <- rownames(metadata)

# Phyloseq object
samples <-  sample_data(metadata)
functions <- otu_table(as.matrix(t(func_table)), taxa_are_rows = F)
ps <- phyloseq(functions,samples)

## FUNGuild output

# Load data
fun_inv <- read.delim("../FUNGuild/tax_4_fun_inv.guilds.txt", sep = "\t", header = TRUE)
fun_nat <- read.delim("../FUNGuild/tax_4_fun_nat.guilds.txt", sep = "\t", header = TRUE)

fun_inv$status <- "Invaded"
fun_nat$status <- "Native"
sample_cols <- grep("^[CSR]\\d+", names(df), value = TRUE)

summarize_guilds <- function(df) {
  sample_cols <- grep("^[CSR]\\d+", names(df), value = TRUE)
  
  df %>%
    filter(!is.na(Guild)) %>%
    mutate(Total = rowSums(across(all_of(sample_cols)))) %>%
    mutate(
      Guild = strsplit(as.character(Guild), "-"),
      Trophic.Mode = strsplit(as.character(Trophic.Mode), "-")
    ) %>%
    unnest(Guild) %>%
    unnest(Trophic.Mode) %>%
    mutate(
      Guild = str_trim(gsub("\\|", "", Guild)),
      Trophic.Mode = str_trim(gsub("\\|", "", Trophic.Mode))
    ) %>%
    group_by(Guild, status) %>%
    summarise(Abundance = sum(Total), .groups = "drop") %>%
    group_by(status) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance))
}

guild_native <- summarize_guilds(fun_nat)
guild_invaded <- summarize_guilds(fun_inv)

guilds_combined <- bind_rows(guild_native, guild_invaded)


# Replace small or blank guilds with "Other", except for key ones
guilds_cleaned <- guilds_combined %>%
  mutate(Guild = ifelse(
    (RelativeAbundance < 0.05 | Guild == "" | is.na(Guild)) & 
      !Guild %in% c("Ectomycorrhizal", "Arbuscular Mycorrhizal"),
    "Other",
    Guild
  )) %>%
  group_by(Guild, status) %>%
  summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop")

library(scales) 

ggplot(guilds_cleaned, aes(x = 2, y = RelativeAbundance, fill = Guild)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ status) +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(fill = "Guild") +
  geom_text(
    aes(
      label = ifelse(
        RelativeAbundance > 0.05,
        paste0(Guild, "\n(", percent(RelativeAbundance, accuracy = 0.001), ")"),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  )


wilcox.test(Abundance ~ status, data = guilds_combined)

results <- guilds%>%
  group_by(Guild) %>%
  summarise(p_value = wilcox.test(RelativeAbundance ~ Environment)$p.value)

