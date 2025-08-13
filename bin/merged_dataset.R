library(phyloseq)
library(microbiome)
library(lmerTest)
library(dplyr)
library(rcompanion)
library(reshape2)
library(vegan)
library(tidyverse)

setwd("~/Documents/seqs/bin")
#Read all phyloseq objects and merge them ----

ps_Gi = readRDS("../data/ps_Gi_16s.RDS")
ps_RC = readRDS("../data/ps_RC_16s.RDS")
ps_Ca = readRDS("../data/ps_Ca_16s.RDS")
ps_SS = readRDS("../data/ps_SS_16s.RDS")

sort(sample_sums(ps_RC))

ps_Gi_f = subset_samples (ps_Gi, status != "Eradicated")
ps_Gi_g = tax_glom(ps_Gi_f, taxrank = "Genus")
ps_RC_g = tax_glom(ps_RC, taxrank ="Genus")
ps_Ca_g = tax_glom(ps_Ca, taxrank ="Genus")

taxa_names(ps_Gi_g) <- as.character(tax_table(ps_Gi_g)[, "Genus"])
taxa_names(ps_RC_g) <- as.character(tax_table(ps_RC_g)[, "Genus"])
taxa_names(ps_Ca_g) <- as.character(tax_table(ps_Ca_g)[, "Genus"])

ps_merge_1 <- merge_phyloseq(ps_Gi_g, ps_RC_g)
ps_merge_genus <- merge_phyloseq(ps_merge_1, ps_Ca_g)

##Rarefaction curve ----
library(iNEXT)
otu_mat <- as.matrix(otu_table(ps_merge_genus))
otu_mat = t(otu_mat)

otu_list <- apply(otu_mat, 2, function(x) {
  x <- x[x > 0]                     # remove 0s
  if (length(x) >= 2 && sum(x) > 0) return(x) else return(NULL)  # only keep valid samples
})

otu_list <- otu_list[!sapply(otu_list, is.null)]
inext_out <- iNEXT(otu_list, q = 0, datatype = "abundance")
?merge_taxa
ggiNEXT(inext_out, type = 1) + 
  theme_minimal() +
  labs(title = "Rarefaction Curves", x = "Sample Size", y = "Species Richness")

sort(sample_sums(ps_merge_genus))
min_sample_size = 14286
ps_rare <- rarefy_even_depth(ps_merge_genus, sample.size = min_sample_size,
                             rngseed = 9999)


##Diversity indexes linear mixed models ----

index<- estimate_richness(ps_rare, measures=c("Observed", "Chao1", "Shannon", "InvSimpson", "ACE"))
sample_df <- sample_data(ps_rare)
index$Sample <- rownames(index)
index <-left_join(index, sample_df, by = "Sample") 
rownames(index)  <- index$Sample

library(bestNormalize)

bn_SH <- bestNormalize(index$Shannon, allow_lambert_s = TRUE)
SH_bn  <- predict(bn_SH)
index$Shannon_st <- SH_bn

bn_SH$chosen_transform


bn_IS <- bestNormalize(index$InvSimpson, allow_lambert_s = FALSE)
IS_bn  <- predict(bn_IS)
index$InvSimpson_st <- IS_bn

bn_IS$chosen_transform


obs <- lmerTest:::lmer(Observed ~ status + (1 |location), data = index)
summary(obs)  


plot(obs) ##Check for linearity 

# check for homocedasticity on the residuals
df_diag <- data.frame(
  status       = sample_data(obs@frame)$status,  # or your original data$status
  resid_pearson = resid(obs, type = "pearson")
)

# Boxplot of residuals for each status
boxplot(resid_pearson ~ status, data = df_diag,
        xlab = "Invasion Status",
        ylab = "Pearson residuals",
        main = "Residuals by Status")
abline(h = 0, lty = 2, col = "gray50")
library(car)
leveneTest(resid_pearson ~ status, data = df_diag)

##Check for normality of residuals
shapiro.test(resid(obs))
qqnorm(resid(obs)); qqline(resid(obs))

library(lattice) 
re_loc <- lme4::ranef(obs)$location[, "(Intercept)"]  ##Check for normality of random effects

# Q–Q plot of the location random effects
qqnorm(re_loc,
       main = "Q–Q Plot of Location Random Effects",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(re_loc, col = "steelblue", lwd = 2)

hist(re_loc, 
     breaks = 10, 
     col = "lightgray", 
     border = "white",
     main = "Distribution of Location Random Effects",
     xlab = "Random intercept (location)")
# optionally overlay a density
lines(density(re_loc), col = "darkred", lwd = 2)

library(influence.ME) #check for influential observations by comparing cook distances
infl <- influence(obs, group = "location")
plot(infl, which = "cook")


acf(resid(obs)) ##check for autocorrelation


sample_data(ps_rare)
table(sample_df$location)

index_obs <- mutate(index,
                    prediction = fitted(obs),
                    resid = Observed - prediction,
                    resid2 = resid^2)
ggplot(index_obs, aes(status, resid)) + 
  stat_summary() + labs(y = 'Mean Squared Error')

VarCorr(obs)$location[1]

obs_plot <-  ggplot(index, aes(x=status, y=Observed, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c( "#E41A1C", "#377EB8"))+
  ylab("Observed Richness")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.text.x = element_blank()) 

obs_plot


##Chao

chao1 <- lmerTest::lmer(Observed ~ status + (1 | location), data = index)
summary(chao1)  
plot(chao1) 

# check for homocedasticity on the residuals
df_diag <- data.frame(
  status       = sample_data(chao1@frame)$status, 
  resid_pearson = resid(chao1, type = "pearson")
)

# Boxplot of residuals for each status
boxplot(resid_pearson ~ status, data = df_diag,
        xlab = "Invasion Status",
        ylab = "Pearson residuals",
        main = "Residuals by Status")
abline(h = 0, lty = 2, col = "gray50")
library(car)
leveneTest(resid_pearson ~ status, data = df_diag)

##Check for normality of residuals
shapiro.test(resid(chao1))
qqnorm(resid(chao1)); qqline(resid(chao1))

library(lattice) 
re_loc <- lme4::ranef(chao1)$location[, "(Intercept)"]  ##Check for normality of random effects

# Q–Q plot of the location random effects
qqnorm(re_loc,
       main = "Q–Q Plot of Location Random Effects",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(re_loc, col = "steelblue", lwd = 2)

hist(re_loc, 
     breaks = 10, 
     col = "lightgray", 
     border = "white",
     main = "Distribution of Location Random Effects",
     xlab = "Random intercept (location)")
# optionally overlay a density
lines(density(re_loc), col = "darkred", lwd = 2)


library(influence.ME) #check for influential observations by comparing cook distances
infl <- influence(chao1, group = "location")
plot(chao1, which = "cook")


acf(resid(chao1)) ##check for autocorrelation


sample_data(ps_rare)
table(sample_df$location)

index_chao1 <- mutate(index,
                    prediction = fitted(chao1),
                    resid = Observed - prediction,
                    resid2 = resid^2)
ggplot(index_chao1, aes(status, resid)) + 
  stat_summary() + labs(y = 'Mean Squared Error')

VarCorr(obs)$location[1]



shannon <- lmerTest::lmer(SH_bn ~ status + (1 | location), data = index)
summary(shannon)  

sh_plot <-  ggplot(index, aes(x=status, y=Shannon, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c("#E41A1C", "#377EB8"))+
  ylab("Shannon Index")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.text.x = element_blank())

sh_plot


plot(shannon) ##Check for linearity 

# check for homocedasticity on the residuals
df_diag <- data.frame(
  status       = sample_data(shannon@frame)$status,  # or your original data$status
  resid_pearson = resid(shannon)
)

# Boxplot of residuals for each status
boxplot(resid_pearson ~ status, data = df_diag,
        xlab = "Invasion Status",
        ylab = "Pearson residuals",
        main = "Residuals by Status")
abline(h = 0, lty = 2, col = "gray50")
library(car)
leveneTest(resid_pearson ~ status, data = df_diag)

##Check for normality of residuals
shapiro.test(resid(shannon))
qqnorm(resid(shannon)); qqline(resid(shannon))

df_diag <- data.frame(
  fitted = fitted(shannon),
  resid  = resid(shannon)
)
ggplot(df_diag, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted",
       x = "Fitted values",
       y = "Residuals") +
  theme_minimal()
library(lattice) 
re_loc <- lme4::ranef(obs)$location[, "(Intercept)"]  ##Check for normality of random effects

# Q–Q plot of the location random effects
qqnorm(re_loc,
       main = "Q–Q Plot of Location Random Effects",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(re_loc, col = "steelblue", lwd = 2)

hist(re_loc, 
     breaks = 10, 
     col = "lightgray", 
     border = "white",
     main = "Distribution of Location Random Effects",
     xlab = "Random intercept (location)")
# optionally overlay a density
lines(density(re_loc), col = "darkred", lwd = 2)

library(influence.ME) #check for influential observations by comparing cook distances
infl <- influence(obs, group = "location")
plot(infl, which = "cook")


acf(resid(obs)) ##check for autocorrelation


sample_data(ps_rare)
table(sample_df$location)

index_obs <- mutate(index,
                    prediction = fitted(obs),
                    resid = Observed - prediction,
                    resid2 = resid^2)
ggplot(index_obs, aes(status, resid)) + 
  stat_summary() + labs(y = 'Mean Squared Error')

VarCorr(obs)$location[1]

simpson <- lmerTest::lmer(InvSimpson_bn ~ status + (1 | location), data = index)
summary(simpson)  

simp_plot <-  ggplot(index, aes(x=status, y=InvSimpson, fill=status))+ 
  geom_violin (width = 0.3) +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  scale_fill_manual(values=c( "#E41A1C", "#377EB8"))+
  ylab("InvSimpson Index")+ xlab("")+ 
  theme_minimal()+ 
  theme(legend.position = "none") 
simp_plot


plot(simpson)


(obs_plot / sh_plot  / simp_plot)


# check for homocedasticity on the residuals
df_diag <- data.frame(
  status       = sample_data(shannon@frame)$status,  # or your original data$status
  resid_pearson = resid(shannon)
)

# Boxplot of residuals for each status
boxplot(resid_pearson ~ status, data = df_diag,
        xlab = "Invasion Status",
        ylab = "Pearson residuals",
        main = "Residuals by Status")
abline(h = 0, lty = 2, col = "gray50")
library(car)
leveneTest(resid_pearson ~ status, data = df_diag)

##Check for normality of residuals
shapiro.test(resid(shannon))
qqnorm(resid(shannon)); qqline(resid(shannon))

df_diag <- data.frame(
  fitted = fitted(shannon),
  resid  = resid(shannon)
)
ggplot(df_diag, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted",
       x = "Fitted values",
       y = "Residuals") +
  theme_minimal()
library(lattice) 
re_loc <- lme4::ranef(obs)$location[, "(Intercept)"]  ##Check for normality of random effects

# Q–Q plot of the location random effects
qqnorm(re_loc,
       main = "Q–Q Plot of Location Random Effects",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(re_loc, col = "steelblue", lwd = 2)

hist(re_loc, 
     breaks = 10, 
     col = "lightgray", 
     border = "white",
     main = "Distribution of Location Random Effects",
     xlab = "Random intercept (location)")
# optionally overlay a density
lines(density(re_loc), col = "darkred", lwd = 2)

library(influence.ME) #check for influential observations by comparing cook distances
infl <- influence(obs, group = "location")
plot(infl, which = "cook")


acf(resid(obs)) ##check for autocorrelation


sample_data(ps_rare)
table(sample_df$location)

index_obs <- mutate(index,
                    prediction = fitted(obs),
                    resid = Observed - prediction,
                    resid2 = resid^2)
ggplot(index_obs, aes(status, resid)) + 
  stat_summary() + labs(y = 'Mean Squared Error')

VarCorr(obs)$location[1]

library(MuMIn)
library(effsize)
library(parameters)
library(tidyverse)

div_long <- index %>%
  pivot_longer(cols = c(Shannon_st, InvSimpson_st, Observed),
               names_to = "index_name", values_to = "value")

cohen_results <- div_long %>%
  group_by(index_name) %>%
  summarise(
    cohen_d = cohen.d(value ~ status, hedges.correction = TRUE)$estimate,
    conf_low = cohen.d(value ~ status, hedges.correction = TRUE)$conf.int[1],
    conf_high = cohen.d(value ~ status, hedges.correction = TRUE)$conf.int[2]
  )

r.squaredGLMM(shannon)
r2_results <- div_long %>%
  group_by(index_name) %>%
  group_map(~{
    model <- lmer(value ~ status + (1 | location), data = .x)
    r2 <- r.squaredGLMM(model)
    tibble(m_r2 = r2[1], c_r2 = r2[2])
  }) %>%
  bind_rows()

effect_df <- bind_cols(cohen_results, r2_results)


ggplot(effect_df, aes(y = index_name)) +
  geom_point(aes(x = cohen_d), color = "tomato", size = 3) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), color = "tomato", height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(x = m_r2), color = "steelblue", size = 3, shape = 17) +
  geom_point(aes(x = c_r2), color = "lightblue", size = 3, shape = 17) +
  labs(
    subtitle = "Cohen's d (red dots) vs. R² (blue triangles)",
    x = "Effect Size",
    y = "Diversity Index"
  ) +
  theme_minimal()


##NMDS Beta diversity ----
ps_rare_comp <- ps_rare %>% microbiome::transform("compositional")
ps_rare_filt <- filter_taxa(ps_rare, function(x) sum(x > 0) > 5, prune = TRUE)
ps_rf_comp <- ps_rare_filt %>% microbiome::transform("compositional")
ps_mg_comp <- ps_merge_genus %>% microbiome::transform("compositional")


ps_r_comp  <- subset_samples(ps_rare_comp, location != "G")

nmds_ord <- ordinate(ps_rare_comp, method = "NMDS", distance = "bray")
stress_val <- round(nmds_ord$stress, 3)


colours <- c("black", "magenta", "red", "green", "darkgreen", "skyblue", "purple",
              "darkblue", "red4", "pink", "orange", "yellow") 

nmds_plot <- plot_ordination(ps_rare_comp, nmds_ord, 
                             color = "location", shape = "status") +
  geom_point(size = 3, alpha = 0.8) +
  annotate("text", x = Inf, y = Inf, label = paste("Stress =", round(stress_val, 3)),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  scale_color_manual(values = colours) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(color = "Location", shape = "Status")

nmds_plot


bray_dist <- phyloseq::distance(ps_rare_comp, method = "bray")
bray_mat <- as.matrix(bray_dist)

bray_df <-bray_df <- reshape2::melt(bray_mat, varnames = c("Sample1", "Sample2"), value.name = "BrayCurtis") %>%
  filter(as.character(Sample1) < as.character(Sample2))

##NMDS with env arrows

library(vegan)

metadata <- read.csv("../data/metadata.csv",  sep=';')
metadata <- sapply(metadata, gsub, pattern = ",", replacement= ".")
metadata$TN <- as.numeric(as.character(metadata$TN))
metadata$TOC <- as.numeric(as.character(metadata$TOC))
metadata$MAT <- as.numeric(as.character(metadata$MAT))
metadata$MAP <- as.numeric(as.character(metadata$MAP))
metadata$CN <- as.numeric(as.character(metadata$CN))
metadata$Lat <- as.numeric(as.character(metadata$Lat))
metadata$Long <- as.numeric(as.character(metadata$Long))
metadata <- as.data.frame(metadata)
metadata$location <- sample_data(ps_merge_genus)$location

#Fit env vectors

meta_rare <- as.data.frame(sample_data(ps_rare_comp))

metadata_nmds <- filter(metadata, sample %in% meta_rare$Sample)

ef <- envfit(
  nmds_ord, 
  metadata_nmds[, c("Lat", "Long", "MAP", "MAT", "TOC", "TN")],  
  permutations = 999
)
?envfit

# Extract the significant vectors (p < 0.05)
sig_vecs <- as.data.frame(scores(ef, display="vectors"))
sig_vecs$pval <- ef$vectors$pvals
sig_vecs <- sig_vecs %>% filter(pval <= 0.05)
sig_vecs$var <- rownames(sig_vecs)



# Get your NMDS site scores + metadata
sites <- as.data.frame(scores(nmds_ord, display="sites")) %>%
  rownames_to_column("sample") %>%
  left_join(metadata_nmds , by="sample")


nmds_ord <- ordinate(ps_rare_comp, method = "NMDS", distance = "bray")
stress_val <- round(nmds_ord$stress, 3)


colours <- c("black", "magenta", "red", "green", "darkgreen", "skyblue", "purple",
             "darkblue", "red4", "pink", "orange", "yellow") 

nmds_plot <- plot_ordination(ps_rare_comp, nmds_ord, 
                             color = "location", shape = "status") +
  geom_point(size = 3, alpha = 0.8) +
  annotate("text", x = Inf, y = Inf, label = paste("Stress =", round(stress_val, 3)),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  scale_color_manual(values = colours) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(color = "Location", shape = "Status")

nmds_plot +
  geom_segment(
    data = sig_vecs, inherit.aes = FALSE,
    aes(x = 0, y = 0, xend = NMDS1 * 1.5, yend = NMDS2 * 1.5),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "black"
  ) +
  geom_text(
    data = sig_vecs, inherit.aes = FALSE,
    aes(x = NMDS1 * 1.6, y = NMDS2 * 1.6, label = var),
    size = 4, fontface = "italic", color = "black"
  )

# Add metadata for both samples in each pair
bray_df <- bray_df %>%
  left_join(meta, by = c("Sample1" = "Sample")) %>%
  rename(status1 = status) %>%
  left_join(meta, by = c("Sample2" = "Sample")) %>%
  rename(status2 = status)

bray_within <- bray_df %>%
  filter(status1 == status2) %>%
  mutate(group = status1)

ggplot(bray_within, aes(x = group, y = BrayCurtis, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  geom_jitter(aes(color = group),
              width = 0.15, height = 0, size = 1.5, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Group",
    y = "Bray–Curtis Dissimilarity") +
  theme_minimal()

meta <- as.data.frame(as.matrix(meta))
meta$status <- as.factor(meta$status)
meta$location <- as.factor(meta$location)
table(meta$status, meta$status)
class(meta)



perma_bray <- adonis2(bray_mat ~ meta$status, strata = meta$location, permutations = 9999, method = bray)
print(perma_bray)



##Core microbiome ----
ps_SS = readRDS("../data/ps_SS_16s.RDS")
ps_SS_g <- tax_glom(ps_SS, taxrank ="Genus")
taxa_names(ps_SS_g) <- as.character(tax_table(ps_SS_g)[, "Genus"])

ps_inv <- subset_samples(ps_merge_genus, status == "Invaded")

ps_merge_core <- merge_phyloseq(ps_inv, ps_SS_g)
ps_core_fam <- tax_glom(ps_merge_core, taxrank ="Family")
fam_names <- as.character(tax_table(ps_core_fam)[, "Family"])
unique_fam_names <- make.unique(fam_names)
taxa_names(ps_core_fam) <- unique_fam_names
ps_coref_comp <- ps_core_fam %>% microbiome::transform("compositional")

# Define combinations
prevalence_vals <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60,0.70)
detection_vals <- c(0.001, 0.003, 0.005, 0.007, 0.01)

# Step 4: Loop through combinations
results <- list()

for (prev in prevalence_vals) {
  for (det in detection_vals) {
    combo_name <- paste0("prev", prev*100, "_abun", det*1000)
    
    core_taxa <- core_members(ps_core_comp, 
                              detection = det,
                              prevalence = prev)
    
    results[[combo_name]] <- length(core_taxa)
  }
}

library(tidyverse)
# Step 5: Convert to data frame
core_summary <- enframe(results, name = "Combo", value = "Num_Fam") %>%
  separate(Combo, into = c("Prevalence", "Abundance"), sep = "_") %>%
  mutate(Prevalence = as.numeric(gsub("prev", "", Prevalence)),
         Abundance = as.numeric(gsub("abun", "", Abundance)) / 10)

core_summary$Num_Fam<- as.numeric(core_summary$Num_Fam)


ggplot(core_summary, aes(x = Abundance, y = Prevalence, fill = Num_Fam)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Num_Fam), color = "black", size = 4) +  # <- add this line
  scale_fill_viridis_c(option = "plasma") +
  scale_x_continuous(breaks = sort(unique(core_summary$Abundance))) +
  scale_y_continuous(breaks = sort(unique(core_summary$Prevalence))) +
  labs(title = "Core Microbiome Size by Threshold",
       x = "Detection threshold (% rel. abundance)",
       y = "Prevalence threshold (%)",
       fill = "Number of Families") +
  theme_minimal()

ps_df <- psmelt(ps_coref_comp)

# Filter only present observations
present_df <- ps_df %>% filter(Abundance > 0)

# Count unique locations per Family
gen_location_counts <- present_df %>%
  group_by(Family) %>%
  summarise(LocationCount = n_distinct(Sample)) %>%
  filter(LocationCount >= 2)

# Extract only those families
valid_fam <- gen_location_counts$Family

# Prune phyloseq object# Prune phFamilyyloseq object
ps_famf_core <- prune_taxa(taxa_names(ps_coref_comp)[tax_table(ps_coref_comp)[, "Family"] %in% valid_fam], ps_coref_comp)
library(pheatmap)
?pheatmap

# Run core analysis on filtered object
core_fam <- core_members(ps_famf_core,
                              detection = 0.005,    # 0.5% abundance
                              prevalence = 0.35)

?core_members

ps_corefam_final <- prune_taxa(core_fam, ps_famf_core)


##Core microbiome per location

otu_df <- as.data.frame(otu_table(ps_corefam_final))
otu_df <- as.data.frame(t(otu_df))

taxa_df <- tax_table(ps_corefam_final) %>%
  as.data.frame() %>%
  rownames_to_column("TaxaID")

otu_df$TaxaID <- rownames(otu_df)
otu_df <- left_join(otu_df, taxa_df, by = "TaxaID")

# Pivot longer to get sample-wise presence
sample_cols <- otu_df %>%
  select(where(is.numeric)) %>%
  colnames()

#Pivot only the numeric columns
otu_long <- otu_df %>%
  pivot_longer(cols = all_of(sample_cols), 
               names_to = "Sample", 
               values_to = "Abundance")%>%
  mutate(Present = ifelse(Abundance > 0, 1, 0))


sample_core <- sample_data(ps_corefam_final)

otu_long <- left_join(otu_long, sample_core, by = "Sample")

core_count_by_location <- otu_long %>%
  filter(Present == 1) %>%
  group_by(Family, location) %>%
  summarise(Present = 1, .groups = "drop") %>%
  count(location, name = "Num_Core_Families")

fam_location_matrix <- otu_long %>%
  filter(Present == 1) %>%
  distinct(Family, location) %>%
  mutate(Presence = 1)

ggplot(fam_location_matrix, aes(x = location, y = Family)) +
  geom_tile(aes(fill = Presence), color = "white") +
  scale_fill_gradient(low = "white", high = "darkgreen", guide = "none") +
  labs(title = "Presence of Core Families Across Locations",
       x = "Location", y = "Core Familiies") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(pheatmap)

abund_loc <- otu_long %>%
  group_by(Family, Phylum, location) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

abund_mat <- abund_loc %>%
  select(Family, location, mean_abundance) %>%
  pivot_wider(names_from  = location,
              values_from = mean_abundance,
              values_fill = 0) %>%
  column_to_rownames("Family") %>%
  as.matrix()

anno_row <- abund_loc %>%
  select(Family, Phylum) %>%
  distinct() %>%
  column_to_rownames("Family")

phyla <- unique(anno_row$Phylum)
unique(anno_row$Phylum)
phy_colors <- c(
  "Acidobacteriota" = "#0072B2",  # blue
  "Actinobacteriota" = "#D55E00",  # vermilion
  "Bacteroidota"= "#009E73",  # green
  "Firmicutes"= "#F0E442",  # yellow
  "Myxococcota"= "#56B4E9",  # sky blue
  "Patescibacteria" ="#8C564B",  # brown
  "Proteobacteria" = "#E69F00",  # orange 
  "Verrucomicrobiota" = "#999999",  # gray
  "Planctomycetota" = "#000000" #black
)
library(viridis)
n_colors <- 100
color_vec <- viridis(n_colors, option = "D")
breaks_vec <- seq(0, 0.1, length.out = n_colors + 1)

soil_info <- c(
  AT = "Eutric Regosol",
  A  = "Albic Arenosol",
  CA = "Albic Arenosol", 
  CI = "Andic Cambisol",
  G  = "Vertic Cambisol",
  L  = "Albic Arenosol",
  P  = "Albic Arenosol",
  S  = "Albic Arenosol",
  SA = "Chromic Luvisol",
  SS = "Chromic Luvisol",
  T  = "Albic Arenosol",
  TA = "Chromic Luvisol",
  V  = "Albic Arenosol"
)
annotation_col <- data.frame(
  Soil = factor(soil_info[colnames(abund_mat)]),
  row.names = colnames(abund_mat)
)
soil_types  <- unique(annotation_col$Soil)
soil_colors <- setNames(
  c("#D95F02",
    "#890",
    "#1B9E77", 
    "#7570B3", 
    "#E7298A"

  ),
  soil_types
)

library(pheatmap)
#Plot the heatmap
pheatmap(
  abund_mat,
  color = viridis(100),       
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row   = anno_row,  
  annotation_col   = annotation_col,
  annotation_colors= list(
    Phylum = phy_colors,
    Soil   = soil_colors
  ),
  border_color = "grey80",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  scale = "row",  
)

#core family

core_frac_df <- data.frame(
  SampleID      = sample_names(ps_famf_core  ),
  TotalCounts   = sample_sums(ps_famf_core   ),
  CoreCounts    = sample_sums(prune_samples(sample_names(ps_famf_core   ), ps_corefam_final ))
) %>%
  mutate(CoreFraction = CoreCounts / TotalCounts)

# add metadata back
core_frac_df <- core_frac_df %>%
  left_join(
    as(sample_data(ps_famf_core   ), "data.frame") %>% 
      rownames_to_column("SampleID"),
    by = "SampleID"
  )

core_frac_df %>% 
  summarise(
    mean_core   = mean(CoreFraction),
    median_core = median(CoreFraction),
    min_core    = min(CoreFraction),
    max_core    = max(CoreFraction)
  )


###Network analysis ----
##Networking
metadata_draft <- sample_data(ps_merge_genus)
metadata_draft <- metadata_draft[,-6]
class(metadata_draft)
metadata_draft <- as.matrix(metadata_draft)

write.csv(metadata_draft, "metadata.csv")
?write.csv

##Invaded 
ps_inv_net <- subset_samples(ps_merge_genus, status == "Invaded")

ps_inv_filt <- filter_taxa(ps_inv_net, function(x) sum(x > 1) > 0.1 * length(x), TRUE)
otu_inv_005 <- as.data.frame(otu_table(ps_inv_filt))

write.table(otu_inv_005, "../data/networking/otu_table_inv_005.csv", sep = ",")


tax_inv <- as.data.frame(tax_table(ps_inv_filt))
write.table(tax_inv, "../data/networking/tax_table_inv_005.csv", sep = ",")

metadata <- read.csv("../data/metadata.csv",  sep=';')
metadata <- sapply(metadata, gsub, pattern = ",", replacement= ".")


meta_inv <- as.data.frame(metadata) %>% filter(status == "Invaded")
rownames(meta_inv) <- meta_inv$sample
meta_inv <- meta_inv [,-1]
meta_inv <- meta_inv [,-1]


write.table(meta_inv, "../data/networking/metadata_inv.csv", sep = ",")

table(colnames(otu_inv_005))

###Native
ps_nat_net <- subset_samples(ps_merge_genus, status == "Native")

ps_nat_filt <- filter_taxa(ps_nat_net, function(x) sum(x > 1) > 0.1 * length(x), TRUE)
otu_nat_005 <- as.data.frame(otu_table(ps_nat_filt))

write.table(otu_nat_005, "../data/networking/otu_table_nat_005.csv", sep = ",")

meta_nat <- as.data.frame(metadata) %>% filter(status == "Native")
rownames(meta_nat) <- meta_nat$sample
meta_nat <- meta_nat [,-1]
meta_nat <- meta_nat [,-1]

write.table (meta_nat, "../data/networking/metadata_nat.csv", sep = ",")


library(igraph)
library(tidyverse)
library(ggraph)
library(tidygraph)

# Invaded set
edges_inv <- read.table("../data/networking/net_inv.edgelist", header=TRUE, sep="\t")
tax_inv <- as.data.frame(tax_table(ps_inv_filt))

graph_inv <- graph_from_data_frame(d = edges_inv, directed = FALSE)

tax_inv$id <- rownames(tax_inv)
tax_inv <- tax_inv %>% rename(name = id)  
nodes_inv <- left_join(data.frame(name = V(graph_inv)$name), tax_inv, by = "name")

# Assign a phylum or class for coloring
V(graph_inv)$phylum <- nodes_inv$Phylum

##Convert weights to absolute value
E(graph_inv)$weight_dist <- 1 / abs(E(graph_inv)$weight)


##Calculate network stats
deg <- degree(graph_inv, mode = "all")
V(graph_inv)$degree <- deg
V(graph_inv)$betweenness <-betweenness(graph_inv, weights = E(graph_inv)$weight_dist, directed = FALSE)
V(graph_inv)$closeness <-closeness(graph_inv, weights = E(graph_inv)$weight_dist, mode = "all")

#Potential hub nodes
hub_nodes <- names(sort(deg, decreasing = TRUE))[1:30]
x <- sort(V(graph_inv)$degree, decreasing = TRUE)[1:10]
hub_nodes


# Extract their degree values
hub_data <- data.frame(
  OTU = hub_nodes,
  Degree = sort(deg, decreasing = TRUE)[1:10]
)


# Create stats dataframe
df_stats <- data.frame(
  Betweenness = V(graph_inv)$betweenness,
  Closeness = V(graph_inv)$closeness,
  Degree = V(graph_inv)$degree,
  Node = V(graph_inv)$name
)

# Z-score transformation to center and standardize values
df_scaled <- df_stats %>%
  mutate(across(c(Betweenness, Closeness, Degree), scale)) %>%
  pivot_longer(cols = c(Betweenness, Closeness, Degree),
               names_to = "Measure", values_to = "Value") %>%
  group_by(Measure) %>%
  mutate(Index = row_number())

# If you have a column like `sign` or `type`:
num_nodes <- gorder(graph_inv)  # Number of vertices
num_edges <- gsize(graph_inv)   # Number of edges

num_cooccur <- sum(E(graph_inv)$weight > 0 ) 
num_exclusion <- sum(E(graph_inv)$weight <0 )


network_label <- paste0("Total interactions: ", num_edges, "\n",
                        "Nodes: ", num_nodes)

network_label <- paste0("Total interaction: ", num_edges, "\n",
                        "Co-occurrence: ", num_cooccur, " (", round(100 * num_cooccur / num_edges), "%)", "\n",
                        "Mutual exclusion: ", num_exclusion, " (", round(100 * num_exclusion / num_edges), "%)")


library(cowplot)

V(graph_inv)$phylum <- ifelse(is.na(V(graph_inv)$phylum),
                              "Environmental node",
                              V(graph_inv)$phylum)
V(graph_inv)$phylum


table(nodes_inv$Phylum)

keep_graph_phyla <- c(
  "Environmental node",
  "Acidobacteriota",
  "Actinobacteriota",
  "Bacteroidota",
  "Firmicutes",
  "Myxococcota",
  "Proteobacteria",
  "Verrucomicrobiota",
  "Planctomycetota"
)

V(graph_inv)$phylum <- ifelse(
  V(graph_inv)$phylum %in% keep_graph_phyla,
  V(graph_inv)$phylum,
  "Other"
)

phy_colors <- c(
  "Environmental node" = "#CC79A7",  # magenta
  "Acidobacteriota" = "#0072B2",  # blue
  "Actinobacteriota" = "#D55E00",  # vermilion
  "Bacteroidota"= "#009E73",  # green
  "Firmicutes"= "#F0E442",  # yellow
  "Myxococcota"= "#56B4E9",  # sky blue
  "Proteobacteria" = "#E69F00",  # orange 
  "Verrucomicrobiota" = "#999999",  # gray
  "Planctomycetota" = "#000000", #black
  "Other" = "purple")


network_plot <- ggraph(graph_inv, layout = "stress") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(aes(color = phylum, size = degree(graph_inv))) +
  scale_color_manual(
    values = phy_colors
  )+
  theme_void() +
  draw_label(network_label, x = 0.8, y = 8.2 , hjust = 0, vjust = 1,
             fontface = "italic", size = 10)
  
network_plot

stats_plot <- # Plot styled like your figure
  ggplot(df_scaled, aes(x = Value, y = Index)) +
  geom_line(linewidth = 0.2, color = "gray50") +   # line connecting points
  geom_point(size = 0.6, color = "black") +        # small dots
  facet_wrap(~Measure, scales = "free_x", nrow = 1) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold", angle = 0, vjust = 0.5),
    strip.background = element_rect(fill = "gray90", color = NA)
  ) +
  labs(y = "All nodes")

hubs_plot<-ggplot(hub_data, aes(x = reorder(OTU, Degree), y = Degree)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(y = "Degree", title = "Top hub nodes")
hubs_plot

inv_plot <- (network_plot /stats_plot/ hubs_plot) + plot_layout(nrow = 3)
inv_plot

ggsave("../results/network_plot.pdf", plot = inv_plot, width = 16, height = 14)


# Zi-Pi plot, within and between module connectivity
communities <- cluster_fast_greedy(graph_inv, weights = E(graph_inv)$weight_dist)
mod_membership <- communities$membership
V(graph_inv)$module <- mod_membership
mod_inv <- modularity(communities)

calculate_zp <- function(graph, membership) {
  df <- data.frame(
    OTU = V(graph)$name,
    Module = membership,
    Degree = degree(graph)
  )
  
  # Calculate Zi (within-module degree z-score)
  df <- df %>%
    group_by(Module) %>%
    mutate(Zi = (Degree - mean(Degree)) / sd(Degree)) %>%
    ungroup()
  
  # Calculate Pi (among-module connectivity)
  A <- as_adjacency_matrix(graph, sparse = FALSE)
  mod_list <- split(df$OTU, df$Module)
  
  Pi <- sapply(df$OTU, function(node) {
    k_total <- sum(A[node, ] > 0)
    if (k_total == 0) return(0)
    k_i_m <- sapply(mod_list, function(members) sum(A[node, members] > 0))
    p_i <- k_i_m / k_total
    1 - sum(p_i^2)
  })
  
  df$Pi <- Pi
  return(df)
}

zp_inv <- calculate_zp(graph_inv, mod_membership)

zp_inv <- zp_inv %>%
  mutate(Role = case_when(
    Zi >= 2.5 & Pi >= 0.62 ~ "Network hub",
    Zi >= 2.5 & Pi < 0.62 ~ "Module hub",
    Zi < 2.5 & Pi >= 0.62 ~ "Connector",
    TRUE ~ "Peripheral"
  ))

zp_inv_plot <- ggplot(zp_inv, aes(x = Pi, y = Zi, color = Role)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0.62, linetype = "dashed") +
  geom_hline(yintercept = 2.5, linetype = "dashed") +
  scale_color_manual(values = c(
    "Peripheral" = "brown",
    "Connector" = "gold",
    "Module hub" = "purple",
    "Network hub" = "darkgreen" 
  )) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Among-module connectivity (Pi)",
    y = "Within-module connectivity (Zi)",
    title = "Zi-Pi Plot — Invasion Present"
  )


# Native set
edges_nat <- read.table("../data/networking/net_nat.edgelist", header=TRUE, sep="\t")
tax_nat <- as.data.frame(tax_table(ps_nat_filt))

graph_nat <- graph_from_data_frame(d = edges_nat, directed = FALSE)


# Join taxonomy
tax_nat$id <- rownames(tax_nat)
tax_nat <- tax_nat %>% rename(name = id) 
nodes_nat <- left_join(data.frame(name = V(graph_nat)$name), tax_nat, by = "name")

# Assign a phylum or class for coloring
V(graph_nat)$phylum <- nodes_nat$Phylum

##Convert weights to absolute value
E(graph_nat)$weight_dist <- 1 / abs(E(graph_nat)$weight)


##Calculate network stats
deg_nat <- degree(graph_nat, mode = "all")
V(graph_nat)$degree <- deg_nat
V(graph_nat)$betweenness <-betweenness(graph_nat, weights = E(graph_nat)$weight_dist, directed = TRUE)
V(graph_nat)$closeness <-closeness(graph_nat, weights = E(graph_nat)$weight_dist, mode = "all")

#Potential hub nodes
hub_nodes_nat <- names(sort(deg_nat, decreasing = TRUE))[1:30]
hub_nodes_nat

# Extract their degree values
hub_data_nat <- data.frame(
  OTU = hub_nodes_nat,
  Degree = sort(deg_nat, decreasing = TRUE)[1:10]
)

hubs_plot_nat<-ggplot(hub_data, aes(x = reorder(OTU, Degree), y = Degree)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(y = "Degree", title = "Top hub nodes")
hubs_plot_nat

# Create stats dataframe
df_stats_nat <- data.frame(
  Betweenness = V(graph_nat)$betweenness,
  Closeness = V(graph_nat)$closeness,
  Degree = V(graph_nat)$degree,
  Node = V(graph_nat)$name
)

# Z-score transformation to center and standardize values
df_scaled_nat <- df_stats_nat %>%
  mutate(across(c(Betweenness, Closeness, Degree), scale)) %>%
  pivot_longer(cols = c(Betweenness, Closeness, Degree),
               names_to = "Measure", values_to = "Value") %>%
  group_by(Measure) %>%
  mutate(Index = row_number())

num_nodes_n <- gorder(graph_nat)  # Number of vertices
num_edges_n <- gsize(graph_nat)   # Number of edges

num_cooccur_n <- sum(E(graph_nat)$weight > 0 )  
num_exclusion_n <- sum(E(graph_nat)$weight <0 )


network_label_n <- paste0("Total interactions: ", num_edges_n, "\n",
                        "Nodes: ", num_nodes_n)

network_label_n <- paste0("Total interaction: ", num_edges_n, "\n",
                        "Co-occurrence: ", num_cooccur_n, " (", round(100 * num_cooccur_n / num_edges_n), "%)", "\n",
                        "Mutual exclusion: ", num_exclusion_n, " (", round(100 * num_exclusion_n / num_edges_n), "%)")


E(graph_nat)$layout_weight <- abs(E(graph_nat)$weight)

V(graph_nat)$phylum <- ifelse(is.na(V(graph_nat)$phylum),
                              "Environmental node",
                              V(graph_inv)$phylum)

table(nodes_nat$Phylum)


V(graph_nat)$phylum <- ifelse(
  V(graph_nat)$phylum %in% keep_graph_phyla,
  V(graph_nat)$phylum,
  "Other"
)


network_plot_nat <- ggraph(graph_nat, layout = "stress",weight = layout_weight) +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(aes(color = phylum, size = degree(graph_nat))) +
  scale_color_manual(
    values = phy_colors
  )+
  theme_void() +
  draw_label(network_label_n, x = -18.2, y = -10.8 , hjust = 0, vjust = 1,
             fontface = "italic", size = 10)

network_plot_nat




communities_nat <- cluster_fast_greedy(graph_nat, weights = E(graph_nat)$weight_dist)
mod_membership_nat <- communities_nat$membership
V(graph_nat)$module <- mod_membership_nat
mod_nat <- modularity(communities_nat)

zp_nat <- calculate_zp(graph_nat, mod_membership_nat)


zp_nat <- zp_nat %>%
  mutate(Role = case_when(
    Zi >= 2.5 & Pi >= 0.62 ~ "Network hub",
    Zi >= 2.5 & Pi < 0.62 ~ "Module hub",
    Zi < 2.5 & Pi >= 0.62 ~ "Connector",
    TRUE ~ "Peripheral"
  ))

zp_nat_plot <- ggplot(zp_nat, aes(x = Pi, y = Zi, color = Role)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0.62, linetype = "dashed") +
  geom_hline(yintercept = 2.5, linetype = "dashed") +
  scale_color_manual(values = c(
    "Peripheral" = "brown",
    "Connector" = "gold",
    "Module hub" = "purple",
    "Network hub" = "darkgreen"
  )) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Among-module connectivity (Pi)",
    y = "Within-module connectivity (Zi)",
    title = "Zi-Pi Plot — Native"
  )

zp_inv <- zp_inv %>% mutate(Status = "Invaded")
zp_nat <- zp_nat %>% mutate(Status = "Native")

zp_all <- bind_rows(zp_inv, zp_nat)

ggplot(zp_all, aes(x = Pi, y = Zi, color = Role)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0.62, linetype = "dashed") +
  geom_hline(yintercept = 2.5, linetype = "dashed") +
  facet_wrap(~Status) +
  scale_color_manual(values = c(
    "Peripheral" = "brown",
    "Connector" = "gold",
    "Module hub" = "purple",
    "Network hub" = "darkgreen"
  )) +
  labs(
    x = "Among-module connectivity (Pi)",
    y = "Within-module connectivity (Zi)",
    color = "Topological role"
  )

summary_df <- tibble(
  Network = c("Invaded", "Native"),
  Nodes = c(vcount(graph_inv), vcount(graph_nat)),
  Edges = c(ecount(graph_inv), ecount(graph_nat)),
  AvgClustering = c(
    transitivity(graph_inv, type = "average"),
    transitivity(graph_nat, type = "average") ),
 AvgPathLength = c(
    mean_distance(graph_inv, directed = FALSE, weights = E(graph_inv)$weight_dist, unconnected = TRUE),
    mean_distance(graph_nat, directed = FALSE, weights = E(graph_nat)$weight_dist, unconnected = TRUE)),
  AvgDegree = c(mean(deg), mean(deg_nat)),
  Keystones = c(
    length(which(zp_inv$Zi > 2.5 | zp_inv$Pi > 0.62)),
    length(which(zp_nat$Zi > 2.5 | zp_nat$Pi > 0.62))
    ),
  Modules = c(length(unique(zp_inv$Module)), length(unique(zp_nat$Module))),
  PctCooccur = c(
    mean(E(graph_inv)$weight > 0 ) ,
    mean(E(graph_nat)$weight > 0 )
  )
)


deg_plot = ggplot(summary_df, aes(x = Network, y = AvgDegree)) +
  geom_col(fill = c( "#E41A1C", "#377EB8")) +
  labs(y = "Average degree") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())

key_plot =ggplot(summary_df, aes(x = Network, y = Keystones)) +
  geom_col(fill = c("#E41A1C", "#377EB8")) +
  labs(y = "Number of Keystones") +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())

mod_plot = ggplot(summary_df, aes(x = Network, y = Modules)) +
  geom_col(fill = c("#E41A1C", "#377EB8")) +
  labs(y = "Number of Modules") +
  theme_minimal()

avg_clust = ggplot(summary_df, aes(x = Network, y = AvgClustering)) +
  geom_col(fill = c("#E41A1C", "#377EB8")) +
  labs(y = "Average clustering coefficient") +
  theme_minimal()

  

# Create sign column for invasion
edges_inv <- edges_inv %>%
  mutate(sign = ifelse(weight > 0, "positive", "negative"))

# Create sign column for native
edges_nat <- edges_nat %>%
  mutate(sign = ifelse(weight > 0, "positive", "negative"))

edge_summary <- bind_rows(
  edges_inv %>% mutate(Network = "Native"),
  edges_nat %>% mutate(Network = "Invaded")
) %>%
  group_by(Network, sign) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Network) %>%
  mutate(Proportion = Count / sum(Count))

coocc_plot = ggplot(edge_summary, aes(x = Network, y = Proportion, fill = sign)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("positive" = "red", "negative" = "blue"),
                    labels = c("co-exclude", "co-occur")) +
  labs(y = "Proportion of link types", fill = "Link type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )
library(patchwork)

(deg_plot + key_plot + mod_plot + avg_clust)

simulate_robustness <- function(graph) {
  nodes <- V(graph)$name
  results <- tibble()

    
  for (i in 1:length(nodes)) {
    sampled <- sample(nodes, i)
    subgraph <- delete_vertices(graph, sampled)
    avg_deg <- if (gorder(subgraph) > 0) mean(degree(subgraph)) else NA
    nat_conn <- if (gorder(subgraph) > 1) eigen(as_adjacency_matrix(subgraph))$values[2] else NA
    results <- bind_rows(results, tibble(Removed = i, AvgDegree = avg_deg, NatConn = nat_conn))
  }
  
  return(results)
}

robust_inv <- simulate_robustness(graph_inv)
robust_inv$Network <- "Invaded"

robust_nat <- simulate_robustness(graph_nat)
robust_nat$Network <- "Native"

robust_deg <- bind_rows(robust_inv, robust_nat) 
  
robust_deg_plot <- ggplot(robust_deg, aes(x = Removed, y = AvgDegree, color = Network)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  scale_color_manual(
    values = c(
      "Invaded" = "#E41A1C", 
      "Native" = "#377EB8")) +
  labs(y = "Average degree") +
  theme_minimal()

robust_deg_plot


simulate_natural_connectivity <- function(graph, max_remove = 200, reps = 1) {
  results <- list()
  
  for (r in 1:reps) {
    g <- graph
    nodes <- sample(V(g)$name)  # random removal order
    conn_vals <- numeric()
    
    for (i in seq_len(min(max_remove, gorder(g)))) {
      g <- delete_vertices(g, nodes[i])
      if (gorder(g) > 1) {
        A <- as_adj(g, sparse = FALSE)
        eigen_vals <- eigen(A, only.values = TRUE)$values
        nat_conn <- log(mean(exp(eigen_vals)))
      } else {
        nat_conn <- NA
      }
      conn_vals[i] <- nat_conn
    }
    
    results[[r]] <- tibble(
      Step = 1:length(conn_vals),
      NaturalConnectivity = conn_vals,
      Replicate = r
    )
  }
  
  bind_rows(results)
}
robust_inv_c <- simulate_natural_connectivity(graph_inv, max_remove = 200, reps = 1) 
robust_inv_c$Network <- "Invaded"

robust_nat_c <- simulate_natural_connectivity(graph_nat, max_remove = 200, reps = 1)
robust_nat_c$Network <- "Native"

robust_con <- bind_rows(robust_inv_c, robust_nat_c)

robust_con_plot <- ggplot(robust_con, aes(x = Step, y = NaturalConnectivity, color = Network))  +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  scale_color_manual(
    values = c(
      "Invaded" = "#E41A1C", 
      "Native" = "#377EB8")) +
  labs(x = "Nodes removed", y = "Natural connectivity") +
  theme_minimal(base_size = 13) 

robust_con_plot


keystones_inv <- zp_inv %>%
  filter(Zi > 2.5 | Pi > 0.62)

keystone_otus_inv <- otu_inv_005[,colnames(otu_inv_005) %in% keystones_inv$OTU ]
keystone_tax_inv <- tax_inv[rownames(tax_inv) %in% keystones_inv$OTU, ]


keystone_tax_inv <- keystone_tax_inv[,-8]

write.csv(keystone_tax_inv, "../data/networking/keystone_tax_inv.csv")
write.csv(keystone_otus_inv, "../data/networking/keystone_otu_inv.csv")


keystones_nat <- zp_nat %>%
  filter(Zi > 2.5 | Pi > 0.62)

keystone_otus_nat <- otu_nat_005[,colnames(otu_nat_005) %in% keystones_nat$OTU]
keystone_tax_nat <- tax_nat[rownames(tax_nat) %in% keystones_nat$OTU, ]

keystone_tax_nat <- keystone_tax_nat[,-8]

write.csv(keystone_tax_nat, "../data/networking/keystone_tax_nat.csv")
write.csv(keystone_otus_nat, "../data/networking/keystone_otu_nat.csv")

env_tables <- function(g, weight_col = "weight") {
  # ensure vertex names exist
  if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
  
  # vertex attributes
  vdf <- data.frame(name = V(g)$name,
                    phylum = V(g)$phylum,
                    stringsAsFactors = FALSE)
  
  # edge list with endpoints by name (plus weights if present)
  edf <- as_data_frame(g, what = "edges")  # columns: from, to, weight (if exists)
  
  # join phyla of endpoints
  edf2 <- edf %>%
    left_join(vdf, by = c("from" = "name")) %>% rename(from_phylum = phylum) %>%
    left_join(vdf, by = c("to"   = "name")) %>% rename(to_phylum   = phylum)
  
  # keep only edges touching an environmental node
  env_edges <- edf2 %>%
    filter(from_phylum == "Environmental node" | to_phylum == "Environmental node") %>%
    mutate(
      EnvironmentalNode = ifelse(from_phylum == "Environmental node", from, to),
      Neighbor          = ifelse(from_phylum == "Environmental node", to,   from),
      NeighborPhylum    = ifelse(from_phylum == "Environmental node", to_phylum, from_phylum),
      Weight            = if (weight_col %in% names(edf2)) .data[[weight_col]] else NA_real_,
      Sign              = case_when(
        !is.na(Weight) & Weight > 0  ~ "positive",
        !is.na(Weight) & Weight < 0  ~ "negative",
        TRUE                         ~ "zero/NA")
    ) %>%
    select(EnvironmentalNode, Neighbor, NeighborPhylum, Weight, Sign)
  
  # degree per environmental node (counts of incident edges in this subtable)
  env_summary <- env_edges %>%
    group_by(EnvironmentalNode) %>%
    summarise(
      Degree   = n(),
      PosEdges = sum(Weight > 0, na.rm = TRUE),
      NegEdges = sum(Weight < 0, na.rm = TRUE),
      ZeroNA   = sum(is.na(Weight) | Weight == 0),
      .groups  = "drop"
    )
  
  list(long = env_edges, summary = env_summary)
}

inv_env <- env_tables(graph_inv)
env_long_inv    <- inv_env$long
env_summary_inv <- inv_env$summary

##Biomarkers random forest ----
ps_fam <- tax_glom(ps_merge_genus, taxrank="Family")
fam_names <- as.character(tax_table(ps_fam)[, "Family"])
unique_fam_names <- make.unique(fam_names)
taxa_names(ps_fam) <- unique_fam_names
ps_fam_comp <- ps_fam %>% microbiome::transform("compositional")
otu_ind <- as.data.frame(otu_table(ps_fam_comp))
keep_otus <- colSums(otu_ind[, -ncol(otu_ind)] > 0) >= 0.0005 * nrow(otu_ind)
otu_filt <- otu_ind[, c(keep_otus, TRUE)]

otu_filt$status <- as.factor(sample_data(ps_fam_comp)$status)
colnames(otu_filt) <- make.names(colnames(otu_filt))


library(randomForest)
set.seed(43)  # for reproducibility
rf_model <- randomForest(
  status ~ ., 
  data = otu_filt, 
  importance = TRUE, 
  ntree = 10000
)
print(rf_model)



importance_df <- as.data.frame(importance(rf_model))
importance_df$otu <- rownames(importance_df)

# Sort by MeanDecreaseGini
top_otu <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(16)

tax_ind <- as.data.frame(tax_table(ps_fam_comp))
tax_ind$otu <- rownames(tax_ind)
top_annotated <- left_join(top_otu, tax_ind, by = "otu")


library(forcats)
library(tidyverse)
# Extract data

otu_filt$sample <- rownames(otu_filt)
otu_long <- otu_filt %>%
  select(where(is.numeric), sample) %>%   # Keep only numeric + SampleID
  pivot_longer(cols = -sample, names_to = "Taxon", values_to = "Abundance")

colnames(otu_long) <- c("Sample", "Family", "Abundance")
# Get sample metadata
meta_ind <- as.data.frame(sample_data(ps_fam_comp))
# Merge metadata and abundance
merged_ind <- left_join(otu_long, meta_ind, by = "Sample")

# Group and summarize mean abundance
abund_summary <- merged_ind%>%
  group_by(Family, status) %>%
  summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = status,
              values_from = MeanAbundance)

colnames(abund_summary) <- c("Family", "Inv_Abundance", "Nat_Abundance")


# Join RF importance results with abundance summaries
rf_data <- left_join(top_annotated, abund_summary, by = "Family")

rf_data <- rf_data %>%
  arrange(MeanDecreaseAccuracy) %>%
  mutate(Family = factor(Family, levels = Family))

abundance_long <- rf_data %>%
  select(Family, Inv_Abundance, Nat_Abundance) %>%
  pivot_longer(cols = c(Inv_Abundance, Nat_Abundance),
               names_to = "status", values_to = "Abundance")


abundance_long$Family <- factor(abundance_long$Family, levels = levels(rf_data$Family))

abundance_long$status <- gsub("Inv_Abundance", "Invaded", abundance_long$status)
abundance_long$status <- gsub("Nat_Abundance", "Native", abundance_long$status)


inv_colors <- c("Invaded" = "#E41A1C", "Native" = "#377EB8" )
phy_colors <- c(
  "Actinobacteriota" = "#D55E00",  # vermilion
  "Bacteroidota"= "#009E73",  # green
  "Firmicutes"= "#F0E442",  # yellow
  "Myxococcota"= "#56B4E9",  # sky blue
  "Proteobacteria" = "#E69F00",  # orange 
  "Verrucomicrobiota" = "#999999",  # gray
  "Nitrospirota" = "pink" #black
)

# Left plot: Mean Decrease Accuracy
p1 <- ggplot(rf_data, aes(x = MeanDecreaseAccuracy, y = Family, color = Phylum)) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, yend = Family), size = 0.7) +
  geom_point(size = 4) +
  scale_color_manual (values =phy_colors)+
  labs(x = "Mean decrease accuracy", y = NULL, color = "Phylum") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "bottom")

p1

p2 <- ggplot(abundance_long, aes(x = Abundance, y = Family, fill = status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = inv_colors) +
  labs(x = "Abundance", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom")

p2
combined_plot <- p1 + p2 
combined_plot


rf_model$confusion[, 'class.error']


library(randomForest)
library(ggplot2)

# Get the families in descending importance
ordered_fams <- rf_data %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  pull(Family)

# 2) Prepare results storage
set.seed(42) 
max_k <- length(ordered_fams)         # or set to 50/100 if you only care about the top X
res <- tibble(
  k         = 1:max_k,
  OOB_Error = NA_real_
)

# 3) Loop over k using your original data (otu_filt),
#    whose response column is InvasionStatus
for(i in seq_len(max_k)) {
  fams_i <- ordered_fams[1:i]
  
  # build formula: InvasionStatus ~ fam1 + fam2 + ... + fam_i
  fml <- as.formula(
    paste("status ~", paste(fams_i, collapse = " + "))
  )
  
  # subset the training data
  df_i <- otu_filt %>%
    select(status, all_of(fams_i))
  
  # fit RF
  rf_mod <- randomForest(
    formula    = fml,
    data       = df_i,
    ntree      = 10000,
    importance = FALSE
  )
  
  # record the final OOB error
  res$OOB_Error[i] <- rf_mod$err.rate[rf_mod$ntree, "OOB"]
}

# 4) Plot the curve
ggplot(res, aes(x = k, y = OOB_Error * 100)) +
  geom_line(size = 1, color = "purple") +
  geom_point(size = 2) +
  geom_vline(xintercept = which.min(res$OOB_Error), 
             linetype = "dashed", color = "grey40") +

  scale_x_continuous(breaks = seq(0, max_k, by = 5)) +
  labs(
    x     = "Number of families in model",
    y     = "OOB error rate (%)",
  ) +
  theme_minimal(base_size = 14)

