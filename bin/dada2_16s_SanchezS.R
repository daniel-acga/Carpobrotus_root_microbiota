library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)

path <- "~/Documents/seqs/SanchezS2023"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(basename(fnFs), function(x) sub("_1\\..*$", "", x))
sample.names


##Identify primers
FWD <- "CCTACGGGNGGCWGCAG" 
REV <- "GACTACHVGGGTATCTAATCC"

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


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = 
        sapply(FWD.orients,
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


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                        primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                   fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(path.cut, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,  minLen = 50,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
tail(out)

?plotQualityProfile
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


save.image("~/Documents/seqs/Canini2024/dada2_Canini_16s.RData")



sample_df <- data.frame(
  Sample = c("C1","C11","C12","C13","C15","C16","C18","C19","C21",
             "C22","C24","C25","C27","C28","C3","C30","C31","C32",
             "C33","C34","C35","C37","C39","C4","C40","C41","C43",
             "C44","C46","C47","C49","C50","C52","C6","C7","C9",
             "SS29", "SS30", "SS31", "SS32", "SS33", "SS34", "SS35",
             "SS36", "SS37", "SS38", "SS39", "SS40"),
  status = c("Native","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated","Eradicated",
             "Native","Native","Native","Native","Native","Native","Native","Native","Native",
             "Eradicated","Eradicated","Invaded","Invaded","Invaded","Native","Invaded","Invaded","Invaded",
             "Invaded","Invaded","Invaded","Invaded","Invaded","Invaded","Native","Eradicated","Eradicated",
             "Invaded", "Invaded", "Invaded", "Invaded", "Native", "Native", "Native",
             "Native", "Invaded", "Invaded", "Native", "Native"),
  location = c("G", "G", "G","G","G","G","G","G","G","G","G","G",
               "G","G","G","G","G","G","G","G","G","G","G","G",
               "G","G","G","G","G","G","G","G","G","G","G","G",
               "AT", "AT", "AT", "AT", "AT", "AT", "AT",
               "AT", "AT", "AT", "AT", "AT") ,
  code = c("GN", "GE", "GE", "GE", "GE", "GE", "GE", "GE", "GE",
           "GN", "GN", "GN", "GN", "GN", "GN", "GN", "GN", "GN",
           "GE", "GE", "GI", "GI", "GI", "GN", "GI","GI","GI",
           "GI","GI","GI","GI","GI","GI","GN","GE","GE",
           "AT", "AT", "AT", "AT", "AT", "AT", "AT",
           "AT", "AT", "AT", "AT", "AT"),
  name =  c("T130_new","T129_new", "T129_new", "T131M", "T131M", "T136_new", "T136_new", "T135", "T135",
            "T147_new", "T147_new", "T134", "T134", "T134", "T130_new", "T134", "T150", "T150",
            "T132_new", "T132_new", "82", "82", "80", "T133_new", "80","81_new","81_new",
            "INV_1_new","INV_1_new","INV_2_new","INV_2_new","84_new","84_new","T133_new","T148_new","T148_new",
            "AT", "AT", "AT", "AT", "AT", "AT", "AT",
            "AT", "AT", "AT", "AT", "AT"),
  year  = c("2023", "2023", "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023","2023",
            "2022", "2022", "2022","2022","2022","2022","2022","2022","2022","2022","2022","2022"))

rownames(sample_df) <- c("C1","C11","C12","C13","C15","C16","C18","C19","C21","C22","C24","C25","C27","C28","C3","C30","C31","C32","C33","C34","C35","C37","C39","C4","C40","C41","C43","C44","C46","C47","C49","C50","C52","C6","C7","C9",
                         "SS29", "SS30", "SS31", "SS32", "SS33", "SS34", "SS35",
                         "SS36", "SS37", "SS38", "SS39", "SS40")

rownames(seqtab.nochim) <- c("C1","C11","C12","C13","C15","C16","C18","C19","C21","C22","C24","C25","C27","C28","C3","C30","C31","C32","C33","C34","C35","C37","C39","C4","C40","C41","C43","C44","C46","C47","C49","C50","C52","C6","C7","C9",
                             "SS29", "SS30", "SS31", "SS32", "SS33", "SS34", "SS35",
                             "SS36", "SS37", "SS38", "SS39", "SS40")


#BUilding a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa.species))


#rename asvs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps_inv <- subset_samples(ps, status=='Invaded')
ps_SS <- subset_samples(ps_inv, location=='AT')

sample_data(ps_SS)

sort(sample_sums(ps_SS))
saveRDS(ps_SS, "ps_SS_16s.RDS")

save.image("~/Documents/seqs/SanchezS2023/dada2_16s_SanchezS.RData")

