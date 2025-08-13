library(dada2)
library(ShortRead)
library(Biostrings)

path <- "~/Documents/seqs/RodriguezC2019"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(basename(fnFs), function(x) sub("_1\\..*$", "", x))

sample.names

plotQualityProfile(fnRs[1:8])

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

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, truncLen=c(240,200),
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)


#Learn the Error Rates
errF <- learnErrors(fnFs, multithread=TRUE)
errR <- learnErrors(fnRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(fnFs, err=errF, multithread=TRUE)
dadaRs <- dada(fnRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, fnFs, dadaRs, fnRs, verbose=TRUE)
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



sample_df <- data.frame(
  Sample = c("S150", "S151", "S152","S153","S154","S155","S156","S157",
             "S158", "S159", "S160", "S161","S162","S163","S164","S165",
             "S166","S167","S168","S169","S178","S179", "S180", "S181",
             "S182", "S183","S184","S185","S186","S187","S190","S191",
             "S192","S193","S228","S229","S230","S231","S232","S233",
             "S234","S235","S236","S237","S248","S249","S250","S251",
             "S252","S253","S254","S255","S256","S257", "S268", "S269",
             "S270","S271","S272","S273","S274","S275","S276","S277"),
  status = c("Invaded", "Invaded",'Invaded','Invaded','Native','Native','Native','Native',
           'Invaded','Invaded','Invaded','Invaded','Native','Native','Native','Native',
           'Invaded','Invaded','Invaded','Invaded','Invaded','Invaded','Native','Native',
           'Invaded','Invaded','Invaded','Invaded','Native','Native','Native','Native',
           'Native','Native','Invaded','Invaded','Invaded','Invaded','Native','Native',
           'Native','Native','Invaded','Invaded','Invaded','Invaded','Invaded','Invaded',
           'Native','Native','Native','Native','Native','Native','Native','Native',
           'Invaded','Invaded','Invaded','Invaded','Native','Native','Native','Native'),
  code = c("TC", "TC",'TC','TC','TN','TN','TN','TN',
           'CIC','CIC','AC','AC','VN','VN','VN','VN',
           'VC','VC','VC','VC','CAC','CAC','CIN','CIN',
           'CAC','CAC','CIC','CIC','CIN','CIN','CAN','CAN',
           'CAN','CAN','LC','LC','AC','AC','AN','AN',
           'AN','AN','LC','LC','PC','PC','PC','PC',
           'LN','LN','LN','LN','PN','PN','SN','SN',
           'SC','SC','SC','SC','PN','PN','SN','SN'),
  location = c("T", "T",'T','T','T','T','T','T',
               'CI','CI','A','A','V','V','V','V',
               'V','V','V','V','CA','CA','CI','CI',
               'CA','CA','CI','CI','CI','CI','CA','CA',
               'CA','CA','L','L','A','A','A','A',
               'A','A','L','L','P','P','P','P',
               'L','L','L','L','P','P','S','S',
               'S','S','S','S','P','P','S','S') ,
  year  = c("2017", "2017",'2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017'))

rownames(sample_df) <- c("S150", "S151", "S152","S153","S154","S155","S156","S157",
                        "S158", "S159", "S160", "S161","S162","S163","S164","S165",
                        "S166","S167","S168","S169","S178","S179", "S180", "S181",
                        "S182", "S183","S184","S185","S186","S187","S190","S191",
                        "S192","S193","S228","S229","S230","S231","S232","S233",
                        "S234","S235","S236","S237","S248","S249","S250","S251",
                        "S252","S253","S254","S255","S256","S257", "S268", "S269",
                        "S270","S271","S272","S273","S274","S275","S276","S277")

rownames(seqtab.nochim) <- c("S150", "S151", "S152","S153","S154","S155","S156","S157",
                         "S158", "S159", "S160", "S161","S162","S163","S164","S165",
                         "S166","S167","S168","S169","S178","S179", "S180", "S181",
                         "S182", "S183","S184","S185","S186","S187","S190","S191",
                         "S192","S193","S228","S229","S230","S231","S232","S233",
                         "S234","S235","S236","S237","S248","S249","S250","S251",
                         "S252","S253","S254","S255","S256","S257", "S268", "S269",
                         "S270","S271","S272","S273","S274","S275","S276","S277")

##Build phylosq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa.species))


#rename asvs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

sample_data(ps)
saveRDS(ps, "ps_RC_16s.RDS")

save.image("~/Documents/seqs/RodriguezC2019/dada2_RodriguezC_16s.RData")
