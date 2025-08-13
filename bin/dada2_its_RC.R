library(dada2)
library(ShortRead)
library(Biostrings)

path <- "~/Documents/seqs/RodriguezC2019/its"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(basename(fnFs), function(x) sub("_1\\..*$", "", x))

sample.names

##Identify primers
FWD <- "GTGAATCATCGAATCTTTGAA"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"

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
#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, '/Users/danielacosta/Documents/seqs/tax/sh_general_release_dynamic_04.04.2024_dev.fasta', 
                       multithread = TRUE, tryRC = TRUE, verbose=TRUE)
save.image("~/Documents/seqs/RodriguezC2019/its/dada2_RodriguezC_its.RData")


sample_df <- data.frame(
  Sample = c("S170", "S171", "S172", "S173", "S174", "S175","S176", "S177",
           "S188", "S189","S194","S195","S196","S197","S198","S199",
           "S200","S201","S202","S203","S204","S205","S206","S207",
           "S208","S209","S210","S211","S212","S213","S214","S215",
           "S216","S217", "S218", "S219","S220","S221","S222","S223",
           "S224","S225","S226","S227","S238","S239","S240","S241",
           "S242","S243","S244","S245","S246","S247","S258","S259",
           "S260","S261","S262","S263","S264","S265","S266","S267"),
  status = c("Invaded", "Invaded",'Invaded','Invaded','Native','Native','Native','Native',
           'Native','Native','Invaded','Invaded','Invaded','Invaded','Native','Native',
           'Invaded','Invaded',"Native", "Native", "Native", "Native", "Invaded", "Invaded",
           'Native','Native','Native','Native','Invaded','Invaded','Invaded','Invaded',
           'Invaded','Invaded','Invaded','Invaded','Native','Native','Native','Native',
           'Invaded','Invaded','Invaded','Invaded','Native','Native','Invaded','Invaded',
           'Invaded','Invaded','Native','Native','Native','Native','Native','Native',
           'Invaded','Invaded','Invaded','Invaded','Native','Native', 'Native','Native'),
  code = c("CAC", "CAC",'CAC','CAC','CAN','CAN','CAN','CAN',
           'VN','VN','VC','VC','VC','VC','AN','AN',
           'AC','AC',"VN", "VN", "AN", "AN", "AC", "AC",
           'LN','LN','LN','LN','LC','LC','LC','LC',
           'PC','PC','PC','PC','PN','PN','PN','PN',
           'SC','SC','SC','SC','TN','TN','CIC','CIC',
           'CIC','CIC','CIN','CIN','CIN','CIN','TN','TN',
           'TC','TC','TC','TC','SN','SN', 'SN','SN'),
  location = c("CA", "CA",'CA','CA','CA','CA','CA','CA',
           'V','V','V','V','V','V','A','A',
           'A','A',"V", "V","A","A", "A","A",
           'L','L','L','L','L','L','L','L',
           'P','P','P','P','P','P','P','P',
           'S','S','S','S','T','T','CI','CI',
           'CI','CI','CI','CI','CI','CI','T','T',
           'T','T','T','T','S','S', 'S','S'),
  year  = c("2017", "2017",'2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017',
            '2017','2017','2017','2017','2017','2017','2017','2017'))

rownames(sample_df) <- c("S170", "S171", "S172", "S173", "S174", "S175","S176", "S177",
                          "S188", "S189","S194","S195","S196","S197","S198","S199",
                          "S200","S201","S202","S203","S204","S205","S206","S207",
                          "S208","S209","S210","S211","S212","S213","S214","S215",
                          "S216","S217", "S218", "S219","S220","S221","S222","S223",
                          "S224","S225","S226","S227","S238","S239","S240","S241",
                          "S242","S243","S244","S245","S246","S247","S258","S259",
                          "S260","S261","S262","S263","S264","S265","S266","S267")

rownames(seqtab.nochim) <- c("S170", "S171", "S172", "S173", "S174", "S175","S176", "S177",
                              "S188", "S189","S194","S195","S196","S197","S198","S199",
                              "S200","S201","S202","S203","S204","S205","S206","S207",
                              "S208","S209","S210","S211","S212","S213","S214","S215",
                              "S216","S217", "S218", "S219","S220","S221","S222","S223",
                              "S224","S225","S226","S227","S238","S239","S240","S241",
                              "S242","S243","S244","S245","S246","S247","S258","S259",
                              "S260","S261","S262","S263","S264","S265","S266","S267")

##Build phylosq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa))

#rename asvs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))



sort(sample_sums(ps))

saveRDS(ps, "ps_RC_its.RDS")

save.image("~/Documents/seqs/RodriguezC2019/its/dada2_RodriguezC_its.RData")



