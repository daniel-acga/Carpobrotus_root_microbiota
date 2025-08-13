library(dada2)
library(ShortRead)
library(Biostrings)

path <- "~/Documents/seqs/Canini2024/its"
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(basename(fnFs), function(x) sub("_1\\..*$", "", x))
sample.names



#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


save.image("~/Documents/seqs/Canini2024/its/dada2_its_Canini.RData")

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
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

#Filter and trim
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path.cut, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
  filtRs <- file.path(path.cut, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                       minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
  head(out)
  
  save.image("~/Documents/seqs/Canini2024/its/dada2_its_Canini.RData")
  
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
  taxa <- assignTaxonomy(seqtab.nochim, '/Users/danielacosta/Documents/seqs/tax/sh_general_release_dynamic_04.04.2024_dev.fasta', 
                         multithread = TRUE, tryRC = TRUE, verbose=TRUE)
  
  save.image("~/Documents/seqs/Canini2024/its/dada2_its_Canini.RData")
  
  



sample_df <- data.frame(
  Sample = c('R432','R433','R434','R435','R436','R437','R438','R439','R440','R441',
             'R443','R444','R445','R446','R447','R448','R449','R450','R451','R452',
             'R456','R457','R458','R459','R460','R461','R462','R463', 'R464','R465',
             'R467','R468','R469','R470','R471','R472','R473','R474','R475', 'R476',
             'R481','R482','R483','R484','R485','R486','R487','R488','R489', 'R490',
             'R492','R493','R494','R495','R496','R497', 'R498','R499','R500','R501'),
  status = c('Invaded','Native','Native','Native','Invaded','Invaded','Invaded','Native','Native','Native',
             'Invaded','Invaded','Invaded','Native','Native','Native','Invaded','Invaded','Invaded','Invaded',
             'Invaded','Invaded','Invaded','Native','Native','Native','Invaded','Invaded','Invaded','Native',
             'Native','Native','Invaded','Invaded','Invaded','Native','Native','Native','Invaded','Invaded',
             'Invaded','Invaded','Invaded','Invaded','Invaded','Invaded','Invaded','Invaded','Native','Native',
             'Native','Native','Native','Native','Native','Native','Native','Native','Native','Native'),
  code = c('SAC','SAN','SAN','SAN','SAC','SAC','SAC','SAN','SAN','SAN',
           'SAC','SAC','SAC','SAN','SAN','SAN','TAC','TAC','TAC','TAC',
           'SSC','SSC','SSC','SSN','SSN','SSN','SSC','SSC','SSC','SSN',
           'SSN','SSN','SSC','SSC','SSC','SSN','SSN','SSN','SAC','SAC',
           'TAC','TAC','TAC','TAC','TAC','TAC','TAC','TAC','TAN','TAN',
           'TAN','TAN','TAN','TAN','TAN','TAN','TAN','TAN','TAN','TAN'),
  location =c('SA','SA','SA','SA','SA','SA','SA','SA','SA','SA',
              'SA','SA','SA','SA','SA','SA','TA','TA','TA','TA',
              'SS','SS','SS','SS','SS','SS','SS','SS','SS','SS',
              'SS','SS','SS','SS','SS','SS','SS','SS','SA','SA',
              'TA','TA','TA','TA','TA','TA','TA','TA','TA','TA',
              'TA','TA','TA','TA','TA','TA','TA','TA','TA','TA'),
  year  = c('2022','2022','2022','2022','2022','2022','2022','2022','2022','2022',
            '2022','2022','2022','2022','2022','2022','2022','2022','2022','2022',
            '2022','2022','2022','2022','2022','2022','2022','2022', '2022','2022',
            '2022','2022','2022','2022', '2022','2022','2022','2022','2022','2022',
            '2022','2022','2022','2022','2022','2022','2022','2022', '2022','2022',
            '2022','2022', '2022','2022','2022','2022','2022','2022','2022','2022'))

rownames(sample_df) <- c('R432','R433','R434','R435','R436','R437','R438','R439','R440','R441',
                         'R443','R444','R445','R446','R447','R448','R449','R450','R451','R452',
                         'R456','R457','R458','R459','R460','R461','R462','R463', 'R464','R465',
                         'R467','R468','R469','R470','R471','R472','R473','R474','R475', 'R476',
                         'R481','R482','R483','R484','R485','R486','R487','R488','R489', 'R490',
                         'R492','R493','R494','R495','R496','R497', 'R498','R499','R500','R501')

rownames(seqtab.nochim) <- c('R432','R433','R434','R435','R436','R437','R438','R439','R440','R441',
                              'R443','R444','R445','R446','R447','R448','R449','R450','R451','R452',
                              'R456','R457','R458','R459','R460','R461','R462','R463', 'R464','R465',
                              'R467','R468','R469','R470','R471','R472','R473','R474','R475', 'R476',
                              'R481','R482','R483','R484','R485','R486','R487','R488','R489', 'R490',
                              'R492','R493','R494','R495','R496','R497', 'R498','R499','R500','R501')

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



saveRDS(ps, "ps_Ca_its.RDS")

save.image("~/Documents/seqs/Canini2024/its/dada2_its_Canini.RData")
