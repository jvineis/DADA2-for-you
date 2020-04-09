# DADA2-for-you
Reinventing the wheel of DADA2 analysis for the Vineis mind and anyone else who thinks like me or is looking for another way to process amplicon data

### Run your data through the DADA2 pipeline - Much of the pipeline in the R scripts that I have here are specific to your samples and I stress that you get to know your data and not rely soely on DADA2 to make everything awesome.  Think about the amplicon size and what you expect for the merging of reads (how much overlap) what should be the size distribution of your merged reads etc... This tutorial assumes that you have demultiplexed Illumina paired end sequencing data. All the work below is carried out using R or in a bash shell terminal.

### So the first step is to look at your filenames.  Sometimes they might look like pairs of read files with a prefix filename followed by R1.fastq and R2.fastq.  Here are read pairs for two samples that we can focus on.  sample 1 = 15cm_S_Alt_2 and sample 2 = 35cm_S_Alt_10_3 

    15cm_S_Alt_2-R1.fastq
    15cm_S_Alt_2-R2.fastq
    35cm_S_Alt_10_3-R1.fastq
    35cm_S_Alt_10_3-R2.fastq

### Now begins the first steps of processing using R. load the required package

    library('dada2');packageVersion('dada2')

### sort files into R1 (forward) and R2 (reverse)

    fnFs <- sort(list.files(path, pattern="-R1.fastq", full.names = TRUE))
    fnRs <- sort(list.files(path, pattern="-R2.fastq", full.names = TRUE))

### You can now derive the sample names in your project based on the read file names. This will be very specific to your sample names.  The "-" character works for me because of the sample names I have thoughtfull named so that I can easily parse things

    sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

### Take a look at the quality plots for read 1 and read 2

    plotQualityProfile(fnFs[1:2])
    plotQualityProfile(fnRs[1:2])

### create the file names and paths that you can assign to the "filterAndTrim" command

    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names

### Now you can filter athe reads using the DADA2 command "filterAndTrim".  You should trim based on the quality plots that you have carefully examined above and also based on the lenght that you need to acheive merging of the paired reads.  If you trim too much, longer amplicons will fail to merge

    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

### Run the magical DADA2 error profiling - you might want to do this for each sequencing run that you have, especially if there are poor looking quality profiles.

    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
    dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

### Now you can build an ASV table
    
    seqtab <- makeSequenceTable(mergers)

### If you want to check the number of samples and ASVs you have, you can check that here.. but you don't have to
      
    dim(seqtab)
    table(nchar(getSequences(seqtab)))

### Remove the chimeric sequences and check the amount of seqs that were chimeric.  

    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    dim(seqtab.nochim)
    sum(seqtab.nochim)/sum(seqtab)

### Summary of the data that you can write to a file if you like to keep track of that sort of thing.. I find it very handy

    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

### If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    head(track)

    
