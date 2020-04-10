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
    
### An example of the output
    
    input filtered denoisedF denoisedR merged nonchim
    J3-3-control 36313    35145     34345     34575  27631    9487
    J3-3-Red     38858    37508     35930     36803  27930   10851

### In the case that you have distinct error structures between runs or for whatever reason, you want to run DADA2 separately for a separate set of data, start by creating a table of sequences.  skip this step if you don't need to run separate DADA2 analyses
    
    seqtab = makeSequenceTable(mergers)

### save the table as an .rds file - Then return to the top and restart a DADA2 analysis for another set of samples.

    saveRDS(seqtab, "~/Dropbox/TRAP-COMPETITION/COMPETITION-DILUTION-EXPERIMENTS/seqtab-comp-dilution.rds")

##### then I would go ahead and save the second DADA2 run after the analysis is complete for the second run and

    saveRDS(seqtab, "~/Dropbox/TRAP-COMPETITION/ORIGINAL-CULTIVATION/seqtab-original.rds")

### Here is an example of how to merge the two runs
    st1 = readRDS("~/Dropbox/TRAP-COMPETITION/COMPETITION-DILUTION-EXPERIMENTS/seqtab-comp-dilution.rds")
    st2 = readRDS("~/Dropbox/TRAP-COMPETITION/ORIGINAL-CULTIVATION/seqtab-original.rds")
    st.all = mergeSequenceTables(st1, st2)
    seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
    
### Whether or not you had a single or multiple DADA2 runs, now its time to write the ASV matrix containing only non-chimeric quality ASVs to a file

    write.table(seqtab, "~/Dropbox/TRAP-COMPETITION/seqtab-original-comp-dil-combined-runs.txt", sep = "\t")

### Then you need to fix the output using some type of txt editor so that the rows are sequences and the colums are samples..   Its just a matter of transpoosing the matrix and adding the "sample" as the column1 header. like this.

    sample 15cm_S_Alt_2_1  15cm_S_Alt_2_2  15cm_S_Alt_2_3
    ACTTGA 0  26  33
    ACTTAA 0  100 30
    AACTGG 25 33  26

### Save it as a text file called dada2-t-count-table.txt, then use this handy script
    
    python ~/scripts/create-fasta-from-seqs.py -n dada2-t-count-table.txt -fa dada2-fasta.fa -o DADA2-MATRIX.txt

##### which produces a count table with simple names instead of the asv as the identifier and a fasta file that you can use to identify the taxonomy

## I ran vsearch on this file to generate taxonmy like this in the terminal

    vsearch --usearch_global dada2-fasta.fa --db ~/scripts/databas/silva119.fa --blast6out NODE-HITS.txt --id 0.6
### Then I ran the wonderful script below in the terminal

    python ~/scripts/mu-dada2-phyloseq-creator.py -hits NODE-HITS.txt -tax_ref ~/scripts/databas/silva_fix.tax -dada2 DADA2-MATRIX.txt -fa dada2-fasta.fa

## which creates the PHYLOSEQ-TAX.txt and PHYLOSEQ-MATRIX.txt that you will use below
library(phyloseq)

mat = read.table("~/Dropbox/TRAP-COMPETITION/PHYLOSEQ-MATRIX.txt", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("~/Dropbox/TRAP-COMPETITION/PHYLOSEQ-TAX.txt", header = TRUE, sep = ";", row.names = 1)

mat = as.matrix(mat)
tax = as.matrix(tax)

OTU = otu_table(mat, taxa_are_rows = TRUE)
TAX = tax_table(tax)

physeq = phyloseq(OTU,TAX)

per = transform_sample_counts(physeq, function (x) x/sum(x)*100)

lowabundnames = filter_taxa(per, function(x) mean(x) > 0.1) ## This filters out anything that has less than a mean of 0.2% relative abundance acrocss all samples
per_abund = prune_taxa(lowabundnames, per)

### Do any other filtering that you think should be required here.  

### Write the taxa to a file, just because its handy to have
    
    write.table(tax_table(per_abund), "~/Dropbox/TRAP-COMPETITION/DADA2-TAX-meangreaterthanpoint1.txt", sep = '\t')

### Write the matrix for further processing - e.g. build a heirarchical tree anvi'o display etc..

    write.table(otu_table(per_abund), "~/Dropbox/TRAP-COMPETITION/DADA2-MATRIX-phyloseq-meangreaterthanpoint1.txt", sep = '\t')

### READ in a table that contains the percent relative abundance of the crosses experiment and export a heirarchical tree that can be used by Anvi'o display

    library(ape)
    library(vegan)

    dat = read.table("~/Dropbox/TRAP-COMPETITION/CROSSES-perabund-table.txt", header = TRUE, row.names = 1)
    dist = vegdist(dat, method = "bray")
    dclust = hclust(dist)
    write.tree(as.phylo(dclust), "~/Dropbox/TRAP-COMPETITION/CROSSES-perabund-table.tre")

    
