Prepare species log2 fold changes
================

This is a first stage to independently determine log2 fold changes for
each species.

First load in the sample information and provide the salmon counts files
directory. For this I am using files that I filtered for only the
*Candida* features (the original alignment was to a concatenated
mouse-fungal transcriptome):

``` r
library("tximport")
library("DESeq2")
library("ggplot2")
library("dplyr")

samples <- read.table("candida_sample_information.txt",sep="\t",header=T)
rownames(samples) <- samples$sample_ID
dir <- "salmon_outputs"
```

### Getting fold changes for *C. albicans*

I’m going to start with *Candida mealbicans* to show the full workflow
and because it’s done slightly different from the others due to the
SC5314 reference having a diploid genome.

First build a sample table with input file names:

``` r
species <- "Ca"
sample.IDs <- samples$sample_ID[samples$species == species] # Extracts sample IDs for C. albicans

files <- file.path(dir,paste(sample.IDs,"quant.no_mouse.sf",sep="/")) # Salmon output with all mouse features removed
names(files) <- sample.IDs
sampleTable <- samples[samples$species == species, c("species", "condition", "replicate")]
sampleTable$condition <- relevel(sampleTable$condition, "only")
sampleTable$file <- files
print(sampleTable)
```

    ##            species condition replicate
    ## Ca_only_1       Ca      only         1
    ## Ca_only_2       Ca      only         2
    ## Ca_mouse_1      Ca     mouse         1
    ## Ca_mouse_2      Ca     mouse         2
    ##                                                   file
    ## Ca_only_1   salmon_outputs/Ca_only_1/quant.no_mouse.sf
    ## Ca_only_2   salmon_outputs/Ca_only_2/quant.no_mouse.sf
    ## Ca_mouse_1 salmon_outputs/Ca_mouse_1/quant.no_mouse.sf
    ## Ca_mouse_2 salmon_outputs/Ca_mouse_2/quant.no_mouse.sf

Now that we have a list of files, we need to create a tx2gene data
frame. This is something that is required for the tximport-DESeq2
analysis - it is a list of transcripts with the corresponding parent
gene, allowing gene-level analysis from transcript-level counts. For
most reference genomes used here, there is only one annotated transcript
per gene, so these are the same thing. For *C. albicans*, however, the
diploid genome means there are two alleles for every gene. We treat
these simply as different transcripts of the same gene, allowing the
improved mapping associated with using a diploid genome while collapsing
the alleles into one gene to allow comparison to haploid genomes.

First, we get a list of gene IDs from the first salmon sample file.
Then, we take advantage of the fact that all alleles are designated
either "\_A" or "\_B", depending on their chromosome of origin (the
exception is mitochondrial genes, but these only have one allele):

``` r
# Get a list of gene IDs
tx_IDs <- read.table(sampleTable$file[1],sep="\t",header=T)$Name

# Make the tx2gene object
tx2gene <- data.frame(tx=tx_IDs, gene=gsub("-T", "", gsub("_A|_B", "", tx_IDs)))
print(head(tx2gene))
```

    ##            tx      gene
    ## 1 C1_00010W_A C1_00010W
    ## 2 C1_00020C_A C1_00020C
    ## 3 C1_00030C_A C1_00030C
    ## 4 C1_00040W_A C1_00040W
    ## 5 C1_00050C_A C1_00050C
    ## 6 C1_00060W_A C1_00060W

Now we are ready to build the DESeq data object for analysis. First we
create a tximport object using the tx2gene data frame, then run DESeq2
analysis to get fold changes. These are finally saved in a text file:

``` r
txi <- tximport(sampleTable$file,type="salmon",tx2gene=tx2gene)
```

    ## reading in files with read_tsv

    ## 1 2 3 4 
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
```

    ## using counts and average transcript lengths from tximport

``` r
dds <- DESeq(dds)
```

    ## estimating size factors
    ## using 'avgTxLength' from assays(dds), correcting for library size
    ## estimating dispersions
    ## gene-wise dispersion estimates
    ## mean-dispersion relationship
    ## final dispersion estimates
    ## fitting model and testing

``` r
res <- results(dds, contrast=c("condition", "mouse", "only"))

write.table(as.data.frame(res), paste("fold_changes/", species, "_pairwise.txt", sep=""),
            quote=F, sep="\t")
```

### Getting fold changes for other species

Now that this is done for C. albicans, let’s do it for the other
species, this time creating a function to do so:

``` r
species_analysis <- function(species) {
  # Creating a sample table file
  sample.IDs <- samples$sample_ID[samples$species == species]
  files <- file.path(dir,paste(sample.IDs,"quant.no_mouse.sf",sep="/"))
  names(files) <- sample.IDs
  sampleTable <- samples[samples$species == species, c("species", "condition", "replicate")]
  sampleTable$condition <- relevel(sampleTable$condition, "only")
  sampleTable$file <- files
  
  # Create the tx2gene data frame - note that for other species, gene and transcript IDs are the same
  tx_IDs <- read.table(sampleTable$file[1],sep="\t",header=T)$Name
  tx2gene <- data.frame(tx=tx_IDs, gene=tx_IDs)
  
  # Do analysis and save the result
  txi <- tximport(sampleTable$file,type="salmon",tx2gene=tx2gene)
  dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition", "mouse", "only"))

  write.table(as.data.frame(res), paste("fold_changes/", species, "_pairwise.txt", sep=""),
            quote=F, sep="\t")
}

for (species in c("Cd", "Ct", "Cp", "Le", "Cl")) {
  species_analysis(species)
}
```

This should have generated a series of text files for each species
comparing DESeq2 results for mouse vs only (i.e. *Candida* = mouse
macrophages vs *Candida* only).
