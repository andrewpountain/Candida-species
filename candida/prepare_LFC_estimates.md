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

I’m going to start with *Candida albicans* to show the full workflow and
because it’s done slightly different from the others due to the SC5314
reference having a diploid genome.

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

Finally, we want to combine all of the LFC estimates into a single file.
To do this, we need to use information from the CGOB (see
“estimating\_counts\_across\_species.md” for more information). First,
prepare the CGOB table:

``` r
species.list <- c("Ca","Cd","Cl","Cp","Ct","Le")

CGOB <- read.table("CGOB_headed.txt", sep="\t",header=T, na.strings="---", quote = "")
CGOB <- CGOB[, species.list]
CGOB <- na.omit(CGOB)
names(CGOB)[names(CGOB) == "Ca"] <- "Ca_orf19" # This is because I'm going to use the header Ca for Assembly 22 IDs

# Provide assembly 22 to ORF19 conversion 
A22_to_orf19 <- read.table("A22_to_orf19.txt", sep="\t",header=T)
CGOB <- merge(CGOB,A22_to_orf19[,c("ORF19_ID","ASSEMBLY22_ID")],by.x="Ca_orf19",by.y="ORF19_ID",all.x=F,all.y=F)
names(CGOB)[names(CGOB) == "ASSEMBLY22_ID"] <- "Ca"

# Strip off the "_A" suffix - note that there are no "_B"-suffixed alleles in the A22 to orf19 conversion set.
# The orf19 IDs refer to a haploid genome, so the conversion uses only the "_A" allele of the diploid Assembly 22.
CGOB$Ca <- gsub("_A", "", CGOB$Ca)

# Alter the naming of C. dubliniensis gene IDs to fit what you get in the data
CGOB$Cd <- gsub("CD36", "Cd36", CGOB$Cd)

# Now intersect this with only genes that made it into the shared quant.sf files
included_A22_IDs <- unique(read.table("tx2gene.txt", header=T, sep="\t")$gene)
CGOB <- CGOB[CGOB$Ca %in% included_A22_IDs,]
```

Then use this file to combine the various LFC estimates (along with
adjusted p-values and mean expression) by sequentially merging the
files:

``` r
combined.df <- CGOB
for (species in species.list) {
  species.pairwise <- read.table(
    paste("fold_changes/", species, "_pairwise.txt", sep=""), header=T, sep="\t")[, c("baseMean", "log2FoldChange", "padj")]
  names(species.pairwise) <- c(paste(species, "_baseMean", sep=""), paste(species, "_l2fc", sep=""), paste(species, "_padj", sep=""))
  species.pairwise$gene <- rownames(species.pairwise)
  combined.df <- merge(combined.df, species.pairwise, by.x=species, by.y="gene")
}

print(head(combined.df))
```

    ##           Le         Ct           Cp         Cl         Cd        Ca   Ca_orf19
    ## 1 LELG_00002 CTRG_01071 CPAR2_800090 CLUG_03997 Cd36_25140 CR_00120C orf19.7523
    ## 2 LELG_00003 CTRG_01072 CPAR2_800100 CLUG_03996 Cd36_25130 CR_00110W orf19.7527
    ## 3 LELG_00004 CTRG_00804 CPAR2_800110 CLUG_05204 Cd36_28020 CR_02950C orf19.2852
    ## 4 LELG_00005 CTRG_00805 CPAR2_800130 CLUG_01030 Cd36_28010 CR_02940C orf19.2851
    ## 5 LELG_00006 CTRG_00806 CPAR2_800140 CLUG_01031 Cd36_28000 CR_02930W orf19.2850
    ## 6 LELG_00009 CTRG_00809 CPAR2_800160 CLUG_01037 Cd36_27980 CR_02910W orf19.2848
    ##   Ca_baseMean      Ca_l2fc      Ca_padj Cd_baseMean     Cd_l2fc   Cd_padj
    ## 1   1891.0062  0.004777736 9.627793e-01   1345.3027 -0.03922999 0.9986782
    ## 2    577.4502  0.279841114 3.620065e-02    163.5682  0.24700093 0.9986782
    ## 3   1898.2891 -0.048546970 6.041740e-01   2269.2494 -0.64553183 0.9986782
    ## 4    263.2109 -0.759273361 2.163460e-05    479.9360 -0.30993712 0.9986782
    ## 5   1763.5277 -0.644140103 7.719753e-13   1285.6038 -0.13316578 0.9986782
    ## 6    575.0480  0.690173655 2.539763e-07    413.2560  0.12910351 0.9986782
    ##   Cl_baseMean     Cl_l2fc      Cl_padj Cp_baseMean     Cp_l2fc   Cp_padj
    ## 1  1915.80210 -0.93772582 2.882713e-08   1132.5526 -0.04832429 0.8822040
    ## 2   216.01044  0.26473482 5.475465e-01   1941.0503 -0.57536364 0.0575047
    ## 3   733.83194  0.33894449 2.161315e-01   1386.5771 -0.42954350 0.1435426
    ## 4    67.32928 -0.50820494 4.538545e-01    336.6131  0.01509386 0.9719963
    ## 5  1157.58481  0.78621558 3.164471e-06   1368.5260 -0.18948399 0.5726805
    ## 6   352.05427 -0.07990118 8.594494e-01    647.4683 -0.46274865 0.1371498
    ##   Ct_baseMean    Ct_l2fc   Ct_padj Le_baseMean     Le_l2fc   Le_padj
    ## 1   1430.2750  0.9125476 0.1221938   2936.5665  0.82674034 0.2877517
    ## 2    435.3381  0.6773683 0.3319344    582.8477 -0.22214553 0.9224055
    ## 3    513.2956 -0.5907639 0.5022867   2178.5545 -0.34353920 0.7945029
    ## 4    306.6037 -0.2688522 0.7381940    568.0312 -0.07428936 0.9594938
    ## 5    571.5184 -1.5619877 0.2406075   1649.1555 -0.58911052 0.5499439
    ## 6    492.0220 -0.1281868 0.9331312    111.3409 -0.19292794 0.9627521

``` r
write.table(combined.df, "combined_LFCs_padj.txt", row.names=F, quote=F, sep="\t")
```
