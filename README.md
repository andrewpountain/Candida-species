# Candida-species

Code used for comparing *Candida* and macrophage responses across different *Candida* species. This is to accompany a manuscript in preparation.

Files are as follows:

- candida: Markdown documents and associated files for *Candida* RNA-seq
- macrophage: Markdown documents and associated files for macrophage RNA-seq
- quantify_species_mappying.py: a Python script for quantifying reads mapping to different species in a BAM file mapped to a concatenated reference

In addition, in each directory, I include the Snakefile for running the Snakemake pipelines for data pre-processing.

I have not included all the counts files used. These can be found from GEO with reference numbers GSE151288 (fungal datasets) and GSE152700 (macrophage datasets). To use these scripts with Salmon outputs found on GEO (i.e. for the *Candida* data), you will need to convert from a format of all counts files in the same directory <sample name>.quant.sf to separate directories, as <sample name>/quant.no_mouse.sf.
  
 Other files required for working with the data are included in the relevant directories.
