# Washer_et_al_microglia_media

This is the code for the paper "Single-cell transcriptomics defines an improved, validated monoculture protocol for differentiation of human iPSCs to microglia". https://www.biorxiv.org/content/10.1101/2022.08.02.502447v1

The pipeline for the **initial processing of the data** (CRAMs are stored on ENA ERP140337) is described in the **Snakefile**, which has the following steps:

- 1. Convert CRAM files to FASTQ with [samtools](http://www.htslib.org/doc/samtools-fasta.html).
- 2. Map the data with [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
- 3. Prepare the genotype file for the donor deconvolution with [tabix](http://www.htslib.org/doc/tabix.html).
- 4. Run the deconvolution with [cellSNP-lite](https://cellsnp-lite.readthedocs.io/en/latest/) and [Vireo](https://vireosnp.readthedocs.io/en/latest/).

For more info on how to run the Snakefile with snakemake see [here](https://snakemake.readthedocs.io/en/stable/). The bash script in the deconvolution step (**deconvoluteCells.sh**) and the Snakefile have some hardcoded paths, so make sure you change those to your own paths and directory structure. The folder "for_deconvolution" contains the sample metadata ("sampleMetadata.txt") with the list of donors needed to run Vireo and replicate the results of the single cell analysis in the paper. An additional genotype file is needed, "merged.genotype.MAF0.01.vcf.gz", which can be downloaded from [zenodo](https://zenodo.org/record/7010323#.Yv-9nC8w1-U).

Once you have performed these steps, you can run the data analysis scripts:
- **1.1.Determine_doublets.R** to gather information on the doublets from Vireo.
- **1.2.Determine_cluster_parameters.R** to decide which parameters to choose for clustering in Seurat.
- **1.3.QC.R** to do the quality control and processing of the cellranger data to Seurat object, with filtering, normalization, clustering, etc.
- **2.Markers and DEA.R** contains the analysis of microlglial and perivascular macrophage markers and differential expression analysis between the different media samples. The perivascular macrophage markers were selected after following the pipeline described in the separate "select_perivascular_markers" folder.
- **3.SingleR_labelTransfer.R** describes the process of transfering labels from two different datasets onto our data.

There is also a separate folder for the **data analysis of qPCR results** ("for_qPCR), with the input data "Micro001_All_Markers.csv" and the script **Differentiation_1_qPCR_Linear_Regression_Analysis.R** to run the analysis.

At the end of these scripts you have the product of the `sessionInfo()` command, with details of the version of R and used packages.
