#!/usr/bin/env python3
# Snakefile for processing media samples
import os
SAMPLE = ["IGBN","IM","IMBN","ITGBN","ITM_ADMEM","ITMG"]
RSAMPLE,= glob_wildcards("../../data/crams/{sample}.cram")

rule all:
	input:
		fastq_R2=expand("../../data/fastqs/{rsample}_R2_001.fastq.gz",rsample=RSAMPLE),
		ids=expand("../../data/{sample}/Vireo_donor_deconvolution/vireoOutput/{sample}/donor_ids.tsv",sample=SAMPLE)

rule cram_to_fastq:
        input:
                cram="../../data/crams/{rsample}.cram"
        output:
                fastq_R1="../../data/fastqs/{rsample}_R1_001.fastq.gz",
                fastq_R2="../../data/fastqs/{rsample}_R2_001.fastq.gz"
        message: "Convert cram to fastq"
        params:
                group="-G teamtrynka",
                queue="-q normal",
                threads="-n 16",
                memory="-M20000 -R'span[hosts=1] select[mem>20000] rusage[mem=20000]'",
                jobname= "-o ../../logs/log_cram_to_fastq.{rsample}.%J.%I",
                error="-e ../../errors/error_cram_to_fastq.{rsample}.%J.%I"
        shell:
                """
                mkdir -p ../../logs
                mkdir -p ../../errors
                echo "Working on ../../data/crams/{wildcards.rsample}.cram"
                /software/teamtrynka/conda/trynka-base/bin/samtools fastq -1 {output.fastq_R1} -2 {output.fastq_R2} -i {input.cram}
                """

rule cellranger:
	input:
		fastq_R1="../../data/fastqs/{sample}_S1_L002_R1_001.fastq.gz",
		fastq_R2="../../data/fastqs/{sample}_S1_L002_R2_001.fastq.gz"
	output:
		dummy="../../data/{sample}_finished"
	message: "Cellranger to process 10X data (no gRNAs). Place all fastqs together in fastqs folder first."
	params:
		group= "-G teamtrynka",
		queue="-q normal",
		threads="-n 32",
		memory="-M300000 -R'span[hosts=1] select[mem>300000] rusage[mem=300000]'",
		jobname= "-o ../../logs/log_{sample}_cellranger.%J.%I",
		error="-e ../../errors/error_{sample}_cellranger.%J.%I",
		transcriptome="/software/teamtrynka/cellranger/refdata-gex-GRCh38-2020-A"
	shell:
		"""
		# Need to run in directory that contains the fastq folder to avoid errors
		cd ../../data/
        samples={wildcards.sample}
        echo "processing ${{samples}}"

		# Run cellranger
		/software/teamtrynka/cellranger/cellranger-6.0.1/cellranger count \
        --id={wildcards.sample} \
        --sample=${{samples}} \
        --fastqs=./fastqs \
        --transcriptome={params.transcriptome} \
		--localmem=290

		touch {output.dummy}
		"""
rule preparing_genotype_file:
    input:
        genotype="../../data/genotypes/merged.genotype.MAF0.01.vcf.gz"
    output:
        index="../../data/genotypes/merged.genotype.MAF0.01.vcf.gz.tbi"
    message: "Indexing genotype file with tabix."
    params:
        group= "-G teamtrynka",
        queue="-q normal",
        threads="-n 20",
        memory="-M4000 -R'span[hosts=1] select[mem>4000] rusage[mem=4000]'",
        jobname= "-o ../../logs/log_preparing_genotype_file.%J.%I",
        error="-e ../../errors/error_preparing_genotype_file.%J.%I"
    shell:
        """
        echo "Generating index file with tabix"
        # Generating index file: also needed for Vireo to work
        /software/teamtrynka/conda/trynka-base/bin/tabix -fp vcf {input.genotype}
        """

rule deconvolution:
    input:
        bam="../../data/{sample}/outs/possorted_genome_bam.bam",
        genotype="../../data/genotypes/merged.genotype.MAF0.01.vcf.gz",
        metadata="../../data/sampleMetadata.txt",
        snps="../../data/genotypes/merged.genotype.MAF0.01.vcf.gz"
    output:
        ids = "../../data/{sample}/Vireo_donor_deconvolution/vireoOutput/{sample}/donor_ids.tsv"
    message: "Deconvoluting pooled donors using cellSNP-lite + vireo. "
    params:
        group= "-G teamtrynka",
        queue="-q normal",
        cores="-n 20",
        memory="-M4000 -R'span[hosts=1] select[mem>4000] rusage[mem=4000]'",
        jobname= "-o ../../logs/log_{sample}_deconv.%J.%I",
        error="-e ../../errors/error_{sample}_deconv.%J.%I",
        output_dir="../../data/{sample}/Vireo_donor_deconvolution/",
    shell:
        """

        echo "Running Vireo"
        ./deconvoluteCells.sh --out {params.output_dir} --genotypes {input.genotype} -s {input.snps} {wildcards.sample} {input.metadata}

        """
