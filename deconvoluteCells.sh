#!/bin/bash
### BEGINNING OF GENOTYPE DECONVOLUTION SCRIPT ###
# AUTHOR: Eddie Cano-Gamez (ecg@sanger.ac.uk, https://github.com/eddiecg)
# Modified by Marta Perez Alcantara (ma23@sanger.ac.uk) to follow her directory structure and change cellSNP to cellSNP-lite
# Defining messages
printHelpMessage() {
	echo ""
	echo ":::::::::::::::::::::::::::: GENOTYPE DECONVOLUTION WITH CELLSNP AND VIREO ::::::::::::::::::::::::::::"
	echo ""
	echo "This script performs deconvolution of single cells by genotypes using the cellSNP and Vireo algorithms"
	echo "This script can either be run on a single sample or parallelized as a job array"
	echo "The script takes as input the name(s) of the sample(s) and a metadata table with donor numbers and donor IDs for each sample."
	echo "For more information, see usage or refer to the README.md file"
	echo ""
	echo "usage: deconvoluteCells.sh [OPTIONAL ARGUMENTS] sampleNames sampleMetadata"
	echo ""
	echo "Optional arguments:"
	echo "[-h | --help]          Display this help message"
	echo "[-s | --snplist]       List of common SNPs used for genotype calling with cellSNP (defaults to a file with hg38 SNPs from 1K genomes)"
	echo "[-g | --genotypes]     VCF file with reference genotypes, if available (defaults to no reference genotypes)"
	echo "[-t | --genotypeTag]   Format of genotype tags, if reference genotypes are available (i.e. GT, DP, PL; defaults to GT)"
	echo "[-o | --out]			 Output directory where Vireo results should be stored (defaults to the current working directory)"
	echo "[-p | --parallelized]  If this flag is added, the script will run for a group of samples in parallel as a job array"
	echo ""
	echo "Positional arguments (REQUIRED):"
	echo "sampleNames         Name of the sample to process (e.g. I1484). If running in --parallelized mode, then a text file with a list of sample names (one per line)"
	echo "sampleMetadata	  Text file containing three columns in the following order: sample IDs, number of donors per sample and sample composition (i.e. donor IDs present in each sample, separated by semicolons)"
	echo ""
	echo "IMPORTANT: Always provide optional arguments before positional arguments"
	echo ""
}

printErrorMessage() {
	echo ""
  	echo "[deconvoluteCells]: ERROR: $1"
  	echo "Please see program usage (--help)"
  	echo ""
}

# Setting help message as default program call
if [ "$#" -eq 0 ]
then
    printHelpMessage
    exit 1
fi

# Setting arguments to default values
snplist="/lustre/scratch117/cellgen/teamtrynka/eddie/resources/1Kgenomes/grch38/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
genotypes=""
genotypeTag="GT"
outputPath="$PWD"
parallel=false

# Fetching optional arguments
# TO DO: Add an option to allow keeping intermediary files
for arg in "$@"
do
    case $arg in
        -h|--help)
        printHelpMessage
        exit 1
        ;;
        -s|--snplist)
        snplist="$2"
        shift # Remove --snplist from argument list
        shift # Remove argument value from argument list
        ;;
        -g|--genotypes)
        genotypes="$2"
        shift # Remove --genotypes from argument list
        shift # Remove argument value from argument list
        ;;
        -t|--genotypeTag)
        genotypeTag="$2"
        shift # Remove --genotypes from argument list
        shift # Remove argument value from argument list
        ;;
        -o|--out)
        outputPath="$2"
        shift # Remove --out from argument list
        shift # Remove argument value from argument list
        ;;
        -p|--parallelized)
        parallel=true
        shift # Remove --parallelized from argument list
        ;;
    esac
done

# Fetching positional arguments
echo ""
echo "[deconvoluteCells]: Running deconvoluteCells script..."

if [ "$#" -lt 2 ]
then
  printErrorMessage "Some arguments are missing"
  exit 2
elif [ "$#" -gt 2 ]
then
  printErrorMessage "Too many input arguments"
  exit 2
fi


if [[ $parallel == true ]]
then
  echo "[deconvoluteCells]: Running script in parallelized mode..."
  echo "[deconvoluteCells]: Fetching list of input samples..."
  sampleList=($(<"$1"))
  echo "[deconvoluteCells]: ${#sampleList[@]} samples found on list..."

  sampleID="${sampleList[$((LSB_JOBINDEX-1))]}"
  echo "[deconvoluteCells]: Processing sample ${sampleID}..."
  sampleMetadata="$2"
elif [[ $parallel == false ]]
then
  echo "[deconvoluteCells]: Running script in standard mode..."
  sampleID="$1"
  echo "[deconvoluteCells]: Processing sample ${sampleID}..."
  sampleMetadata="$2"
fi


# Verifying that input files for cellSNP exist
bamFile="../../data/${sampleID}/outs/possorted_genome_bam.bam"
barcodeFile="../../data/${sampleID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

echo "[deconvoluteCells]: Fetching BAM file..."
if [[ -f $bamFile ]]
then
	:
else
	printErrorMessage "BAM file does not exist [${sampleID}]"
	exit 1
fi

echo "[deconvoluteCells]: Fetching 10X cell barcodes..."
if [[ -f $barcodeFile ]]
then
	:
else
	printErrorMessage "Barcode file does not exist [${sampleID}]"
	exit 1
fi

echo "[deconvoluteCells]: Fetching SNP list..."
if [[ -f $snplist ]]
then
	:
else
	printErrorMessage "SNP list file does not exist"
	exit 1
fi

# Creating temporary directory for intermediary files
temporaryDirectory="_temp_${sampleID}"
mkdir -p $temporaryDirectory

# Running cellSNP
cellsnpOutput="${temporaryDirectory}"

echo "[deconvoluteCells]: Using cellSNP-lite to call genotypes from 10X reads..."
echo "[deconvoluteCells]: Non default parameters on this session: all default (minMAF 0.01, minCount 20, minLen 30)"

/software/teamtrynka/conda/cellSNP-lite/bin/cellsnp-lite \
    -s $bamFile \
    -b $barcodeFile \
    -O $cellsnpOutput \
    -R $snplist \
    -p 20 \
    --minMAF 0.01 \
    --minCOUNT 20 \
		--genotype \
		--minLEN 30 \
		--gzip

# Verifying that input files for Vireo exist
cellsnpVCF="${cellsnpOutput}/cellSNP.cells.vcf.gz"

echo "[deconvoluteCells]: Fetching cellSNP output..."
if [[ -f $cellsnpVCF ]]
then
	:
else
	printErrorMessage "cellSNP VCF file does not exist [${sampleID}]"
	echo ""
	echo "[deconvoluteCells]: Stopping now..."
	exit 1
fi

echo "[deconvoluteCells]: Fetching metadata file..."
if [[ -f $sampleMetadata ]]
then
	:
else
	printErrorMessage "Metadata file does not exist"
	echo ""
	echo "[deconvoluteCells]: Stopping now..."
	exit 1
fi

# Creating output directory for Vireo
echo "[deconvoluteCells]: Creating Vireo output directory if it does not exist..."
if [[ -d $outputPath ]]
then
	:
else
	printErrorMessage "The specified output directory does not exist"
	echo ""
	echo "[deconvoluteCells]: Stopping now..."
	exit 1
fi

outputPathWithoutSlash=$(echo "$outputPath" | sed 's/\/$//g')
vireoDirectory="${outputPathWithoutSlash}/vireoOutput"

mkdir -p $vireoDirectory

outputDirectory="${vireoDirectory}/${sampleID}"

# Running deconvolution without reference genotypes
if [[ -z "$genotypes" ]]
then
	echo "[deconvoluteCells]: Running deconvolution without reference genotypes..."

	# Identifying donors in sample pool
	echo "[deconvoluteCells]: Fetching number of donors in input sample..."
	numberOfDonors=$(cat "$sampleMetadata" | grep -w "$sampleID" | awk '{print $2}')
	echo "[deconvoluteCells]: ${numberOfDonors} donors found in sample..."

	# Running Vireo
	echo "[deconvoluteCells]: Using Vireo to deconvolute cells..."
	/software/teamtrynka/conda/trynka-base/bin/vireo \
		-c $cellsnpVCF \
		-N $numberOfDonors \
		-o $outputDirectory

# Running with reference genotypes
else
	echo "[deconvoluteCells]: Running deconvolution with reference genotypes..."

	echo "[deconvoluteCells]: Fetching reference genotypes..."
	if [[ -f $genotypes ]]
	then
		:
	else
		printErrorMessage "Reference genotypes file does not exist"
		echo ""
		echo "[deconvoluteCells]: Stopping now..."
		exit 1
	fi

	# Identifying donors in sample pool
	detectedDonors="${temporaryDirectory}/${sampleID}_donorNames.txt"

	echo "[deconvoluteCells]: Finding donors present in the input sample..."
	cat "$sampleMetadata" | grep -w "$sampleID" | awk '{print $3}' | tr ';' '\n' > $detectedDonors

	# Running bcftools
	echo "[deconvoluteCells]: Using bcftools to subset reference genotypes..."

	subsetFile="${temporaryDirectory}/${sampleID}_subset.vcf.gz"

	/software/teamtrynka/conda/trynka-base/bin/bcftools view \
		$genotypes \
		-R $cellsnpVCF \
		-S $detectedDonors \
		-Oz \
		-o $subsetFile

	# Running Vireo
	echo "[deconvoluteCells]: Using Vireo to deconvolute cells..."

	/software/teamtrynka/conda/trynka-base/bin/vireo \
		-c $cellsnpVCF \
		-d $subsetFile \
		-t $genotypeTag \
		-o $outputDirectory

fi

# Removing intermediary files
echo "[deconvoluteCells]: NOT cleaning up..."
#rm -r $temporaryDirectory


echo "[deconvoluteCells]: ...DONE!"
echo ""

### END OF GENOTYPE DECONVOLUTION SCRIPT ###
