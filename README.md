# FinalTask-EPIGENOMICS
EN‐TEx ATAC‐seq data: downstream analyses AND Distal regulatory activity


#### 4. EN‐TEx ATAC‐seq data: downstream analyses ####


1. Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq. 

1.1. Run the docker
$ sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course


1.2. Clone the epigenomics repository and move to the ATAC-seq folder.
$ git clone https://github.com/bborsari/epigenomics_uvic

$ cd epigenomics_uvic/ATAC-seq

1.3. Go to ENCODE (ENTEX) to look for the metadata and copy it to the ATAC-seq folder (where we will do the analysis).
$ cp /mnt/c/Users/asole/Downloads/files.txt /home/solcos/chipseq/epigenomics_uvic/ATAC-seq/files.txt


1.4. Move to folder ATAC-seq, download the files from the experiment and check metadata.tsv (may need to exit the docker or not).
$ ./bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assay_slims=DNA+accessibility&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assembly=GRCh38&assay_title=ATAC-seq" 

# Check metadata file
$ head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'

1.5. Then create folders to store the different files. 

# Create the analyses and data folders
$ mkdir analyses data

# Create bigBED files folder (only bigBED, because in bigBED we have coordinates peaks (that is what we need) and bigWigs we have coordinates, pvalues, logFC...)
$ mkdir data/bigBed.files

#Create peaks analyses files folder
$ mkdir analyses/peak.analysis


2. Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.

2.1. Extract the necessary information from the metadata.tsv to create bigBED file with peaks from ATAC-seq experiments. 
# bigBed peak calling files (bigBed narrow, pseudoreplicated peaks, assembly GRCh38, 
most recent file for each tissue): (we do not have the target (H3K4me3) to search, we are in ATAC-seq data).
$ grep -F "bigBed_narrowPeak" metadata.tsv |grep -F "pseudoreplicated_peaks" |grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |sort -k2,2 -k1,1r |sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

2.2. Download them
$ cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

2.3. Check the integrity of the downloaded files.
$ for file_type in bigBed; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done

2.4. Convert the files to BED files.
$ mkdir data/bed.files

$ cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done


3. For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions). 

3.1. Create an annotation folder 
$ mkdir annotation

3.3.1. Download and move from a link the list of promoters ([-2 kb, +2 Kb] from TSS) of protein-coding genes. Store this file inside the annotation folder.

# List of promoters
$ mv /mnt/c/Users/asole/Downloads/gencode.v24.protein.coding.non.redundant.TSS.bed /home/solcos/atacseq/epigenomics_uvic/ATAC-seq/annotation

3.3.2. Download and move the full primary assembly annotation. Store this file inside the annotation folder.
# Assembly annotation
$ mv /home/asole/Downloads/gencode.v24.primary_assembly.annotation.gtf.gz /home/solcos/atacseq/epigenomics_uvic/ATAC-seq/annotation
# Unzip the file
$ gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz

3.4. The number of peaks that intersect promoter regions
$ cut -f-2 analyses/bigBed.peaks.ids.txt |while read filename tissue; do    bedtools intersect -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -a data/bed.files/"$filename".bed -u |  cut -f7 |  sort -u > analyses/peak.analysis/genes.with.peaks."$tissue".promoter.txt; done

# Number of peaks for sigmoid colon in promoter regions:
$ wc -l analyses/peak.analysis/genes.with.peaks.sigmoid_colon.promoter.txt

	## 38071 

# Number of peaks for stomach in promoter regions:
$ wc -l analyses/peak.analysis/genes.with.peaks.stomach.promoter.txt

	## 33169

3.5. The number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions).

3.5.1. Retrieve the protein coding genes for getting the gene coordinates.
$ awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed

3.5.2. Retrieve peaks that fall outside gene coordinates.
$ cut -f-2 analyses/bigBed.peaks.ids.txt |while read filename tissue; do    bedtools intersect -b annotation/gencode.v24.protein.coding.gene.body.bed -a data/bed.files/"$filename".bed -u |  cut -f7 |  sort -u > analyses/peak.analysis/genes.with.peaks."$tissue".outside.txt; done

# Number of peaks for sigmoid colon outside gene coordinates:
$ wc -l analyses/peak.analysis/genes.with.peaks.sigmoid_colon.outside.txt

	## 45568

# Number of peaks for stomach outside gene coordinates:
$ wc -l analyses/peak.analysis/genes.with.peaks.stomach.outside.txt
	
	## 38526

##########################################################################################

#### 5. Distal regulatory activity ####

Task 1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.
$ mkdir regulatory_elements
$ cd regulatory_elements

# Create all the folders for further steps
$ mkdir analyses
$ mkdir analyses/peaks.analyses
$ mkdir annotation
$ mkdir data
$ mkdir data/bed.files
$ mkdir data/bigBed.files

Task 2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

2.1. Download the files (bigBed files) containing the peaks for H3K4me1 and H3K27ac using the metadata obtained in the previous tasks. 

# Copy the metadata.tsv file to the folder regulatory_elements to work better.
$ cp ../ChIP-seq/metadata.tsv ../regulatory_elements/

# Grep the data filename that we want to download (H3K4me1)
$ grep -F "H3K4me1" metadata.tsv | grep -F "bigBed_narrowPeak" |grep -F "pseudoreplicated_peaks" | grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1,$11}' | sort -k2,2 -k1,1r |sort -k2,2 -u > analyses/bigBed.peaks.H3K4me1.txt

# Download it
$ cut -f1 analyses/bigBed.peaks.H3K4me1.txt |while read filename; do   wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; done

# Grep the data filename that we want to download (H3K27ac)
$ grep -F "H3K27ac" metadata.tsv | grep -F "bigBed_narrowPeak" |grep -F "pseudoreplicated_peaks" | grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1,$11}' | sort -k2,2 -k1,1r |sort -k2,2 -u > analyses/bigBed.peaks.H3K27ac.txt

# Download it
$ cut -f1 analyses/bigBed.peaks.H3K27ac.txt |while read filename; do   wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; done

2.2. Check the integrity of the downloaded files in order to verify their MD5 hash.

# For H3K4me1 file
$ for file_type in bigBed; do   ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".peaks.H3K4me1.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt;    cat data/"$file_type".files/md5sum.txt |  while read filename original_md5sum; do      md5sum data/"$file_type".files/"$filename"."$file_type" |    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' ;   done > tmp ;   mv tmp data/"$file_type".files/md5sum.txt;    awk '$2!=$3' data/"$file_type".files/md5sum.txt;  done

# For H3K27ac file
$ for file_type in bigBed; do   ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".peaks.H3K27ac.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt;    cat data/"$file_type".files/md5sum.txt |  while read filename original_md5sum; do      md5sum data/"$file_type".files/"$filename"."$file_type" |    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' ;   done > tmp ;   mv tmp data/"$file_type".files/md5sum.txt;    awk '$2!=$3' data/"$file_type".files/md5sum.txt;  done

2.3. Convert bigBed files to bed files.

# For H3K4me1 file
$ cut -f1 analyses/bigBed.peaks.H3K4me1.txt |while read filename; do   bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed; done

# For H3K27ac file
$ cut -f1 analyses/bigBed.peaks.H3K27ac.txt |while read filename; do   bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed; done

2.4. Retrieve those peaks (from H3K4me1 or H3K27ac) that are outside the promoter region.

# We need to perform an intersection between the peaks outside the coordinates of the genes (obtained in the ATAC-seq exercise (Task 4)) and the peaks of each mark in both tissues.

# Knowing which filename corresponds to each tissue:

# For Stomach tissue
$ bedtools intersect -a ../ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed -b data/bed.files/ENCFF977LBD.bed data/bed.files/ENCFF844XRN.bed -u |sort -u > analyses/peaks.analyses/genes.peaks.stomach.txt

$ wc -l analyses/peaks.analyses/genes.peaks.stomach.txt

	# 15905

# For Sigmoid Colon tissue
$ bedtools intersect -a ../ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed -b data/bed.files/ENCFF724ZOF.bed data/bed.files/ENCFF872UHN.bed -u |sort -u > analyses/peaks.analyses/genes.peaks.sigmoid_colon.txt

$ wc -l analyses/peaks.analyses/genes.peaks.sigmoid_colon.txt

	# 16036


Task 3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

# For Stomach tissue
$ grep -F "chr1" analyses/peaks.analyses/genes.peaks.stomach.txt | awk 'BEGIN{FS=OFS="\t"} {print $4,$2}' > analyses/peaks.analyses/starts.stomach.tsv

$ wc -l analyses/peaks.analyses/starts.stomach.tsv

	# 8439

# For Sigmoid colon tissue
$ grep -F "chr1" analyses/peaks.analyses/genes.peaks.sigmoid_colon.txt | awk 'BEGIN{FS=OFS="\t"} {print $4,$2}' > analyses/peaks.analyses/starts.sigmoid_colon.tsv

$ wc -l analyses/peaks.analyses/starts.sigmoid_colon.tsv

	# 8465

Task 4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point:

$ awk 'BEGIN{FS=OFS="\t"} $1=="chr1" { if ($6=="+"){start=$2} else {start=$3}; print $4, start}' ../ChIP-seq/gencode.v24.protein.coding.gene.body.bed > analyses/genes.starts.tsv

$ wc -l analyses/genes.starts.tsv

	# 2047

Task 5: Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

# Scripts modified

#!/usr/bin/env python


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)


#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
        gene, y = line.strip().split('\t') # split the line into two columns based on a tab
        # define a variable called position that correspond to the integer of the start of the gene
        position= int(y)
        # compute the absolute value of the difference between position and enhancer_start
        absolut_pos = abs(position - enhancer_start)
	# if this absolute value is lower than x
        if absolut_pos < x:
                # this value will now be your current x
                x = absolut_pos
                # save gene as selectedGene
                selectedGene = gene
                # save position as selectedGeneStart
                selectedGeneStart = position

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])


# Make sure that the Script is working
$ python ../bin/get.distances.py --input analyses/genes.starts.tsv --start 980000

	# ENSG00000187642.9       982093  2093


Task 6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:

# For stomach tissue
$ cat analyses/peaks.analyses/starts.stomach.tsv | while read element start; do python ../bin/get.distances.py --input analyses/genes.starts.tsv --start $start > analyses/peaks.analyses/regulatoryElements.genes.distances.stomach.tsv done

# For sigmoid colon tissue
$ cat analyses/peaks.analyses/starts.sigmoid_colon.tsv | while read element start; do python ../bin/get.distances.py --input analyses/genes.starts.tsv --start $start > analyses/peaks.analyses/regulatoryElements.genes.distances.sigmoid_colon.tsv done


Task 7: Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.

# First we open R and load the files

# For Stomach tissue
$ stomach <- read.table("analyses/peaks.analyses/regulatoryElements.genes.distances.stomach.tsv", sep="\t")

# Mean
$ mean(stomach[,2])

	# 155141.2

# Median
$ median(stomach[,2])

	# 51573

# For Sigmoid colon tissue
$ sigmoid_colon <- read.table("analyses/peaks.analyses/regulatoryElements.genes.distances.sigmoid_colon.tsv", sep="\t")

# Mean
$ mean(sigmoid_colon[,2])

	# 142698.2

# Median
$ median(sigmoid_colon[,2])

	# 47592
