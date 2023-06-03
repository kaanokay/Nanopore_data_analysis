#!/bin/bash

# To count total number of reads with indels greater than 2 bases in target site of CRISPR, we need to use samtools mpileup output as input in script below: --- start

#!/bin/bash

filename="Kmt2a_1_aligned_mapping_quality_30_filtered_400bp_read_length_filtered_mpileup_output_2.txt"
count=0

while IFS=$'\t' read -r -a columns; do
    fifth_column="${columns[4]}"
    cleaned_column=$(sed 's/[^0-9+2GA]//g' <<< "$fifth_column")

    for ((i=0; i<${#cleaned_column}; i++)); do
        char="${cleaned_column:i:1}"
        if [[ "$char" =~ [3-9] ]]; then
            count=$((count + 1))
        fi
    done
done < "$filename"

echo "Total number of reads with indels greater than 2 bases: $count"

# With this script we can count number of reads with indels greater than 2 bases using samtools mpileup output.
# if you wanna count number of reads with indels greater than 4 bases "if [[ "$char" =~ [5-9] ]]; then" part in script should be written.
# if you wanna count number of reads with indels greater than 2 bases "if [[ "$char" =~ [3-9] ]]; then" part in script should be written.
# if you wanna count number of reads with indels greater than 5 bases "if [[ "$char" =~ [6-9] ]]; then" part in script should be written.
# This script takes samtools mpileup output as input to count number of reads with small indels!
# rationale: go 5th column of mpileup output and count integer values which stand for deletion and insertion in reads.

To count reads with indels greater than 2 bases in target site of CRISPR, we need to use samtools mpileup output as input in script below: --- end

# Some important sources for knockout efficiency calculation:
# 1. https://bioconductor.org/packages/devel/bioc/vignettes/CrispRVariants/inst/doc/user_guide.pdf
# 2. https://cmdcolin.github.io/posts/2022-02-06-sv-sam

# Kasper's suggestion is that including supplementary alignments in coverage calculation, that is, do not get rid of them in coverage calculation and counting CRISPR-editing reads!
# reads with small indels larger than 4 or 5 bases should be counted to find CRISPR-editing reads because reads with small indels less than 4 or 5 bases might be sequencing errors!

# To calculate knockout efficiency in ONT data

# --- Start

# 1. Go to target site which is 5bp up and downstream of gRNA coordinates and exclude supplementary reads which do not belong to target site of interest

# samtools mpileup --region chr9:44881220-44881249 --excl-flags 0x800 Ctr6_3_exon1.bam | cut -f 5 | tr '[a-z]' '[A-Z]' | fold -w 1 | sort | uniq -c > Ctr6_3_target_site_mpileup.txt

# --region chr9:44881220-44881249 is 5bp up and downstream of gRNA location for Kmt2a exon1
# --excl-flags excludes reads with 0x800 flags which means supplementary reads
# In Ctr6_3_target_site_mpileup.txt output, number of $ and ^ characters gives us number of reads with large indels, whereas number of + and - characters gives us number of reads with small indels

# 2. Calculation of coverage in target site of interest. In this step, supplementary reads should be exluceded because those reads coming from second target site which is targeted by second gRNA!

# samtools coverage --region chr9:44881220-44881249 --ff 0x800 Ctr6_3_exon1.bam

# 3. Knockout efficiency calculation

sum of number of reads with large and small indels divide by total number of reads in target site (coverage)

# Note: * character in samtools mpileup output represents deletion of reference genome instead representing a deletion in reads. Therefore, in KO efficiency calculation, + and - characters should be used to count reads with small indels instead of * character.
# because * character is not representative for bases in reads, its representative for bases in reference genome

# --- End

# Extra notes about CIGAR string and samtools mpileup

# --- Start

# When you want to have read ids in output of samtools mpileup, you need to add --output-QNAME argument

samtools mpileup --region chr9:44881220-44881249 --output-QNAME --excl-flags 0x800 Ctr6_3_exon1.bam > Ctr6_3_exon1.txt

# When you want to exclude supplementary reads in a particular region of genome and get output as bam file, you can do this with command below

samtools view -h -b -F 0x800 input.bam chr9:44881220-44881249 > output.bam

# Subset a bam file by coordinates of target sites (for example, Kmt2a exon 1 coordinates)

samtools view -@ 128 -b -h Kmt2a_1_merged_and_aligned.bam chr9:44880849-44881296 > Kmt2a_1_exon1.bam

# Get CIGAR strings for each reads in a bam file. 6th column of bam file stands for CIGAR strings of each read.

samtools view -@ 128 Kmt2a_1_exon1.bam | less | cut -f6 >  Kmt2a_1_exon1_CIGARs.txt

# When you want to see only read ids in a bam file

samtools view Ctr6_3_exon1.bam | less | cut -f1 | head

# Burada CIGAR dosyasında 1base deletion iceren kac line var onu hesaplıyoruz! Herbir line’ da read’lerin oldugu düsünülürse aslında read sayısını hesaplamıs oluyoruz asagıdaki command’ lar ile. Yani kac tane read 1base deletion
# veya 20base’lik deletion iceriyor onun sayısını bulmus oluyoruz.

# To find how many lines (reads) have deletion (each line in Kmt2a_1_exon1_CIGARs.txt file has CIGAR string of a read). D letter in CIGAR string stands for deletion in reference.

grep -c "D" Kmt2a_1_exon1_CIGARs.txt

# To find how many lines (reads) have insertion. I letter in CIGAR string stands for insertion.

grep -c "I" Kmt2a_1_exon1_CIGARs.txt

# Find how many reads with 1 base length insertion:

grep -c "1I" Kmt2a_1_exon1_CIGARs.txt

# Find largest deletion in target site:

grep -o "[0-9]*D" Kmt2a_1_exon1_CIGARs.txt | sort -nr | head -1

# For example above command gave us 29D which means the largest deletion is 29 base length.

# The important note: # When you limit coordinates in bam file with gRNA coordinates, that is, target site coordinates, CIGAR strings give me whole read size, it does not give me target site limited CIGAR strings! Thus, when you want to count indels in target site
# to getting CIGAR strings is not good option, instead you can use samtools mpileup to get target site specific indel information!

# Command below allow us to find largest deletion in each read from CIGAR string

#!/bin/bash

while read line; do
    largest_del=$(echo "$line" | grep -o "[0-9]*D" | sort -nr | head -1)
    echo "Line: $line, Largest Deletion: $largest_del"
done < Kmt2a_1_exon1_CIGARs.txt

# To find how many reads have 1 base, 2 bases, 3 bases, 4 bases, 5 bases, 6 bases, 7 bases, 8 bases, 9 bases, and 10 bases length deletion. Note/ each line in CIGAR file is a read.

#!/bin/bash

for length in {1..10}; do
    echo "Number of lines with ${length} base deletion: $(grep -c "${length}D" Kmt2a_1_exon1_CIGARs.txt)"
done

# Note: If you write for length in {1..20}; do, which allows to get how many reads with deletion ranging from 1 base to 20 bases

# When you want to subset a bam file into a bam file with one read, command below can be used

samtools view -@ 128 -h input.bam | awk -v READ_ID="9af91409-808d-4f77-bafc-ca6c23d0ee05" '$1 ~ "@" || $1 == READ_ID' | samtools view -bS - > 9af91409-808d-4f77-bafc-ca6c23d0ee05.bam

# check whether this bam file has only this read or not

samtools view 9af91409-808d-4f77-bafc-ca6c23d0ee05.bam | less | cut -f1 | head

# --- End
