#!/bin/bash

# Increasing how many file can be read at the same time

-Sn 50000

# samtools merge to merge bam files into one larger bam file

samtools merge -n -c -@ 128 -o Ctr6_1_merged.bam bam_pass/*.bam
exit
