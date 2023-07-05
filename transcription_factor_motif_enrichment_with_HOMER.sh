# How to install HOMER and use it for motif enrichment analysis of differentially methylated regions?

# Start -------

# 1. create a directory with name of homer
# 2. create another directory with name of bin in this homer directory
# 3. Downloading "configureHomer.pl" perl script from link below
# "http://homer.ucsd.edu/homer/introduction/install.html"
# 4. Place this perl script to homer directory and then call it like:
/usr/bin/perl configureHomer.pl -install
# Note/ without creating bin directory, we're not able to install HOMER!
# 5. To call any script of HOMER in bin directory, we can add path to bin directory like:
export PATH=/path/to/homer/bin:$PATH
# 6. Lets call "findMotifsGenome.pl" script like:
findMotifsGenome.pl allKOs_vs_all_controls_DMRs.txt /home/ko/homer/data/genomes/mm10 24_DMRs_results_cpg_normalization -fdr -seqlogo -cpg -p 8
# allKOs_vs_all_controls_DMRs.txt is bed format file containing coordinates of differentially methylated regions
# /home/ko/homer/data/genomes/mm10 is a path to mm10 UCSC mouse genome.
# To create this mouse genome we can run: "configureHomer.pl -list" to see genomes can be called and then install mm10 UCSC mouse genome "configureHomer.pl -install mm10"
# 24_DMRs_results_cpg_normalization is output directory where all results will be saved
# -size argument is fragment size to use for motif finding, default=200
# -len is motif length, default=8
# -fdr allows calculation of empirical FDR for de novo discovery #=number of randomizations
# -seqlogo is use weblogo/seqlogo/ghostscript to generate logos, default uses SVG now.
# -cpg is (use CpG% instead of GC% for sequence content normalization)
# Note/ for CpG-level methylation data, I used CpG normalization for motif search and this made difference than GC% normalization in terms of reported number of transcription factors and 
# type of transcription factors as well!

# Homer can be updated like: "perl configureHomer.pl -update"

# To learn rationale behind HOMER read following paper: "https://www.sciencedirect.com/science/article/pii/S1097276510003667?via%3Dihub"

# After detection of enriched motifs, how to find which input region (in this case differentially methylated regions) has identified enriched DNA motifs?
# How to find coordinates of enriched DNA motifs (transcription factor binding sites)?

$ annotatePeaks.pl \
allKOs_vs_all_controls_DMRs.txt \ # (differentially methylated regions)
/home/ko/homer/data/genomes/mm10 \ # (path to mm10 mouse genome)
-m /home/ko/homer/24_DMRs_results_cpg_normalization/homerMotifs.all.motifs \ # (all identified motifs from homer output)
> Finding_instance_of_specific_motifs_24_DMRs_results_cpg_normalization_2.txt # (output)

# After getting Finding_instance_of_specific_motifs_24_DMRs_results_cpg_normalization_2.txt output, go UCSC genome browser and select conservation track
# to see the DNA sequence of enriched motifs.

# http://homer.ucsd.edu/homer/ngs/peakMotifs.html (the link to find location of enriched motifs).

# How to find methylated motifs, that is, transcription factor binding sites which are methylated and this is related to binding affinity!
# Check following links: "https://hemtools.readthedocs.io/en/latest/content/Motif_analysis/methylmotifs.html" & https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748772/
# & http://wanglab.ucsd.edu/star/mEpigram/

# End -------
