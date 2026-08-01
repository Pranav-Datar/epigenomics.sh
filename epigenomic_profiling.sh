#Epigenomic profiling would require an aligned bam file. So the raw bam file from pacbio output should be aligned alongthe assembled genome.

conda activate pbmm2_env
pbmm2 align genome_assembly.fasta m84144_260212_101037_s2.hifi_reads.bc2052-001.bam aligned.bam --sort
samtools index aligned.bam

conda create -n pb-cpg-tools_env -c bioconda pb-cpg-tools
conda activate pb-cpg-tools_env
aligned_bam_to_cpg_scores --bam aligned.bam --output-prefix methylation --threads 20

#to calculate % methylation at CpG sites
zcat methylation.combined.bed.gz | awk '
BEGIN{
meth=0;
total=0
}
!/^#/ {
meth += $7;
total += ($7 + $8)
}
END {
print "Weighted CpG methylation (%) =", (meth/total)*100
}'

##homozygous and heterozygous methylation
