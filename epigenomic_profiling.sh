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

zcat methylation/methylation.combined.bed.gz | awk '
BEGIN{
    low=0;
    intermediate=0;
    high=0;
    total=0;
}
!/^#/ && $6>=10 {

    total++;

    if($9 >= 70)
        high++;
    else if($9 >= 30)
        intermediate++;
    else
        low++;
}
END{
    printf "CpGs analyzed (coverage >=10): %d\n", total;
    printf "High methylation (>=70%%): %d (%.2f%%)\n", high, 100*high/total;
    printf "Intermediate (30-70%%): %d (%.2f%%)\n", intermediate, 100*intermediate/total;
    printf "Low methylation (<30%%): %d (%.2f%%)\n", low, 100*low/total;
}'

##homozygous and heterozygous methylation
#variant calling

podman run --rm \
  -v /home/pranav/genome_assemblies/primary_data/V_brevicauda/brevicauda_combined/epigenomics:/work \
  -v /home/pranav/genome_assemblies/primary_data/V_brevicauda/brevicauda_combined_reanalysis2/assembly_files/genome_assembly:/ref \
  docker.io/google/deepvariant:1.10.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref /ref/primary.fasta \
  --reads /work/alignment/aligned.bam \
  --output_vcf /work/variants/variants.vcf.gz \
  --output_gvcf /work/variants/variants.g.vcf.gz \
  --num_shards 30
