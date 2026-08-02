# Genome-wide CpG Methylation Profiling from PacBio HiFi Reads

## 1. Align HiFi reads to the reference genome

#Epigenomic profiling requires an aligned BAM file. Align the PacBio HiFi reads to the assembled reference genome using **pbmm2**, then index the resulting BAM file.

conda activate pbmm2_env

pbmm2 align genome_assembly.fasta \
    m84144_260212_101037_s2.hifi_reads.bc2052-001.bam \
    aligned.bam \
    --sort

samtools index aligned.bam
```

---

##2. Call CpG methylation

#Install and activate **pb-cpg-tools**, then generate genome-wide CpG methylation calls.


conda create -n pb-cpg-tools_env -c bioconda pb-cpg-tools
conda activate pb-cpg-tools_env

aligned_bam_to_cpg_scores \
    --bam aligned.bam \
    --output-prefix methylation \
    --threads 20
```

#This generates:

#methylation.combined.bed.gz — per-CpG methylation calls
#methylation.combined.bw — genome browser track
#methylation.log — run log

---

## 3. Calculate genome-wide weighted CpG methylation

This computes the weighted average methylation across all CpG sites using the estimated methylated and unmethylated read counts.

```bash
zcat methylation.combined.bed.gz | awk '
BEGIN{
    meth=0;
    total=0;
}
!/^#/{
    meth += $7;
    total += ($7 + $8);
}
END{
    printf("Weighted CpG methylation (%%) = %.2f\n",
           (meth/total)*100);
}'
```

---

## 4. Classify CpGs by methylation level

The following script classifies CpGs into three methylation categories after filtering for a minimum coverage of **10×**.

* **High methylation:** ≥70%
* **Intermediate methylation:** 30–70%
* **Low methylation:** <30%

```bash
zcat methylation.combined.bed.gz | awk '
BEGIN{
    low=0;
    intermediate=0;
    high=0;
    total=0;
}
!/^#/ && $6>=10{

    total++;

    if($9>=70)
        high++;
    else if($9>=30)
        intermediate++;
    else
        low++;
}
END{
    printf("CpGs analyzed (coverage ≥10×): %d\n", total);
    printf("High methylation (≥70%%): %d (%.2f%%)\n",
           high,100*high/total);
    printf("Intermediate methylation (30–70%%): %d (%.2f%%)\n",
           intermediate,100*intermediate/total);
    printf("Low methylation (<30%%): %d (%.2f%%)\n",
           low,100*low/total);
}'
```

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

#check the variants
bcftools stats variants.vcf.gz > variants.stats

#variant filtering

#keep PASS variants
bcftools view -f PASS variants.vcf.gz -Oz -o variants_PASS.vcf.gz
bcftools index variants_PASS.vcf.gz

#GQ and DP filtering

#to set the threshold, examine the GQ and DP distribution
bcftools query -f '[%GQ\n]' variants_PASS.vcf.gz > GQ.txt
bcftools query -f '[%DP\n]' variants_PASS.vcf.gz > DP.txt

awk '
{
sum+=$1
n++
if(min=="" || $1<min) min=$1
if($1>max) max=$1
}
END{
print "Variants:",n
print "Min:",min
print "Mean:",sum/n
print "Max:",max
}' GQ.txt

awk '
{
sum+=$1
n++
if(min=="" || $1<min) min=$1
if($1>max) max=$1
}
END{
print "Variants:",n
print "Min:",min
print "Mean:",sum/n
print "Max:",max
}' DP.txt

conda deactivate
conda activate r_env
R
GQ <- scan("GQ.txt")
DP <- scan("DP.txt")
length(GQ)
length(DP)
summary(GQ)
summary(DP)

quantile(GQ,
         probs = c(0,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,1))

quantile(DP,
         probs = c(0,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,1))hist(GQ,
     breaks = 100,
     main = "Genotype Quality (GQ)",
     xlab = "Genotype Quality",
     ylab = "Number of variants")

     hist(GQ[GQ <= 50],
     breaks = 50,
     main = "Genotype Quality (0-50)",
     xlab = "Genotype Quality",
     ylab = "Number of variants")

     hist(DP,
     breaks = 100,
     main = "Read Depth (DP)",
     xlab = "Read Depth",
     ylab = "Number of variants")

     hist(DP[DP <= 60],
     breaks = 60,
     main = "Read Depth (0-60)",
     xlab = "Read Depth",
     ylab = "Number of variants")

     pdf("Variant_QC_histograms.pdf", width = 8, height = 6)

hist(GQ[GQ <= 50],
     breaks = 50,
     main = "Genotype Quality (0-50)",
     xlab = "GQ",
     ylab = "Number of variants")

hist(DP[DP <= 60],
     breaks = 60,
     main = "Read Depth (0-60)",
     xlab = "DP",
     ylab = "Number of variants")

dev.off()

 bcftools filter -i 'FORMAT/GQ>=20 && FORMAT/DP>=10 && FORMAT/DP <=50' variants_PASS.vcf.gz -Oz -o variants_GQ20_DP10_50_PASS.vcf.gz

#total number of variants left
bcftools view -H variants_GQ20_DP10_50_PASS.vcf.gz | wc -l

#keep only heterozygous variants
bcftools view -g het variants_GQ20_DP10_50_PASS.vcf.gz -Oz -o variants_het_GQ20_DP10_50_PASS.vcf.gz

bcftools index variants_GQ20_DP10_50_PASS.vcf.gz
bcftools index variants_het_GQ20_DP10_50_PASS.vcf.gz
 
#WhatsHap for phasing, it uses variant calling data and aligned.bam file. Assigns variants to paternal and maternal haplotypes

conda create -n whatshap_env -c bioconda whatshap
