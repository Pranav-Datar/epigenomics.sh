# Allele-Specific DNA Methylation Analysis Using PacBio HiFi Reads

#This workflow performs genome-wide and allele-specific DNA methylation analysis from PacBio HiFi sequencing data.

#The pipeline consists of:

#1. Align HiFi reads to the reference genome.
#2. Call genome-wide CpG methylation.
#3. Perform variant calling.
#4. Filter high-confidence variants.
#5. Phase heterozygous variants.
#6. Haplotag sequencing reads.
#7. Generate haplotype-specific methylomes.
#8. Quantify epigenetic heterozygosity by comparing methylation levels between homologous chromosomes.



## Directory Structure

#epigenomics/
#│
#├── alignment/
#├── haplotypes/
#├── methylation/
#│   ├── whole_genome/
#│   ├── assembled_haplotypes/
#│   ├── haplotype_specific/
#│   └── allele_specific_methylation/
#│
#├── variants/
#├── phasing/
#├── statistics/
#├── figures/
#├── scripts/
#├── logs/
#└── reference/


#---

## Required Software

# pbmm2, samtools, DeepVariant (v1.10), bcftools, WhatsHap, pb-cpg-tools, R (base R is sufficient)


## Input Files

#Reference genome: reference/primary.fasta
#PacBio HiFi reads:reads.bam


##################################################################################

# 2. Read Alignment

#Align PacBio HiFi reads against the assembled reference genome.

## Activate pbmm2

conda activate pbmm2_env

## Align reads

pbmm2 align \
    reference/primary.fasta \
    reads.bam \
    alignment/aligned.bam \
    --sort


## Index the BAM file

samtools index alignment/aligned.bam


### Output

#alignment/
#├── aligned.bam
#└── aligned.bam.bai

#####################################################

# 3. Whole-Genome CpG Methylation Calling

#Generate genome-wide CpG methylation estimates directly from PacBio HiFi reads.

## Activate pb-cpg-tools

conda activate pb-cpg-tools_env

## Run methylation calling

aligned_bam_to_cpg_scores \
    --bam alignment/aligned.bam \
    --output-prefix methylation/whole_genome/methylation \
    --threads 40


### Output


#methylation/
#└── whole_genome/
#    ├── methylation.combined.bed.gz
#    ├── methylation.combined.bed.gz.tbi
#    └── methylation.combined.bw


##########################################

# 4. Whole-Genome Methylation Statistics

## Calculate weighted genome-wide CpG methylation

#The weighted methylation level accounts for sequencing depth by summing methylated and unmethylated read counts across all CpG sites.

zcat methylation/whole_genome/methylation.combined.bed.gz | \
awk '
BEGIN{
    methylated=0
    total=0
}

!/^#/{
    methylated += $7
    total += ($7+$8)
}

END{
    printf("Weighted CpG methylation = %.2f%%\n",
           100*methylated/total)
}'


## Classify CpGs by methylation level

#CpGs can be classified into three categories based on discretized methylation percentage.

# Low (<30%), Intermediate (30–70%), High (>70%)

zcat methylation/whole_genome/methylation.combined.bed.gz | \
awk '
BEGIN{
    low=0
    intermediate=0
    high=0
    total=0
}

!/^#/ && $6>=10{

    total++

    if($9>=70)
        high++
    else if($9>=30)
        intermediate++
    else
        low++
}

END{

    printf("CpGs analysed (coverage ≥10): %d\n\n",total)

    printf("High methylation (≥70%%)        : %d (%.2f%%)\n",
            high,100*high/total)

    printf("Intermediate methylation (30–70%%): %d (%.2f%%)\n",
            intermediate,100*intermediate/total)

    printf("Low methylation (<30%%)         : %d (%.2f%%)\n",
            low,100*low/total)

}'


## Methylation Distribution

#Generate a genome-wide methylation distribution table.

zcat methylation/whole_genome/methylation.combined.bed.gz | \
awk '
BEGIN{
    OFS="\t"
}

!/^#/ && $6>=10{

    bin=int($9/10)

    if(bin==10)
        bin=9

    bins[bin]++
    total++

}

END{

    print "Bin","Count","Percent"

    for(i=0;i<10;i++){

        start=i*10
        end=start+10

        printf("%d-%d\t%d\t%.2f\n",
               start,
               end,
               bins[i],
               100*bins[i]/total)

    }

}'


### Outputs

#Weighted genome-wide methylation
#Number of CpGs analysed
#High / Intermediate / Low methylated CpGs
#Genome-wide methylation distribution


### Biological interpretation

#This section summarizes the global methylation landscape of the genome before allele-specific analyses.
#The weighted methylation percentage provides an overall estimate of CpG methylation, while the methylation distribution describes the proportion of CpGs that are constitutively methylated, unmethylated, or partially methylated across the genome.

###################################################################

# 5. Variant Calling

## Activate DeepVariant

REF=reference/primary.fasta
BAM=alignment/aligned.bam

podman run --rm \
    -v $(pwd):/work \
    -v $(dirname "$REF"):/ref \
    docker.io/google/deepvariant:1.10.0 \
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref /ref/$(basename "$REF") \
        --reads /work/$BAM \
        --output_vcf /work/variants/variants.vcf.gz \
        --output_gvcf /work/variants/variants.g.vcf.gz \
        --num_shards 40

### Output

#variants/
#├── variants.vcf.gz#
#├── variants.vcf.gz.tbi
#├── variants.g.vcf.gz
#└── variants.g.vcf.gz.tbi

####################################################

# 6. Variant Quality Assessment


conda activate bcftools_env

#Generate variant statistics.

bcftools stats variants/variants.vcf.gz \
> variants/variants.stats

#Inspect summary statistics.


grep "^SN" variants/variants.stats

#Inspect transition/transversion ratio.

grep "^TSTV" variants/variants.stats

############################################################

# 7. Initial Variant Filtering

#Retain only variants that passed all DeepVariant filters (PASS filter).


bcftools view \
-f PASS \
-Oz \
-o variants/variants_PASS.vcf.gz \
variants/variants.vcf.gz


Index the filtered VCF.


bcftools index variants/variants_PASS.vcf.gz


######################################################################

# 8. Determine Genotype Quality (GQ) and Read Depth (DP) Thresholds

#Rather than choosing arbitrary thresholds, inspect the empirical distributions of genotype quality and read depth.

#Extract genotype quality.

bcftools query \
-f '[%GQ\n]' \
variants/variants_PASS.vcf.gz \
> variants/GQ.txt

#Extract read depth.

bcftools query \
-f '[%DP\n]' \
variants/variants_PASS.vcf.gz \
> variants/DP.txt


#---

## Visualize GQ and DP distributions


conda activate r_env

R

GQ <- scan("variants/GQ.txt")
DP <- scan("variants/DP.txt")

summary(GQ)
summary(DP)

quantile(GQ,
         probs=c(0,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,1))

quantile(DP,
         probs=c(0,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,1))

pdf("figures/Variant_QC_histograms.pdf",
    width=10,
    height=5)

par(mfrow=c(1,2))

hist(GQ,
     breaks=100,
     main="Genotype Quality",
     xlab="GQ")

abline(v=20,
       lwd=2,
       lty=2)

hist(DP,
     breaks=100,
     main="Read Depth",
     xlab="DP")

abline(v=10,
       lwd=2,
       lty=2)

dev.off()


#---

## Select Filtering Thresholds

#Filtering thresholds should be chosen based on the empirical distributions above.

#Typical values are
#GQ ≥ 20
#DP ≥ 10
#Upper DP threshold determined from the depth distribution

#The upper depth threshold helps remove regions with abnormally high coverage, which often correspond to repetitive sequences, collapsed duplications, or mapping artifacts.

###########################################################

# 9. Final Variant Filtering

#Apply the selected quality thresholds.

#Example

bcftools filter \
-i 'FORMAT/GQ>=20 && FORMAT/DP>=10 && FORMAT/DP<=50' \
variants/variants_PASS.vcf.gz \
-Oz \
-o variants/variants_filtered.vcf.gz


#Index the filtered VCF.

bcftools index variants/variants_filtered.vcf.gz


#---

## Evaluate the filtered variant set

#Generate updated statistics.


bcftools stats \
variants/variants_filtered.vcf.gz \
> variants/variants_filtered.stats


#Inspect summary statistics.


grep "^SN" variants/variants_filtered.stats


#Inspect Ti/Tv ratio.


grep "^TSTV" variants/variants_filtered.stats


3Count heterozygous variants.

bcftools view \
-g het \
variants/variants_filtered.vcf.gz \
| bcftools view -H \
| wc -l


#---

### Biological interpretation

This section generates a high-confidence set of variants suitable for downstream phasing.

Filtering low-quality variants minimizes phasing errors, while removing excessively high-depth variants reduces the influence of repetitive regions and mapping artifacts.


#########################################################################################################


# 10. Read-Based Haplotype Phasing

#High-confidence heterozygous variants are phased using WhatsHap. Phasing reconstructs 
#the two homologous chromosome sequences by determining which alleles are inherited together.

conda activate whatshap_env

#Run phasing.


whatshap phase \
    --reference reference/primary.fasta \
    --output variants/phased.vcf.gz \
    variants/variants_filtered.vcf.gz \
    alignment/aligned.bam


#Index the phased VCF.


tabix -p vcf variants/phased.vcf.gz


### Output

#variants/
#├── phased.vcf.gz
#└── phased.vcf.gz.tbi


#---

### Biological interpretation

#Variant calling identifies heterozygous variants but does not determine which alleles belong to the same chromosome.

#For example,

#Position 1000 : A/G
#Position 2000 : C/T

#could represent either

#Haplotype 1 : A C
#Haplotype 2 : G T


#or


#Haplotype 1 : A T
#Haplotype 2 : G C


#WhatsHap reconstructs these haplotypes by leveraging long PacBio HiFi reads spanning multiple heterozygous variants.

#---

##############################################################

# 11. Haplotagging PacBio Reads

#Assign each sequencing read to its corresponding haplotype.

whatshap haplotag \
    --reference reference/primary.fasta \
    --output alignment/haplotagged.bam \
    variants/phased.vcf.gz \
    alignment/aligned.bam

#Index the haplotagged BAM.


samtools index alignment/haplotagged.bam


### Output

#alignment/
#├── haplotagged.bam
#└── haplotagged.bam.bai


#---

### Biological interpretation

#After phasing, the haplotype structure is known.
#Haplotagging labels every PacBio HiFi read according to the haplotype from which it originated.
#Instead of analysing all reads together,

#Read1
#Read2
#Read3
#Read4


#the reads become

#Haplotype 1

#Read1
#Read2
#Read5

#Haplotype 2

#Read3
#Read4
#Read6

#This enables methylation to be estimated independently for each homologous chromosome.

######################################################

# 12. Allele-Specific CpG Methylation Calling

conda activate pb-cpg-tools_env


aligned_bam_to_cpg_scores \
    --bam alignment/haplotagged.bam \
    --output-prefix methylation/haplotype_specific/allele_specific_methylation \
    --threads 40


### Output

#methylation/
#└── haplotype_specific/

    allele_specific_methylation.combined.bed.gz

    allele_specific_methylation.hap1.bed.gz

    allele_specific_methylation.hap2.bed.gz

    allele_specific_methylation.combined.bw

    allele_specific_methylation.hap1.bw

    allele_specific_methylation.hap2.bw


#---

### Biological interpretation

#Unlike the whole-genome methylation analysis, which estimates methylation by combining reads from both homologous chromosomes, 
#haplotagging enables pb-cpg-tools to estimate methylation independently for each chromosome copy.

#Each CpG therefore has two methylation estimates:
#for example
#Haplotype 1 : 82%
#Haplotype 2 : 18%

#rather than a single genome-wide estimate of

#50%

#This separation enables the identification of **allele-specific methylation (ASM)** throughout the genome.



#############################################################################

# 13. Coverage Quality Control

#Coverage is evaluated independently for each haplotype before comparing methylation levels.

#Extract coverage values.

zcat methylation/haplotype_specific/allele_specific_methylation.hap1.bed.gz \
| awk '!/^#/ {print $6}' \
> statistics/hap1_coverage.txt

zcat methylation/haplotype_specific/allele_specific_methylation.hap2.bed.gz \
| awk '!/^#/ {print $6}' \
> statistics/hap2_coverage.txt


#Coverage distributions are then inspected in R to determine an appropriate minimum coverage threshold for downstream analyses.

#Typically,
#Coverage ≥10 provides a good balance between reliability and CpG retention, 
#although the exact threshold should be determined empirically for each dataset.

conda activate r_env
R
hap1 <- scan("hap1_coverage.txt")
hap2 <- scan("hap2_coverage.txt")

summary(hap1)
summary(hap2)

quantile(hap1, probs=c(0,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,1))
quantile(hap2, probs=c(0,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,1))

pdf("Coverage_distribution.pdf", width=10, height=8)

par(mfrow=c(2,2))

hist(hap1, breaks=100, xlim=c(0,100),
     main="Haplotype 1 Coverage",
     xlab="Coverage", ylab="CpG count")
abline(v=10, lty=2, lwd=2)
abline(v=15, lty=3, lwd=2)

hist(hap2, breaks=100, xlim=c(0,100),
     main="Haplotype 2 Coverage",
     xlab="Coverage", ylab="CpG count")
abline(v=10, lty=2, lwd=2)
abline(v=15, lty=3, lwd=2)

boxplot(hap1, horizontal=TRUE, main="Haplotype 1")
boxplot(hap2, horizontal=TRUE, main="Haplotype 2")

cat("\nRetention (%) at different coverage thresholds\n")
for(i in c(4,6,8,10,12,15,20)){
  cat(sprintf("Coverage >= %-2d : Hap1 = %6.2f%%   Hap2 = %6.2f%%\n",
              i,
              100*mean(hap1>=i),
              100*mean(hap2>=i)))
}

dev.off()


############################################################################################

# 14. Coverage Filtering

#To ensure reliable methylation estimates, CpGs with insufficient read support are removed independently from each haplotype.
#The minimum coverage threshold should be selected based on the empirical coverage distribution (Section 13).

#Example (coverage ≥10):

zcat methylation/haplotype_specific/allele_specific_methylation.hap1.bed.gz \
| awk '!/^#/ && $6>=10 {print $1"\t"$2"\t"$6"\t"$9}' \
> statistics/hap1_methylation_cov10.tsv

zcat methylation/haplotype_specific/allele_specific_methylation.hap2.bed.gz \
| awk '!/^#/ && $6>=10 {print $1"\t"$2"\t"$6"\t"$9}' \
> statistics/hap2_methylation_cov10.tsv


#Sort both files by genomic coordinate.

sort -k1,1 -k2,2n \
statistics/hap1_methylation_cov10.tsv \
> statistics/hap1.cov10.sorted.tsv

sort -k1,1 -k2,2n \
statistics/hap2_methylation_cov10.tsv \
> statistics/hap2.cov10.sorted.tsv



#######################################################################


# 15. Identify Shared CpG Sites

#Only CpGs observed in both haplotypes are retained for comparison.

#For each shared CpG, calculate

# Δ methylation = | Hap1 − Hap2 |

awk '
BEGIN{
    FS=OFS="\t"
}

FNR==NR{

    key=$1":"$2

    hap1[key]=$4

    next

}

{

    key=$1":"$2

    if(key in hap1){

        diff=hap1[key]-$4

        if(diff<0)
            diff=-diff

        print $1,$2,hap1[key],$4,diff

    }

}
' \
statistics/hap1.cov10.sorted.tsv \
statistics/hap2.cov10.sorted.tsv \
> statistics/haplotype_methylation_difference.tsv


#Output columns

#Chromosome

#Position

#Haplotype 1 methylation (%)

#Haplotype 2 methylation (%)

#Absolute methylation difference (Δ)


#---

######################################################################

# 16. Summary Statistics

conda activate r_env
R
delta <- read.table(
    "statistics/haplotype_methylation_difference.tsv"
)$V5

summary(delta)

quantile(
    delta,
    probs=c(
        0,
        0.01,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99,
        1
    )
)

sd(delta)

var(delta)

IQR(delta)


#These statistics summarize the genome-wide distribution of allele-specific methylation differences.

###################################################################################################################################

# 17. Visualizing Epigenetic Heterozygosity

R
d <- read.table(
    "statistics/haplotype_methylation_difference.tsv",
    header=FALSE
)

delta <- d$V5

## Density and ECDF are generated from a random subset
## to reduce figure size without affecting interpretation.

set.seed(123)

delta.sample <- sample(
    delta,
    min(length(delta),100000)
)

png(
    "figures/Delta_methylation_statistics.png",
    width=3200,
    height=2400,
    res=300
)

par(
    mfrow=c(2,2),
    mar=c(4,4,2,1),
    las=1
)

hist(
    delta,
    breaks=seq(0,100,2),
    col="grey85",
    border="grey40",
    main="Histogram",
    xlab="Absolute methylation difference (%)"
)

plot(
    density(delta.sample),
    lwd=2,
    main="Density",
    xlab="Absolute methylation difference (%)",
    ylab="Probability density"
)

boxplot(
    delta,
    horizontal=TRUE,
    col="grey85",
    main="Boxplot",
    xlab="Absolute methylation difference (%)"
)

plot(
    ecdf(delta.sample),
    main="Empirical cumulative distribution",
    xlab="Absolute methylation difference (%)",
    ylab="Cumulative proportion"
)

dev.off()


########################################

# 18. Distribution of Allelic Methylation Differences

#Summarize Δ methylation into biologically meaningful classes.


awk '
{

    if($5<20)
        low++

    else if($5<50)
        moderate++

    else if($5<80)
        high++

    else
        extreme++

    total++

}

END{

    printf("Low (<20%%)\t%d\t%.2f%%\n",
        low,
        100*low/total)

    printf("Moderate (20–50%%)\t%d\t%.2f%%\n",
        moderate,
        100*moderate/total)

    printf("High (50–80%%)\t%d\t%.2f%%\n",
        high,
        100*high/total)

    printf("Extreme (≥80%%)\t%d\t%.2f%%\n",
        extreme,
        100*extreme/total)

}
' statistics/haplotype_methylation_difference.tsv


#---

# Biological Interpretation

#Whole-genome methylation combines reads originating from both homologous chromosomes and therefore represents an average methylation level.

#By separating reads into haplotypes before methylation calling, methylation levels can be estimated independently for each chromosome copy.

#The absolute methylation difference
# Δ = | Hap1 − Hap2 |

#provides a quantitative measure of allele-specific methylation (ASM).

#Small Δ values indicate similar methylation states between homologous chromosomes, whereas large Δ values indicate strong allele-specific methylation.

#The genome-wide distribution of Δ methylation therefore represents a measure of **epigenetic heterozygosity**, allowing direct comparison among individuals or species.

#---

# Final Outputs
Whole-genome methylome
Haplotype-specific methylomes
High-confidence phased variants
Haplotagged BAM
Genome-wide Δ methylation table
Summary statistics
Coverage QC figures
Variant QC figures
Genome-wide Δ methylation figure
Epigenetic heterozygosity estimates
