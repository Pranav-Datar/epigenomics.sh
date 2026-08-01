# Genome-wide CpG Methylation Profiling from PacBio HiFi Reads

## 1. Align HiFi reads to the reference genome

Epigenomic profiling requires an aligned BAM file. Align the PacBio HiFi reads to the assembled reference genome using **pbmm2**, then index the resulting BAM file.

```bash
conda activate pbmm2_env

pbmm2 align genome_assembly.fasta \
    m84144_260212_101037_s2.hifi_reads.bc2052-001.bam \
    aligned.bam \
    --sort

samtools index aligned.bam
```

---

## 2. Call CpG methylation

Install and activate **pb-cpg-tools**, then generate genome-wide CpG methylation calls.

```bash
conda create -n pb-cpg-tools_env -c bioconda pb-cpg-tools
conda activate pb-cpg-tools_env

aligned_bam_to_cpg_scores \
    --bam aligned.bam \
    --output-prefix methylation \
    --threads 20
```

This generates:

* `methylation.combined.bed.gz` — per-CpG methylation calls
* `methylation.combined.bw` — genome browser track
* `methylation.log` — run log

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
