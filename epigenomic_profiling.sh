conda activate pbmm2_env
pbmm2 align /home/pranav/genome_assemblies/primary_data/V_panoptes/panoptes_reanalysis2/assembly_files/genome_assembly/primary.fasta /home/pranav/genome_assemblies/primary_data/V_panoptes/B_01_L27/m84144_260212_101037_s2.hifi_reads.bc2052-001.bam aligned.bam --sort
samtools index aligned.bam

conda create -n pb-cpg-tools_env -c bioconda pb-cpg-tools
conda activate pb-cpg-tools_env
aligned_bam_to_cpg_scores --bam aligned.bam --output-prefix methylation --threads 20

