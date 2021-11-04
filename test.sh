ref=/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna
fq1=/mnt/d/breakpoints/HGT/uhgg_snp/species10_snp0.01_depth20_reads150_sample_1_HGT.1.fq
fq2=/mnt/d/breakpoints/HGT/uhgg_snp/species10_snp0.01_depth20_reads150_sample_1_HGT.2.fq
#ID=species10_snp0.01_depth20_reads150_sample_0
outdir=/mnt/d/breakpoints/HGT/uhgg_snp_results/
#fq1=/mnt/d/breakpoints/HGT/cami_fq/species10_snp0.08_depth20_reads150_sample_0_high_HGT.1.fq
#fq2=/mnt/d/breakpoints/HGT/cami_fq/species10_snp0.08_depth20_reads150_sample_0_high_HGT.2.fq
ID=species10_snp0.08_depth20_reads150_sample_0_high

bash pipeline.sh $ref $fq1 $fq2 $ID $outdir 
