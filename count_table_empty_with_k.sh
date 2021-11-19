for k in {20..32}
do
for r in {1..9}
do

./count_diff_kmer /mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.01_depth50_reads150_sample_0_high.1.fq /mnt/d/breakpoints/HGT/uhgg_snp//species20_snp0.01_depth50_reads150_sample_0_high.2.fq $k ${r}0
done
done
