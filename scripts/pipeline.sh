#!/bin/bash

original_ref=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
accept_hit_ratio=$6
accept_perfect_hit_ratio=$7
thread=$8
k=$9
max_peak=${10}
coder_num=${11}
seed=${12}

interval_file=$outdir/$ID.interval.txt
sample=$outdir/$ID
extracted_ref=$outdir/$ID.specific.ref.fasta
start=$(date +%s)
dir=$(cd `dirname $0`; pwd)


if [ ! -d $outdir ]; then
  mkdir $outdir
fi

# :<<!
# echo "$dir/extract_ref $fq1 $fq2 $original_ref $interval_file $accept_hit_ratio $accept_perfect_hit_ratio $thread $k $max_peak $coder_num $seed"
$dir/extract_ref $fq1 $fq2 $original_ref $interval_file $accept_hit_ratio $accept_perfect_hit_ratio $thread $k $max_peak $coder_num $seed
python $dir/get_bed_file.py $original_ref $interval_file > ${sample}.log

samtools faidx -r ${interval_file}.bed $original_ref > $extracted_ref

bwa index $extracted_ref
samtools faidx $extracted_ref
end=$(date +%s)
take=$(( end - start ))
echo Time taken to prepare ref is ${take} seconds. >> ${sample}.log

##################skip bam-sorting#################
bwa mem -M -t $thread -R "@RG\tID:id\tSM:sample\tLB:lib" $extracted_ref $fq1 $fq2 | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.bam $sample.unsort.bam

mv $sample.bam $sample.unique.bam
samtools view -h $sample.unique.bam \
| python3 $dir/extractSplitReads_BwaMem.py -i stdin \
| samtools view -Sb > $sample.unsort.splitters.bam

samtools sort -o $sample.splitters.bam $sample.unsort.splitters.bam



samtools index $sample.splitters.bam
samtools index $sample.unique.bam

end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds. >> ${sample}.log

python $dir/get_raw_bkp.py -t $thread -u $sample.unique.bam -o $sample.raw.csv
# !
python $dir/accurate_bkp.py -g $original_ref -u $sample.unique.bam -b ${interval_file}.bed \
-s $sample.splitters.bam -a $sample.raw.csv -o $sample.repeat.acc.csv -t $thread

python $dir/remove_repeat.py $sample.repeat.acc.csv $sample.acc.csv
rm $sample.repeat.acc.csv

if [ ! -f "$sample.acc.csv" ]; then
    echo "Error: Final HGT breakpoint file is not generated."
else
    echo "Final HGT breakpoint file is generated."
    wc -l $sample.acc.csv
fi

rm $extracted_ref*
rm $sample.unsort.splitters.bam
rm $sample.unsort.bam
# rm $sample.splitters.bam

end=$(date +%s)
take=$(( end - start ))
echo Time taken to execute commands is ${take} seconds. >> ${sample}.log
echo "Final result is in $sample.acc.csv"
echo "--------------------------"
echo "Finished!"














##################original bam-sorting version################
# bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $extracted_ref $fq1 $fq2 \
#   | samtools view -bhS -> $sample.unsort.bam
# # Sort bam file
# samtools sort -o $sample.bam $sample.unsort.bam

# end=$(date +%s)
# take=$(( end - start ))
# echo Time taken to map reads is ${take} seconds. >> ${sample}.log
# rm $sample.unsort.bam

# # Extract split reads
# samtools view -h $sample.bam \
#   | python /mnt/d/breakpoints/script/extractSplitReads_BwaMem.py -i stdin \
#   | samtools view -Sb > $sample.unsort.splitters.bam
# # Sort split reads bam file
# samtools sort -o $sample.splitters.bam $sample.unsort.splitters.bam
# # Extract unique reads bam file
# samtools view -q 20 -b $sample.bam > $sample.unique.bam
# samtools index $sample.unique.bam
# samtools index $sample.splitters.bam

# python /mnt/d/breakpoints/script/get_raw_bkp_v2.py -u $sample.unique.bam \
# -o $sample.raw.txt
# python /mnt/d/breakpoints/script/accurate_bkp.py -g $original_ref -u $sample.unique.bam \
# -s $sample.splitters.bam -a $sample.raw.txt -o $sample.acc.txt



