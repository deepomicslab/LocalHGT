original_ref=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5

interval_file=$outdir/$ID.interval.txt
sample=$outdir/$ID
extracted_ref=$outdir/$ID.specific.ref.fasta
start=$(date +%s)

samtools faidx $original_ref
/mnt/d/breakpoints/script/extract_ref $fq1 $fq2 $original_ref $interval_file $6 $7
cat ${interval_file}_tmp_* >$interval_file
rm ${interval_file}_tmp_*
python /mnt/d/breakpoints/script/extract_ref_seq.py $original_ref $extracted_ref $interval_file
:<<!
bwa index $extracted_ref
samtools faidx $extracted_ref



bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $extracted_ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam

# Sort bam file
samtools sort -o $sample.bam $sample.unsort.bam

end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.
rm $sample.unsort.bam

# Extract split reads
samtools view -h $sample.bam \
  | python /mnt/d/breakpoints/script/extractSplitReads_BwaMem.py -i stdin \
  | samtools view -Sb > $sample.unsort.splitters.bam
# Sort split reads bam file
samtools sort -o $sample.splitters.bam $sample.unsort.splitters.bam
# Extract unique reads bam file
samtools view -q 20 -b $sample.bam > $sample.unique.bam
samtools index $sample.unique.bam
samtools index $sample.splitters.bam
# Calculate coverage
# bedtools genomecov -ibam $sample.bam -bg > $sample.coverage.txt
#rm $sample.bam



python /mnt/d/breakpoints/script/get_raw_bkp_v2.py -u $sample.unique.bam \
-o $sample.raw.txt
python /mnt/d/breakpoints/script/accurate_bkp.py -g $original_ref -u $sample.unique.bam \
-s $sample.splitters.bam -a $sample.raw.txt -o $sample.acc.txt
!

end=$(date +%s)
take=$(( end - start ))
echo Time taken to execute commands is ${take} seconds.
