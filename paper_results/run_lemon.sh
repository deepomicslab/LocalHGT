
ref=$1
fq1=$2
fq2=$3
outdir=$5
ID=$4
sample=$outdir/$ID



start=$(date +%s)
#:<<!
if [ ! -f $ref.pac ];then
samtools faidx $ref
bwa index $ref
fi

bwa mem -M -t 10 -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam

# Sort bam file
samtools sort -o $sample.bam $sample.unsort.bam

end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.  >> ${sample}.log
rm $sample.unsort.bam

# Extract split reads
samtools view -h $sample.bam \
  | python3 /mnt/d/breakpoints/script/scripts/extractSplitReads_BwaMem.py -i stdin \
  | samtools view -Sb > $sample.unsort.splitters.bam
# Sort split reads bam file
samtools sort -o $sample.splitters.bam $sample.unsort.splitters.bam
# Extract unique reads bam file
samtools view -q 20 -b $sample.bam > $sample.unique.bam
samtools index $sample.unique.bam
# Calculate coverage
# bedtools genomecov -ibam $sample.bam -bg > $sample.coverage.txt
#rm $sample.bam
!
python /mnt/d/breakpoints/lemon/get_raw_bkp.py -t 10 -r $ref -u $sample.unique.bam -o $sample.raw.csv

/mnt/d/breakpoints/lemon//get_acc_bkp -t 10 -r $ref -u $sample.unique.bam -s $sample.splitters.bam -b $sample.raw.csv -o $sample.acc.csv

end=$(date +%s)
take=$(( end - start ))
echo Time taken to run with LEMON is ${take} seconds.  >> ${sample}.log
