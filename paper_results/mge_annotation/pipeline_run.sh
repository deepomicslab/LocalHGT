
perl 0.merge_bp.pl > identified_event.merge_bp.csv

perl 1.note_inter_intra_phylum.pl > identified_event.merge_bp.anno_type

# Download gff.gz files of representive genomes from UHGG database to genome/ folder

for i in `cat identified_event.merge_bp.anno_type |awk -F ',' '{print $3"\n"$5}'|sort|uniq`
do
	perl 2.extract_fa.pl $i
done

perl 3.stat_seq_len.pl > scaffold_len 


perl 4.mk_donor_seq.pl inter > inter_donor_seq.fa
perl 4.mk_donor_seq.pl intra > intra_donor_seq.fa

perl 5.mk_recipent_seq.pl inter > inter_recip_seq.fa
perl 5.mk_recipent_seq.pl intra > intra_recip_seq.fa



# ImmeDB database
for i in `ls *fa|sed -e 's/.fa//'`
do 
	blastn -db ImmeDB_Data1_MGE_sequences.rename.blastdb -query $i.fa -outfmt 6 > $i.MEG.m6
done
for i in `ls *fa|sed -e 's/.fa//'`
do 
	perl a.filter_res.pl $i MGE > $i.MGE.m6.anno
done

for i in `ls *fa|sed -e 's/.anno//'`
do
	perl b.best_hit.pl $i > $i.select
done
