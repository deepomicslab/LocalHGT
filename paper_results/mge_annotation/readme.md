### MGEs annotation
| Inputs  | Description |
| :------------- | :------------- |
|identified_event.csv|HGT events output from localHGT|
|genomes/|*.gff.gz from UHGG database can obtained from ```build_UHGG_reference.py```|
|representive_genome.metadata.tsv|[UHGG metadata](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv)|
|ImmeDB_Data1_MGE_sequences.blastdb|blastndb of [ImmeDB](https://jianglabnlm.com/immedb/data/Data1_MGE_sequences.fasta)|


| Scripts  | Description |
| :------------- | :------------- |
| 0.merge_bp.pl | Merge breakpoints of all HGT events|
| 1.note_inter_intra_phylum.pl | Annotated HGT as inter or intra-phylum HGTs |
| 2.extract_fa.pl, 3.stat_seq_len.pl | Generate the scaffold length file |
| 4.mk_donor_seq.pl, 5.mk_recipent_seq.pl | Generate the fasta files of ±5kb flanking sequences around donor/recipient breakpoints |
| a.filter_res.pl, b.best_hit.pl | Select the best hit of every event donor/recipient sequence |
| pipeline_run.sh | Shell scripts to generate all intermedia and final results |