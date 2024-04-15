## The scripts for the results generation in the paper


### Method validation 
| Files  | Description |
| :------------- | :------------- |
| simulation.py| simulate benchmark data to evaluate LocalHGT|
|generate_run_scripts.py| generate batch running scripts for LocalHGT in evaluation|
|evaluation.py| evaluate LocalHGT's accuracy in benchmark data|
|count_table_empty_with_k.py| evaluate hash collision after kmer counting in the complex sample|
|validate_bkp_match.py|Validate the identified HGT events using matched long-read data|
| batch_validate_match.py | Perform the event validation in a batch manner |
|run_lemon.sh| Run LEMON|



### Real-data analyses
| Files  | Description |
| :------------- | :------------- |
|build_UHGG_reference.py| construct the UHGG reference|
| basic_statistics.py | Characterize HKP distribution among samples |
| analyze_transfer_gene.py | Analyze the function of HGT-related genes and analyze the HGT transfer patterns  |
| microhomology.py | Test the enrichment of microhomology in HGT breakpoint junctions |
|  mechanism.py| Assign the mutational mechanism for HGT events |
|ana_time_lines.py| Examine whether HGT can serve as fingerprints in time-series cohort| 
| kegg_enrichment.py | Given KO list, perform KEGG enrichment analysis |
| association_study.py | Analyze the functional association between HGT and diseases |
| HGT_classifier.py | Identify differential HGTs, construct the classifier for each disease |
| CRC_LODO_Analysis_v2.py | Evaluate the integration of HGT and abundance biomarkers in eight CRC cohorts |
| additional_validation.py  | Get information on the independent CRC cohort|
| HGT_network.py | Perform individual HGT network analyses |
| network_analysis/* | Perform populational HGT network analyses |
| mge_annotation/* | Perform mobile gene element (MGE) analyses in HGT events |
|[Detected HGTs in real data](https://doi.org/10.5281/zenodo.10906354)|HGT breakpoints and events detected by LocalHGT from all real samples. |




#### MGEs annotation
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




