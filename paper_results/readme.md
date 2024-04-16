## The scripts for the results generation in the paper


### Method validation 
| Files  | Description |
| :------------- | :------------- |
| simulation.py| simulate benchmark data to evaluate LocalHGT|
|generate_run_scripts.py| generate batch running scripts for LocalHGT in evaluation|
|evaluation.py| evaluate LocalHGT's accuracy in benchmark data|
|count_table_empty_with_k.py| evaluate hash collision after kmer counting in the complex sample|
|validate_bkp_match.py|Validate the identified HGT events using matched long-read data|
|run_lemon.sh| Run LEMON|



### Real-data analyses
| Files  | Description |
| :------------- | :------------- |
|build_UHGG_reference.py| construct the UHGG reference|
| basic_statistics.py | Characterize HGT breakpoint distribution among samples |
| analyze_transfer_gene.py | Analyze the function of HGT-related genes and analyze the HGT transfer patterns  |
| microhomology.py | Test the enrichment of microhomology in HGT breakpoint junctions |
| mechanism.py| Assign the mutational mechanism for HGT events |
|ana_time_lines.py| Examine whether HGT can serve as fingerprints in time-series cohort| 
| kegg_enrichment.py | Given KO list, perform KEGG enrichment analysis |
| association_study.py | Analyze the functional association between HGT and diseases |
| HGT_classifier.py | Identify differential HGTs, construct the classifier for each disease |
| CRC_LODO_Analysis_v2.py | Evaluate the integration of HGT and abundance biomarkers in eight CRC cohorts |
| additional_validation.py  | Get information on the independent CRC cohort|
| HGT_network.py | Perform individual HGT network analyses |
| extract_phenotype.py | Extract the phenotype information for real data |
| network_analysis/* | Perform populational HGT network analyses |
| mge_annotation/* | Perform mobile gene element (MGE) analyses in HGT events |
|[Detected HGTs in real data](https://doi.org/10.5281/zenodo.10906354)|HGT breakpoints and events detected by LocalHGT from all real samples. |









