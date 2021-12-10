"""
The script to download UHGG genomes.

Sometimes, we fail to download some species, we can use this script iteratively, 
until all species are downloaded.

After downloading, please ungzip all the genome files, and cat them to a single reference file.

Sometimes the fasta formate might have error, we can use the tool "seqkit seq" to reformate it.

"""


import os

downloaded_num = 0
for i in range(1, 4645):
    pre = int(i / 100)
    pre = str(pre).zfill(3)
    i = str(i).zfill(5)
    link = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue/MGYG-HGUT-%s/MGYG-HGUT-%s/genome/MGYG-HGUT-%s.fna"%(pre, i, i)
    file = 'genomes/MGYG-HGUT-%s.fna'%(i)
    #print (i, pre, link)
    if not os.path.isfile(file):
        print ('%s is downloaded now!'%(file))
        down = "wget %s -P genomes/"%(link)
        os.system(down)
        #print (down)
    else:
        downloaded_num += 1
        print ('%s was downloaded already!'%(file))
print (downloaded_num, 'were downloaded.')
