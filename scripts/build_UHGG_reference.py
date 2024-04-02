#!/usr/bin/env python3

"""
The script to download UHGG genomes.

Sometimes, we fail to download some species, we download iteratively until all species are downloaded.

After downloading, cat all the genomes to a single reference file.

Sometimes the fasta formate might have error, we use the tool "seqkit seq" to refine the reference file.

"""


import os
import argparse
import sys
import shutil
from shutil import which


def check_software_availability(software_name):
    if shutil.which(software_name) is not None:
        print(f"{software_name} is available in the system environment.")
    else:
        print(f"{software_name} is not available in the system environment.")
        print ("please construct the conda environment using the given *.yml file.")
        sys.exit()

def get_link(i):
    """
    Given an index, infer the download link for the corresponding genome, and infer the genome name.
    """
    pre = int(i / 100)
    pre = str(pre).zfill(3)
    i = str(i).zfill(5)
    link = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue/MGYG-HGUT-%s/MGYG-HGUT-%s/genome/MGYG-HGUT-%s.fna"%(pre, i, i)
    file = '%s/MGYG-HGUT-%s.fna'%(args['b'] , i)
    return link, file

def download(args):
    if not os.path.isdir(args['b']):
        os.system("mkdir %s"%(args['b']))

    downloaded_num = 0
    for i in range(1, 4645):
        link, file = get_link(i)
        #print (i, pre, link)
        if not os.path.isfile(file):
            print ('download %s...'%(file))
            down = "wget %s -P %s/"%(link, args['b'])
            os.system(down)
        else:
            downloaded_num += 1
            # print ('%s was downloaded already!'%(file))

        # if i > 99:  # just for test
        #     print ("Just for testing, only download at most 100 genomes.")
        #     break 
    # print (downloaded_num, 'genomes were downloaded.')

def iterate(args):
    """
    Iteratively download genomes until all the genomes are downloaded.
    """
    iter_times = 1
    while iter_times <= args['m']:
        finished_flag = True
        # test if all the genomes are downloaded
        downloaded_num = 0
        for i in range(1, 4645):
            link, file = get_link(i)
            #print (i, pre, link)
            if not os.path.isfile(file):
                finished_flag = False
            else: 
                downloaded_num += 1
        print ("Iteration times: %s, N.O. of genomes: %s."%(iter_times - 1, downloaded_num))
        if finished_flag:
            print ("All the genomes are downloaded.")
            break
        download(args)
        iter_times += 1

def merge(args):
    """
    Merge all the downloaded genomes into a single fasta file as the reference required for LocalHGT, and refine the fasta format
    """
    command = f"""
    echo Merge the genomes...
    cat {args['b']}/MGYG-HGUT-*.fna > {args['r']}.raw.fasta
    # echo Ensure seqkit is installed.
    seqkit seq {args['r']}.raw.fasta > {args['r']}
    rm {args['r']}.raw.fasta
    samtools faidx {args['r']}
    echo Reference is generated.
    """
    os.system(command)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct the reference from gut-specific UHGG V1 database.", add_help=False, \
    usage="python %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    optional.add_argument("-r", type=str, default="uhgg_v1.rep.fasta", help="<str> Generated reference file.", metavar="\b")
    optional.add_argument("-b", type=str, default="genomes/", help="<str> Folder saves all the downloaded assemblies.", metavar="\b")
    optional.add_argument("-m", type=int, default=10, help="<int> Try this number of times until all the genomes are downloaded.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    args = vars(parser.parse_args())


    if len(sys.argv)==1:
        # print (f"see python {sys.argv[0]} -h")
        os.system(f"python {sys.argv[0]} -h")
    else:
        print ("start building UHGG database...")
        print ("result will be stored in %s"%(args['r']))
        check_software_availability("wget")
        check_software_availability("samtools")
        check_software_availability("seqkit")

        iterate(args)
        merge(args)
