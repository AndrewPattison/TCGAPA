#!/usr/bin/python
"""This pipeline is desinged to donwlaod TCGA bam files and count features mapping to known APADB sites. Also 
generates coverage plots"""

# To use TCGAPA GDC_manifest GDC_token APA_site_file cores (default = 4)

import sys
import subprocess  
from subprocess import call
import glob
import os
import re
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool


# Do one file at a time. 
# Make each entry its own manifest so I can go one at a time. 
def make_manifest(manifest):
    with open(manifest, 'rw') as man:
        man2 = man.readlines()
        head = man2[0]
        count = 1
        for line in man2[1:]:
            with open("Manifest/mans/" + "manifest" + str(count) + ".mans", 'w') as man_write:
                man_write.write(man2[0] + line)
            count = count + 1    

# Downlaod and rename the bam and bai files
def setup_bams(custom_manifest):
    number = re.sub ("[^0-9]", "", custom_manifest)

    print "downloading " + number

    call(["mkdir", "downloaded_bams/"+ str(number)])

    code = os.system ("GDC_scripts/gdc-client download --debug -m " + custom_manifest + " -t " + token + " -d downloaded_bams/" + str(number) + " > log" + str(number) + " 2> " + str(number) + "log_error.txt" )
    print code
    bam_file = glob.glob("downloaded_bams/" + str(number) + "/*/" + "*.bam")[0]
    
    p = subprocess.Popen("samtools view -H " + bam_file + " | grep SM:TCGA | cut -f 3", shell = True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = p.stdout.readlines()[0].rstrip(' \t\n\r')[3:]
    print out
    index = glob.glob("downloaded_bams/" + str(number) + "/*/" + "*.bai")[0].rstrip(' \t\n\r')
    print index

    call(["mkdir", "new_bams/"+ str(number)])

    os.rename(bam_file, "new_bams/" + str(number) + "/" + out + ".bam")
    os.rename(index, "new_bams/" + str(number) + "/" + out + ".bam.bai")

    # Run feature counts
    run_r_feat_counts(number)
    # Delete the bam files
    cleanup(number)

# Extract the coverage from the bam files

def run_r_feat_counts(number):

    command = 'Rscript'

    r_sc_path = 'GDC_scripts/run_fc.R'

    bam_file = glob.glob("new_bams/" + str(number) + "/" + "*.bam")[0]

    print "bam for R is " + bam_file

    args = bam_file

    cmd = [command, r_sc_path] 

    cmd.append(args)

    cmd.append(str(number))

    x = subprocess.check_output(cmd, universal_newlines = True)

    print x

def cleanup(number):
    
    call(["rm", "-r", "downloaded_bams/" + str(number)])
    call(["rm", "-r", "new_bams/" + str(number)])
    call(["rm", "Manifest/mans/" + "manifest" + str(number) + ".mans"])
    call(["rm","log" + str(number)])
    call(["rm",str(number) + "log_error.txt"])
     

def run_it_all():
    files = glob.glob("Manifest/mans/*.mans")
    pool = ThreadPool(4)
    pool.map(setup_bams,files)


manifest = sys.argv[1]

# Comment this to start again without remaking the manifests. 
make_manifest(manifest)

token = sys.argv[2]

run_it_all()
#threads = sys.argv[3]

#setup_bams(manifest, token, number)

# view -H 907c31f8-881e-47e3-bf45-a9588f232c3e_gdc_realn_rehead.bam | grep SM:TCGA

# ~/bioinformatics_tools/bedtools2/bin/genomeCoverageBed -split -bg  -ibam $i -g ~/genome_from_bam.txt > $(basename ${i%%-*})





# What I'm running to generate some of these files

#samtools view -H 25def54c-3234-43c3-8eed-88bb27b898db_gdc_realn_rehead.bam | grep SM:TCGA | cut -f 3

#bedtools sort -i HG_38_utr_regions.gff > utr_regions_sorted.gff

#bedtools merge -i utr_regions_sorted.gff > utr_regions_sorted_merged.bed

#~/bioinformatics_tools/bedGraphToBigWig coverage.bedGraph ~/bioinformatics/genome_locations/genome_from_bam.txt cov.bw

# Files that did not work:
# 244