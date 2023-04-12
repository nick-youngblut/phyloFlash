#!/usr/bin/env python
# import
## batteries
import os
import re
import sys
import gzip
import pickle
import random
import shutil
import logging
import urllib.request
from subprocess import Popen, PIPE


# Dict to map IUPAC ambiguous bases to [ATGC]
IUPAC_DECODE = {
    'X' : 'ACGT',
    'R' : 'AG',
    'Y' : 'CT',
    'M' : 'CA',
    'K' : 'TG',
    'W' : 'TA',
    'S' : 'CG',
    'B' : 'CTG',
    'D' : 'ATG',
    'H' : 'ATC',
    'V' : 'ACG',
    'N' : 'ACTG'
}

def run_job(cmd: str) -> None:
    """
    Run a shell command and check for errors.
    cmd: str, shell command
    """ 
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    ## check for errors
    rc = p.returncode
    if rc != 0:
        raise ValueError(f'Error running {cmd}: {err.decode()}')
    
def which(exe: str) -> None:
    """
    Check if an executable is in PATH.
    exe: str, executable name
    """
    if shutil.which(exe) is None:
        raise ValueError(f'{exe} not found in PATH')
    
def univec_download(univec_url: str, outdir: str, debug=False):
    """
    Download the latest version of the univec database from ncbi.
    univec_url: str, URL to download the univec database
    outdir: str, output directory
    debug: bool, if True, do not download the database if it already exists
    """
    logging.info('Downloading UniVec database from NCBI...')
    univec_file = os.path.join(outdir, 'UniVec')
    if debug is True and os.path.isfile(univec_file):
        return univec_file
    urllib.request.urlretrieve(univec_url, univec_file)
    return univec_file
    
def silva_download(silva_url: str, outdir: str, debug=False):
    """
    Download the latest version of the SILVA SSU RefNR database from www.arb-silva.de
    """
    logging.info('Downloading SILVA database from www.arb-silva.de...')
    silva_file = os.path.join(outdir, 'SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz')
    if debug is True and os.path.isfile(silva_file):
        return silva_file
    urllib.request.urlretrieve(silva_url, silva_file)
    return silva_file

def silva_uncompress(silva_file: str, outdir: str, num_lines=None) -> str:
    """
    uncompress the SILVA database (fasta) file
    """
    logging.info('Uncompressing the SILVA database...')
    out_file = os.path.join(outdir, 'SILVA_SSU.fasta')
    with gzip.open(silva_file, 'rb') as inF:
        with open(out_file, 'w') as outF:
            for i,line in enumerate(inF):
                outF.write(line.decode())
                if num_lines is not None and i >= num_lines:
                    break
    return out_file

def run_barrnap(cmd: str, barrnap_results: list, domain: str) -> None: 
    """
    Run barrnap_HGV and extract sequences with potential LSU contamination
    """
    logging.info(f'Running barrnap_HGV for the "{domain}" domain...')
    # run command
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    ## check for errors
    rc = p.returncode
    if rc != 0:
        raise ValueError(f'Error running {cmd}: {err.decode()}')
    # extract LSU contamination in SSU RefNR
    regex = re.compile(r'23S_rRNA|28S_rRNA')
    for line in output.decode().split('\n'):
        line = line.split('\t')
        if line[0] == '' or len(line) < 8:
            continue
        if regex.search(line[8]):
            barrnap_results.add(line[0])

def subset_fasta(silva_file: str, n=1000) -> str:
    """
    Debug: subset the SILVA database to n sequences
    """
    logging.warning('NOTE: Using a subset of the SILVA database for debugging purposes...')
    out_file,ext = os.path.splitext(silva_file)
    out_file = out_file + '.subset' + ext
    with open(silva_file, 'r') as inF, open(out_file, 'w') as outF:
        for i, line in enumerate(inF):
            outF.write(line)
            if i == n:
                break
    return out_file

def remove_LSU_contamination(silva_file: str, threads=1) -> str:
    """
    Remove sequences with potential LSU contamination
    """
    logging.info('Removing sequences with potential LSU contamination...')
    # check if barrnap_HGV is in PATH
    base_dir = os.path.split(os.path.realpath(__file__))[0]
    exe = os.path.join(base_dir, 'barrnap-HGV', 'bin', 'barrnap_HGV')
    which(exe)
    # run barrnap_HGV on each domain
    barrnap_results = set()
    for domain in ['bac', 'arch', 'euk']:
        cmd = f'{exe} --kingdom {domain} --threads {threads} --evalue 1e-10 --gene lsu --reject 0.01 {silva_file}'
        run_barrnap(cmd, barrnap_results, domain)
    # remove SILVA sequences with potential LSU contamination
    out_file,ext = os.path.splitext(silva_file)
    out_file = out_file + '.noLSU' + ext
    to_keep = False
    with open(silva_file, 'r') as inF, open(out_file, 'w') as outF:
        for line in inF:
            if line.startswith('>'):
                if line.lstrip('>') in barrnap_results:
                    to_keep = False
                else:
                    to_keep = True
            if to_keep is True:
                outF.write(line)
    # remove the original SILVA file
    #os.remove(silva_file)
    return out_file  

def mask_repeats(silva_file: str, threads: int, memory: int) -> str:
    """
    Mask repetitive regions in the SILVA database
    """
    logging.info('Masking repetitive regions in the SILVA database...')
    # check if bbmask is in PATH
    exe = 'bbmask.sh'
    which(exe)
    out_file = 'SILVA_SSU.noLSU.masked.fasta'
    # run bbmask
    cmd = f'{exe} -Xmx{memory}g overwrite=t threads={threads} in={silva_file} out={out_file}'
    cmd += ' minkr=4 maxkr=8 mr=t minlen=20 minke=4 maxke=8 fastawrap=0'
    ## run command
    run_job(cmd)
    return out_file

def univec_trim(univec_file: str, silva_file: str, threads: int, memory: int, min_length=800) -> str:
    """
    Run bbduk to trim sequences with UniVec contamination
    """
    logging.info('Trimming sequences with UniVec contamination...')
    # check if bbduk is in PATH
    exe = 'bbduk.sh'
    which(exe)
    # run bbduk
    out_file = 'SILVA_SSU.noLSU.masked.trimmed.fasta'
    cmd = f'{exe} -Xmx{memory}g threads={threads} overwrite=t ref={univec_file}'
    cmd += f' fastawrap=0 overwrite=t ktrim=r ow=t minlength={min_length} mink=11 hdist=1'
    cmd += f' in={silva_file} out={out_file} stats={out_file}.UniVec_contamination_stats.txt'
    ## run command
    run_job(cmd)
    return out_file

def make_vsearch_udb(silva_file: str, threads=1) -> str:
    """
    Make a vsearch database from the SILVA database
    """
    logging.info('Making a vsearch database from the SILVA database...')
    exe = 'vsearch'
    which(exe)
    # create shell command
    out_file = 'SILVA_SSU.noLSU.masked.trimmed.udb'
    cmd = f'{exe} --threads {threads} --notrunclabels --makeudb_usearch {silva_file} --output {out_file}'
    ## run command
    run_job(cmd)
    return out_file

def cluster(silva_file: str, seqid=0.99, threads=1) -> str:
    """
    Cluster sequences in the SILVA database with vsearch
    """
    logging.info('Clustering sequences in the SILVA database...')
    # check if vsearch is in PATH
    exe = 'vsearch'
    which(exe)
    # set output file
    out_file,ext = os.path.splitext(silva_file)
    id = 'NR' + str(int(round(seqid* 100, 0)))
    out_file = f'{out_file}.{id}.fasta'
    # create shell command
    cmd = f'{exe} --threads {threads} --cluster_fast {silva_file} --id {seqid} --centroids {out_file} --notrunclabels'
    ## run command
    run_job(cmd)
    return out_file

def fasta_copy_iupac_randomize(silva_file: str) -> str:
    """
    Creates a normalized FASTA file from FASTA file $source.
    * removes alignment characters ("." and "-")
    * converts all bases to uppercase
    * convert RNA bases to DNA bases
    * replace IUPAC coded ambiguous base with a base randomly chosen from the set of options
      * i.e. replaces B with C, T or G
    """
    logging.info('Formatting SILVA database...')
    # output file
    out_file,ext = os.path.splitext(silva_file)
    out_file = out_file + '.fixed' + ext
    # format fasta
    regex1 = re.compile(r'[.-]')
    regex2 = re.compile(r'([^AGCT])')
    with open(silva_file) as inF, open(out_file, 'w') as outF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            # maintain header as-is
            if line.startswith('>'):
                outF.write(line + '\n')
                continue
            # remove alignment, uppercase, turn U into T
            line = regex1.sub('', line.upper().replace('U', 'T'))
            # replace any character (x) in line with a random selection from IUPAC_DECODE[x]
            line = regex2.sub(lambda x: random.choice(IUPAC_DECODE[x.group(0)]), line)
            # write to file
            outF.write(line + '\n')
    return out_file
    
def bbmap_db(silva_file: str, out_dir: str, threads=1, memory=4) -> None:
    """
    Create a bbmap database from the SILVA database
    """
    logging.info('Creating a bbmap database from the SILVA database...')
    # check executable
    exe = 'bbmap.sh'
    which(exe)
    # run bbmap
    cmd = f'{exe} -Xmx{memory}g threads={threads} overwrite=t ref={silva_file} path={out_dir}'
    ## run command
    run_job(cmd)
    
def sortmerna_index(silva_file: str, memory=4) -> None:
    """
    Create a sortmerna index from the SILVA database
    """
    logging.info('Creating a sortmerna database from the SILVA database...')
    # check executable
    exe = 'indexdb_rna'
    which(exe)
    # create shell command
    memory = int(round(memory * 1000,0))
    cmd = f'{exe} -m {memory} --ref {silva_file},SILVA_SSU'
    ## run command
    run_job(cmd)  
    
def hash_SILVA_acc_taxstrings_from_fasta(silva_file: str) -> str:
    """
    * Hash of accession numbers and taxonomy strings from SILVA fasta headers
    * Store as a perl hash image with Storable
    * For later when wrangling sortmerna SAM output to bbmap-like format
    """
    logging.info('Hashing accession numbers and taxonomy strings from SILVA fasta headers...')
    prefix = os.path.splitext(silva_file)[0]
    hash = dict()
    with open(silva_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line.startswith('>'):
                id,taxsplit = line.split(' ',1)
                hash[id] = taxsplit
    # pickle the hash
    out_file = prefix + '.acc2taxstring.hashimage'
    with open(out_file, 'wb') as out_file:
        pickle.dump(hash, out_file)
    return out_file


def main(args: dict) -> None:
    # Debug status
    if args.debug:
        logging.info('NOTE: running in DEBUG mode')
    else:
        args.num_lines = None
        
    # Create database directory
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    
    # Download the latest version of the univec database from ncbi
    univec_file = univec_download(args.univec_url, args.outdir, args.debug)

    # Download latest SSU RefNR from www.arb-silva.de
    silva_file = silva_download(args.silva_url, args.outdir, args.debug)
    
    # Uncompress the SILVA database file
    silva_file = silva_uncompress(silva_file, args.outdir, num_lines = args.num_lines)
    
    # Remove sequences with potential LSU contamination
    silva_file = remove_LSU_contamination(silva_file, threads = args.threads)
    
    # Mask repeats in SILVA SSU sequences
    silva_file = mask_repeats(silva_file, threads=args.threads, memory = args.memory)
    
    # Screen SILVA db against UniVec and trim matching sequences with bbduk
    silva_file = univec_trim(univec_file, silva_file, threads = args.threads, memory = args.memory)
    
    # Use Vsearch to index the SILVA database and create a UDB file
    silva_ubd_file = make_vsearch_udb(silva_file, threads = args.threads)
    
    # Cluster the SILVA database at 99% identity using Vsearch
    silva99_file = cluster(silva_file, seqid = 0.99)
    
    # Format the sequence data
    silva99_file = fasta_copy_iupac_randomize(silva99_file)
    
    # Create bbmap index from SILVA database
    bbmap_db(silva99_file, args.outdir, threads = args.threads, memory = args.memory)
    
    # Cluster the SILVA database at 96% identity using Vsearch
    silva96_file = cluster(silva_file, seqid = 0.96)

    # Format the sequence data
    silva96_file = fasta_copy_iupac_randomize(silva96_file)
    
    # Create sortmerna index from SILVA database
    if not args.skip_sortmerna:
        # sortmerna_index
        sortmerna_index(silva96_file)
        # Create dict of taxonomy strings from SILVA fasta headers
        hash_SILVA_acc_taxstrings_from_fasta(silva96_file)
    

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
    
