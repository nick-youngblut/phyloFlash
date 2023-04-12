#!/usr/bin/env python
# import
## batteries
import os
import re
import shutil
import pathlib
import platform
import importlib.resources


def which(exe):
    """
    Check that executable is in PATH
    """
    if shutil.which(exe) is None:
        raise OSError(f"Executable '{exe}' not found in PATH")

def which_barrnap():
    """
    Check if barrnap nhmmer executables ares in package data path
    """
    data_path = importlib.resources.path('phyloflash.data', 'blank')
    data_dir = None
    with data_path as path:
        data_dir = str(pathlib.Path(path).parent)
        for exe in ['nhmmer-darwin', 'nhmmer-linux']:
            exe_path = os.path.join(data_dir, 'barrnap-HGV', exe)
            if not os.path.exists(exe_path):
                raise OSError(f"Executable '{exe_path}' not found in package data path")
    return data_dir

def check_database(db_home: str, use_sortmerna=False):
    """
    Check whether the required database files are present in "db_home".
    """
    emirge_db   = "SILVA_SSU.noLSU.masked.trimmed.NR96.fixed"
    vsearch_db  = "SILVA_SSU.noLSU.masked.trimmed"
    sortmerna_db = emirge_db
    required = ['ref/genome/1/summary.txt',
                #f'{emirge_db}.fasta',
                'SILVA_SSU.noLSU.masked.trimmed']
    if use_sortmerna:
        required += [f'${sortmerna_db}.bursttrie_0.dat',
                     f'${sortmerna_db}.acc2taxstring.hashimage']
    required = [os.path.join(db_home, x) for x in required]
    return required

def main(args):
    # environment
    ## check for executables needed by phyloFlash
    for exe in ['bbmap.sh', 'reformat.sh', 'vsearch', 'mafft', 
                'fastaFromBed', 'sed', 'grep', 'awk', 'cat']:
        which(exe)
    ## check that nhmmer is in package data path
    data_dir = which_barrnap()
    
    ## Check operating system to decide which nhmmer to use
    os_name = platform.system()
    if os_name.lower() not in ('linux', 'darwin'):
        raise OSError(f'Operating system not supported: {os_name}')
    
    # database
    ## verify precence of database
    check_database(args.db_home, use_sortmerna=args.sortmerna)
    
    