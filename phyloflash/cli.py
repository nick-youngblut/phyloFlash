#!/usr/bin/env python
# import
## batteries
import logging
import argparse
## package
from phyloflash import core
from phyloflash import make_db

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# parser description formatting
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def cmd_make_db(subparsers):
     # subcommand: make-db
    desc = 'Create database for phyloFlash'
    epi = """DESCRIPTION:
    Download and format phyloFlash database files.
    All output files are stored in the output directory.
    """
    parser_make_db = subparsers.add_parser("make-db", formatter_class=CustomFormatter,
                                           description = desc, epilog = epi)
    parser_make_db.set_defaults(func=make_db.main)
    ## add arguments
    parser_make_db.add_argument("outdir", type=str, 
                                help = "Output directory path")
    parser_make_db.add_argument("-s", "--skip-sortmerna", action='store_true', default=False,
                                help = "Skip sortmerna database")
    parser_make_db.add_argument("-t", "--threads", type=int, default=4,
                                help = "Number of threads to use")
    parser_make_db.add_argument("-m", "--memory", type=int, default=12,
                                help = "Memory limit in GB")
    parser_make_db.add_argument("-U", "--univec-url", type=str, 
                                default='https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec',
                                help = "UniVec database URL")
    parser_make_db.add_argument("-S", "--silva-url", type=str, 
                                default='https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz',
                                help = "SILVA database URL")
    parser_make_db.add_argument("-d", "--debug", action='store_true', default=False,
                                help = "Debug mode")
    parser_make_db.add_argument("-N", "--num-lines", type=int, default=1000,
                                help = "Max number of lines to used from the SILVA database, if debug=True")

def cmd_run(subparsers):
       # subcommand: run
    desc = 'Run phyloFlash pipeline'
    epi = """DESCRIPTION:
    Run phyloFlash pipeline.
    """
    parser_run = subparsers.add_parser("run", formatter_class=CustomFormatter,
                                       description = desc, epilog = epi)
    parser_run.set_defaults(func=core.main)
    ## Add arguments specific to cmd2, e.g.
    parser_run.add_argument('--read1', type=str, help='Forward read file')
    parser_run.add_argument('--read2', type=str, help='Reverse read file')
    parser_run.add_argument('--lib', type=str, help='Output file basename')
    parser_run.add_argument('--db-home', type=str, help='phyloFlash DB folder')
    parser_run.add_argument('--interleaved', action='store_true', help='Interleaved reads')
    parser_run.add_argument('--read-length', type=int, help='Read length')
    parser_run.add_argument('--read-limit', type=int, help='Read limit')
    parser_run.add_argument('--amp-limit', type=int, help='Amplimit')
    parser_run.add_argument('--max-insert', type=int, help='Maxinsert')
    parser_run.add_argument('--id', type=int, help='Read mapping identity')
    parser_run.add_argument('--cluster-id', type=int, help='Clustering identity')
    parser_run.add_argument('--tax-level', type=int, default=4, help='Taxon report level')
    parser_run.add_argument('--threads', type=int, default=1, help='threads')
    parser_run.add_argument('--html', action='store_true', help='HTML flag')
    parser_run.add_argument('--treemap', action='store_true', help='Treemap flag')
    parser_run.add_argument('--crlf', action='store_true', help='CRLF flag')
    parser_run.add_argument('--decimal-comma', action='store_true', help='Decimal comma flag')
    parser_run.add_argument('--sortmerna', action='store_true', help='Use sortmerna')
    parser_run.add_argument('--evalue-sortmerna', type=str, help='E-value sortmerna')
    parser_run.add_argument('--emirge', action='store_true', help='Use emirge')
    parser_run.add_argument('--skip-emirge', action='store_true', help='Skip emirge')
    parser_run.add_argument('--skip-spades', action='store_true', help='Skip spades')
    parser_run.add_argument('--trusted', type=str, help='Trusted contigs')
    parser_run.add_argument('--poscov', action='store_true', help='Positional coverage flag')
    parser_run.add_argument('--sc', action='store_true', help='SC flag')
    parser_run.add_argument('--zip', action='store_true', help='Zip flag')
    parser_run.add_argument('--log', action='store_true', help='Save log flag')
    parser_run.add_argument('--keeptmp', action='store_true', help='Keep temporary files')
    parser_run.add_argument('--everything', action='store_true', help='Everything flag')
    parser_run.add_argument('--almost-everything', action='store_true', help='Almost everything flag')
    parser_run.add_argument('--tophit', action='store_true', help='Top hit flag')
    parser_run.add_argument('--check-env', action='store_true', help='Check environment flag')
    parser_run.add_argument('--outfiles', action='store_true', help='Output description flag')

def main():
    parser = argparse.ArgumentParser(
        description="A pipeline to rapidly reconstruct the SSU rRNAs",
        formatter_class=CustomFormatter
    )
    subparsers = parser.add_subparsers()

    # subcommands
    cmd_make_db(subparsers)
    cmd_run(subparsers)

    # parse arguments
    args = parser.parse_args()
    ## call subcommand function
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
