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

# parser = argparse.ArgumentParser(description=desc, epilog=epi,
#                                  formatter_class=CustomFormatter)
# def make_db(args):
#     # Implement your 'cmd1' subcommand here
#     make_db.main(args)

def run(args):
    # Implement your 'cmd2' subcommand here
    print("Command 2 executed with arguments:", args)

def main():
    parser = argparse.ArgumentParser(
        description="A pipeline to rapidly reconstruct the SSU rRNAs",
        formatter_class=CustomFormatter
    )
    subparsers = parser.add_subparsers()

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

    # subcommand: run
    desc = 'Run phyloFlash pipeline'
    epi = """DESCRIPTION:
    Run phyloFlash pipeline.
    """
    parser_run = subparsers.add_parser("run", formatter_class=CustomFormatter,
                                       description = desc, epilog = epi)
    parser_run.set_defaults(func=run)
    ## Add arguments specific to cmd2, e.g.
    # parser_cmd2.add_argument("-o", "--output", help="Output file", default="output.txt")

    # parse arguments
    args = parser.parse_args()
    ## call subcommand function
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
