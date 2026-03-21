"""
scptseq: Computational suite for scPT-seq data

Usage:
  scptseq count        <run_name> <gene_name> <target_info_file> <genome_file> <control_bam_dir> <perturbed_bam_dir> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  scptseq splitfastqs  <fastq_files> --output-dir=<output_dir> [-v | -vv | -vvv]
  scptseq refsplice
  scptseq call

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  count         Count the haplotyped mutation information per cell 
  splitfastqs   Split fastq files by barcode, one file per barcode
"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .preprocessing import haplotyped_mutation_preprocessing
from .split_fastq_by_bc import splitfastqs


def main(**kwargs):
    docopt_args = docopt(__doc__, version=__version__)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
        'count': haplotyped_mutation_preprocessing,
        'splitfastqs': splitfastqs,
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
