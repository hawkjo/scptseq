"""
scptseq: Computational suite for scPT-seq data

Usage:
  scptseq splitfastqs  <fastq_files> --output-dir=<output_dir> [-v | -vv | -vvv]
  scptseq refsplice    <run_name> <gene_name> <target_info_file> <control_bam_dir> [--results-dir=<results_dir>] [-v | -vv | -vvv]
  scptseq count        <run_name> <gene_name> <target_info_file> <genome_file> <perturbed_bam_dir> [--results-dir=<results_dir>] [-v | -vv | -vvv]
  scptseq call         <run_name> <gene_name> <target_info_file> [--results-dir=<results_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  splitfastqs   Split fastq files by barcode, one file per barcode
  count         Count the haplotyped mutation information per cell 
  call          Call mutations per cell 
"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .preprocessing import haplotyped_mutation_preprocessing, find_all_ref_splice_junctions
from .split_fastq_by_bc import splitfastqs
from .call_mutations import call_mutations


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
        'splitfastqs': splitfastqs,
        'refsplice': find_all_ref_splice_junctions,
        'count': haplotyped_mutation_preprocessing,
        'call': call_mutations,
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
