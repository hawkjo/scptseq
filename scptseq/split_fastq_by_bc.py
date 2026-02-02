import logging
import sys
import os
from Bio import SeqIO
from .misc import gzip_friendly_open

log = logging.getLogger(__name__)


class FileHandler(dict):
    def __init__(self):
        pass

    def __getitem__(self, out_fpath):
        if out_fpath not in self.__dict__:
            self.__dict__[out_fpath] = open(out_fpath, 'w')
        return self.__dict__[out_fpath]

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for fh in self.__dict__.values():
            fh.close()

def split_fastq_by_bc(fastq_fpaths, out_dir):
    if os.path.exists(out_dir):
        sys.exit('Output directory already exists')
    os.mkdir(out_dir)

    with FileHandler() as fh:
        for fastq_fpath in fastq_fpaths:
            for rec in SeqIO.parse(gzip_friendly_open(fastq_fpath), 'fastq'):
                bc = str(rec.id).split('_')[0]
                out_fpath = os.path.join(out_dir, f'{bc}.fq')
                SeqIO.write(rec, fh[out_fpath], 'fastq')

def splitfastqs(arguments):
    split_fastq_by_bc(arguments.fastq_files, arguments.output_dir)   
