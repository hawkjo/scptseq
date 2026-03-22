#!/bin/bash

scptseq splitfastqs \
    example_data/TX46_Prosalpha3_samp.fq `# <fastq_files>: comma separated list of fastq files` \
    --output-dir=output_splitfastqs `# <output_dir>: output directory for split fastq files` \
    -vvv `# [-v | -vv | -vvv]:  verbosity`
