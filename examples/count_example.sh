#!/bin/bash

scptseq count \
    TX46_Prosalpha3 \
    Prosalpha3 \
    example_data/targets.yml \
    example_data/Drosophila_melanogaster.BDGP6.28.dna.toplevel.chrm_2R.fa \
    example_data/TX46_Prosalpha3_bams \
    --results-dir=results \
    -vvv
