#!/bin/bash

scptseq refsplice \
    TX46_Prosalpha3 \
    Prosalpha3 \
    example_data/targets.yml \
    example_data/TX46_Prosalpha3_ctrl_bams \
    --results-dir=results \
    -vvv
