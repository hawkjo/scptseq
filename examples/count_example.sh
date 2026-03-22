#!/bin/bash

scptseq count \
    example_data/Drosophila_melanogaster.BDGP6.28.dna.toplevel.chrm_2R.fa \ # <genome_file>
    example_data/TX46_Prosalpha3_bams \ # <perturbed_bam_dir>: dirwith aligned perturbed files by barcode  
    TX46_Prosalpha3 \           # <run_name>: user-selected name for the experiment  
    Prosalpha3 \                # <gene_name>: gene name, matching <target_info_file> entry  
    example_data/targets.yml \  # <target_info_file>: target info, see file for format  
    --results-dir=results \     # <results_dir>: results directory for mutation calling pipeline  
    -vvv                        # [-v | -vv | -vvv]:  verbosity  
