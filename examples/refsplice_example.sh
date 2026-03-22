#!/bin/bash

scptseq refsplice \
    example_data/TX46_Prosalpha3_ctrl_bams `# <control_bam_dir>: dir with aligned control files by barcode` \
    TX46_Prosalpha3             `# <run_name>: user-selected name for the experiment` \
    Prosalpha3                  `# <gene_name>: gene name, matching <target_info_file> entry` \
    example_data/targets.yml    `# <target_info_file>: target info, see file for format` \
    --results-dir=results       `# <results_dir>: results directory for mutation calling pipeline` \
    -vvv                        `# [-v | -vv | -vvv]:  verbosity`
