# scPT-seq: Single-Cell Perturbation and Transcriptome Sequencing

A suite of tools for processing scPT-seq data, as described in

### Direct detection of CRISPR mutations and transcriptional responses at single cell resolution in vivo

**John A. Hawkins, Siamak Redhai, Svenja Leible, Mireia Osuna Lopez, Hilal Ozgur, Tianyu Wang, Michaela Holzem, Michael Boutros, Oliver Stegle**

*bioRxiv*. Dec 24, 2025. https://doi.org/10.64898/2025.12.23.696319


### Installation

For easiest installation, use pip:

```
pip install scptseq
```

The following instructions should also work for manual installation across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/hawkjo/scptseq.git
```

Optionally, you can install into a virtual environment (recommended):

```
cd scptseq
python -m venv envscptseq
. envscptseq/bin/activate
```

Now install required packages listed in `setup.py` and install scptseq with `setup.py`:

```
python -m pip install .
```

### Usage

```
Usage:
  scptseq splitfastqs  <fastq_files> --output-dir=<output_dir> [-v | -vv | -vvv]
  scptseq refsplice    <control_bam_dir> <run_name> <gene_name> <target_info_file> [--results-dir=<results_dir>] [-v | -vv | -vvv]
  scptseq count        <genome_file> <perturbed_bam_dir> <run_name> <gene_name> <target_info_file> [--results-dir=<results_dir>] [-v | -vv | -vvv]
  scptseq call         <run_name> <gene_name> <target_info_file> [--results-dir=<results_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  splitfastqs   Split fastq files by barcode, one file per barcode
  refsplice     Find all splice junctions in control data, including non-canonical
  count         Count the haplotyped mutation information per cell 
  call          Call mutations per cell 
```

### Top level functions
#### splitfastqs
A helper function for separating the output of
[freediv10xbarcodes](https://github.com/hawkjo/freediv10Xcellbcs) into separate files, one per
barcode. Done for both control and perturbed files, which are then aligned using mapping software
of choice.

#### refsplice
Used for identifying all splice junctions in a control sample, both canonical and non-canonical, 
to distinguish perturbation-caused splice junction changes from low-frequency non-canonical
junctions.

#### count
Preprocesses the aligned perturbation reads to collate key statistics per mutation per cell per
haplotype.

#### call
Perform the mutation calling, outputing mutation status per chromosome per cell, as well as initial
overview statistics and figures.

### Examples
Examples demonstrating use of each of the top level functions with further descriptions are given
in the `examples` folder. 
