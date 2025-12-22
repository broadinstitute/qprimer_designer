# Setup
## Prerequisites 
The setup script assumes that you have `conda` installed, and a environment called `snakemake` activated. This also assumes `gcc` is accessible through the `$PATH`.

To install and setup conda, follow [the instructions here](https://www.anaconda.com/docs/getting-started/miniconda/install#linux).
Once installed, create and activate the the environment `snakemake` as follows:

```
conda create -n snakemake
conda activate snakemake
```

Next, run the setup script:

```
./setup.sh
```

This will install the correct dependencies into the environment.
To verify the install has completed successfully, we can run the example Snakefile to test the functionality of the pipeline.
GPU acceleration is assumed to be included. TODO: test functionality without GPU

```
snakemake -s Snakefile.example --resources gpu=1 --cores all
```

If configured correctly, the output should result in a `.csv` file in a newly created `final` directory upon completion. In addition, this includes a cleaned FASTA of human transcript sequences for ease of designing diagnostics with minimal human activity.

## Run
To set up custom runs, the script assumes a directory structure as follows:

```
.
└── target_seqs/
    └── original/
        ├── target1.fa
        ├── target2.fa
        └── offtarget.fa
```

The FASTA files in `./target_seqs/original` are multi-sequence, unaligned FASTAs that are a pool of related sequences which must be named `*.fa`.

To specify which inputs are used in the operation of the pipeline, we specify which sequences we are maximizing activity and coverage for as well as minimizing off-target activity for through the following section in our Snakefile template, `Snakefile.template`: It is important that the files follow the naming convention strictly. For the target file `foo.fa`, the file extension must be `.fa`, and the input string for the Snakefiles be `foo`. TODO: explain CROSS

```
...
TARGET      = [ 'target1' ] 
CROSS       = [ ] 
HOST        = [ 'offtarget' ]
...

```

To run this command, we can use the command `snakemake --resources gpu=1 --cores all`.

### Training

Pretrained models are located in the `src` directory. However, training and data processing scripts are available in the `training` directory.
