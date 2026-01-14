# Setup
## Prerequisites 
The setup script assumes that you have `conda` installed, and a environment called `snakemake` activated.

To install and setup conda, follow [the instructions here](https://www.anaconda.com/docs/getting-started/miniconda/install#linux).
Once installed, create and activate the the environment `snakemake` as follows:

```
conda create -n snakemake
conda activate snakemake
```

## Downloading the repo

To download the repo, we can download by clicking the code button, downloading the zip file and extracting the archive. Alternatively, you can also download the repo using the `git clone` functionality on command line.

## Pre-install
Next, run the setup script:

```
./setup.sh
```

This will install the correct dependencies into the environment. The setup will take ten to fifteen minutes.

We've included an example file with a fasta of H5N1 sequences (`H5.fa`) along with a fasta of human transcripts (`HUMAN.fa`) .

To verify the install has completed successfully, we can run the example Snakefile to test the functionality of the pipeline and generate H5N1 primers that are minimally off-target against human transcripts.

If GPU support is available on the machine, you can include the flag `--resources gpu=1` in the call to snakemake, but the performance gains are not drastic and CPU performance is quite acceptable.

### Single-plex Example
Thus, we can run the following command to execute the pipeline for a single plex reaction:

```
snakemake -s Snakefile.example --cores all
```

If configured correctly, the output should result in a `.csv` named for the input sequence file in a newly created `final` directory in the main directory upon completion. This run takes approximately twenty minutes for completion.


### Multiplex Example
We've also set up an example of a multiplex reaction in the same Snakefile. To run the pipeline for the multiplex reaction, we use the following command:

```
snakemake -s Snakefile.example --config multiplex=1 --cores all
```
This will also result in a file named `multiplex_output.csv` in the `/final` directory, which contains the top 3 candidates of each target, and minimizing the activity of the primer against the other sequences, as well as the probe sequences for the multiplex reaction (PROBE DESIGN TO BE IMPLEMENTED).

# Run

To set up runs on your arbitrary sequences of choice, the script assumes a directory structure as follows:

```
.
└── target_seqs/
    └── original/
        ├── target1.fa
        ├── target2.fa
        └── offtarget.fa
```

The FASTA files in `./target_seqs/original` are multi-sequence, unaligned FASTAs that are a pool of related sequences which must be named `*.fa`. You must put all FASTA files into the `./target_seqs/original` directory in order for the pipeline to access them.

To specify which inputs are used in the operation of the pipeline, we specify which sequences we are maximizing activity and coverage for as well as minimizing off-target activity for through the following section in our Snakefile. It is important that the files follow the naming convention strictly. For the target file `foo.fa`, the file extension must be `.fa`, and the input string for the Snakefiles be `foo`. Cross reactivity (CROSS) is treated the same as HOST, in which we attempt to minimize the activity of any matches of the designed primers.

## Singleplex qPCR

In this scenario, we've setup our Snakemake Snakefile to design primers with high sensitivity for `target1.fa` sequences, but minimal activity on `offtarget.fa`. This is a singleplex run, where we design primers and probes for the single target sequence. We specify them in the singleplex region of the Snakefile.

```
...
TARGETS      = [ 'target1' ] 
CROSS       = [ ] 
...
HOST        = [ 'offtarget' ]
...

```

To run this command, we can use the command `snakemake -s YOUR_SNAKEFILE_HERE --cores all`. Multiplex is disabled by default.

## Multiplex qPCR
In addition, we can also design primers and probe sets which have minimal interference with each other. 

To specify which targets we want to design a panel for, we specify them in the arguments for our Snakefile.

```
PANEL = ['target1', 'target2']
```
We can specify which host this belongs to with the same `HOST` field that we used to specify in the singleplex reaction.

To run this command, we can use the command `snakemake -s YOUR_SNAKEFILE_HERE -config multiplex=1 --cores all`.

# Retraining Models

Pretrained models are located in the `src` directory. However, training and data processing scripts are available in the `training` directory for reference.

To run the training process for these scripts, we can run each script to generate the appropriate data processing/produce an equivalent model to our provided pretrained models (raw dataset available upon request).
