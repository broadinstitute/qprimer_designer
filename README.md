# Setup
## Prerequisites 
The setup script assumes that you have `conda` installed, and a environment called `snakemake` activated. This also assumes `gcc` is accessible through the `$PATH`.

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

This will install the correct dependencies into the environment. The setup will take approximately half an hour.

After installation, we can create a copy of the `Snakefile.template` template file as well as the `params.txt.template` file.

```
cp Snakefile.template Snakefile.example
cp params.txt.template params.txt
```

We then want to edit the file `Snakefile.example` to reflect our run settings. We've included an example file with a fasta of H5N1 sequences (`H5.fa`) along with a fasta of human transcripts (`HUMAN.fa`) .

In your text editor of choice, edit the Snakefile.example file to have the following segment.

```
...
TARGET      = [ 'H5' ] 
CROSS       = [ ] 
HOST        = [ 'HUMAN' ]
...

```

To verify the install has completed successfully, we can run the example Snakefile to test the functionality of the pipeline.
If GPU support is available on the machine, you can set gpu=1, but the performance gains are not drastic and CPU performance is quite acceptable.


```
snakemake -s Snakefile.example --resources gpu=0 --cores all
```

If configured correctly, the output should result in a `.csv` file in a newly created `final` directory in the main directory upon completion. This run takes approximately half an hour for completion.

## Run


To set up runs on your sequences of choice, the script assumes a directory structure as follows:

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

```
...
TARGET      = [ 'target1' ] 
CROSS       = [ ] 
HOST        = [ 'offtarget' ]
...

```

To run this command, we can use the command `snakemake -s YOUR_SNAKEFILE_HERE --resources gpu=0 --cores all`.

### Retraining Models

Pretrained models are located in the `src` directory. However, training and data processing scripts are available in the `training` directory for reference.

To run the training process for these scripts, we can run each script to generate the appropriate data processing/produce an equivalent model to our provided pretrained models (raw dataset available upon request).
