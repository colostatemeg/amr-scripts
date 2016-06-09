# amr-scripts
Scripts used in the MEG pipelines

### amrPipeline.sh
This is the standard pipeline script for the MEG group.  It accepts as input
paired raw fastq reads and outputs antimicrobial resistance information and
optionally microbiome/pathway information.  The workflow is as follows:

1. Run Trimmomatic on the raw fastqs, output to trimmed file directory
2. Align the trimmed files to the bovine genome using BWA
3. Filter out reads that align to the bovine genome using SAMTools
4. Align the remaining reads to the master AMR database
5. Run SamRatio to get relevant information on the resistome
6. (Optional) Profile the microbial community with [Kraken](https://ccb.jhu.edu/software/kraken/)
7. (Optional) Profile the microbial community with [MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan)
8. (Optional) Profile the microbial pathways with [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2)

Two logfiles are kept when running amrPipeline.sh:

1. logfile.txt (or whatever you name it)
2. LabNotebook.txt

The logfile is a direct copy of all standard out and standard error from the
pipline.  The LabNotebook file is a high-level summary of each sample's
status as it runs through the pipeline.  The LabNotebook file also contains
information about the versions of the software used and the paths to the
input/output directories when the script was run.  It also timestamps each run.

The pipeline has built in file checking and error handling, so if the script
exits during a part of the pipeline due to an error, you can simply rerun
the script after fixing the error and it won't repeat samples that have
already been run.  Make sure to delete the output file that threw the
error though, since that one will need to be regenerated.

Usage:

* amrPipline.sh -i "raw\_sequence\_reads/\*" \[options\]

Note: you MUST pass the input files in double quotes if you want to
use a wildcard to run more than one sample through at a time.  The
wildcard character must be at the end of the path, so put the fastq files in
their own directory with no subdirectories, then pass the wildcard after
that directory path, as shown in the above example.

Note: Please pass directory paths with a terminal slash,
e.g. /s/bovine/index/projs/canada/analysis/AMR/

Options:

* \-a | \-\-kraken DIR	perform kraken analysis, output to the indicated DIR
* \-h | \-\-help		display the help menu and exit
* \-i | \-\-input FILES	input fastq.gz files (must use quotes with wildcard)
* \-k | \-\-keep		save intermediate files (ex. sam alignment files)
* \-l | \-\-logfile	file name for logfile output, default is logfile.txt
* \-m | \-\-metaphlan DIR	perform metaphlan analysis, output to the indicated DIR
* \-n | \-\-nonhost DIR	directory for nonhost fastq output
* \-o | \-\-output DIR	directory for AMR file output
* \-t | \-\-threads INT	number of threads to use [1]
* \-tr| \-\-trimmed DIR	directory for trimmed file output
* \-u | \-\-humann2	perform humann2 analysis, output to the indicated DIR

Defaults for options:

* \-\-kraken \[False\]
* \-\-input \[none\]
* \-\-keep \[False\]
* \-\-logfile \[logfile.txt\]
* \-\-metaphlan \[False\]
* \-\-nonhost \[non\_bovine\_fastq/\]
* \-\-output \[analysis/AMR/\]
* \-\-threads \[1\]
* \-\-trimmed \[trimmed\_fastq\_files/\]
* \-\-humann2 \[False\]


### subsampleFastq.sh

This is a script for subsampling fastqs for rarefaction analysis.  This script
will do the following:

1. Call Trimmomatic on the input files
2. Subsample the files using the dist.py script to given depths

Two logfiles will be kept for this run:

1. subsample\_logfile.txt (or whatever you name it)
2. LabNotebook.txt

The logfile is a direct copy of all standard out and standard error from the
pipline.  The LabNotebook file is a high-level summary of each sample's
status as it runs through the pipeline.  The LabNotebook file also contains
information about the versions of the software used and the paths to the
input/output directories when the script was run.  It also timestamps each run.

Usage:

* subsampleFastq.sh -i "/s/bovine/.../fastq/\*" -o DIR -d "INTS" [options]

Note: Input files AND input depths MUST be passed in quotes (e.g. -d "0.2 0.5 0.6")

Options:

* \-d | \-\-depths INTS	depths to subsample the fastq file
* \-h | \-\-help	display help and exit
* \-i | \-\-input FILES	input fastq files (must use quotes with wildcard)
* \-l | \-\-logfile FILE	output file for logging
* \-o | \-\-output DIR	directory for subsampled outputs
* \-t | \-\-threads INT	number of threads to use [1]
* \-tr| \-\-trimmed DIR	directory to output trimmed reads

Defaults for options:

* \-\-depths \[none\]
* \-\-input \[none\]
* \-\-logfile \[subsample\_logfile.txt\]
* \-\-output \[none\]
* \-\-threads \[1\]
* \-\-trimmed \[none\]
