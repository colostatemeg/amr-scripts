#!/bin/bash

# Batch script for the MEG subsampling of fastq files
# Outputs fastqs trimmed and subsampled to the specified depths

# Dependencies: dist.py script and python with the numpy package installed

########
# Help #
########

#Display help menu

displayHelp () {
    echo "
    Usage: subsampleFastq.sh -i \"/s/bovine/.../fastq/*\" -o DIR -d \"INTS\" [options]

    -h | --help             help info

    Input options:

	-i | --input FILES	    input fastq files (e.g. \"/s/bovine/.../fastq/*\")
	-o | --output DIR       directory to output the subsampled files
	-d | --depths INT_ARRAY any number of depths to subsample, (e.g. -d \"0.1 0.5 0.8\")
	-l | --logfile FILE     output file for logging [subsample_logfile.txt]
	-t | --threads INT      number of threads to use [1]
	-tr| --trimmed DIR      directory to output trimmed reads [output_dir]
	
	Note: Input files AND input depths MUST be passed in quotes (e.g. -d \"0.2 0.5 0.8\")
	
    "
}

#############
# Variables #
#############

inputs=""
output_dir=""
depths=""
logfile="subsample_logfile.txt"
threads=1

distPATH="/s/bovine/index/projs/canada/dist.py"
trimmPATH="/s/bovine/e/nobackup/common/tools/Trimmomatic-0-1.32/"

###########
# Methods #
###########

## Check that the program and directory paths exist
validatePaths() {
    echo -e "\n@.@\nFastQ subsampling started on `date`\n" >> LabNotebook.txt

    local missing=""

    echo "
    Paths currently set:
    dist.py:                ${distPATH}
    Trimmomatic:            ${trimmPATH}
    Subsampled output:      ${output_dir}
    Trimmed output:         ${output_dir_trimmed}
    Logfile:                ${logfile}

    Validating....
    
    " | tee -i -a LabNotebook.txt
    
    if [ ! -e "${distPATH}" ]; then
        local missing="$missing:dist.py;"
    fi
    
    if [ ! -d "${trimmPATH}" ]; then
        local missing="$missing:Trimmomatic;"
    fi
    
    if [ ! -d "${output_dir}" ]; then
        local missing="$missing:OutputDirectory;"
    fi
    
    if [ ! -d "${output_dir_trimmed}" ]; then
        local missing="$missing:TrimmedOutputDirectory;"
    fi
    
    ## Check if components are missing
    if [ ! -e $missing ]; then
        echo -e "\nDirectories or Paths are missing, check the following:\n"
        echo -e $missing | sed 's/;/\n/g'
        exit 1
    else
        echo -e "\nPaths and Directories are valid. Proceeding...\n"
    fi
}
    
validateFiles() {
    local shortpath="$3"
    local read1="$1"
    local read2="$2"
    
    if [ ! -e ${read1} ] || [ ! -e ${read2} ]; then
        echo -e "\n${read1} or ${read2} does not exist."
        echo -e "\n${read1} or ${read2} is an invalid file" >> LabNotebook.txt
        exit 1
    else
        echo -e "\nBegin pipeline for ${shortpath}" | sed 's/_$//g' >> LabNotebook.txt
    fi
}
    
## Output the versions of the programs used at the time of running the script
## and other relevant information to the LabNotebook.txt file
getVersions() {
        echo -e "\nTrimmomatic: " >> LabNotebook.txt
        echo -e ${trimmPATH} | grep -o -e "-[0-9]-[0-9A-Za-z\._]*" >> LabNotebook.txt
}

## Use Trimmomatic to trim the raw input reads
trimmomatic() {
    local read1="$1"
    local read2="$2"
    local shortpath="$3"
    
    if [ -e "${read1}" ] && [ -e "${read2}" ]; then
        java -jar "${trimmPATH}trimmomatic-0.32.jar" PE -threads "${threads}" -phred33 "${read1}" "${read2}" "${output_dir_trimmed}${shortpath}"1P.fastq "${output_dir_trimmed}${shortpath}"1U.fastq "${output_dir_trimmed}${shortpath}"2P.fastq "${output_dir_trimmed}${shortpath}"2U.fastq ILLUMINACLIP:"${trimmPATH}adapters/TruSeq3-PE.fa:2:30:10:3:TRUE" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        echo -e "\nerror: please enter a valid fastq file path\n"
        echo -e "\tTrimmomatic: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e ${output_dir_trimmed}${shortpath}1P.fastq ] || [ ! -e ${output_dir_trimmed}${shortpath}1U.fastq ] || [ ! -e ${output_dir_trimmed}${shortpath}2P.fastq ] || [ ! -e ${output_dir_trimmed}${shortpath}2U.fastq ]; then
        echo -e "\nTrimmomatic did not complete properly for sample ${shortpath}\n"
        echo -e "\tTrimmomatic: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tTrimmomatic: Completed" >> LabNotebook.txt
    fi
}

########
# Main #
########

while [[ "${1+defined}"  ]]; do
    case "$1" in
	-i | --input)
	    inputs=($2)
	    shift 2
	    ;;
        -h | --help)
            displayHelp #help function
            exit 0
            ;;
        -d | --depths)
            depths=($2)
            shift 2
            ;;
        -l | --logfile)
            logfile="$2"
            shift 2
            ;;
        -o | --output)
            output_dir="$2" #insert check for existing or create
            shift 2
            ;;
        -t | --threads)
            threads="$2"
            shift 2
            ;;
	-tr | --trimmed)
	    output_dir_trimmed="$2"
	    shift 2
	    ;;
        --) #End of options
            shift 1
            break
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *) #No more options
            break
            ;;
    esac
done


# Copy STDOUT and STDERR to logfile
exec > >(tee -i -a "$logfile")
exec 2> >(tee -i -a "$logfile" >&2)

## Validate that input dirs exist and programs are in paths
validatePaths

## Output the versions of the programs to LabNotebook.txt
getVersions

## Run Trimmomatic on the raw input fastqs
for i in ${!inputs[*]}; do
    if [[ ${inputs[$i]} == *R2*.fastq ]]; then
	    continue
    else
    	echo -e "\n^.^\n"
        prefix=$( echo ${inputs[$i]} | sed -r 's/.fastq[A-Za-z0-9._\*]*//g' ) #Remove trailing .fastq
        shortpath=$( echo ${inputs[$i]} | sed -r 's/[A-Za-z0-9_.:\*-]*\///g' | sed -r 's/.fastq[A-Za-z0-9._\*]*//g' | sed '/s/R1//g' ) #Remove trailing .fastq and leading filepath, leaving only the file name stem
        read1=${inputs[$i]}
        read2=$( echo ${inputs[$i]} | sed -r 's/R1/R2/g' )
        validateFiles $read1 $read2 $shortpath #Make sure the input files are fastq and exist
        if [ ! -e ${output_dir_trimmed}${shortpath}1P.fastq ]; then
	        echo -e "Trimming ${read1} and ${read1}...\n"
	        trimmomatic $read1 $read2 $shortpath
                if [ $? -ne 0 ]; then #Catch errors or SIGINTs
                    echo -e "Ctrl^C"
                    echo -e "\tTrimmomatic: Error or Command Interrupt Ctrl^C" >> LabNotebook.txt
                    exit 1
                fi
	        echo -e "\nTrimming complete!\n"
        else
	        echo -e "${output_dir_trimmed}${shortpath}1P.fastq already exists.  Proceeding with next file...\n"
	        echo -e "\tTrimmomatic: Completed Previously, skipping..." >> LabNotebook.txt
        fi
    fi
    echo -e "\nBeginning Subsampling\n"
    python ${distPATH} -i ${output_dir_trimmed}${shortpath}1P.fastq -o ${output_dir} ${depths[*]}
    echo -e "\nFinished Subsampling\n"
done
