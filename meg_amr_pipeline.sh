#!/usr/bin/env bash

# Version: 1.01
# Batch script for the MEG standard AMR pipeline
# Optional integration of microbiome (Kraken, MetaPhlAn) and pathway (HUMAnN2) data

# Workflow:
# 1. Trim reads with Trimmomatic
# 2. Align reads to bovine genome
# 3. Filter those reads out with SAMTools
# 4. Align the remaining reads to the AMR database
# 5. Use SamRatio to calculate number of AMR gene hits
# 6. Run Kraken on the non-bovine fastqs (Optional)
# 7. Run MetaPhlAn on the non-bovine fastqs (Optional)
# 8. Run HUMAnN2 on t he non-bovine fastqs (Optional), not yet implemented

# Dependencies:

# Contact: Steven Lakin (Steven.Lakin@colostate.edu)
# Last Updated: Sept 8 2015

########
# Help #
########

#Display help menu
displayHelp () {
    echo "
    Usage: amrPipeline.sh -i \"raw_sequence_reads/*\" [options]

        -h | --help             help info

    Input options:

        -i | --input FILES	input fastq.gz files

    Pipeline parameters:

        -t | --threads INT      number of threads to use [1]
        -a | --kraken DIR	perform kraken analysis, output to a directory
        -m | --metaphlan	perform metaphlan analysis, output to a directory

    Output options:

        -k | --keep             save intermediate files (ex. sam/bam alignment files)
        -l | --logfile		file name for logfile output [./logfile.txt]
        -n | --nonhost DIR	directory for nonhost fastqs
        -o | --output DIR       directory for AMR file output
        -tr | --trimmed DIR      directory for output of trimmed fastq files

    Note: You must pass the fastq file paths in double quotes if they contain a wildcard
"
}

#############
# Variables #
#############
shopt -s extglob

## The below are default paths to the tools at the time this script was
## written.  They will need to be changed manually if the paths change.
## The script will check for correct paths before running.
## Be sure to use terminal slashes for directories.
trimmPATH="/s/angus/index/common/tools/Trimmomatic-0-1.32/"
krakenPATH="/usr/local/bin/kraken/kraken"
krakenDatabasePATH="/s/angus/index/databases/kraken_databases/Standard_kraken_10.14.db/"
bamToFastqPATH="/s/angus/index/common/tools/bedtools2/bin/bamToFastq"
hostGenomePATH="/s/angus/index/databases/bwa_indexes/Homo_sapiens_GRCh38.p3/Homo_sapiens_GRCh38.fasta"
AMRdatabasePATH="/s/angus/index/databases/bait_creation_files/databases/2015Oct22_combinedDB.fa"
AMRindexedPATH="/s/angus/index/databases/bait_creation_files/databases/2015Oct22_combinedDB.fa"
samRatioPATH="/s/angus/index/common/tools/samratio.jar"
bowtie2PATH="/usr/bin/bowtie2"
metaphlanPATH="/s/angus/i/nobackup/tools/bin/metaphlan2/"
humann2PATH="/s/angus/i/nobackup/tools/bin/humann2/humann2/humann2.py"
diamondPATH="/s/angus/i/nobackup/tools/bin/diamond"

## These flags determine if the kraken, metaphlan, and humann2 pipelines
## are run.  The keep flag determines if intermediate files are kept.
krakenflag="False"
metaphlanflag="False"
humann2flag="False"
keep="False"

## The below are default values; they will be changed by the inputs if
## the parameter is set by the user
output_dir_trimmed="${PWD}/trimmed_fastq_files/"
output_dir_amr="${PWD}/analysis/AMR/"
output_dir_nonhost="${PWD}/non_bovine_fastq/"
output_kraken="${PWD}/analysis/microbiome/kraken/"
kraken_report="${PWD}/analysis/microbiome/kraken/kraken_reports/"
output_metaphlan="${PWD}/analysis/microbiome/metaphlan/"
output_humann2="${PWD}/analysis/pathways/"
logfile="${PWD}/logfile.txt"
threads=1

###########
# Methods #
###########

## Check that the program and directory paths exist
validatePaths() {
    echo -e "\n@.@\nAMR Pipeline started on `date`\n" >> LabNotebook.txt

    local missing=""

    echo "
    Paths currently set:
    Trimmomatic:            ${trimmPATH}
    Kraken:                 ${krakenPATH}
    Kraken Database:        ${krakenDatabasePATH}
    bamToFastq:             ${bamToFastqPATH}
    Host Genome:            ${hostGenomePATH}
    AMR Database:           ${AMRdatabasePATH}
    Indexed AMR Database:   ${AMRindexedPATH}
    Sam Ratio:              ${samRatioPATH}
    MetaPhlAn:              ${metaphlanPATH}
    bowtie2:                ${bowtie2PATH}
    
    Directories selected:
    Trimmed output:         ${output_dir_trimmed}
    AMR output:             ${output_dir_amr}
    Non-host output:        ${output_dir_nonhost}
    Kraken output:          ${output_kraken}
    MetaPhlAn output:       ${output_metaphlan}
    Logfile:                ${logfile}
    
    Validating....
    
    " | tee -i -a LabNotebook.txt
    
    if [ ! -d "${trimmPATH}" ]; then
        local missing="$missing:Trimmomatic;"
    fi
    
    if [ ! -e "${krakenPATH}" ]; then
        local missing="$missing:Kraken;"
    fi
    
    if [ ! -e "${krakenDatabasePATH}" ]; then
        local missing="$missing:KrakenDatabase;"
    fi
    
    if [ ! -e "${bamToFastqPATH}" ]; then
        local missing="$missing:bamToFastq;"
    fi
    
    if [ ! -e "${hostGenomePATH}" ]; then
        local missing="$missing:BovineGenome;"
    fi
    
    if [ ! -e "${AMRdatabasePATH}" ]; then
        local missing="$missing:AMRdatabase;"
    fi
    
    if [ ! -e "${AMRindexedPATH}.sa" ]; then
        local missing="$missing:AMRindexedDatabase;"
    fi
    
    if [ ! -e "${samRatioPATH}" ]; then
        local missing="$missing:SamRatio;"
    fi
    
    if [ ! -e "${bowtie2PATH}" ]; then
        local missing="$missing:bowtie2;"
    fi
    
    if [ ! -d "${metaphlanPATH}" ]; then
        local missing="$missing:MetaPhlAn;"
    fi
    
    if [ ! -d "${output_dir_trimmed}" ]; then
        local missing="$missing:TrimmedOutputDir;"
    fi
    
    if [ ! -d "${output_dir_amr}" ]; then
        local missing="$missing:AMR_OutputDir;"
    fi
    
    if [ ! -d "${output_dir_nonhost}" ]; then
        local missing="$missing:NonHostOutputDir;"
    fi
    
    if [ ! -d "${output_kraken}" ] && [ $krakenflag == "True" ]; then
        local missing="$missing:KrakenOutputDir;"
    fi
    
    if [ ! -d "${output_metaphlan}" ] && [ $metaphlanflag == "True" ]; then
        local missing="$missing:MetaphlanOutputDir;"
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


## Check that the input files are valid
validateFiles() {
    if [ ! -e $1 ] || [ ! -e $2 ]; then
        echo -e "$1 or $2 does not exist."
        echo -e "$1 or $2 is an invalid file" >> LabNotebook.txt
        exit 1
    else
        echo -e "\nBegin pipeline for ${3}" | sed 's/_$//g' >> LabNotebook.txt
    fi
}


## Output the versions of the programs used at the time of running the script
## and other relevant information to the LabNotebook.txt file
getVersions() {
    echo -e "\nTrimmomatic: " >> LabNotebook.txt
	    echo -e ${trimmPATH} | grep -o -e "-[0-9]-[0-9A-Za-z\._]*" >> LabNotebook.txt
    echo -e "\nKraken: " >> LabNotebook.txt
        ${krakenPATH} -v | grep "version" >> LabNotebook.txt
    echo -e "\nKraken Database: " >> LabNotebook.txt
        echo -e ${krakenDatabasePATH} | sed -r 's/[A-Za-z0-9._]*\///g' | sed 's/.db/\n/g' >> LabNotebook.txt
    echo -e "\nbamToFastq: " >> LabNotebook.txt
        ${bamToFastqPATH} 2> temp.txt; cat temp.txt | grep "Version" >> LabNotebook.txt; rm temp.txt
    echo -e "\nGenome: " >> LabNotebook.txt
        head ${hostGenomePATH} | grep ">" >> LabNotebook.txt
    echo -e "\nSamRatio: " >> LabNotebook.txt
        java -jar ${samRatioPATH} -v >> LabNotebook.txt
    echo -e "\nMetaPhlAn2: " >> LabNotebook.txt
        ${metaphlanPATH}metaphlan2.py -v >> LabNotebook.txt
    echo -e "\nbowtie2: " >> LabNotebook.txt
        ${bowtie2PATH} --version >> LabNotebook.txt
    echo -e "\nDiamond: " >> LabNotebook.txt
        ${diamondPATH} -v >> LabNotebook.txt
   # echo -e "\nHUMAnN2: " >> LabNotebook.txt
    #    ${humann2PATH} --version >> LabNotebook.txt
}


trimmomatic() {
    if [ -e "$1" ] && [ -e "$2" ]; then
        java -jar "${trimmPATH}trimmomatic-0.32.jar" PE -threads "${threads}" -phred33 "$1" "$2" "${output_dir_trimmed}${3}"_1P.fastq "${output_dir_trimmed}${3}"_1U.fastq "${output_dir_trimmed}${3}"_2P.fastq "${output_dir_trimmed}${3}"_2U.fastq ILLUMINACLIP:"${trimmPATH}adapters/TruSeq3-PE.fa:2:30:10:3:TRUE" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        echo -e "\nerror: please enter a valid FASTQ file path\n"
        echo -e "\tTrimmomatic: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e ${output_dir_trimmed}${3}_1P.fastq ] || [ ! -e ${output_dir_trimmed}${3}_1U.fastq ] || [ ! -e ${output_dir_trimmed}${3}_2P.fastq ] || [ ! -e ${output_dir_trimmed}${3}_2U.fastq ]; then
        echo -e "\nTrimmomatic did not complete properly for sample ${3}\n"
        echo -e "\tTrimmomatic: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tTrimmomatic: Completed" >> LabNotebook.txt
    fi
}


## Use BWA aln to filter host dna from trimmed fastqs
hostfilter() {
    if [ ! -d intermediate_files ]; then
        mkdir intermediate_files
    fi
    
    if [ ! -e intermediate_files/${1}_f.sai ] || [ ! -e intermediate_files/${1}_r.sai ]; then
        bwa aln ${hostGenomePATH} ${output_dir_trimmed}${1}_1P.fastq -t ${threads} > intermediate_files/${1}_f.sai
        bwa aln ${hostGenomePATH} ${output_dir_trimmed}${1}_2P.fastq -t ${threads} > intermediate_files/${1}_r.sai
    elif [ -e intermediate_files/${1}_f.sai ] && [ -e intermediate_files/${1}_r.sai ]; then
        echo -e "\n.sai files already present for sample ${1}.  Proceeding with bwa sampe..."
    else
        echo -e "\nNon-specific error in bwa aln for Host Filtering. Check syntax (shell script)\n"
        echo -e "\tHost Filtering, bwa sampe: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e intermediate_files/${1}_f.sai ] || [ ! -e intermediate_files/${1}_r.sai ]; then
        echo -e "\nbwa aln for host filtering did not complete properly for sample ${1}\n"
        echo -e "\tHost Filtering, bwa aln: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tHost Filtering, bwa aln: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e intermediate_files/${1}_cow.sam ]; then
        bwa sampe ${hostGenomePATH} intermediate_files/${1}_f.sai intermediate_files/${1}_r.sai ${output_dir_trimmed}${1}_1P.fastq ${output_dir_trimmed}${1}_2P.fastq > intermediate_files/${1}_cow.sam
    elif [ -e intermediate_files/${1}_cow.sam ]; then
        echo -e "\nSAM file already present for sample ${1}.  Proceeding with conversion to BAM...\n"
    else
        echo -e "\nNon-specific error in bwa sampe for Host Filtering.  Check syntax (shell script)\n"
        echo -e "\tHost Filtering, bwa sampe: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e intermediate_files/${1}_cow.sam ]; then
        echo -e "\nbwa sampe for host filtering did not complete properly for sample ${1}\n"
        echo -e "\tHost Filtering, bwa sampe: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tHost Filtering, bwa sampe: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e intermediate_files/${1}_cow.bam ]; then
        samtools view -hbS intermediate_files/${1}_cow.sam > intermediate_files/${1}_cow.bam
    elif [ -e intermediate_files/${1}_cow.bam ]; then
        echo -e "\nBAM file already present for sample ${1}.  Proceeding with BAM sorting...\n"
    else
        echo -e "\nNon-specific error in SAM to BAM conversion for Host Filtering.  Check syntax (shell script)\n"
        echo -e "\tHost Filtering, samtools view: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e intermediate_files/${1}_cow.bam ]; then
        echo -e "\nSAM to BAM conversion for host filtering  did not complete properly for sample ${1}\n"
        echo -e "\tHost Filtering, samtools view: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tHost Filtering, samtools view: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e intermediate_files/${1}_sorted_cow.bam ]; then
        samtools sort intermediate_files/${1}_cow.bam intermediate_files/${1}_sorted_cow
    elif [ -e intermediate_files/${1}_sorted_cow.bam ]; then
        echo -e "\nBAM file already present for sample ${1}. Proceeding with host DNA removal..."
    else
        echo -e "\nNon-specific error in BAM sorting for Host Filtering. Check syntax (shell script)\n"
        echo -e "\tHost Filtering, samtools sort: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e intermediate_files/${1}_sorted_cow.bam ]; then
        echo -e "\nsamtools sort for host filtering did not complete properly for sample ${1}\n"
        echo -e "\tHost Filtering, samtools sort: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tHost Filtering, samtools sort: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e intermediate_files/${1}_sorted_noncow.bam ]; then
        samtools view -h -f 4 -b intermediate_files/${1}_sorted_cow.bam > intermediate_files/${1}_sorted_noncow.bam
    elif [ -e intermediate_files/${1}_sorted_noncow.bam ]; then
        echo -e "\nFiltered BAM already present for sample ${1}. Proceeding with fastq conversion..."
    else
        echo -e "\nNonspecific error in BAM filtering for Host Filtering. Check syntax (shell script)\n"
        echo -e "\tHost Filtering, samtools host removal: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e intermediate_files/${1}_sorted_noncow.bam ]; then
        echo -e "\nsamtools host DNA removal did not complete properly for sample ${1}\n"
        echo -e "\tHost Filtering, samtools host removal: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tHost Filtering, samtools host removal: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e ${output_dir_nonhost}${1}_noncow_2p.fastq ] || [ ! -e ${output_dir_nonhost}${1}_noncow_1p.fastq ]; then
        ${bamToFastqPATH} -i intermediate_files/${1}_sorted_noncow.bam -fq ${output_dir_nonhost}${1}_noncow_1p.fastq -fq2 ${output_dir_nonhost}${1}_noncow_2p.fastq
    elif [ -e ${output_dir_nonhost}${1}_noncow_2p.fastq ] && [ -e ${output_dir_nonhost}${1}_noncow_1p.fastq ]; then
        echo -e "\nConverted fastqs already present for sample ${1}. Proceeding with AMR pipeline..."
    else
        echo -e "\nNon-specific error in BAM to fastq conversion. Check syntax (shell script)\n"
        echo -e "\tHost Filtering, bamToFastq: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    fi
    
    if [ ! -e ${output_dir_nonhost}${1}_noncow_2p.fastq ] || [ ! -e ${output_dir_nonhost}${1}_noncow_1p.fastq ]; then
        echo -e "\nBAM to fastq conversion during host removal did not complete properly for sample ${1}\n"
        echo -e "\tHost Filtering, bamToFastq: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tHost Filtering, bamToFastq: Completed" >> LabNotebook.txt
    fi
    
    mv intermediate_files/${1}_sorted_noncow.bam bam_files/
    
    if [ ! $keep == "keep" ]; then
        rm intermediate_files/${1}_f.sai
        rm intermediate_files/${1}_r.sai
        rm intermediate_files/${1}_cow.sam
        rm intermediate_files/${1}_cow.bam
        rm intermediate_files/${1}_sorted_cow.bam
    fi
}


## Use BWA aln to align the nonhost reads against the AMR master database
amrAlign() {
    if [ ! -d intermediate_files ]; then
        mkdir intermediate_files
    fi
    
    if [ ! -e intermediate_files/${1}_f_amr.sai ] && [ ! -e intermediate_files/${1}_r_amr.sai ]; then
        bwa aln ${AMRindexedPATH} ${output_dir_nonhost}${1}_noncow_1p.fastq -t ${threads} > intermediate_files/${1}_f_amr.sai
        bwa aln ${AMRindexedPATH} ${output_dir_nonhost}${1}_noncow_2p.fastq -t ${threads} > intermediate_files/${1}_r_amr.sai
    else
        echo -e "\n.sai files already present for sample ${1}, proceeding with bwa sampe..."
    fi
    
    
    if [ ! -e intermediate_files/${1}_f_amr.sai ] || [ ! -e intermediate_files/${1}_r_amr.sai ]; then
        echo -e "\nbwa aln did not complete properly for sample ${1}\n"
        echo -e "\tAMR, bwa aln: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tAMR, bwa aln: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e intermediate_files/${1}_AMR.sam ]; then
        bwa sampe -n 1000 -N 1000 ${AMRindexedPATH} intermediate_files/${1}_f_amr.sai intermediate_files/${1}_r_amr.sai ${output_dir_nonhost}${1}_noncow_1p.fastq ${output_dir_nonhost}${1}_noncow_2p.fastq > intermediate_files/${1}_AMR.sam
    else
        echo -e "\nSAM file already present for sample ${1}, proceeding with BAM generation..."
    fi
    
    
    if [ ! -e intermediate_files/${1}_AMR.sam ]; then
        echo -e "\nbwa sampe did not complete properly for sample ${1}\n"
        echo -e "\tAMR, bwa sampe: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tAMR, bwa sampe: Completed" >> LabNotebook.txt
    fi
    
    if [ ! -e intermediate_files/${1}_AMR.bam ] && [ ! -e bam_files/${1}_AMR.bam ]; then
        samtools view -hbS intermediate_files/${1}_AMR.sam > intermediate_files/${1}_AMR.bam
    else
        echo -e "\n BAM file already present for sample ${1}, proceeding wtih samratio..."
    fi
    
    
    java -jar ${samRatioPATH} -d ${AMRdatabasePATH} -i intermediate_files/${1}_AMR.sam -t 1 -m ${output_dir_amr}${1}_mismatch -o ${output_dir_amr}${1}_parsed
    
    mv intermediate_files/${1}_AMR.bam bam_files/
    
    if [ ! -e ${output_dir_amr}${1}_mismatch ] || [ ! -e ${output_dir_amr}${1}_parsed ]; then
        echo -e "\nsamratio did not complete properly for sample ${1}\n"
        echo -e "\tAMR, SamRatio: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tAMR, SamRatio: Completed" >> LabNotebook.txt
    fi
    
    if [ ! $keep == "keep" ]; then
        rm intermediate_files/${1}_f_amr.sai
        rm intermediate_files/${1}_r_amr.sai
        rm intermediate_files/${1}_AMR.sam
    fi
}


## Use Kraken to profile the mirobiome from the metagenomic reads
krakenProfile() {
    if [ ! -e ${output_kraken}${1}_kraken_output ]; then
        ${krakenPATH} --preload --db ${krakenDatabasePATH} --threads ${threads} --fastq-input --paired ${output_dir_nonhost}${1}_noncow_1p.fastq ${output_dir_nonhost}${1}_noncow_2p.fastq > ${output_kraken}${1}_kraken_output
    else
        echo -e "\nkraken output file already present for sample ${1}, proceeding with report generation...\n"
    fi
    
    
    if [ ! -e ${output_kraken}${1}_kraken_output ]; then
        echo -e "\nkraken did not complete properly for sample ${1}\n"
        echo -e "\tKraken, classify: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tKraken, classify: Completed" >> LabNotebook.txt
    fi
    
    ${krakenPATH}-report -db ${krakenDatabasePATH} ${output_kraken}${1}_kraken_output > ${kraken_report}${1}_kraken_reports
    
    if [ ! -e ${kraken_report}${1}_kraken_reports ]; then
        echo -e "\nkraken reports did not complete properly for sample ${1}\n"
        echo -e "\tKraken, report: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\tKraken, report: Completed" >> LabNotebook.txt
    fi
}


metaphlanProfile() {
    ${metaphlanPATH}metaphlan2.py ${output_dir_nonhost}${1}_noncow_1p.fastq,${output_dir_nonhost}${1}_noncow_2p.fastq ${output_metaphlan}${1}_metaphlan_output.tsv --input_type fastq --bowtie2db ${metaphlanPATH}db_v20/mpa_v20_m200 --mpa_pkl ${metaphlanPATH}db_v20/mpa_v20_m200.pkl --bowtie2out ${output_metaphlan}${1}.bowtie2.bz2 --bowtie2_exe ${bowtie2PATH} -t rel_ab_w_read_stats --nproc ${threads}
    
    if [ ! -e ${output_metaphlan}${1}_metaphlan_output.tsv ]; then
        echo -e "\nmetaphlan did not complete properly for sample ${1}\n"
        echo -e "\tMetaPhlAn: Error. Pipeline stopped" >> LabNotebook.txt
        exit 1
    else
        echo -e "\nMetaPhlAn: Completed" >> LabNotebook.txt
    fi
    
}


combineOutputs() {
    echo -e "\nMerging Files...\n"
    for i in ${infiles[@]}; do
		if [[ $i =~ (_R1_) ]]; then
			sample_name=$( basename $i | sed -r 's/(.fastq|.fastq.gz)//' )
			echo -e "Merging files for sample ${sample_name}" >> LabNotebook.txt
			sed 's/$/\t'"$sample_name"'/g' ${output_kraken}${sample_name}_kraken_output >> ${output_kraken}master_kraken.tsv
			sed 's/$/\t'"$sample_name"'/g'  ${kraken_report}${sample_name}_kraken_reports >> ${kraken_report}master_kraken_reports.tsv
			sed 's/$/\t'"$sample_name"'/g' ${output_dir_amr}${sample_name}_parsed >> ${output_dir_amr}master_amr_parsed.tsv
			sed 's/$/\t'"$sample_name"'/g' ${output_dir_amr}${sample_name}_mismatch >> ${output_dir_amr}master_amr_mismatch.tsv
			sed 's/$/\t'"$sample_name"'/g' ${output_metaphlan}${sample_name}_metaphlan_output.tsv >> ${output_metaphlan}master_metaphlan.tsv

		fi
    done   
    echo -e "Merging Completed:"
    echo "
	${output_kraken}master_kraken.tsv
	${kraken_report}master_kraken_reports.tsv
	${output_dir_amr}master_amr_parsed.tsv
	${output_dir_amr}master_amr_mismatch.tsv
	
	If metaphlan was run:
	${output_metaphlan}master_metaphlan.tsv
    "
    echo -e "\nMerging complete for all files" >> LabNotebook.txt
}


setPermissions() {
    chmod 666 ${logfile} ${output_dir_trimmed}/* ${output_dir_nonhost}/* ${output_dir_amr}/* LabNotebook.txt
    
    if [ $krakenflag == "True" ]; then
        chmod 666 ${output_kraken}/*output ${output_kraken}kraken_reports/*
    fi
    
    if [ $metaphlanflag == "True" ] && [ $humann2flag == "False" ]; then
        chmod 666 ${output_metaphlan}/*
    fi
}


##########
## Main ##
##########
[[ $@ ]] || { display_help; exit 1; }

exec > >(tee -i -a "$logfile")
exec 2> >(tee -i -a "$logfile" >&2)

while [[ "${1+defined}"  ]]; do
    case "$1" in
	    -i | --input)
	        infiles=($2)
	        shift 2
	        ;;
        -a | --kraken)
            krakenflag="True"
            output_kraken="$2"
            kraken_report="$2kraken_reports/"
            if [ ! -d "${kraken_report}" ]; then
                mkdir "${kraken_report}"
            fi
            shift 2
            ;;
        -h | --help)
            displayHelp #help function
            exit 0
            ;;
        -k | --keep)
            keep="keep"
            shift
            ;;
        -l | --logfile)
            logfile="$2"
            shift 2
            ;;
        -m | --metaphlan)
            metaphlanflag="True"
            output_metaphlan="$2"
            shift 2
            ;;
		-n | --nonhost)
			output_dir_nonhost="$2"
			shift 2
			;;
        -o | --output)
            output_dir_amr="$2"
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


for i in ${infiles[@]}; do
	if [[ $i =~ (_R1_) ]]; then
		echo -e "\n=^.^=\n" #Marker for start of sample run. Meow.
		sample_name=$( basename $i | sed -r 's/(.fastq|.fastq.gz)//' )
		prefix=$( echo "$i" | sed -r 's/(.fastq|.fastq.gz)//' )
		revfile=$( echo "$i" | sed 's/_R1_/_R2_/' )
		validateFiles $i $revfile $sample_name  #Make sure the input files are fastq.gz and exist
		if [ ! -e ${output_dir_trimmed}${sample_name}_1P.fastq ]; then
	        echo -e "Trimming ${i} and ${revfile}...\n"
	        trimmomatic $i $revfile $sample_name
            if [ $? -ne 0 ]; then #Catch errors or SIGINTs
                echo -e "Ctrl^C"
                echo -e "\tTrimmomatic: Error or Command Interrupt Ctrl^C" >> LabNotebook.txt
                exit 1
            fi
	    echo -e "\nTrimming complete!\n"
        else
	        echo -e "${output_dir_trimmed}${sample_name}_1P.fastq already exists.  Proceeding with next file...\n"
	        echo -e "\tTrimmomatic: Completed Previously, skipping..." >> LabNotebook.txt
       	fi
    fi
    
    ## Remove Bovine DNA, move to BAM, output fastq and BAM
    if [ ! -e ${output_dir_trimmed}${sample_name}_1P.fastq ] || [ ! -e ${output_dir_trimmed}${sample_name}_2P.fastq ]; then
        echo -e "\nError in sample ${sample_name}, trimmed fastq files were not generated or do not exist"
        echo -e "\tHost Filtering: Error. Pipeline Stopped" >> LabNotebook.txt
        exit 1
    elif [ -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] && [ -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        echo -e "\nNon-host fastqs already exist for ${sample_name} samples.  Proceeding with next file...\n"
        echo -e "\tHost Filtering: Completed Previously, skipping..." >> LabNotebook.txt
    elif [ ! -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] && [ ! -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ] && [ -e ${output_dir_trimmed}${sample_name}_1P.fastq ] && [ -e ${output_dir_trimmed}${sample_name}_2P.fastq ]; then
        echo -e "\nFiltering non-host DNA from ${sample_name}_1P.fastq and ${sample_name}_2P.fastq...\n"
        hostfilter $sample_name
        if [ $? -ne 0 ]; then #Catch errors or SIGINTs
            echo -e "Ctrl^C"
            echo -e "\tHost Filtering: Error or Command Interrupt Ctrl^C" >> LabNotebook.txt
            exit 1
        fi
        echo -e "\n Host DNA filtering complete!\n"
    else
        echo -e "\nNon-specific error in host-DNA removal (shell script).  Check the syntax.\n"
        echo -e "\tHost Filtering: Error. Pipeline Stopped" >> LabNotebook.txt
        exit 1
    fi
    
    ## Align to AMR master database, output mismatch and parsed files
    if [ ! -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] || [ ! -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        echo -e "\nError in sample ${sample_name}, one or more noncow fastq files were not generated or do not exist"
        echo -e "\tAMR: Error, Pipeline Stopped" >> LabNotebook.txt
        exit 1
    elif [ -e ${output_dir_amr}${sample_name}_mismatch ] && [ -e ${output_dir_amr}${sample_name}_parsed ]; then
        echo -e "\nAMR mismatch and parsed files already exist for ${sample_name} samples.  Proceeding with next file...\n"
        echo -e "\tAMR: Completed Previously, skipping..." >> LabNotebook.txt
    elif [ -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] && [ -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        amrAlign $sample_name
        if [ $? -ne 0 ]; then #Catch errors or SIGINTs
            echo -e "Ctrl^C"
            echo -e "\tAMR: Error or Command Interrupt Ctrl^C" >> LabNotebook.txt
            exit 1
        fi
    else
        echo -e "\nNon-specific error in AMR alignment pipeline (shell script).  Check the syntax.\n"
        echo -e "\tAMR: Error. Pipeline Stopped" >> LabNotebook.txt
        exit 1
    fi
    
    ## If the kraken flag is set, then run kraken on the nonhost fastqs, output kraken and kraken report files
    if [ ! -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] || [ ! -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        echo -e "\nError in sample ${sample_name}, one or more noncow fastq files were not generated or do not exist"
        echo -e "\tKraken: Error. Pipeline Stopped" >> LabNotebook.txt
        exit 1
    elif [ ${krakenflag} == "True" ] && [ -e ${output_kraken}${sample_name}_kraken_output ] && [ -e ${kraken_report}${sample_name}_kraken_reports ]; then
        echo -e "\nKraken output and report files already exist for ${sample_name} samples.  Proceeding with next file...\n"
        echo -e "\tKraken: Completed Previously, skipping..." >> LabNotebook.txt
    elif [ ${krakenflag} == "True" ] && [ -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] && [ -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        krakenProfile $sample_name
        if [ $? -ne 0 ]; then #Catch errors or SIGINTs
            echo -e "Ctrl^C"
            echo -e "\tKraken: Error or Command Interrupt Ctrl^C" >> LabNotebook.txt
            exit 1
        fi
    elif [ ${krakenflag} == "False" ]; then
    	echo -e "\nKraken flag is False, skipping kraken pipeline...\n"
    	echo -e "\tKraken: Flag is False, skipping..." >> LabNotebook.txt
    else
        echo -e "\nNon-specific error in Kraken pipeline (shell script).  Check the syntax.\n"
        echo -e "\tKraken: Error. Pipeline Stopped" >> LabNotebook.txt
        exit 1
    fi
    
    ## If the metaphlan flag is set, then run metaphlan on the nonhost fastqs, output to metaphlan folder in analysis/microbiome
    if [ ${metaphlanflag} == "True" ] && [ -e ${output_metaphlan}${sample_name}_metaphlan_output.tsv ]; then
        echo -e "\nMetaPhlAn output file already exists for ${sample_name} samples.  Proceeding with next file...\n"
        echo -e "\tMetaPhlAn: Completed Previously, skipping..." >> LabNotebook.txt
    elif [ ${metaphlanflag} == "True" ] && [ ! -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] || [ ! -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        echo -e "\nError in sample ${sample_name}, one or more noncow fastq files were not generated or do not exist\n"
        echo -e "\tMetaPhlAn: Error. Pipeline Stopped" >> LabNotebook.txt
    elif [ ${metaphlanflag} == "True" ] && [ -e ${output_dir_nonhost}${sample_name}_noncow_1p.fastq ] && [ -e ${output_dir_nonhost}${sample_name}_noncow_2p.fastq ]; then
        metaphlanProfile $sample_name
        if [ $? -ne 0 ]; then #Catch errors or SIGINTs
            echo -e "Ctrl^C"
            echo -e "\tMetaPhlAn: Error or Command Interrupt Ctrl^C" >> LabNotebook.txt
            exit 1
        fi
    elif [ ${metaphlanflag} == "False" ]; then
    	echo -e "\nMetaPhlAn: flag is False, skipping kraken pipeline...\n"
    	echo -e "\tMetaPhlAn: Flag is False, skipping..." >> LabNotebook.txt
    else
        echo -e "\nNon-specific error in MetaPhlAn pipeline (shell script).  Check the syntax.\n"
        echo -e "\tMetaPhlAn: Error. Pipeline Stopped" >> LabNotebook.txt
        exit 1
    fi
done

echo -e "\nPipeline complete!\n"
echo -e "\nPipeline complete!\n" >> LabNotebook.txt

## Combine output reports
if [ ! -e ${output_dir_amr}master_amr_mismatch.tsv ] || [ ! -e ${output_dir_amr}master_amr_parsed.tsv ] || [ ! -e ${output_kraken}master_kraken.tsv ] || [ ! -e ${kraken_report}master_kraken_reports.tsv ] || [ ! -e ${output_metaphlan}master_metaphlan.tsv ]; then
    combineOutputs
fi

## Set permissions for output files to read/write access
setPermissions

exit 0

