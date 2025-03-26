#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script is meant to be run on the contigs and the sequencing reads to generate coverage information
# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.

# Some of the programs this pipeline uses are from the binning.sh file in MetaWRAP.
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: bash gen_coverage_file.sh [options] -a assembly.fa -o output_dir readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]"
	echo "  or : bash gen_coverage_file.sh [options] -a assembly.fa -o output_dir -f forward_reads.fastq -r reverse_reads.fastq"
	echo "  or : bash gen_coverage_file.sh [options] -a assembly.fa -o output_dir basename        # Will look for basename_{1,2}.fastq"
	echo "  or : bash gen_coverage_file.sh [options] -a assembly.fa -o output_dir basename_*.fastq"
	echo "  or : bash gen_coverage_file.sh [options] --text-files input.txt -o output_dir"
	echo ""
	echo "Note1: Make sure to provide all your separately replicate read files, not the joined file."
	echo "Note2: You may provide single end or interleaved reads as well with the use of the correct option"
	echo "Note3: If the output already has the .bam alignments files from previous runs, the module will skip re-aligning the reads"
	echo "Note4: Compressed files (.fastq.gz) are also supported"
	echo ""
	echo "Options:"
	echo ""
	echo "	-a STR    metagenomic assembly file"
	echo "	-o STR    output directory (to save the coverage files)"
	echo "	-b STR    directory for the bam files"
	echo "	-t INT    number of threads (default=1)"
	echo "	-m INT    amount of RAM available (default=4)"
	echo "	-l INT    minimum contig length to bin (default=1000bp)."
	echo "	--single-end	non-paired reads mode (provide *.fastq files)"
	echo "	--interleaved	the input read files contain interleaved paired-end reads"
	echo "	--text-files STR    input CSV file with format: contig_file,read1_file,read2_file"
	echo "	-f STR    Forward read file or suffix for paired reads (default="_1.fastq")"
	echo "	-r STR    Reverse read file or suffix for paired reads (default="_2.fastq")"
	echo "";}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# setting scripts and databases from config file (should be in same folder as main script)
#config_file=$(which config-metawrap)
#source $config_file

SOFT="$( cd "$( dirname "$0"  )" && pwd  )"


#/home/wangzy/tools/test_tool/metabinner/scripts
chmod +x ${SOFT}/print_comment.py
comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }



# default params
threads=1; mem=4; len=1000; out=false; ASSEMBLY=false; bout=false
# long options defaults
read_type=paired

F_reads_suffix=_1.fastq
R_reads_suffix=_2.fastq
FORWARD_READS=false
REVERSE_READS=false
DIRECT_INPUT_MODE=false
TEXT_FILES=false

# load in params

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
		            -m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
                -b) bout=$2; shift 2;;
                -a) ASSEMBLY=$2; shift 2;;
		            -l) len=$2; shift 2;;
		            -f) 
                    if [[ "$2" == *fastq* ]] || [[ "$2" == *fq* ]]; then
                        FORWARD_READS=$2
                        DIRECT_INPUT_MODE=true
                    else
                        F_reads_suffix=$2
                    fi
                    shift 2;;
		            -r) 
                    if [[ "$2" == *fastq* ]] || [[ "$2" == *fq* ]]; then
                        REVERSE_READS=$2
                        DIRECT_INPUT_MODE=true
                    else
                        R_reads_suffix=$2
                    fi
                    shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
		            --single-end) read_type=single; shift 1;;
		            --interleaved) read_type=interleaved; shift 1;;
                --text-files) TEXT_FILES=$2; shift 2;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done

echo ${threads}
echo ${len}
#exit 1

# Function to process a line from the text file
process_text_file_line() {
    local line=$1
    local assembly_file=$(echo "$line" | cut -d',' -f1 | tr -d ' ')
    local forward_reads=$(echo "$line" | cut -d',' -f2 | tr -d ' ')
    local reverse_reads=$(echo "$line" | cut -d',' -f3 | tr -d ' ')
    
    if [[ ! -s "$assembly_file" ]]; then
        error "Assembly file not found: $assembly_file"
    fi
    
    if [[ ! -s "$forward_reads" ]]; then
        error "Forward read file not found: $forward_reads"
    fi
    
    if [[ ! -s "$reverse_reads" ]]; then
        error "Reverse read file not found: $reverse_reads"
    fi
    
    # Set the variables for this iteration
    ASSEMBLY="$assembly_file"
    FORWARD_READS="$forward_reads"
    REVERSE_READS="$reverse_reads"
    DIRECT_INPUT_MODE=true
    
    comm "Processing assembly: $ASSEMBLY with reads: $FORWARD_READS and $REVERSE_READS"
}

# Process the text file if provided
if [[ $TEXT_FILES != false ]]; then
    # Check for conflicts with other parameters
    if [[ $ASSEMBLY != false ]]; then
        error "--text-files option cannot be used together with -a option"
    fi
    
    if [[ $DIRECT_INPUT_MODE = true ]]; then
        error "--text-files option cannot be used together with -f or -r options"
    fi
    
    if [[ ! -s "$TEXT_FILES" ]]; then
        error "Text file $TEXT_FILES not found or empty"
    fi
    
    # Check if out parameter is set
    if [[ $out = false ]]; then
        error "Output directory (-o) must be specified when using --text-files option"
    fi

    # Skip the header line
    is_first_line=true
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $is_first_line = true ]]; then
            is_first_line=false
            continue
        fi
        
        if [[ -z "$line" || "$line" == "#"* ]]; then
            continue  # Skip empty lines and comment lines
        fi
        
        process_text_file_line "$line"
        
        # Process this assembly and read pair
        # We'll continue with the existing workflow from here
        
########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

        # Run validation checks for this line's input
        if [[ $read_type != "paired" ]]; then
            error "When using --text-files, only paired-end mode is supported"
        fi

        # Make sure the output directory exists
        if [ ! -d $out ]; then mkdir $out;
        else
            echo "Warning: $out already exists."
        fi

        # Extract assembly name for subfolder creation
        assembly_filename=$(basename "$ASSEMBLY")
        # Remove common extensions to create clean folder name
        assembly_basename=${assembly_filename}
        assembly_basename=${assembly_basename%.fastq.gz}
        assembly_basename=${assembly_basename%.fq.gz}
        assembly_basename=${assembly_basename%.fa.gz}
        assembly_basename=${assembly_basename%.fasta.gz}
        assembly_basename=${assembly_basename%.fastq}
        assembly_basename=${assembly_basename%.fq}
        assembly_basename=${assembly_basename%.fa}
        assembly_basename=${assembly_basename%.fasta}

        # Create assembly-specific directory
        assembly_dir="${out}/${assembly_basename}"
        if [ ! -d ${assembly_dir} ]; then 
            mkdir ${assembly_dir}
            comm "Created assembly-specific directory: ${assembly_dir}"
        fi

        # Update bam output directory if not specified
        if [[ $bout = false ]]; then
            bout=${assembly_dir}
        fi
        if [ ! -d ${bout} ]; then mkdir ${bout}; fi

        # Handle compressed assembly file
        if [[ $ASSEMBLY == *.gz ]]; then
            comm "Input assembly file is compressed, decompressing to output directory"
            gunzip -c $ASSEMBLY > ${assembly_dir}/assembly.fa
            assembly_copied=true
        elif [ -f ${assembly_dir}/assembly.fa ]; then
            comm "Looks like the assembly file is already copied, but will re-transfer just in case to avoid truncation problems."
            cp $ASSEMBLY ${assembly_dir}/assembly.fa
            assembly_copied=true
        else
            comm "Making copy of assembly file $ASSEMBLY"
            cp $ASSEMBLY ${assembly_dir}/assembly.fa
            assembly_copied=true
        fi

        # Index the assembly
        if [ -f ${assembly_dir}/assembly.fa.bwt ]; then
            comm "Looks like there is a index of the assembly already. Skipping..."
        else
            comm "Indexing assembly file"
            bwa index ${assembly_dir}/assembly.fa
            if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi
        fi

        # Extract sample name from filename
        tmp=${FORWARD_READS##*/}
        if [[ $tmp == *"_1.fastq"* ]]; then
            sample=${tmp%_1.fastq*}
        elif [[ $tmp == *".1.fastq"* ]]; then
            sample=${tmp%.1.fastq*}
        elif [[ $tmp == *"_1.fq"* ]]; then
            sample=${tmp%_1.fq*}
        elif [[ $tmp == *".1.fq"* ]]; then
            sample=${tmp%.1.fq*}
        elif [[ $tmp == *"_1.fastq.gz"* ]]; then
            sample=${tmp%_1.fastq.gz*}
        elif [[ $tmp == *".1.fastq.gz"* ]]; then
            sample=${tmp%.1.fastq.gz*}
        elif [[ $tmp == *"_1.fq.gz"* ]]; then
            sample=${tmp%_1.fq.gz*}
        elif [[ $tmp == *".1.fq.gz"* ]]; then
            sample=${tmp%.1.fq.gz*}
        else
            # Use a more generic approach if pattern doesn't match
            sample=$(basename "$FORWARD_READS" | sed 's/\.[^.]*$//')
        fi
        
        if [[ ! -f ${bout}/${sample}.bam ]]; then
            comm "Aligning $FORWARD_READS and $REVERSE_READS back to assembly"
            
            # Handle compressed files
            if [[ $FORWARD_READS == *.gz && $REVERSE_READS == *.gz ]]; then
                bwa mem -v 1 -t $threads ${assembly_dir}/assembly.fa <(zcat $FORWARD_READS) <(zcat $REVERSE_READS) > ${bout}/${sample}.sam
            else
                bwa mem -v 1 -t $threads ${assembly_dir}/assembly.fa $FORWARD_READS $REVERSE_READS > ${bout}/${sample}.sam
            fi
            
            if [[ $? -ne 0 ]]; then error "Something went wrong with aligning $FORWARD_READS and $REVERSE_READS reads to the assembly. Exiting"; fi
            
            comm "Sorting the $sample alignment file"
            samtools sort -T ${assembly_dir}/tmp-samtools -@ $threads -O BAM -o ${bout}/${sample}.bam ${bout}/${sample}.sam
            if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiting..."; fi
            rm ${bout}/${sample}.sam
        else
            comm "skipping aligning $sample reads to assembly because ${bout}/${sample}.bam already exists."
        fi
        
    done < "$TEXT_FILES"
    
    # Exit after processing the text file
    announcement "The process of generating bam files from the input text file finished!!!"
    exit 0
fi

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################
# Make sure at least one binning method was chosen

# check if all parameters are entered
if [ $out = false ] || [ $ASSEMBLY = false ] ; then
	comm "Non-optional parameters -a and/or -o were not entered"
	help_message; exit 1
fi

#check if the assembly file exists
if [ ! -s $ASSEMBLY ]; then error "$ASSEMBLY does not exist. Exiting..."; fi

#check if parameter for bout dir is entered
if [[ $bout = false ]]; then
    # Extract assembly name for proper path creation
    assembly_filename=$(basename "$ASSEMBLY")
    # Create assembly basename by removing extensions
    assembly_basename=${assembly_filename}
    assembly_basename=${assembly_basename%.fastq.gz}
    assembly_basename=${assembly_basename%.fq.gz}
    assembly_basename=${assembly_basename%.fa.gz}
    assembly_basename=${assembly_basename%.fasta.gz}
    assembly_basename=${assembly_basename%.fastq}
    assembly_basename=${assembly_basename%.fq}
    assembly_basename=${assembly_basename%.fa}
    assembly_basename=${assembly_basename%.fasta}
    # Set bout directly to assembly_basename subfolder in output directory
    bout=${out}/${assembly_basename}
fi

# Function to check if a file exists and handle compressed files
check_file_exists() {
    local file=$1
    if [[ -s "$file" ]]; then
        echo "$file"
        return 0
    elif [[ -s "${file}.gz" ]]; then
        echo "${file}.gz"
        return 0
    else
        return 1
    fi
}

# Process input files based on different patterns
process_input_files() {
    local reads_found=false
    local forward_files=()
    local reverse_files=()
    
    # Case 1: Direct specification through -f and -r
    if [[ $DIRECT_INPUT_MODE = true ]]; then
        comm "Using directly specified read files"
        if [[ $FORWARD_READS != false ]]; then
            # For direct file paths, don't use check_file_exists, just check directly
            if [[ -s "$FORWARD_READS" ]]; then
                forward_files+=("$FORWARD_READS")
                reads_found=true
            else
                error "Forward read file $FORWARD_READS does not exist!"
            fi
        fi
        
        if [[ $REVERSE_READS != false ]]; then
            if [[ -s "$REVERSE_READS" ]]; then
                reverse_files+=("$REVERSE_READS")
                reads_found=true
            else
                error "Reverse read file $REVERSE_READS does not exist!"
            fi
        fi
    
    # Case 2: Just a base name is provided (first positional argument)
    elif [[ $# -eq 1 && ! "$1" == *"*"* && ! "$1" == *"fastq"* && ! "$1" == *"fq"* ]]; then
        base_name=$1
        comm "Processing basename: $base_name"
        
        # Check for {1,2}.fastq and other variations
        for suffix1 in "_1.fastq" "_1.fq" ".1.fastq" ".1.fq"; do
            for suffix2 in "_2.fastq" "_2.fq" ".2.fastq" ".2.fq"; do
                f_file=$(check_file_exists "${base_name}${suffix1}")
                r_file=$(check_file_exists "${base_name}${suffix2}")
                
                if [[ $? -eq 0 ]] && [[ -s "$r_file" ]]; then
                    forward_files+=("$f_file")
                    reverse_files+=("$r_file")
                    reads_found=true
                    break 2
                fi
            done
        done
        
        if [[ $reads_found = false ]]; then
            error "Could not find paired files with basename $base_name"
        fi
    
    # Case 3: Wildcard pattern or explicit file names
    else
        comm "Processing wildcard or explicit file patterns"
        # Original processing logic for multiple files
        for num in "$@"; do
            if [[ $read_type = paired ]]; then
                if [[ $num == *"${F_reads_suffix}"* ]]; then
                    f_file=$(check_file_exists "$num")
                    if [[ $? -eq 0 ]]; then
                        forward_files+=("$f_file")
                        
                        r_name=${num%${F_reads_suffix}}${R_reads_suffix}
                        r_file=$(check_file_exists "$r_name")
                        if [[ $? -eq 0 ]]; then
                            reverse_files+=("$r_file")
                            reads_found=true
                        else
                            error "Matching reverse read for $num (expected: $r_name) does not exist!"
                        fi
                    fi
                fi
            else
                # For single-end or interleaved
                if [[ $num == *".fastq"* || $num == *".fq"* ]]; then
                    file=$(check_file_exists "$num")
                    if [[ $? -eq 0 ]]; then
                        forward_files+=("$file")
                        reads_found=true
                    fi
                fi
            fi
        done
    fi
    
    if [[ $reads_found = false ]]; then
        error "No valid read files were found! Check your input parameters."
    fi
    
    # Return the results through global arrays
    PROCESSED_FORWARD_FILES=("${forward_files[@]}")
    PROCESSED_REVERSE_FILES=("${reverse_files[@]}")
}

# Process input files
PROCESSED_FORWARD_FILES=()
PROCESSED_REVERSE_FILES=()
process_input_files "$@"

comm "Entered read type: $read_type"

if [[ ${#PROCESSED_FORWARD_FILES[@]} -eq 0 ]]; then
    error "No valid forward read files found!"
fi

if [[ $read_type = paired && ${#PROCESSED_REVERSE_FILES[@]} -eq 0 ]]; then
    error "No valid reverse read files found for paired-end mode!"
fi

if [[ $read_type = paired && ${#PROCESSED_FORWARD_FILES[@]} -ne ${#PROCESSED_REVERSE_FILES[@]} ]]; then
    error "Number of forward and reverse read files must be the same!"
fi

comm "${#PROCESSED_FORWARD_FILES[@]} forward read files detected"
if [[ $read_type = paired ]]; then
    comm "${#PROCESSED_REVERSE_FILES[@]} reverse read files detected"
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################


########################################################################################################
########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################
########################################################################################################
announcement "ALIGNING READS TO MAKE COVERAGE FILES"

# setting up the output folder
if [ ! -d $out ]; then mkdir $out;
else
	echo "Warning: $out already exists."
fi

# Extract assembly name for subfolder creation
assembly_filename=$(basename "$ASSEMBLY")
# Remove common extensions to create clean folder name
assembly_basename=${assembly_filename}
assembly_basename=${assembly_basename%.fastq.gz}
assembly_basename=${assembly_basename%.fq.gz}
assembly_basename=${assembly_basename%.fa.gz}
assembly_basename=${assembly_basename%.fasta.gz}
assembly_basename=${assembly_basename%.fastq}
assembly_basename=${assembly_basename%.fq}
assembly_basename=${assembly_basename%.fa}
assembly_basename=${assembly_basename%.fasta}

# Create assembly-specific directory
assembly_dir="${out}/${assembly_basename}"
if [ ! -d ${assembly_dir} ]; then 
    mkdir ${assembly_dir}
    comm "Created assembly-specific directory: ${assembly_dir}"
fi

# Update bam output directory if not specified
if [[ $bout = false ]]; then
    bout=${assembly_dir}
fi
if [ ! -d ${bout} ]; then mkdir ${bout}; fi

# Handle compressed assembly file
if [[ $ASSEMBLY == *.gz ]]; then
	comm "Input assembly file is compressed, decompressing to output directory"
	gunzip -c $ASSEMBLY > ${assembly_dir}/assembly.fa
	assembly_copied=true
elif [ -f ${assembly_dir}/assembly.fa ]; then
	comm "Looks like the assembly file is already copied, but will re-transfer just in case to avoid truncation problems."
	cp $ASSEMBLY ${assembly_dir}/assembly.fa
	assembly_copied=true
else
	comm "Making copy of assembly file $ASSEMBLY"
	cp $ASSEMBLY ${assembly_dir}/assembly.fa
	assembly_copied=true
fi

# Index the assembly
if [ -f ${assembly_dir}/assembly.fa.bwt ]; then
	comm "Looks like there is a index of the assembly already. Skipping..."
else
	comm "Indexing assembly file"
	bwa index ${assembly_dir}/assembly.fa
	if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi
fi

# Process the read files for alignment
for ((i=0; i<${#PROCESSED_FORWARD_FILES[@]}; i++)); do
    if [ $read_type = paired ]; then
        reads_1="${PROCESSED_FORWARD_FILES[$i]}"
        reads_2="${PROCESSED_REVERSE_FILES[$i]}"
        
        # Extract sample name from filename
        tmp=${reads_1##*/}
        if [[ $tmp == *"_1.fastq"* ]]; then
            sample=${tmp%_1.fastq*}
        elif [[ $tmp == *".1.fastq"* ]]; then
            sample=${tmp%.1.fastq*}
        elif [[ $tmp == *"_1.fq"* ]]; then
            sample=${tmp%_1.fq*}
        elif [[ $tmp == *".1.fq"* ]]; then
            sample=${tmp%.1.fq*}
        elif [[ $tmp == *"_1.fastq.gz"* ]]; then
            sample=${tmp%_1.fastq.gz*}
        elif [[ $tmp == *".1.fastq.gz"* ]]; then
            sample=${tmp%.1.fastq.gz*}
        elif [[ $tmp == *"_1.fq.gz"* ]]; then
            sample=${tmp%_1.fq.gz*}
        elif [[ $tmp == *".1.fq.gz"* ]]; then
            sample=${tmp%.1.fq.gz*}
        else
            # Use a more generic approach if pattern doesn't match
            sample=$(basename "$reads_1" | sed 's/\.[^.]*$//')
        fi
        
        if [[ ! -f ${bout}/${sample}.bam ]]; then
            comm "Aligning $reads_1 and $reads_2 back to assembly"
            
            # Handle compressed files
            if [[ $reads_1 == *.gz && $reads_2 == *.gz ]]; then
                bwa mem -v 1 -t $threads ${assembly_dir}/assembly.fa <(zcat $reads_1) <(zcat $reads_2) > ${bout}/${sample}.sam
            else
                bwa mem -v 1 -t $threads ${assembly_dir}/assembly.fa $reads_1 $reads_2 > ${bout}/${sample}.sam
            fi
            
            if [[ $? -ne 0 ]]; then error "Something went wrong with aligning $reads_1 and $reads_2 reads to the assembly. Exiting"; fi
            
            comm "Sorting the $sample alignment file"
            samtools sort -T ${assembly_dir}/tmp-samtools -@ $threads -O BAM -o ${bout}/${sample}.bam ${bout}/${sample}.sam
            if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiting..."; fi
            rm ${bout}/${sample}.sam
        else
            comm "skipping aligning $sample reads to assembly because ${bout}/${sample}.bam already exists."
        fi
    else
        # single end or interleaved reads
        reads="${PROCESSED_FORWARD_FILES[$i]}"
        tmp=${reads##*/}
        sample=$(basename "$reads" | sed 's/\.[^.]*$//')
        
        if [[ ! -f ${bout}/${sample}.bam ]]; then
            comm "Aligning $reads back to assembly, and sorting the alignment"
            
            # Handle compressed files
            if [[ $reads == *.gz ]]; then
                if [ $read_type = single ]; then
                    bwa mem -t $threads ${assembly_dir}/assembly.fa <(zcat $reads) > ${bout}/${sample}.sam
                elif [ $read_type = interleaved ]; then
                    bwa mem -v 1 -p -t $threads ${assembly_dir}/assembly.fa <(zcat $reads) > ${bout}/${sample}.sam
                fi
            else
                if [ $read_type = single ]; then
                    bwa mem -t $threads ${assembly_dir}/assembly.fa $reads > ${bout}/${sample}.sam
                elif [ $read_type = interleaved ]; then
                    bwa mem -v 1 -p -t $threads ${assembly_dir}/assembly.fa $reads > ${bout}/${sample}.sam
                fi
            fi
            
            if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
            
            comm "Sorting the $sample alignment file"
            samtools sort -T ${assembly_dir}/tmp-samtools -@ $threads -O BAM -o ${bout}/${sample}.bam ${bout}/${sample}.sam
            if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiting..."; fi
            rm ${bout}/${sample}.sam
        else
            comm "skipping aligning $sample reads to assembly because ${bout}/${sample}.bam already exists."
        fi
    fi
done

announcement "The process of generating bam files finished!!!"s