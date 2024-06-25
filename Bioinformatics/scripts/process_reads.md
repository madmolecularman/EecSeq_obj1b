# This file analyzes the raw sequenced files and processes them for analysis

## Enmasse copying files

Files are in the /RAID_STORAGE2/Raw_Data/NB_and_EAGER_2024/01.RawData/ directory. There are multiple directories for an individual sample which look like different lanes. I am not sure about the naming convention. but initially I will pull all of these samples in the /home/jgreen/eager_obj1b/01_process_reads/raw directory.

```bash
## Code for bash script
# Define an array with the names of all directories in the current directory
directories=(*/)

# Start a loop over each directory in the array
for dir in "${directories[@]}"; do

  # Check if the directory name does not start with the letter 'D'
  if [[ ! $dir =~ ^D ]]; then

    # Use the find command to search for files in the current directory
    # -type f specifies that we're looking for files
    # -name "*.fq.gz" specifies that we're looking for files that end with .fq.gz
    # -exec cp {} /home/jgreen/eager_obj1b/01_process_reads/raw/ \; specifies that for each file found, 
    # we should execute the cp command to copy the file to the /home/jgreen/eager_obj1b/01_process_reads/raw/ directory
    find "$dir" -type f -name "*.fq.gz" -exec cp {} /home/jgreen/eager_obj1b/01_process_reads/raw/ \;
  fi

# End the loop
done


## Code to run in /RAID_STORAGE2/Raw_Data/NB_and_EAGER_2024/01.RawData/
directories=(*/)
for dir in "${directories[@]}"; do
  if [[ ! $dir =~ ^D ]]; then
    find "$dir" -type f -name "*.fq.gz" -exec cp {} /home/jgreen/eager_obj1b/01_process_reads/raw/ \;
  fi
done
```

## Edits names to remove redudancy

```bash
# Define an array with the strings to be removed from file names
strings=("CKDL2400187-1A_H3VJJDSXC_" "CKDL240018761-1A_HHC5NDSXC_" "CKDL240018760-1A_HHC5NDSXC_" "CKDL240018756-1A_H3VJJDSXC_" "CKDL240018757-1A_H3VJJDSXC_" "CKDL240018758-1A_H3VJJDSXC_" "CKDL240018759-1A_H3VJJDSXC_")

# Start a loop over each string in the array
for str in "${strings[@]}"; do

  # Find files in the /home/jgreen/eager_obj1b/01_process_reads/raw/ directory that contain the current string in their names
  # and rename those files to remove the string
  find /home/jgreen/eager_obj1b/01_process_reads/raw/ -type f -name "*_${str}*" | while read file; do mv "$file" "${file/${str}/}"; done

# End the loop
done
```

## Run fastqc and multiqc to visualize files
```bash
fastqc -f fastq -t 20 *.fastq.gz
multiqc -o raw_multiqc_run .
```

## Results
Multiqc plots show an uneven distribution of reads across samples particularly b2 files. Looking at the "Per Base Sequence Content" in the report shows probes in samples b2, b2_2, b2_3, b6, b6_2, b6_3, and b7. These files will have to have probes filtered out of them. [multiqc report](https://github.com/madmolecularman/EecSeq_obj1b/blob/main/Bioinformatics/01_process_reads/multiqc_report_240620.html)


## Renaming reads

The main issues with the way samples have been labeled is that they are from 4 different lanes, each with a different capture, and each with multiple samples within each one. I am hesitant to treat each files as it is labeled.

## Processing short reads on files that are not demultiplexed

I would like to atleast bin them out into their proper Capture groups. We know that Lane 3 is Capture 3, Lane 1 is Capture 10, Lane 2 is Capture 11, and Lane 4 is Capture 12 due to their 100% adapter content in the multiqc file. We will mv the files in order to add their capture designation to the front of the sequence. I also need to get rid of the double underscore from the filenames.

## Add proper capture name to each file
```bash
# Scary way 
# Define an associative array with the patterns as keys and the prefixes as values
declare -A patterns=( ["*L3*"]="capture_3_" ["*L1*"]="capture_10_" ["*L2*"]="capture_11_" ["*L4*"]="capture_12_" )

# Loop through the associative array
for pattern in "${!patterns[@]}"; do
  # Use find to locate files matching the current pattern and rename them by prepending the corresponding prefix
  find /home/jgreen/eager_obj1b/01_process_reads/raw -type f -name "$pattern" | while read file; do
    mv "$file" "$(dirname "$file")/${patterns[$pattern]}$(basename "$file")"
  done
done

# Less scary but more reassuring way
find /home/jgreen/eager_obj1b/01_process_reads/raw -type f -name '*L3*' | while read file; do mv "$file" "$(dirname "$file")/capture_3_$(basename "$file")"; done
find /home/jgreen/eager_obj1b/01_process_reads/raw -type f -name '*L1*' | while read file; do mv "$file" "$(dirname "$file")/capture_10_$(basename "$file")"; done
find /home/jgreen/eager_obj1b/01_process_reads/raw -type f -name '*L2*' | while read file; do mv "$file" "$(dirname "$file")/capture_11_$(basename "$file")"; done
find /home/jgreen/eager_obj1b/01_process_reads/raw -type f -name '*L4*' | while read file; do mv "$file" "$(dirname "$file")/capture_12_$(basename "$file")"; done
```

## Change _1. and _2. to R1 and R2

```bash
find . -type f \( -name '*_1.*' -o -name '*_2.*' \) | while read file; do 
  mv "$file" "${file//_1./_R1.}"; 
  mv "$file" "${file//_2./_R2.}"; 
done
```
## Demultiplexing capture files

Use process_shortreads to demultiplex files from raw/ to the demux/ folder. One of the major problems is that we are not sure if novogene was able to separate out the captures and samples for each file. We need to create a demultiplexing run that will filter the index and the inline barcode into the proper file. The question here is can I demultiplex just the index and then the adapter inline barcode?

Its important to understand what ![process_shortreads](https://catchenlab.life.illinois.edu/stacks/comp/process_shortreads.php) is doing and how to create the barcode files ![Stacks manual](https://catchenlab.life.illinois.edu/stacks/manual/index.php#clean)


## Using duplicated barcodes and inline_inline to investigate ambigous barcode drop all files

New directory for double inline barcodes. I am wanting to test out how we retain reads with double barcodes.

Make new barcode file

``bash
nano single_barcode.txt
```

Make sure to put tabs as the spaces

```text
TTAGGC  TTAGGC  capture147
TAGCTT  TAGCTT  capture10
GGCTAC  GGCTAC  capture13
CTTGCA  CTTGCA  capture16
```

```bash
mkdir demultiplex
cd demultiplex
```

Make directories for each sample by making a sample list file, then cat the file but add the parenthesis and space needed for the loop.

``bash
nano directory_list.txt
```

```bash
"capture_1_b2"
"capture_1_g1"
"capture_1_g7"
"capture_1_k1"
"capture_1_k3"
"capture_1_m5"
"capture_1_m7"
"capture_1_n1"
"capture_1_n2"
"capture_4_b1"
"capture_4_b2"
"capture_4_g1"
"capture_4_g7"
"capture_4_k1"
"capture_4_k3"
"capture_4_m5"
"capture_4_m7"
"capture_4_n1"
"capture_4_n2"
"capture_7_b1"
"capture_7_b2"
"capture_7_g1"
"capture_7_g7"
"capture_7_k1"
"capture_7_k3"
"capture_7_m5"
"capture_7_m7"
"capture_7_n1"
"capture_7_n2"
"capture_10_b10"
"capture_10_b1"
"capture_10_b2"
"capture_10_b6"
"capture_10_b7"
"capture_10_g1"
"capture_10_g5"
"capture_10_g6"
"capture_10_g7"
"capture_10_g8"
"capture_10_k1"
"capture_10_k2"
"capture_10_k3"
"capture_10_k5"
"capture_10_k8"
"capture_10_m1"
"capture_10_m3"
"capture_10_m5"
"capture_10_m7"
"capture_10_m9"
"capture_10_n1"
"capture_10_n2"
"capture_10_n3"
"capture_10_n4"
"capture_10_n5"
"capture_13_b10"
"capture_13_b1"
"capture_13_b2"
"capture_13_b6"
"capture_13_b7"
"capture_13_g1"
"capture_13_g5"
"capture_13_g6"
"capture_13_g7"
"capture_13_g8"
"capture_13_k1"
"capture_13_k2"
"capture_13_k3"
"capture_13_k5"
"capture_13_k8"
"capture_13_m1"
"capture_13_m3"
"capture_13_m5"
"capture_13_m7"
"capture_13_m9"
"capture_13_n1"
"capture_13_n2"
"capture_13_n3"
"capture_13_n4"
"capture_13_n5"
"capture_16_b10"
"capture_16_b1"
"capture_16_b2"
"capture_16_b6"
"capture_16_b7"
"capture_16_g1"
"capture_16_g5"
"capture_16_g6"
"capture_16_g7"
"capture_16_g8"
"capture_16_k1"
"capture_16_k2"
"capture_16_k3"
"capture_16_k5"
"capture_16_k8"
"capture_16_m1"
"capture_16_m3"
"capture_16_m5"
"capture_16_m7"
"capture_16_m9"
"capture_16_n1"
"capture_16_n2"
"capture_16_n3"
"capture_16_n4"
"capture_16_n5"
```

## Making directories for all Capture files
```bash
declare -a StringArray=("capture_1_b2" "capture_1_g1" "capture_1_g7" "capture_1_k1" "capture_1_k3" "capture_1_m5" "capture_1_m7" "capture_1_n1" "capture_1_n2" "capture_4_b1" "capture_4_b2" "capture_4_g1" "capture_4_g7" "capture_4_k1" "capture_4_k3" "capture_4_m5" "capture_4_m7" "capture_4_n1" "capture_4_n2" "capture_7_b1" "capture_7_b2" "capture_7_g1" "capture_7_g7" "capture_7_k1" "capture_7_k3" "capture_7_m5" "capture_7_m7" "capture_7_n1" "capture_7_n2" "capture_10_b10" "capture_10_b1" "capture_10_b2" "capture_10_b6" "capture_10_b7" "capture_10_g1" "capture_10_g5" "capture_10_g6" "capture_10_g7" "capture_10_g8" "capture_10_k1" "capture_10_k2" "capture_10_k3" "capture_10_k5" "capture_10_k8" "capture_10_m1" "capture_10_m3" "capture_10_m5" "capture_10_m7" "capture_10_m9" "capture_10_n1" "capture_10_n2" "capture_10_n3" "capture_10_n4" "capture_10_n5" "capture_13_b10" "capture_13_b1" "capture_13_b2" "capture_13_b6" "capture_13_b7" "capture_13_g1" "capture_13_g5" "capture_13_g6" "capture_13_g7" "capture_13_g8" "capture_13_k1" "capture_13_k2" "capture_13_k3" "capture_13_k5" "capture_13_k8" "capture_13_m1" "capture_13_m3" "capture_13_m5" "capture_13_m7" "capture_13_m9" "capture_13_n1" "capture_13_n2" "capture_13_n3" "capture_13_n4" "capture_13_n5" "capture_16_b10" "capture_16_b1" "capture_16_b2" "capture_16_b6" "capture_16_b7" "capture_16_g1" "capture_16_g5" "capture_16_g6" "capture_16_g7" "capture_16_g8" "capture_16_k1" "capture_16_k2" "capture_16_k3" "capture_16_k5" "capture_16_k8" "capture_16_m1" "capture_16_m3" "capture_16_m5" "capture_16_m7" "capture_16_m9" "capture_16_n1" "capture_16_n2" "capture_16_n3" "capture_16_n4" "capture_16_n5")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

## Running process_shortreads on all files


Use single barcode file located here: 

> /home/jgreen/eager_obj1b/01_process_reads/barcodes/single_barcode.txt

Code for process_shortreads on all Captures

Test on one sample

```bash
declare -a StringArray=("capture_1_b2")

for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/raw/${i}_R1.fq.gz -2 $WORKING_DIR/01_process_reads/raw/${i}_R2.fq.gz -o $WORKING_DIR/01_process_reads/demultiplex/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/single_barcode.txt --inline_inline -r -c -q -D
done
```

```bash
declare -a StringArray=("capture_1_b2" "capture_1_g1" "capture_1_g7" "capture_1_k1" "capture_1_k3" "capture_1_m5" "capture_1_m7" "capture_1_n1" "capture_1_n2" "capture_4_b1" "capture_4_b2" "capture_4_g1" "capture_4_g7" "capture_4_k1" "capture_4_k3" "capture_4_m5" "capture_4_m7" "capture_4_n1" "capture_4_n2" "capture_7_b1" "capture_7_b2" "capture_7_g1" "capture_7_g7" "capture_7_k1" "capture_7_k3" "capture_7_m5" "capture_7_m7" "capture_7_n1" "capture_7_n2" "capture_10_b10" "capture_10_b1" "capture_10_b2" "capture_10_b6" "capture_10_b7" "capture_10_g1" "capture_10_g5" "capture_10_g6" "capture_10_g7" "capture_10_g8" "capture_10_k1" "capture_10_k2" "capture_10_k3" "capture_10_k5" "capture_10_k8" "capture_10_m1" "capture_10_m3" "capture_10_m5" "capture_10_m7" "capture_10_m9" "capture_10_n1" "capture_10_n2" "capture_10_n3" "capture_10_n4" "capture_10_n5" "capture_13_b10" "capture_13_b1" "capture_13_b2" "capture_13_b6" "capture_13_b7" "capture_13_g1" "capture_13_g5" "capture_13_g6" "capture_13_g7" "capture_13_g8" "capture_13_k1" "capture_13_k2" "capture_13_k3" "capture_13_k5" "capture_13_k8" "capture_13_m1" "capture_13_m3" "capture_13_m5" "capture_13_m7" "capture_13_m9" "capture_13_n1" "capture_13_n2" "capture_13_n3" "capture_13_n4" "capture_13_n5" "capture_16_b10" "capture_16_b1" "capture_16_b2" "capture_16_b6" "capture_16_b7" "capture_16_g1" "capture_16_g5" "capture_16_g6" "capture_16_g7" "capture_16_g8" "capture_16_k1" "capture_16_k2" "capture_16_k3" "capture_16_k5" "capture_16_k8" "capture_16_m1" "capture_16_m3" "capture_16_m5" "capture_16_m7" "capture_16_m9" "capture_16_n1" "capture_16_n2" "capture_16_n3" "capture_16_n4" "capture_16_n5")

for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/raw/${i}_R1.fq.gz -2 $WORKING_DIR/01_process_reads/raw/${i}_R2.fq.gz -o $WORKING_DIR/01_process_reads/demultiplex/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/single_barcode.txt --inline_inline -r -c -q -D
done
```



## Make all of the process_shortreads short logs for review

```bash
declare -a StringArray=("capture_1_b2" "capture_1_g1" "capture_1_g7" "capture_1_k1" "capture_1_k3" "capture_1_m5" "capture_1_m7" "capture_1_n1" "capture_1_n2" "capture_4_b1" "capture_4_b2" "capture_4_g1" "capture_4_g7" "capture_4_k1" "capture_4_k3" "capture_4_m5" "capture_4_m7" "capture_4_n1" "capture_4_n2" "capture_7_b1" "capture_7_b2" "capture_7_g1" "capture_7_g7" "capture_7_k1" "capture_7_k3" "capture_7_m5" "capture_7_m7" "capture_7_n1" "capture_7_n2" "capture_10_b10" "capture_10_b1" "capture_10_b2" "capture_10_b6" "capture_10_b7" "capture_10_g1" "capture_10_g5" "capture_10_g6" "capture_10_g7" "capture_10_g8" "capture_10_k1" "capture_10_k2" "capture_10_k3" "capture_10_k5" "capture_10_k8" "capture_10_m1" "capture_10_m3" "capture_10_m5" "capture_10_m7" "capture_10_m9" "capture_10_n1" "capture_10_n2" "capture_10_n3" "capture_10_n4" "capture_10_n5" "capture_13_b10" "capture_13_b1" "capture_13_b2" "capture_13_b6" "capture_13_b7" "capture_13_g1" "capture_13_g5" "capture_13_g6" "capture_13_g7" "capture_13_g8" "capture_13_k1" "capture_13_k2" "capture_13_k3" "capture_13_k5" "capture_13_k8" "capture_13_m1" "capture_13_m3" "capture_13_m5" "capture_13_m7" "capture_13_m9" "capture_13_n1" "capture_13_n2" "capture_13_n3" "capture_13_n4" "capture_13_n5" "capture_16_b10" "capture_16_b1" "capture_16_b2" "capture_16_b6" "capture_16_b7" "capture_16_g1" "capture_16_g5" "capture_16_g6" "capture_16_g7" "capture_16_g8" "capture_16_k1" "capture_16_k2" "capture_16_k3" "capture_16_k5" "capture_16_k8" "capture_16_m1" "capture_16_m3" "capture_16_m5" "capture_16_m7" "capture_16_m9" "capture_16_n1" "capture_16_n2" "capture_16_n3" "capture_16_n4" "capture_16_n5")

for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}.process_shortreads.short.log
done
```

Need to pull Capture[1-4] reads into the clean directory

```bash
#Capture1
#Note All Capture 1 are read as Capture 2 in the process_shortreads output. Need to pull Capture2.1.fq.gz and Capture2.2.fq.gz and rename them to the proper Capture category. The only ones to omit from this are G3 and N2.
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

declare -a StringArray=("Capture1_G3"  "Capture1_N2")
for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture2
#Note All Capture 2 are read as Capture 2 in the process_shortreads output. Need to pull Capture1.1.fq.gz and Capture1.2.fq.gz and rename them to the proper Capture category. The only ones to omit from this are G3 and N2.
declare -a StringArray=( "Capture2_B3" "Capture2_B4" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

declare -a StringArray=("Capture2_G3"  "Capture2_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture3
declare -a StringArray=("Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture3.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture3.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture4
declare -a StringArray=("Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture4.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture4.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done
```

Using chatGPT we can simplify this code to the following

```bash
# List of prefixes for each capture
declare -a captures=("Capture1" "Capture2")

for capture in "${captures[@]}"
do
    # List of samples for each capture
    declare -a samples=("N1" "N2")
    
    for sample in "${samples[@]}"
    do
        # Construct the sample ID
        sample_id="${capture}_${sample}"
        
        # Copy the files
        cat "$WORKING_DIR/01_process_reads/demux/${sample_id}/${capture}.1.fq.gz" >> "$WORKING_DIR/01_process_reads/clean/Capture2_N2.F.fq.gz"
        cat "$WORKING_DIR/01_process_reads/demux/${sample_id}/${capture}.2.fq.gz" >> "$WORKING_DIR/01_process_reads/clean/Capture2_N2.R.fq.gz"
    done
done
```

## Analyze clean files

## Run fastqc and multiqc to visualize files

```bash
fastqc -f fastq -t 20 *.fq.gz
multiqc . -o post_demux_multiqc_report
```

## Merging PE files that have close insert sizes

We will be merging PE FASTQ files that are overlapping into SE files. Using Geneious we found that Capture 2,3, and 4 have read intervals less than 2x their read length which can lead to more overlaps between PE read. 

![Alt text](B3_subsample_insertsize-1.png)

The result will be a PE and SE fastq file which will be mapped to the reference. From here we will merged the PE and SE bam files together and proceed with normalization.

## Merging PE reads from a single individual sample

To merge PE reads we will be using the program [FASTP](https://github.com/OpenGene/fastp)

```bash
#Single file
fastp -m --in1 $WORKING_DIR/01_process_reads/clean/Capture1_B3.F.fq.gz --out1 Capture1_B3.F.fastp.fq.gz --in2 $WORKING_DIR/01_process_reads/clean/Capture1_B3.R.fq.gz --out2 Capture1_B3.R.fastp.fq.gz --merged_out Capture1_B3_merged.fq.gz -D -p

#All files
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do
fastp -m --in1 ../clean/${i}.F.fq.gz --out1 ${i}.R1.fq.gz --in2 ../clean/${i}.R.fq.gz --out2 ${i}.R2.fq.gz --merged_out ${i}.R1.fq.gz -D -p -R "${i}" -j ${i}.json -h ${i}.html
done
```

## Run fastqc and multiqc to visualize files

```bash
fastqc -f fastq -t 20 *.fq.gz
multiqc . -o post_merge_multiqc_report
```

## Create directory in 02_ddocent for single and PE fastp filtered reads
```bash
mkdir 01_PE
mkdir 02_SE
mkdir 04_probes
```

Note: Naming convention for dDocent files is very important
* PE must contain F.fq.gz and R.fq.gz for their identifier with correct prefix names listed in the popmap file
* SE must contain F.fq.gz for their identifier with correct prefix names listed in popmap file.

If you need to rename the filenames

```bash
find . -type f -name "*.R1*" -exec bash -c 'mv "$0" "${0/R1/F}"' {} \;
find . -type f -name "*.R2*" -exec bash -c 'mv "$0" "${0/R2/R}"' {} \;
find . -type f -name "*R1*" -exec bash -c 'mv "$0" "${0/R1/F}"' {} \;
```

## Create a reference.file from reference.fasta that is originally the masked genome fasta

```bash
#Copy genome
cp /RAID_STORAGE2/Shared_Data/Oyster_Genome/masked/masked.cvir.genome.fasta .
#Index genome
samtools faidx masked.cvir.genome.fasta
#Make genome.file
mawk -v OFS='\t' {'print $1,$2'} masked.cvir.genome.fasta.fai > masked.genome.file
```

## Download new dDocent script in scripts directory
```bash
wget "https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/blob/main/scripts/dDocent_ngs"
```

I receive an error when downloading this as most of the file is in an html format. Create a file using nano. Copy and paste the whole script from the github into this new file

```bash 
nano dDocent_ngs.sh
```

change permissions on the files to allow it to be executable

```bash
chmod +x dDocent_ngs.sh
```

## Run the dDocent script in the directory containing our newly named files

```bash
/../scripts/dDocent_ngs.sh config.file
```

## Masking old bed files with haplotig bed file

Copy old bed files

```bash
cp /home/Genomic_Resources/C_virginica/bed_and_GFF/*bed .
```

Pull haplotig bed from github repo https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/blob/main/Haplotig_Masking/Haplotig_Masking.md

```bash
wget https://raw.githubusercontent.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/main/Haplotig_Masking/Output/haplotigs.bed
```

Use bedtools subtract to remove haplotigs present in the bed file from the other bed files

```bash
declare -a StringArray=("sorted.ref3.0.CDS" "sorted.ref3.0.exon" "sorted.ref3.0.gene" "sorted.ref3.0.UTR" )
for i in "${StringArray[@]}"
do
bedtools subtract -a ${i}.bed -b haplotigs.bed > ${i}.hmask.bed
done

declare -a StringArray=("sorted.ref3.0.CDS.sc" "sorted.ref3.0.exon.sc" "sorted.ref3.0.gene.sc" "sorted.ref3.0.UTR.sc" )
for i in "${StringArray[@]}"
do
bedtools subtract -a ${i}.bed -b haplotigs.bed > ${i}.hmask.bed
done
```

## Link all the .F.bam files from ddocent analysis to the 03_mapping folder

```bash
ln -s $WORKING_DIR/02_ddocent/*.F.bam .
```