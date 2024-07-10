# This file shows the work in analyzing individual fastq, fasta, and bam files.

In process_reads.md we were able to rename, demultiplex, trim, and map reads. process_short_reads program was able to discard probes that were present in our gDNA reads. Using fastp we found that nearly 1/3 rd of reads in each sample were converted into single end (SE) reads. For ddocent mapping we separated the PE and SE reads into two different runs. As a result we have bam files for both the PE and SE reads.

## Run ddocent on all individual files

Make config.file

```bash
nano config.file

Number of Processors
12
Maximum Memory
0
Trimming
no
Assembly?
no
Type_of_Assembly
PE
Clustering_SimilarityPercent
0.90
Minimum within individual coverage level to include a read for assembly (K1)
2
Minimum number of individuals a read must be present in to include for assembly (K2)
2
Mapping_Reads?
yes
Mapping_Match_Value
1
Mapping_MisMatch_Value
3
Mapping_GapOpen_Penalty
5
Calling_SNPs?
no
Email
gree9242@uri.edu
```

## Run the dDocent script in the directory containing our newly named files

```bash
../../scripts/dDocent_ngs.sh config.file
```

# Merge SE and PE bam files into a single bam file

```bash
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_b1.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_b1.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_b1.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_b2.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_b2.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_b2.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_g1.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_g1.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_g1.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_g7.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_g7.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_g7.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_k1.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_k1.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_k1.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_k3.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_k3.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_k3.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_m5.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_m5.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_m5.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_m7.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_m7.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_m7.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_n1.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_n1.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_n1.F.bam
done

declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_n2.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_n2.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_n2.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_b6.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_b6.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_b6.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_b7.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_b7.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_b7.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_b10.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_b10.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_b10.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_g5.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_g5.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_g5.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_g6.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_g6.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_g6.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_g8.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_g8.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_g8.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_k2.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_k2.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_k2.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_k5.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_k5.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_k5.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_k8.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_k8.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_k8.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_m1.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_m1.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_m1.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_m3.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_m3.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_m3.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_m9.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_m9.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_m9.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_n3.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_n3.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_n3.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_n4.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_n4.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_n4.F.bam
done

declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
samtools merge --threads 20 -o /home/jgreen/eager_obj1b/02_ddocent/03_mergedPESE/${i}_n5.F.bam /home/jgreen/eager_obj1b/03_mapping/02_PE/${i}_n5.F.bam /home/jgreen/eager_obj1b/03_mapping/03_SE/${i}_n5.F.bam
done

## Simplified code for samples in all captures

```bash
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
declare -a SampleSuffixes=("b1" "b2" "g1" "g7" "k1" "k3" "m5" "m7" "n1" "n2")

for capture in "${StringArray[@]}"; do
  for suffix in "${SampleSuffixes[@]}"; do
    samtools merge --threads 20 -o /home/jgreen/eager_obj1b/03_mapping/01_merged/${capture}_${suffix}.F.bam /home/jgreen/eager_obj1b/02_ddocent/02_PE/${capture}_${suffix}.F.bam /home/jgreen/eager_obj1b/02_ddocent/03_SE/${capture}_${suffix}.F.bam
  done
done

## Simplified code for samples only in captures 10, 13, and 16

```bash
declare -a StringArray=("capture10" "capture13" "capture16")
declare -a SampleSuffixes=("b6" "b7" "b10" "g5" "g6" "g8" "k2" "k5" "k8" "m1" "m3" "m9" "n3" "n4" "n5")

for capture in "${StringArray[@]}"; do
  for suffix in "${SampleSuffixes[@]}"; do
    samtools merge --threads 20 -o /home/jgreen/eager_obj1b/03_mapping/01_merged/${capture}_${suffix}.F.bam /home/jgreen/eager_obj1b/02_ddocent/02_PE/${capture}_${suffix}.F.bam /home/jgreen/eager_obj1b/02_ddocent/03_SE/${capture}_${suffix}.F.bam
  done
done

```

## Normalizing F.bam files for each individual subsample

These lines of code were ran in each of the directories 03_mergedPESE

Figure out a normal value for either each individual groups or across all groups

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
alignment_reads=$(samtools view -c ${i}.F.bam)
echo " ${i} , Number of alignment reads: $alignment_reads" >> bamalignedreads.txt
echo " ${i},$alignment_reads" >> bamalignedreads.csv 
done

for i in  "${StringArray[@]}"
do
alignment_reads=$(samtools view -c ${i}.F.bam)
echo "${i},$alignment_reads" >> bamalignedreads.csv 
done

echo {Sample,Reads} | tr ' ' , > header
cat header bamalignedreads.csv > bamalignedread.csv

declare -a StringArray=("b1" "b2" "b6" "b7" "b10" "g1" "g5" "g6" "g7" "g8" "k1" "k2" "k3" "k5" "k8" "m1" "m3" "m5" "m7" "m9" "n1" "n2" "n3" "n4" "n5")
for i in  "${StringArray[@]}"
do
grep "${i}" bamalignedreads.txt
done
```

Use the previous bamalignedreads.txt file to set the value in $frac for subsampling bam files


# Normalizing for the single end files 

```bash
#b1
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_b1.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=815950/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_b1.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_b1.norm.F.bam
done
#b2
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_b2.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=996100/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_b2.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_b2.norm.F.bam
done
#b6
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_b6.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=2093515/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_b6.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_b6.norm.F.bam
done
#b7
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_b7.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1100453/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_b7.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_b7.norm.F.bam
done
#b10
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_b10.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1342956/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_b10.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_b10.norm.F.bam
done

#g1
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_g1.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1050485/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_g1.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_g1.norm.F.bam
done
#g7
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_g7.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=857529/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_g7.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_g7.norm.F.bam
done
#g5
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_g5.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1169907/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_g5.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_g5.norm.F.bam
done
#g6
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_g6.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=933016/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_g6.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_g6.norm.F.bam
done
#g8
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_g8.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=934304/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_g8.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_g8.norm.F.bam
done

#k1
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_k1.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=870346/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_k1.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_k1.norm.F.bam
done
#k3
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_k3.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=916597/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_k3.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_k3.norm.F.bam
done
#k2
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_k2.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1617777/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_k2.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_k2.norm.F.bam
done
#k5
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_k5.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1176563/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_k5.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_k5.norm.F.bam
done
#k8
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_k8.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1970251/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_k8.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_k8.norm.F.bam
done

#m5
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_m5.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=936435/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_m5.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_m5.norm.F.bam
done
#m7
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_m7.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=807617/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_m7.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_m7.norm.F.bam
done
#m1
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_m1.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=901508/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_m1.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_m1.norm.F.bam
done
#m3
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_m3.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=898476/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_m3.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_m3.norm.F.bam
done
#m9
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_m9.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1002751/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_m9.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_m9.norm.F.bam
done


#n1
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_n1.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=771212/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_n1.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_n1.norm.F.bam
done
#n2
declare -a StringArray=("capture1" "capture4" "capture7" "capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_n2.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=936427/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_n2.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_n2.norm.F.bam
done
#n3
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_n3.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=894749/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_n3.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_n3.norm.F.bam
done
#n4
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_n4.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1025477/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_n4.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_n4.norm.F.bam
done
#n5
declare -a StringArray=("capture10" "capture13" "capture16")
for i in "${StringArray[@]}"
do
frac=$( samtools idxstats ${i}_n5.F.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=981831/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs $frac ${i}_n5.F.bam > /home/jgreen/eager_obj1b/03_mapping/01_merged/${i}_n5.norm.F.bam
done
```

# Analysis for 04_coverage_analysis /01_genome region run from 
```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")

#CDS
#Looping for all samples CDS
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.CDS.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/01_merged/${i}.AllCDS.all.split.txt
done
#Exon
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/01_merged/${i}.AllExon.all.split.txt
done
#Gene
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.gene.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/01_merged/${i}.AllGene.all.split.txt
done
#UTR
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.UTR.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/01_merged/${i}.AllUTR.all.split.txt
done
```

R code for calculating PercentBases > Depth for all Captures and Samples

```R all Capture B3 subsample
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)


s# Set your working directory
setwd("/home/jgreen/eager_obj1b/04_coverage_analysis/01_genome_region/01_merged/")

# Define the capture you want to isolate
capture <- "Capture1"

# Initialize an empty list to store ggplot objects
plot_list <- list()

# List of all samples
samples <- c("B3", "B4", "G3", "G5", "K3", "K4", "M3", "M4", "N2")

# Loop through each sample
for (sample in samples) {
  # Create the file name based on the sample and capture
  file <- paste0(capture, sample)
  
  # Print the current file being processed
  print(file)
  
  # Get a list of files for the current capture and sample
  files_list <- list.files(pattern = file)
  print(files_list)
  
  # Specify the files you want to process for each sample
  files <- c(
    paste0(file, ".AllCDS.all.split.txt"),
    paste0(file, ".AllExon.all.split.txt"),
    paste0(file, ".AllGene.all.split.txt"),
    paste0(file, ".AllUTR.all.split.txt")
  )
  print(files)
  
  # Labels for different data types
  labs <- c("CDS", "Exon", "Gene", "UTR")
  
  # Initialize a list to store data for each data type
  cov <- list()
  
  # Loop through each file and data type
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])[, c(2, 5)]
    cov_cumul <- 1 - cumsum(cov[[i]][, 2])
    cov[[i]]$cov_cumul <- c(1, cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample <- labs[i]
  }
  
  # Combine data for all data types
  cov_df <- do.call("rbind", cov)
  names(cov_df)[1:2] <- c("depth", "fraction")
  
  # Color palette
  pcbPalette <- c("#009E73", "#D55E00", "#CC79A7", "#56B4E9")
  
  # Create a ggplot
  p1 <- ggplot(cov_df, aes(x = depth, y = cov_cumul, color = sample, linetype = sample)) +
    xlim(0, 10) +
    scale_alpha(guide = 'none') +
    geom_line(size = 0.5) +
    scale_color_manual(values = pcbPalette) +
    scale_fill_manual(values = pcbPalette) +
    ggtitle(paste("Read Depth Distribution for", capture, "(", sample, ")")) +
    ylab("Percent of Bases > Depth") +
    xlab("Depth") +
    theme_bw() +
    theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5,face = "bold")) +
    theme(plot.title = element_text(size = 5, face = "bold", hjust = 0.5)) +
    theme(legend.position = "none")
  
  # Append the plot to the list
  plot_list[[sample]] <- p1
  
  # Save the individual plot as a PNG file
  png(
    filename = paste0("Figure2_", capture, "_", sample, ".png"),
    type = "cairo",
    units = "px",
    width = 5600,
    height = 3000,
    res = 600,
    bg = "transparent"
  )
  
  print(p1)
  dev.off()
}

# Combine all the ggplot objects in the list using grid.arrange()
pcombined <- grid.arrange(
  grobs = plot_list,
  ncol = 3,  # Adjust the number of columns as needed
  top = textGrob(paste("Read Depth Distribution for", capture, "across Samples"))
)

# Save the combined plot as a PNG file
ggsave(file = paste0("Figure2_", capture, "_Combined.png"), pcombined)
```

# Generate data for Table with mapped read stats

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")

for i in "${StringArray[@]}"
do 
nom=$(samtools view -@32 $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -c -L $WORKING_DIR/Genome/mtDNA.bed); denom=$(samtools view -@32 $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -c); dupPE=$(mawk '/Unknown/' $WORKING_DIR/02_ddocent/02_PE/logfiles/${i}_dup_metrics.txt | cut -f9); paste <(echo $i) <(echo $(( `zcat $WORKING_DIR/01_process_reads/raw/${i}_R1.fq.gz | wc -l` /4 ))) <(echo $(( `zcat $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz | wc -l` /4 ))) <(samtools view -@ 32 $WORKING_DIR/02_ddocent/02_PE/${i}-RGmd.bam -c) <(python -c "print(round("$dupPE" * 100,2))") <(echo $denom) <(python -c "print(round("$nom"/"$denom" *100,2))") 
done > data.table2
echo -e "Pool\tRaw_Reads\tFiltered_Reads\tMapped_Reads\tPercent_Duplicate\tFiltered_Mapped_Reads\tPercent_mapping_to_mitochondrial_genome" > header
cat header data.table2 > Table2.txt
```

## Generate data for Figure 3

Calculate Exon percentiles

Get total coverage counts for all merged

```bash
bedtools coverage -b ~/eager_obj1a/03_mapping/01_merged/all.filter.merged.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -mean -split > $WORKING_DIR/04_coverage_analysis/02_exon_stats/all.merged.cov.mean.exon.stats
```

Get total coverage counts per exon per capture
```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")

for i in "${StringArray[@]}"
do
bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/${i}.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -mean -split > $WORKING_DIR/04_coverage_analysis/02_exon_stats/01_merged/${i}.indiv.cov.mean.filtered.exon.stats
done
```

Remove mtDNA from stat files

```bash

#Merged files
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")

for i in  "${StringArray[@]}"
do
mawk '!/NC_007175.2/' ${i}.indiv.cov.mean.filtered.exon.stats > ${i}.indiv.DNA.mean.exon.stats
done

#All files
mawk '!/NC_007175.2/' all.merged.cov.mean.exon.stats > all.DNA.merged.mean.exon.stats
```

Calculate lower 10th percentile of exon sizes


```bash
for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}.indiv.DNA.mean.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
done
```

Result: 58 for all samples


Calculate upper 10th percentile of exon sizes


```bash
for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}.indiv.DNA.mean.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
done
```

Result: 524 for all samples


Mark exons into size classes based on size distribution and create data table
```bash
for i in  "${StringArray[@]}"
do
mawk '{if ( $3 -$2 > 524 ) print $0 "\tUpper"; else if ( $3 - $2 < 58 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 524) print $0 "\tMiddle" }' ${i}.indiv.DNA.mean.exon.stats > ${i}.indiv.DNA.mean.cov.exon.stats.class
done

echo -e "Chrom\tStart\tEnd\tDNA_Coverage\tExon_Size_Class" > header

for i in  "${StringArray[@]}"
do
cat header ${i}.indiv.DNA.mean.cov.exon.stats.class > ${i}.indiv.ExonMeanCoverage.txt
done
```

```R Set Libraries and Working Directory
library(MASS)
library(fields)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)
```

```R
{r All density analysis combined with pairwise T-test and ANOVA}
setwd("/home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged/")
# Define a list of capture names and sample names
# Define a list of capture names and sample names
captures <- c("Capture1", "Capture2", "Capture3", "Capture4")
samples <- c("B3", "B4", "G3", "G5", "K3", "K4", "M3", "M4", "N2")

# Initialize an empty list to store data frames
data_frames <- list()

# Directory to store result files
output_dir <- "/home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged/stat_results/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through captures and samples to read data files
for (capture in captures) {
  for (sample in samples) {
    file_name <- paste(capture, "_", sample, ".insert.merged.ExonMeanCoverage.txt", sep = "")
    file_path <- file.path("/home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged/", file_name)
    df <- read.table(file_path, header = TRUE)
    df <- as.data.frame(df)
    # Add capture and sample information as columns
    df$Capture <- capture
    df$Sample <- sample
    data_frames[[paste(capture, sample, sep = "_")]] <- df
  }
}

# Perform pairwise comparisons between captures for each sample
for (sample in samples) {
  # Create an empty list to store comparison results
  comparison_results <- list()
  
  for (capture1 in captures) {
    for (capture2 in captures) {
      if (capture1 != capture2) {
        # Get the data frames for the two captures and the same sample
        df1 <- data_frames[[paste(capture1, sample, sep = "_")]]
        df2 <- data_frames[[paste(capture2, sample, sep = "_")]]
        
        # Perform a t-test for the DNA_Coverage values between captures
        t_test_result <- t.test(df1$DNA_Coverage, df2$DNA_Coverage)
        
        # Store the comparison result in a list
        comparison_results[[paste(capture1, "_vs_", capture2, sep = "")]] <- t_test_result
      }
    }
  }
  
  # Create a file to store t-test results for this sample
  t_test_output_file <- file.path(output_dir, paste("t_test_results_", sample, ".txt", sep = ""))
  
  # Print or further process the comparison results as needed
  # For example, print p-values for each pairwise comparison
  cat("Sample:", sample, "\n")
  for (capture_comparison in names(comparison_results)) {
    p_value <- comparison_results[[capture_comparison]]$p.value
    cat(capture_comparison, "t-test p-value:", p_value, "\n")
    cat(capture_comparison, "t-test p-value:", p_value, "\n", file = t_test_output_file, append = TRUE)
  }
  
  # Perform a one-way ANOVA test for the same sample across captures
  sample_data <- lapply(captures, function(capture) {
    data_frames[[paste(capture, sample, sep = "_")]]$DNA_Coverage
  })
  anova_result <- aov(DNA_Coverage ~ as.factor(Capture), data = data.frame(Capture = captures, DNA_Coverage = unlist(sample_data)))
  anova_p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
  
  # Create a file to store ANOVA results for this sample
  anova_output_file <- file.path(output_dir, paste("anova_results_", sample, ".txt", sep = ""))
  
  # Open the file for writing
  cat("ANOVA p-value for", sample, "across captures:", anova_p_value, "\n\n", file = anova_output_file)
  cat("ANOVA p-value for", sample, "across captures:", anova_p_value, "\n\n")
}
```

```R
{r Make plotting function}
perform_comparison_and_plot <- function(sample, capture1, capture2) {
  # Get the data frames for the two captures and the same sample
  df1 <- data_frames[[paste(capture1, sample, sep = "_")]]
  df2 <- data_frames[[paste(capture2, sample, sep = "_")]]
  
  # Merge the data frames on "Start" and select relevant columns
  merged_df <- merge(df1, df2[, c("Start", "DNA_Coverage", "Exon_Size_Class")], by = "Start")
  
  # Filter rows where both DNA_Coverage values are not zero
  TotalExon <- merged_df[merged_df$DNA_Coverage.x != 0 & merged_df$DNA_Coverage.y != 0, ]
  TotalExon <- TotalExon[, -5]  # Remove unnecessary column
  
  # Modify Exon_Size_Class as needed
  TotalExon$Exon_Size_Class <- factor(TotalExon$Exon_Size_Class.y, levels = c("Lower", "Middle", "Upper"))
  TotalExon$Exon_Size_Class <- revalue(TotalExon$Exon_Size_Class.y, c("Lower" = "Lower 10Percent", "Upper" = "Upper 10Percent", "Middle" = "Middle 80Percent"))
  
  # Calculate density
  get_density <- function(x, y, n = 100) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  TotalExon$density <- get_density(TotalExon$DNA_Coverage.x, TotalExon$DNA_Coverage.y)
  
  # Set color palette
  cbPalette <- c("#009E73", "#D55E00", "#56B4E9", "#0072B2", "#E69F00", "#F0E442", "#999999", "#CC79A7")
  t <- cbPalette[1]
  cbPalette[1] <- cbPalette[2]
  cbPalette[2] <- t
  
  # Scale DNA_Coverage values
  TotalExon$DNA_Coverage.y <- TotalExon$DNA_Coverage.y / 6
  TotalExon$DNA_Coverage.x <- TotalExon$DNA_Coverage.x / 6
  
  # Create the plot
  b1 <- ggplot(TotalExon, aes(x = DNA_Coverage.x + 1, y = DNA_Coverage.y + 1, alpha = 1 / density),
               fill = Exon_Size_Class, color = Exon_Size_Class) +
    geom_point(aes(color = TotalExon$Exon_Size_Class, fill = TotalExon$Exon_Size_Class, shape = TotalExon$Exon_Size_Class)) +
    geom_smooth(method = "auto", alpha = 0.5, size = 0, se = TRUE) +
    geom_abline(linetype = "dashed", color = "red") +
    stat_smooth(geom = "line", alpha = 0.75, size = 0.5, linetype = "dashed") +
    scale_alpha_continuous(guide = "none", range = c(.2, .95)) +
    scale_shape_manual(values = c(15, 16, 17), name = "Exon Size Percentile") +
    scale_fill_manual(values = cbPalette, name = "Exon Size Percentile") +
    scale_color_manual(values = cbPalette, name = "Exon Size Percentile") +
    xlab(paste("Mean", capture1, "DNA Reads per Exon Base Pair")) +
    ylab(paste("Mean", capture2, "DNA Reads per Exon Base Pair")) +
    theme_bw() +
    theme(legend.position = c(0.9, 0.15))
  
  # Define the output file name
  output_file <- file.path(output_dir, paste("Figure_indiv_", sample, "_", capture1, "_vs_", capture2, ".png", sep = ""))
  
  # Save the plot as an image file
  png(filename = output_file, type = "cairo", units = "px", width = 5600, 
      height = 3000, res = 600, bg = "transparent")
  print(b1)
  dev.off()
}

# Loop through captures and samples to perform comparisons and
```

```R
{r Run perform comparison function on all samples}
# Define a list of capture names and sample names
captures <- c("Capture1", "Capture2", "Capture3", "Capture4")
samples <- c("B3", "B4", "G3", "G5", "K3", "K4", "M3", "M4", "N2")

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through captures and samples to read data files
for (capture in captures) {
  for (sample in samples) {
    file_name <- paste(capture, "_", sample, ".insert.merged.ExonMeanCoverage.txt", sep = "")
    file_path <- file.path("/home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged/", file_name)
    df <- read.table(file_path, header = TRUE)
    df <- as.data.frame(df)
    # Add capture and sample information as columns
    df$Capture <- capture
    df$Sample <- sample
    data_frames[[paste(capture, sample, sep = "_")]] <- df
  }
}

# Directory to store result files and plots
output_dir <- "/home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged/results/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through captures and samples to perform comparisons and generate plots
for (sample in samples) {
  for (capture1 in captures) {
    for (capture2 in captures) {
      if (capture1 != capture2) {
        # Call the function to perform comparison and generate plot
        perform_comparison_and_plot(sample, capture1, capture2)
      }
    }
  }
}
```

Next step is to find exons with minimum thresholds of gDNA.  These will be our "target" sets along with confidence intervals.  Based on overall DNA coverage, we chose 35X as our "target" set and choose 2X boundaries around that number for confidence intervals.  We will create multiple `bed` files from our exon coverage stats.

```bash
#Run in /home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
mawk 'BEGIN { FS = "\t" } ; $4 > 1' ${i}.indiv.cov.mean.filtered.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/${i}.indiv.EiRc2.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 6' ${i}.indiv.cov.mean.filtered.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/${i}.indiv.EiRc7.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 11' ${i}.indiv.cov.mean.filtered.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/${i}.indiv.EiRc12.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 17' ${i}.indiv.cov.mean.filtered.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/${i}.indiv.EiRc18.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 24' ${i}.indiv.cov.mean.filtered.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/${i}.indiv.EiRc25.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 34' ${i}.indiv.cov.mean.filtered.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/${i}.indiv.EiRc35.bed
done

#Run in /home/jgreen/eager_obj1b/04_coverage_analysis/02_exon_stats/01_merged
mawk 'BEGIN { FS = "\t" } ; $4 > 1' all.DNA.merged.mean.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/all.filter.merged.EiRc2.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 6' all.DNA.merged.mean.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/all.filter.merged.EiRc7.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 11' all.DNA.merged.mean.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/all.filter.merged.EiRc12.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 17' all.DNA.merged.mean.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/all.filter.merged.EiRc18.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 24' all.DNA.merged.mean.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/all.filter.merged.EiRc25.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 34' all.DNA.merged.mean.exon.stats > $WORKING_DIR/04_coverage_analysis/03_target_interval/01_merged/all.filter.merged.EiRc35.bed
```

### Calculating data for table 3

Next step is to find exons with minimum thresholds of DNA coverage. We will use a BASH function to automate this for us:

# This function uses the individual bed files to compare coverage between the Bam file and Bed file

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 2X targets with more than 1X coverage
X2XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $1.indiv.EiRc2.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 7X targets with more than 1X coverage
X7XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $1.indiv.EiRc7.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 12X targets with more than 1X coverage
X12XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $1.indiv.EiRc12.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 18X targets with more than 1X coverage
X18XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $1.indiv.EiRc18.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 25X targets with more than 1X coverage
X25XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $1.indiv.EiRc25.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $1.indiv.EiRc35.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed | wc -l )
X2X=$(cat all.filter.merged.EiRc2.bed | wc -l)
X7X=$(cat all.filter.merged.EiRc7.bed | wc -l)
X12X=$(cat all.filter.merged.EiRc12.bed | wc -l)
X18X=$(cat all.filter.merged.EiRc18.bed | wc -l)
X25X=$(cat all.filter.merged.EiRc25.bed | wc -l)
X35X=$(cat all.filter.merged.EiRc35.bed | wc -l)


#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"Percent"
echo `python -c "print(round("$X2XC"/"$X2X" * 100,1))"`"Percent"
echo `python -c "print(round("$X7XC"/"$X7X" * 100,1))"`"Percent"
echo `python -c "print(round("$X12XC"/"$X12X" * 100,1))"`"Percent"
echo `python -c "print(round("$X18XC"/"$X18X" * 100,1))"`"Percent"
echo `python -c "print(round("$X25XC"/"$X25X" * 100,1))"`"Percent"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"Percent"
  
}

export -f counts_per_target
```

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
echo "<(counts_per_target ${i})"
done
```
Now with this function we can use `paste` and subshells to produce the table
```bash

#Initial test with one sample
paste <(echo -e "Targets\nAll Exons\n2XR Exons\n7XR Exons\n12XR Exons\n18XR Exons\n25XR Exons\n35XR Exons") <(counts_per_target capture1_b1) > Table3_indiv_test.txt

#All samples
paste <(echo -e "Targets\nAll Exons\n2XR Exons\n7XR Exons\n12XR Exons\n18XR Exons\n25XR Exons\n35XR Exons") <(counts_per_target capture1_b1) <(counts_per_target capture1_b2) <(counts_per_target capture1_g1) <(counts_per_target capture1_g7) <(counts_per_target capture1_k1) <(counts_per_target capture1_k3) <(counts_per_target capture1_m5) <(counts_per_target capture1_m7) <(counts_per_target capture1_n1) <(counts_per_target capture1_n2) <(counts_per_target capture4_b1) <(counts_per_target capture4_b2) <(counts_per_target capture4_g1) <(counts_per_target capture4_g7) <(counts_per_target capture4_k1) <(counts_per_target capture4_k3) <(counts_per_target capture4_m5) <(counts_per_target capture4_m7) <(counts_per_target capture4_n1) <(counts_per_target capture4_n2) <(counts_per_target capture7_b1) <(counts_per_target capture7_b2) <(counts_per_target capture7_g1) <(counts_per_target capture7_g7) <(counts_per_target capture7_k1) <(counts_per_target capture7_k3) <(counts_per_target capture7_m5) <(counts_per_target capture7_m7) <(counts_per_target capture7_n1) <(counts_per_target capture7_n2) <(counts_per_target capture10_b10) <(counts_per_target capture10_b1) <(counts_per_target capture10_b2) <(counts_per_target capture10_b6) <(counts_per_target capture10_b7) <(counts_per_target capture10_g1) <(counts_per_target capture10_g5) <(counts_per_target capture10_g6) <(counts_per_target capture10_g7) <(counts_per_target capture10_g8) <(counts_per_target capture10_k1) <(counts_per_target capture10_k2) <(counts_per_target capture10_k3) <(counts_per_target capture10_k5) <(counts_per_target capture10_k8) <(counts_per_target capture10_m1) <(counts_per_target capture10_m3) <(counts_per_target capture10_m5) <(counts_per_target capture10_m7) <(counts_per_target capture10_m9) <(counts_per_target capture10_n1) <(counts_per_target capture10_n2) <(counts_per_target capture10_n3) <(counts_per_target capture10_n4) <(counts_per_target capture10_n5) <(counts_per_target capture13_b10) <(counts_per_target capture13_b1) <(counts_per_target capture13_b2) <(counts_per_target capture13_b6) <(counts_per_target capture13_b7) <(counts_per_target capture13_g1) <(counts_per_target capture13_g5) <(counts_per_target capture13_g6) <(counts_per_target capture13_g7) <(counts_per_target capture13_g8) <(counts_per_target capture13_k1) <(counts_per_target capture13_k2) <(counts_per_target capture13_k3) <(counts_per_target capture13_k5) <(counts_per_target capture13_k8) <(counts_per_target capture13_m1) <(counts_per_target capture13_m3) <(counts_per_target capture13_m5) <(counts_per_target capture13_m7) <(counts_per_target capture13_m9) <(counts_per_target capture13_n1) <(counts_per_target capture13_n2) <(counts_per_target capture13_n3) <(counts_per_target capture13_n4) <(counts_per_target capture13_n5) <(counts_per_target capture16_b10) <(counts_per_target capture16_b1) <(counts_per_target capture16_b2) <(counts_per_target capture16_b6) <(counts_per_target capture16_b7) <(counts_per_target capture16_g1) <(counts_per_target capture16_g5) <(counts_per_target capture16_g6) <(counts_per_target capture16_g7) <(counts_per_target capture16_g8) <(counts_per_target capture16_k1) <(counts_per_target capture16_k2) <(counts_per_target capture16_k3) <(counts_per_target capture16_k5) <(counts_per_target capture16_k8) <(counts_per_target capture16_m1) <(counts_per_target capture16_m3) <(counts_per_target capture16_m5) <(counts_per_target capture16_m7) <(counts_per_target capture16_m9) <(counts_per_target capture16_n1) <(counts_per_target capture16_n2) <(counts_per_target capture16_n3) <(counts_per_target capture16_n4) <(counts_per_target capture16_n5) > Table3_indiv_bed.txt


# This function is differen tas it uses the merged bed files from all captures to compare coverage between the Bam file and Bed file

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 2X targets with more than 1X coverage
X2XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a all.filter.merged.EiRc2.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 1' | wc -l) 
#Calculate number of 7X targets with more than 1X coverage
X7XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a all.filter.merged.EiRc7.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 1' | wc -l) 
#Calculate number of 12X targets with more than 1X coverage
X12XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a all.filter.merged.EiRc12.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 1' | wc -l) 
#Calculate number of 18X targets with more than 1X coverage
X18XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a all.filter.merged.EiRc18.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 1' | wc -l) 
#Calculate number of 25X targets with more than 1X coverage
X25XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a all.filter.merged.EiRc25.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 1' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/01_merged/$1.norm.F.bam -a all.filter.merged.EiRc35.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 1' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed | wc -l )
X2X=$(cat all.filter.merged.EiRc2.bed | wc -l)
X7X=$(cat all.filter.merged.EiRc7.bed | wc -l)
X12X=$(cat all.filter.merged.EiRc12.bed | wc -l)
X18X=$(cat all.filter.merged.EiRc18.bed | wc -l)
X25X=$(cat all.filter.merged.EiRc25.bed | wc -l)
X35X=$(cat all.filter.merged.EiRc35.bed | wc -l)


#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"Percent"
echo `python -c "print(round("$X2XC"/"$X2X" * 100,1))"`"Percent"
echo `python -c "print(round("$X7XC"/"$X7X" * 100,1))"`"Percent"
echo `python -c "print(round("$X12XC"/"$X12X" * 100,1))"`"Percent"
echo `python -c "print(round("$X18XC"/"$X18X" * 100,1))"`"Percent"
echo `python -c "print(round("$X25XC"/"$X25X" * 100,1))"`"Percent"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"Percent"
  
}

export -f counts_per_target
```

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
echo "<(counts_per_target ${i})"
done
```
Now with this function we can use `paste` and subshells to produce the table
```bash

#Initial test with one sample
paste <(echo -e "Targets\nAll Exons\n2XR Exons\n7XR Exons\n12XR Exons\n18XR Exons\n25XR Exons\n35XR Exons") <(counts_per_target Capture1_B3) > Table3_indiv_test.txt

#All samples
paste <(echo -e "Targets\nAll Exons\n2XR Exons\n7XR Exons\n12XR Exons\n18XR Exons\n25XR Exons\n35XR Exons") <(counts_per_target capture1_b1) <(counts_per_target capture1_b2) <(counts_per_target capture1_g1) <(counts_per_target capture1_g7) <(counts_per_target capture1_k1) <(counts_per_target capture1_k3) <(counts_per_target capture1_m5) <(counts_per_target capture1_m7) <(counts_per_target capture1_n1) <(counts_per_target capture1_n2) <(counts_per_target capture4_b1) <(counts_per_target capture4_b2) <(counts_per_target capture4_g1) <(counts_per_target capture4_g7) <(counts_per_target capture4_k1) <(counts_per_target capture4_k3) <(counts_per_target capture4_m5) <(counts_per_target capture4_m7) <(counts_per_target capture4_n1) <(counts_per_target capture4_n2) <(counts_per_target capture7_b1) <(counts_per_target capture7_b2) <(counts_per_target capture7_g1) <(counts_per_target capture7_g7) <(counts_per_target capture7_k1) <(counts_per_target capture7_k3) <(counts_per_target capture7_m5) <(counts_per_target capture7_m7) <(counts_per_target capture7_n1) <(counts_per_target capture7_n2) <(counts_per_target capture10_b10) <(counts_per_target capture10_b1) <(counts_per_target capture10_b2) <(counts_per_target capture10_b6) <(counts_per_target capture10_b7) <(counts_per_target capture10_g1) <(counts_per_target capture10_g5) <(counts_per_target capture10_g6) <(counts_per_target capture10_g7) <(counts_per_target capture10_g8) <(counts_per_target capture10_k1) <(counts_per_target capture10_k2) <(counts_per_target capture10_k3) <(counts_per_target capture10_k5) <(counts_per_target capture10_k8) <(counts_per_target capture10_m1) <(counts_per_target capture10_m3) <(counts_per_target capture10_m5) <(counts_per_target capture10_m7) <(counts_per_target capture10_m9) <(counts_per_target capture10_n1) <(counts_per_target capture10_n2) <(counts_per_target capture10_n3) <(counts_per_target capture10_n4) <(counts_per_target capture10_n5) <(counts_per_target capture13_b10) <(counts_per_target capture13_b1) <(counts_per_target capture13_b2) <(counts_per_target capture13_b6) <(counts_per_target capture13_b7) <(counts_per_target capture13_g1) <(counts_per_target capture13_g5) <(counts_per_target capture13_g6) <(counts_per_target capture13_g7) <(counts_per_target capture13_g8) <(counts_per_target capture13_k1) <(counts_per_target capture13_k2) <(counts_per_target capture13_k3) <(counts_per_target capture13_k5) <(counts_per_target capture13_k8) <(counts_per_target capture13_m1) <(counts_per_target capture13_m3) <(counts_per_target capture13_m5) <(counts_per_target capture13_m7) <(counts_per_target capture13_m9) <(counts_per_target capture13_n1) <(counts_per_target capture13_n2) <(counts_per_target capture13_n3) <(counts_per_target capture13_n4) <(counts_per_target capture13_n5) <(counts_per_target capture16_b10) <(counts_per_target capture16_b1) <(counts_per_target capture16_b2) <(counts_per_target capture16_b6) <(counts_per_target capture16_b7) <(counts_per_target capture16_g1) <(counts_per_target capture16_g5) <(counts_per_target capture16_g6) <(counts_per_target capture16_g7) <(counts_per_target capture16_g8) <(counts_per_target capture16_k1) <(counts_per_target capture16_k2) <(counts_per_target capture16_k3) <(counts_per_target capture16_k5) <(counts_per_target capture16_k8) <(counts_per_target capture16_m1) <(counts_per_target capture16_m3) <(counts_per_target capture16_m5) <(counts_per_target capture16_m7) <(counts_per_target capture16_m9) <(counts_per_target capture16_n1) <(counts_per_target capture16_n2) <(counts_per_target capture16_n3) <(counts_per_target capture16_n4) <(counts_per_target capture16_n5) > Table3_indiv_merged.txt
```

## Generate data for figure 4

Figure 4 is per bp sensitivity looking at coverage across our target sets, near targets (definied as 150 bp around the edge of targets, and off target (everything that is not near or on target).

First steps involve creating our different interval sets using bedtools.

Lets link all the norm.F.bam files here

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
ln -s /home/jgreen/eager_obj1b/03_mapping/01_merged/$i.norm.F.bam .
done
```


```bash
bedtools flank -i all.filter.merged.EiRc2.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc2.150.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc2.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc2.150.slop.bed 
bedtools complement -i all.filter.merged.EiRc2.150.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc2.150.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc7.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc7.150.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc7.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc7.150.slop.bed 
bedtools complement -i all.filter.merged.EiRc7.150.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc7.150.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc12.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc12.150.neartarget.bed
bedtools slop -i all.filter.merged.EiRc12.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc12.150.slop.bed
bedtools complement -i all.filter.merged.EiRc12.150.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc12.150.offtarget.bed

bedtools flank -i all.filter.merged.EiRc2.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc2.300.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc2.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc2.300.slop.bed 
bedtools complement -i all.filter.merged.EiRc2.300.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc2.300.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc7.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc7.300.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc7.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc7.300.slop.bed 
bedtools complement -i all.filter.merged.EiRc7.300.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc7.300.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc12.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc12.300.neartarget.bed
bedtools slop -i all.filter.merged.EiRc12.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc12.300.slop.bed
bedtools complement -i all.filter.merged.EiRc12.300.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc12.300.offtarget.bed
```

With the target sets defined we again use bedtools to calculate coverage levels across these various genomic regions, and below we use GNU-parallel to speed things up.

```bash
#test with one file
parallel 'bedtools coverage -hist -b /home/jgreen/eager_obj1b/03_mapping/01_merged/capture1_b1.norm.F.bam -a all.filter.merged.EiRc2.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > capture1_b1.hist.EiRc2.150.neartarget.all.txt'

#Fixed ls piping and running on all samples
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc2.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc2.150.neartarget.all.txt'
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc2.150.offtarget.bed  -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc2.150.offtarget.all.txt' 
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc2.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc2.150.all.txt' 

ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc7.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc7.150.neartarget.all.txt' 
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc7.150.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc7.150.offtarget.all.txt' 


ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc12.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc12.150.neartarget.all.txt' 
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc12.150.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc12.150.offtarget.all.txt' 
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc12.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc12.150.all.txt' 
```

Do the same for the 300 bp segments

```bash
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc2.300.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc2.300.neartarget.all.txt'
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc2.300.offtarget.bed  -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc2.300.offtarget.all.txt'
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc2.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc2.300.all.txt'

ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc7.300.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc7.300.neartarget.all.txt'
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc7.300.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc7.300.offtarget.all.txt'
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc7.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc7.300.all.txt'

ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc12.300.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc12.300.neartarget.all.txt' 
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc12.300.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc12.300.offtarget.all.txt' 
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel 'bedtools coverage -hist -b {}.norm.F.bam -a all.filter.merged.EiRc12.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc12.300.all.txt'
```

```R
{r Set functions for figure4 150 bp intervals}
setwd("/home/jgreen/eager_obj1b/04_coverage_analysis/03_target_interval/01_merged/")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

make_graph <- function(j){
  print(files <- list.files(pattern=paste(j,".hist.*EiRc2.150*", sep = "")))
  labs <- c("On target","Near target", "Off target")
  
  cov <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul=1-cumsum(cov[[i]][,2])
    cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample=labs[i]
  }
  
  cov_df=do.call("rbind",cov)
  names(cov_df)[1:2]=c("depth","fraction")
  
  print(files <- list.files(pattern=paste(j,".hist.*EiRc7.150*", sep = "")))
  cov2 <- list()
  for (i in 1:length(files)) {
    cov2[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul2=1-cumsum(cov2[[i]][,2])
    cov2[[i]]$cov_cumul <- c(1,cov_cumul2[-length(cov_cumul2)])
    cov2[[i]]$sample=labs[i]
  }
  cov2_df=do.call("rbind",cov2)
  names(cov2_df)[1:2]=c("depth","fraction")
  
  print(files <- list.files(pattern=paste(j,".hist.*EiRc12.150*", sep = "")))
  cov3 <- list()
  for (i in 1:length(files)) {
    cov3[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul3=1-cumsum(cov3[[i]][,2])
    cov3[[i]]$cov_cumul <- c(1,cov_cumul3[-length(cov_cumul3)])
    cov3[[i]]$sample=labs[i]
  }
  cov3_df=do.call("rbind",cov3)
  names(cov3_df)[1:2]=c("depth","fraction")
  
  cov_df <- subset(cov_df, depth <51)
  cov2_df <- subset(cov2_df, depth <51)
  cov3_df <- subset(cov3_df, depth <51)
  
  cov2_df$high <-cov3_df$cov_cumul
  
  cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7")
  cbPalettep <- c("#0072B2", "#009E73","#D55E00" )
  
  cov_df$sample <- factor(cov_df$sample,levels=c("On target","Near target", "Off target"))
  cov2_df$sample <- factor(cov2_df$sample,levels=c("On target","Near target", "Off target"))
  
  p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  +
    xlim(0,8)+
    geom_ribbon(data=cov2_df,aes(ymin=cov_cumul,ymax=high, color=sample, fill=sample, alpha=0.4)) +
    scale_alpha(guide = 'none') +
    geom_line(size=1.5)+ 
    scale_color_manual(values=cbPalettep) +
    scale_fill_manual(values=cbPalettep) +
    ylab("Percent of Target Bases > Depth")+
    #scale_x_log10()+
    xlab("Depth")+
    ggtitle(eval(j)) +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.position="none")
  #theme(legend.position=c(0.92,0.88))
  
  return(p)
}
```

```R
{r Apply function to all 150 bp samples, warning=FALSE}
sample_names <- c("B3", "B4", "G3", "G5", "K3", "K4", "M3", "M4", "N2")

# Set the working directory
setwd("/home/jgreen/eager_obj1b/04_coverage_analysis/03_target_interval/01_merged/")

# Create an output directory if it doesn't exist
output_directory <- "image_results"
if (!file.exists(output_directory)) {
  dir.create(output_directory)
}

# Create a function to generate and save the comparison plot for a single sample
generate_and_save_sample_plot <- function(sample_name) {
  # Initialize an empty list to store the plots for each capture
  capture_plots <- list()
  
  # Generate the plots for each capture
  for (capture in c("Capture1", "Capture2", "Capture3", "Capture4")) {
    current_plot <- make_graph(paste(capture, sample_name, sep = "_"))
    capture_plots[[capture]] <- current_plot
  }
  
  # Set the PNG file path for saving the plot
  png_file <- file.path(output_directory, paste("target_interval_150",sample_name, ".png", sep = "_"))
  
  # Save the combined plot for the sample in a PNG file
  png(file =png_file, width = 1400, height = 650, bg = "transparent")  # 2x2 layout
  multiplot(plotlist = capture_plots, cols = 2)
  dev.off()
}

# Generate and save the combined plot for each sample
for (sample_name in sample_names) {
  generate_and_save_sample_plot(sample_name)
}
```

## Calculating Specificity

First let's calculate near and off-target intervals for all exons

```bash
# link hmash.bed files in Genome directory to /home/jgreen/eager_obj1b/04_coverage_analysis/04_specificity/01_merged directory
bedtools flank -i sorted.ref3.0.exon.sc.hmask.bed -b 150 -g masked.genome.file | bedtools sort -faidx masked.genome.file >  sorted.ref3.0.exon.sc.hmask.neartarget.bed
bedtools slop -i sorted.ref3.0.exon.sc.hmask.bed -b 150 -g masked.genome.file > sorted.ref3.0.exon.sc.hmask.slop.bed
bedtools complement -i sorted.ref3.0.exon.sc.hmask.slop.bed -g masked.genome.file > sorted.ref3.0.exon.sc.hmask.offtarget.bed
```

# Specificty table for the 150 bp segments

Lets link all the norm.F.bam files here

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
ln -s /home/jgreen/eager_obj1b/03_mapping/01_merged/$i.norm.F.bam .
done
```

Now we can create a specificity table for all exons and for expressed targets at 150 bp flanking regions using a few more BASH functions

```bash
specExon(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L sorted.ref3.0.exon.sc.hmask.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L sorted.ref3.0.exon.sc.hmask.neartarget.bed)
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L sorted.ref3.0.exon.sc.hmask.bed | samtools view - -@10 -c -L sorted.ref3.0.exon.sc.hmask.neartarget.bed)
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L sorted.ref3.0.exon.sc.hmask.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L sorted.ref3.0.exon.sc.hmask.slop.bed | samtools view - -@10 -c -L sorted.ref3.0.exon.sc.hmask.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e $1"\t"`python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f specExon
```

```bash
spec2X150(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc2.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc2.150.neartarget.bed )
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc2.bed | samtools view - -@10 -c -L all.filter.merged.EiRc2.150.neartarget.bed )
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc2.150.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc2.150.neartarget.bed  | samtools view - -@10 -c -L all.filter.merged.EiRc2.150.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f spec2X150
```

```bash
spec7X150(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc7.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc7.150.neartarget.bed )
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc7.bed | samtools view - -@10 -c -L all.filter.merged.EiRc7.150.neartarget.bed )
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc7.150.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc7.150.neartarget.bed  | samtools view - -@10 -c -L all.filter.merged.EiRc7.150.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f spec7X150
```

```bash
spec12X150(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc12.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc12.150.neartarget.bed )
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc12.bed | samtools view - -@10 -c -L all.filter.merged.EiRc12.150.neartarget.bed )
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc12.150.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc12.150.neartarget.bed  | samtools view - -@10 -c -L all.filter.merged.EiRc12.150.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f spec12X150
```

Now we use all the functions to create a table
```bash
echo -e "Pool\tPercent_in_Exons\tPercent_Near_Exons\tPercent_Off_Target_Exons\tPercent_on_Target2X\tPErcent_Near_Target2X\tPercentOff_Target2X\tPercent_on_Target7X\tPercent_Near_Target7X\tPercentOff_Target7X\tPercent_on_Target12X\tPercent_Near_Target12X\tPercentOff_Target12X" > Spec.150.Table

# Only use once to generate the list of commands below
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
paste <(specExon ${i}) <(spec2X150 ${i}) <(spec7X150 ${i}) <(spec12X150 ${i}) >> Spec.150.Table
done

########################################################### Run Commands Below #################################################

paste <(specExon capture1_b1) <(spec2X150 capture1_b1) <(spec7X150 capture1_b1) <(spec12X150 capture1_b1) >> Spec.150.Table
paste <(specExon capture1_b2) <(spec2X150 capture1_b2) <(spec7X150 capture1_b2) <(spec12X150 capture1_b2) >> Spec.150.Table
paste <(specExon capture1_g1) <(spec2X150 capture1_g1) <(spec7X150 capture1_g1) <(spec12X150 capture1_g1) >> Spec.150.Table
paste <(specExon capture1_g7) <(spec2X150 capture1_g7) <(spec7X150 capture1_g7) <(spec12X150 capture1_g7) >> Spec.150.Table
paste <(specExon capture1_k1) <(spec2X150 capture1_k1) <(spec7X150 capture1_k1) <(spec12X150 capture1_k1) >> Spec.150.Table
paste <(specExon capture1_k3) <(spec2X150 capture1_k3) <(spec7X150 capture1_k3) <(spec12X150 capture1_k3) >> Spec.150.Table
paste <(specExon capture1_m5) <(spec2X150 capture1_m5) <(spec7X150 capture1_m5) <(spec12X150 capture1_m5) >> Spec.150.Table
paste <(specExon capture1_m7) <(spec2X150 capture1_m7) <(spec7X150 capture1_m7) <(spec12X150 capture1_m7) >> Spec.150.Table
paste <(specExon capture1_n1) <(spec2X150 capture1_n1) <(spec7X150 capture1_n1) <(spec12X150 capture1_n1) >> Spec.150.Table
paste <(specExon capture1_n2) <(spec2X150 capture1_n2) <(spec7X150 capture1_n2) <(spec12X150 capture1_n2) >> Spec.150.Table
paste <(specExon capture4_b1) <(spec2X150 capture4_b1) <(spec7X150 capture4_b1) <(spec12X150 capture4_b1) >> Spec.150.Table
paste <(specExon capture4_b2) <(spec2X150 capture4_b2) <(spec7X150 capture4_b2) <(spec12X150 capture4_b2) >> Spec.150.Table
paste <(specExon capture4_g1) <(spec2X150 capture4_g1) <(spec7X150 capture4_g1) <(spec12X150 capture4_g1) >> Spec.150.Table
paste <(specExon capture4_g7) <(spec2X150 capture4_g7) <(spec7X150 capture4_g7) <(spec12X150 capture4_g7) >> Spec.150.Table
paste <(specExon capture4_k1) <(spec2X150 capture4_k1) <(spec7X150 capture4_k1) <(spec12X150 capture4_k1) >> Spec.150.Table
paste <(specExon capture4_k3) <(spec2X150 capture4_k3) <(spec7X150 capture4_k3) <(spec12X150 capture4_k3) >> Spec.150.Table
paste <(specExon capture4_m5) <(spec2X150 capture4_m5) <(spec7X150 capture4_m5) <(spec12X150 capture4_m5) >> Spec.150.Table
paste <(specExon capture4_m7) <(spec2X150 capture4_m7) <(spec7X150 capture4_m7) <(spec12X150 capture4_m7) >> Spec.150.Table
paste <(specExon capture4_n1) <(spec2X150 capture4_n1) <(spec7X150 capture4_n1) <(spec12X150 capture4_n1) >> Spec.150.Table
paste <(specExon capture4_n2) <(spec2X150 capture4_n2) <(spec7X150 capture4_n2) <(spec12X150 capture4_n2) >> Spec.150.Table
paste <(specExon capture7_b1) <(spec2X150 capture7_b1) <(spec7X150 capture7_b1) <(spec12X150 capture7_b1) >> Spec.150.Table
paste <(specExon capture7_b2) <(spec2X150 capture7_b2) <(spec7X150 capture7_b2) <(spec12X150 capture7_b2) >> Spec.150.Table
paste <(specExon capture7_g1) <(spec2X150 capture7_g1) <(spec7X150 capture7_g1) <(spec12X150 capture7_g1) >> Spec.150.Table
paste <(specExon capture7_g7) <(spec2X150 capture7_g7) <(spec7X150 capture7_g7) <(spec12X150 capture7_g7) >> Spec.150.Table
paste <(specExon capture7_k1) <(spec2X150 capture7_k1) <(spec7X150 capture7_k1) <(spec12X150 capture7_k1) >> Spec.150.Table
paste <(specExon capture7_k3) <(spec2X150 capture7_k3) <(spec7X150 capture7_k3) <(spec12X150 capture7_k3) >> Spec.150.Table
paste <(specExon capture7_m5) <(spec2X150 capture7_m5) <(spec7X150 capture7_m5) <(spec12X150 capture7_m5) >> Spec.150.Table
paste <(specExon capture7_m7) <(spec2X150 capture7_m7) <(spec7X150 capture7_m7) <(spec12X150 capture7_m7) >> Spec.150.Table
paste <(specExon capture7_n1) <(spec2X150 capture7_n1) <(spec7X150 capture7_n1) <(spec12X150 capture7_n1) >> Spec.150.Table
paste <(specExon capture7_n2) <(spec2X150 capture7_n2) <(spec7X150 capture7_n2) <(spec12X150 capture7_n2) >> Spec.150.Table
paste <(specExon capture10_b10) <(spec2X150 capture10_b10) <(spec7X150 capture10_b10) <(spec12X150 capture10_b10) >> Spec.150.Table
paste <(specExon capture10_b1) <(spec2X150 capture10_b1) <(spec7X150 capture10_b1) <(spec12X150 capture10_b1) >> Spec.150.Table
paste <(specExon capture10_b2) <(spec2X150 capture10_b2) <(spec7X150 capture10_b2) <(spec12X150 capture10_b2) >> Spec.150.Table
paste <(specExon capture10_b6) <(spec2X150 capture10_b6) <(spec7X150 capture10_b6) <(spec12X150 capture10_b6) >> Spec.150.Table
paste <(specExon capture10_b7) <(spec2X150 capture10_b7) <(spec7X150 capture10_b7) <(spec12X150 capture10_b7) >> Spec.150.Table
paste <(specExon capture10_g1) <(spec2X150 capture10_g1) <(spec7X150 capture10_g1) <(spec12X150 capture10_g1) >> Spec.150.Table
paste <(specExon capture10_g5) <(spec2X150 capture10_g5) <(spec7X150 capture10_g5) <(spec12X150 capture10_g5) >> Spec.150.Table
paste <(specExon capture10_g6) <(spec2X150 capture10_g6) <(spec7X150 capture10_g6) <(spec12X150 capture10_g6) >> Spec.150.Table
paste <(specExon capture10_g7) <(spec2X150 capture10_g7) <(spec7X150 capture10_g7) <(spec12X150 capture10_g7) >> Spec.150.Table
paste <(specExon capture10_g8) <(spec2X150 capture10_g8) <(spec7X150 capture10_g8) <(spec12X150 capture10_g8) >> Spec.150.Table
paste <(specExon capture10_k1) <(spec2X150 capture10_k1) <(spec7X150 capture10_k1) <(spec12X150 capture10_k1) >> Spec.150.Table
paste <(specExon capture10_k2) <(spec2X150 capture10_k2) <(spec7X150 capture10_k2) <(spec12X150 capture10_k2) >> Spec.150.Table
paste <(specExon capture10_k3) <(spec2X150 capture10_k3) <(spec7X150 capture10_k3) <(spec12X150 capture10_k3) >> Spec.150.Table
paste <(specExon capture10_k5) <(spec2X150 capture10_k5) <(spec7X150 capture10_k5) <(spec12X150 capture10_k5) >> Spec.150.Table
paste <(specExon capture10_k8) <(spec2X150 capture10_k8) <(spec7X150 capture10_k8) <(spec12X150 capture10_k8) >> Spec.150.Table
paste <(specExon capture10_m1) <(spec2X150 capture10_m1) <(spec7X150 capture10_m1) <(spec12X150 capture10_m1) >> Spec.150.Table
paste <(specExon capture10_m3) <(spec2X150 capture10_m3) <(spec7X150 capture10_m3) <(spec12X150 capture10_m3) >> Spec.150.Table
paste <(specExon capture10_m5) <(spec2X150 capture10_m5) <(spec7X150 capture10_m5) <(spec12X150 capture10_m5) >> Spec.150.Table
paste <(specExon capture10_m7) <(spec2X150 capture10_m7) <(spec7X150 capture10_m7) <(spec12X150 capture10_m7) >> Spec.150.Table
paste <(specExon capture10_m9) <(spec2X150 capture10_m9) <(spec7X150 capture10_m9) <(spec12X150 capture10_m9) >> Spec.150.Table
paste <(specExon capture10_n1) <(spec2X150 capture10_n1) <(spec7X150 capture10_n1) <(spec12X150 capture10_n1) >> Spec.150.Table
paste <(specExon capture10_n2) <(spec2X150 capture10_n2) <(spec7X150 capture10_n2) <(spec12X150 capture10_n2) >> Spec.150.Table
paste <(specExon capture10_n3) <(spec2X150 capture10_n3) <(spec7X150 capture10_n3) <(spec12X150 capture10_n3) >> Spec.150.Table
paste <(specExon capture10_n4) <(spec2X150 capture10_n4) <(spec7X150 capture10_n4) <(spec12X150 capture10_n4) >> Spec.150.Table
paste <(specExon capture10_n5) <(spec2X150 capture10_n5) <(spec7X150 capture10_n5) <(spec12X150 capture10_n5) >> Spec.150.Table
paste <(specExon capture13_b10) <(spec2X150 capture13_b10) <(spec7X150 capture13_b10) <(spec12X150 capture13_b10) >> Spec.150.Table
paste <(specExon capture13_b1) <(spec2X150 capture13_b1) <(spec7X150 capture13_b1) <(spec12X150 capture13_b1) >> Spec.150.Table
paste <(specExon capture13_b2) <(spec2X150 capture13_b2) <(spec7X150 capture13_b2) <(spec12X150 capture13_b2) >> Spec.150.Table
paste <(specExon capture13_b6) <(spec2X150 capture13_b6) <(spec7X150 capture13_b6) <(spec12X150 capture13_b6) >> Spec.150.Table
paste <(specExon capture13_b7) <(spec2X150 capture13_b7) <(spec7X150 capture13_b7) <(spec12X150 capture13_b7) >> Spec.150.Table
paste <(specExon capture13_g1) <(spec2X150 capture13_g1) <(spec7X150 capture13_g1) <(spec12X150 capture13_g1) >> Spec.150.Table
paste <(specExon capture13_g5) <(spec2X150 capture13_g5) <(spec7X150 capture13_g5) <(spec12X150 capture13_g5) >> Spec.150.Table
paste <(specExon capture13_g6) <(spec2X150 capture13_g6) <(spec7X150 capture13_g6) <(spec12X150 capture13_g6) >> Spec.150.Table
paste <(specExon capture13_g7) <(spec2X150 capture13_g7) <(spec7X150 capture13_g7) <(spec12X150 capture13_g7) >> Spec.150.Table
paste <(specExon capture13_g8) <(spec2X150 capture13_g8) <(spec7X150 capture13_g8) <(spec12X150 capture13_g8) >> Spec.150.Table
paste <(specExon capture13_k1) <(spec2X150 capture13_k1) <(spec7X150 capture13_k1) <(spec12X150 capture13_k1) >> Spec.150.Table
paste <(specExon capture13_k2) <(spec2X150 capture13_k2) <(spec7X150 capture13_k2) <(spec12X150 capture13_k2) >> Spec.150.Table
paste <(specExon capture13_k3) <(spec2X150 capture13_k3) <(spec7X150 capture13_k3) <(spec12X150 capture13_k3) >> Spec.150.Table
paste <(specExon capture13_k5) <(spec2X150 capture13_k5) <(spec7X150 capture13_k5) <(spec12X150 capture13_k5) >> Spec.150.Table
paste <(specExon capture13_k8) <(spec2X150 capture13_k8) <(spec7X150 capture13_k8) <(spec12X150 capture13_k8) >> Spec.150.Table
paste <(specExon capture13_m1) <(spec2X150 capture13_m1) <(spec7X150 capture13_m1) <(spec12X150 capture13_m1) >> Spec.150.Table
paste <(specExon capture13_m3) <(spec2X150 capture13_m3) <(spec7X150 capture13_m3) <(spec12X150 capture13_m3) >> Spec.150.Table
paste <(specExon capture13_m5) <(spec2X150 capture13_m5) <(spec7X150 capture13_m5) <(spec12X150 capture13_m5) >> Spec.150.Table
paste <(specExon capture13_m7) <(spec2X150 capture13_m7) <(spec7X150 capture13_m7) <(spec12X150 capture13_m7) >> Spec.150.Table
paste <(specExon capture13_m9) <(spec2X150 capture13_m9) <(spec7X150 capture13_m9) <(spec12X150 capture13_m9) >> Spec.150.Table
paste <(specExon capture13_n1) <(spec2X150 capture13_n1) <(spec7X150 capture13_n1) <(spec12X150 capture13_n1) >> Spec.150.Table
paste <(specExon capture13_n2) <(spec2X150 capture13_n2) <(spec7X150 capture13_n2) <(spec12X150 capture13_n2) >> Spec.150.Table
paste <(specExon capture13_n3) <(spec2X150 capture13_n3) <(spec7X150 capture13_n3) <(spec12X150 capture13_n3) >> Spec.150.Table
paste <(specExon capture13_n4) <(spec2X150 capture13_n4) <(spec7X150 capture13_n4) <(spec12X150 capture13_n4) >> Spec.150.Table
paste <(specExon capture13_n5) <(spec2X150 capture13_n5) <(spec7X150 capture13_n5) <(spec12X150 capture13_n5) >> Spec.150.Table
paste <(specExon capture16_b10) <(spec2X150 capture16_b10) <(spec7X150 capture16_b10) <(spec12X150 capture16_b10) >> Spec.150.Table
paste <(specExon capture16_b1) <(spec2X150 capture16_b1) <(spec7X150 capture16_b1) <(spec12X150 capture16_b1) >> Spec.150.Table
paste <(specExon capture16_b2) <(spec2X150 capture16_b2) <(spec7X150 capture16_b2) <(spec12X150 capture16_b2) >> Spec.150.Table
paste <(specExon capture16_b6) <(spec2X150 capture16_b6) <(spec7X150 capture16_b6) <(spec12X150 capture16_b6) >> Spec.150.Table
paste <(specExon capture16_b7) <(spec2X150 capture16_b7) <(spec7X150 capture16_b7) <(spec12X150 capture16_b7) >> Spec.150.Table
paste <(specExon capture16_g1) <(spec2X150 capture16_g1) <(spec7X150 capture16_g1) <(spec12X150 capture16_g1) >> Spec.150.Table
paste <(specExon capture16_g5) <(spec2X150 capture16_g5) <(spec7X150 capture16_g5) <(spec12X150 capture16_g5) >> Spec.150.Table
paste <(specExon capture16_g6) <(spec2X150 capture16_g6) <(spec7X150 capture16_g6) <(spec12X150 capture16_g6) >> Spec.150.Table
paste <(specExon capture16_g7) <(spec2X150 capture16_g7) <(spec7X150 capture16_g7) <(spec12X150 capture16_g7) >> Spec.150.Table
paste <(specExon capture16_g8) <(spec2X150 capture16_g8) <(spec7X150 capture16_g8) <(spec12X150 capture16_g8) >> Spec.150.Table
paste <(specExon capture16_k1) <(spec2X150 capture16_k1) <(spec7X150 capture16_k1) <(spec12X150 capture16_k1) >> Spec.150.Table
paste <(specExon capture16_k2) <(spec2X150 capture16_k2) <(spec7X150 capture16_k2) <(spec12X150 capture16_k2) >> Spec.150.Table
paste <(specExon capture16_k3) <(spec2X150 capture16_k3) <(spec7X150 capture16_k3) <(spec12X150 capture16_k3) >> Spec.150.Table
paste <(specExon capture16_k5) <(spec2X150 capture16_k5) <(spec7X150 capture16_k5) <(spec12X150 capture16_k5) >> Spec.150.Table
paste <(specExon capture16_k8) <(spec2X150 capture16_k8) <(spec7X150 capture16_k8) <(spec12X150 capture16_k8) >> Spec.150.Table
paste <(specExon capture16_m1) <(spec2X150 capture16_m1) <(spec7X150 capture16_m1) <(spec12X150 capture16_m1) >> Spec.150.Table
paste <(specExon capture16_m3) <(spec2X150 capture16_m3) <(spec7X150 capture16_m3) <(spec12X150 capture16_m3) >> Spec.150.Table
paste <(specExon capture16_m5) <(spec2X150 capture16_m5) <(spec7X150 capture16_m5) <(spec12X150 capture16_m5) >> Spec.150.Table
paste <(specExon capture16_m7) <(spec2X150 capture16_m7) <(spec7X150 capture16_m7) <(spec12X150 capture16_m7) >> Spec.150.Table
paste <(specExon capture16_m9) <(spec2X150 capture16_m9) <(spec7X150 capture16_m9) <(spec12X150 capture16_m9) >> Spec.150.Table
paste <(specExon capture16_n1) <(spec2X150 capture16_n1) <(spec7X150 capture16_n1) <(spec12X150 capture16_n1) >> Spec.150.Table
paste <(specExon capture16_n2) <(spec2X150 capture16_n2) <(spec7X150 capture16_n2) <(spec12X150 capture16_n2) >> Spec.150.Table
paste <(specExon capture16_n3) <(spec2X150 capture16_n3) <(spec7X150 capture16_n3) <(spec12X150 capture16_n3) >> Spec.150.Table
paste <(specExon capture16_n4) <(spec2X150 capture16_n4) <(spec7X150 capture16_n4) <(spec12X150 capture16_n4) >> Spec.150.Table
paste <(specExon capture16_n5) <(spec2X150 capture16_n5) <(spec7X150 capture16_n5) <(spec12X150 capture16_n5) >> Spec.150.Table
```

# Specificity tables for 300 bp segements

Now we can create a specificity table for all exons and for expressed targets at 300 bp flanking regions using a few more BASH functions

```bash
specExon(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L sorted.ref3.0.exon.sc.hmask.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L sorted.ref3.0.exon.sc.hmask.neartarget.bed)
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L sorted.ref3.0.exon.sc.hmask.bed | samtools view - -@10 -c -L sorted.ref3.0.exon.sc.hmask.neartarget.bed)
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L sorted.ref3.0.exon.sc.hmask.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L sorted.ref3.0.exon.sc.hmask.slop.bed | samtools view - -@10 -c -L sorted.ref3.0.exon.sc.hmask.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e $1"\t"`python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f specExon
```

```bash
spec2X300(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc2.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc2.300.neartarget.bed )
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc2.bed | samtools view - -@10 -c -L all.filter.merged.EiRc2.300.neartarget.bed )
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc2.300.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc2.300.neartarget.bed  | samtools view - -@10 -c -L all.filter.merged.EiRc2.300.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f spec2X300
```

```bash
spec7X300(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc7.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc7.300.neartarget.bed )
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc7.bed | samtools view - -@10 -c -L all.filter.merged.EiRc7.300.neartarget.bed )
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc7.300.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc7.300.neartarget.bed  | samtools view - -@10 -c -L all.filter.merged.EiRc7.300.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f spec7X300
```

```bash
spec12X300(){

exon_reads=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc12.bed)
exon_nearr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc12.300.neartarget.bed )
exon_nearo=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc12.bed | samtools view - -@10 -c -L all.filter.merged.EiRc12.300.neartarget.bed )
exon_offtr=$(samtools view -@10 $1.norm.F.bam -c -L all.filter.merged.EiRc12.300.offtarget.bed)
exon_nearO=$(samtools view $1.norm.F.bam  -h -@10 -L all.filter.merged.EiRc12.300.neartarget.bed  | samtools view - -@10 -c -L all.filter.merged.EiRc12.300.offtarget.bed)
total=$(samtools view -@10 $1.norm.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"Percent\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"Percent\t"
  
}

export -f spec12X300
```

Now we use all the functions to create a table
```bash
echo -e "Sample\tPercent_in_exons\tPercent_near_exons\tPercent_off_target_exons\tPercent_on_target2X\tPercent_near_target2X\tPercent_off_target2X\tPercent_on_target7X\tPercent_tear_target7X\tPercent_off_target7X\tPercent_on_target12X\tPercent_near_target12X\tPercent_off_target12X" > Spec.300.Table

# Only use once to generate the list of commands below
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
paste <(specExon ${i}) <(spec2X300 ${i}) <(spec7X300 ${i}) <(spec12X300 ${i}) >> Spec.300.Table
done
```

```R
{r Load Spec.150.table.csv and summarize data}
data <- read.csv("/home/jgreen/eager_obj1b/04_coverage_analysis/04_specificity/Spec.150.Table.csv")

# Split the Sample column into Capture and Sample
data <- data %>%
  separate(Sample, into = c("Capture", "Sample"), sep = "_")

# Summary statistics
summary <- summary(data)

# ANOVA
anova_result_150 <- aov(Percent_in_exons ~ Capture, data = data)

# Print summary statistics
print(summary)

# Print ANOVA results
print(summary(anova_result))
```

```R
{r Create boxplots 150 bp}
# Create a box plot with larger facets
boxplot_data <- data %>%
  gather(key = "Variable", value = "Value", -Capture, -Sample)

# Define the desired order of variables
desired_order <- c(
  "Percent_in_exons", "Percent_near_exons", "Percent_off_target_exons",
  "Percent_on_target2X", "Percent_near_target2X", "Percent_off_target2X",
  "Percent_on_target7X", "Percent_near_target7X", "Percent_off_target7X",
  "Percent_on_target12X", "Percent_near_target12X", "Percent_off_target12X"
)

# Reorder the levels of the "Variable" factor based on the desired order
boxplot_data$Variable <- factor(boxplot_data$Variable, levels = desired_order)

boxplot_plot <- boxplot_data %>%
  ggplot(aes(x = Capture, y = Value, fill = Capture)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free", nrow = 4) +
  labs(title = "Differences Between 150 bp Capture Targets Across Variables",
       x = "Capture Type", y = "Value") +
  scale_x_discrete(limits=c("capture1", "capture4", "capture7", "capture10", "capture13", "capture16")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))

# Print the box plot
print(boxplot_plot)

# Save the plot with the specific order of facets
save_path <- "/home/jgreen/eager_obj1b/04_coverage_analysis/04_specificity/image_results/"
ggsave(filename = paste0(save_path, "Spec_150_boxplot_plot.png"), plot = boxplot_plot, width = 12, height = 8)
```

## Generating data for figure 7


The first step is to generate depth per bp data for each capture pool.

```bash
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel "samtools depth -aa {}.norm.F.bam > {}.genome.depth"
```

Next, we extract the region of interest.

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
for i in  "${StringArray[@]}"
do
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ${i}.genome.depth | mawk '$2 < 32866205' > ${i}.indiv.hsp.graph.depth
done
```

```bash
# Simplified code for above
ls capture*.norm.F.bam | sed 's/.norm.F.bam//g' | parallel "samtools depth -aa {}.norm.F.bam | mawk '\$1 ~ /NC_035780.1/ && \$2 > 32736205 && \$2 < 32866205' > {}.indiv.hsp.graph.depth"
```

```bash
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' Capture4_M4.genome.depth | mawk '$2 < 32866205' > Capture4_M4.indiv.hsp.graph.depth
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' Capture4_N2.genome.depth | mawk '$2 < 32866205' > Capture4_N2.indiv.hsp.graph.depth
```

Next, we add a column for the pool identifier and concatenate into single data table.


```bash
# Make header
echo -e "Contig\tbp\tDepth\tSample" > header
```

Use a for loop of sample names to generate SED code

```bash
declare -a StringArray=("capture1_b1" "capture1_b2" "capture1_g1" "capture1_g7" "capture1_k1" "capture1_k3" "capture1_m5" "capture1_m7" "capture1_n1" "capture1_n2" "capture4_b1" "capture4_b2" "capture4_g1" "capture4_g7" "capture4_k1" "capture4_k3" "capture4_m5" "capture4_m7" "capture4_n1" "capture4_n2" "capture7_b1" "capture7_b2" "capture7_g1" "capture7_g7" "capture7_k1" "capture7_k3" "capture7_m5" "capture7_m7" "capture7_n1" "capture7_n2" "capture10_b10" "capture10_b1" "capture10_b2" "capture10_b6" "capture10_b7" "capture10_g1" "capture10_g5" "capture10_g6" "capture10_g7" "capture10_g8" "capture10_k1" "capture10_k2" "capture10_k3" "capture10_k5" "capture10_k8" "capture10_m1" "capture10_m3" "capture10_m5" "capture10_m7" "capture10_m9" "capture10_n1" "capture10_n2" "capture10_n3" "capture10_n4" "capture10_n5" "capture13_b10" "capture13_b1" "capture13_b2" "capture13_b6" "capture13_b7" "capture13_g1" "capture13_g5" "capture13_g6" "capture13_g7" "capture13_g8" "capture13_k1" "capture13_k2" "capture13_k3" "capture13_k5" "capture13_k8" "capture13_m1" "capture13_m3" "capture13_m5" "capture13_m7" "capture13_m9" "capture13_n1" "capture13_n2" "capture13_n3" "capture13_n4" "capture13_n5" "capture16_b10" "capture16_b1" "capture16_b2" "capture16_b6" "capture16_b7" "capture16_g1" "capture16_g5" "capture16_g6" "capture16_g7" "capture16_g8" "capture16_k1" "capture16_k2" "capture16_k3" "capture16_k5" "capture16_k8" "capture16_m1" "capture16_m3" "capture16_m5" "capture16_m7" "capture16_m9" "capture16_n1" "capture16_n2" "capture16_n3" "capture16_n4" "capture16_n5")
# For loop to print out sed command
for i in "${StringArray[@]}"
do
sed -i "s/$/\t${i}/g" ${i}.indiv.hsp.graph.depth
cat header ${i}.indiv.hsp.graph.depth > TotalCov${i}.indiv.txt
done

# For loop for cat command
for i in "${StringArray[@]}"
do
echo "cat header ${i}.indiv.hsp.graph.depth > TotalCov${i}.indiv.txt"
done
```

For the graph, we all need the annotations for the gene regions, the exons, and CDS

```bash
mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "exon"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/Percent.*//g' > exons
cat <(echo -e "Contig\tStart\tEnd\tTreatment") exons > exon.list

mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "mRNA"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/Percent.*//g' > genes
cat <(echo -e "Contig\tStart\tEnd\tTreatment") genes > genes.list

mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "CDS"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/Percent.*//g' > CDS
cat <(echo -e "Contig\tStart\tEnd\tTreatment") CDS > CDS.list
```


R code for Figure 7

```r
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")
DepCap1 <- read.table("TotalCovCap1.indiv.txt", header = TRUE)
DepCap1 <- read.table("TotalCovCap2.indiv.txt", header = TRUE)
DepCap1 <- read.table("TotalCovCap3.indiv.txt", header = TRUE)
DepCap1 <- read.table("TotalCovCap4.indiv.txt", header = TRUE)


DepC <- as.data.frame(DepCap1)
DepC$Sample <- factor(DepC$Sample,levels=c("capture1"))
DepR <- as.data.frame(DepCap2)
DepR$Sample <- factor(DepR$Sample,levels=c("capture2"))

exons <- read.table("exon.list", header = TRUE, sep = "\t")
exons <- as.data.frame(exons)

genes <- read.table("genes.list", header = TRUE, sep = "\t")
genes <- as.data.frame(genes)

cds <- read.table("cds.list", header = TRUE, sep = "\t")
cds <- as.data.frame(cds)

subDepC <-subset(DepC, bp <32755000 & bp > 32739000)
subDepR <-subset(DepR, bp <32755000 & bp > 32739000)
subexons <-subset(exons, End <32755205 & End > 32740205)
subgenes <-subset(genes, End <32800757 & Start < 32754201)
subcds <-subset(cds, End <32800757 & Start < 32755000)
#subDepR$Depth <- subDepR$Depth / -1
submean.cov <- ddply(subDepC, .(Contig,bp), summarize,  Depth=mean(Depth))
submeanR.cov <- ddply(subDepR, .(Contig,bp), summarize,  Depth=mean(Depth))
subgenes$End[4] <- 32755000


pie(rep(1, length(cbPalette)), labels = sprintf("Percentd (Percents)", seq_along(cbPalette), 
                                                cbPalette), col = cbPalette)
redcol <-"#940000"

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")
cbPalettedd <- c( "#009E73","#D55E00", "#E69F00")


dd <- ggplot(subDepC, aes(x= bp, y=Depth)) +
  geom_area(aes(group=Sample),position = "identity",color=alpha("grey30",0.25),fill=cbPalette[4], alpha=0.1, linetype="dotted")+  
  geom_line(data=submean.cov,aes(y=rollmean(Depth, 100, na.pad=TRUE)),colour=cbPalette[4], size =1.0, alpha=0.9)  +
  geom_line(data=submeanR.cov,aes(y=rollmean(Depth, 100, na.pad=TRUE)),colour=redcol, size =1.0, alpha=0.9)  +
  geom_area(data=subDepR, aes(group=Sample),position = "identity",color=alpha("grey30",0.25),fill=redcol, alpha=0.1, linetype="dotted")+
  scale_color_manual(values=cbPalettedd) +
  geom_segment(data=subgenes, aes(x = Start, y = 715, xend = End, yend = 715), size = 6,color=cbPalette[9], alpha=1)+
  geom_segment(data=subexons,aes(x = Start, y = 715, xend = End, yend = 715, color=Treatment),size = 4, alpha=1) +
  geom_segment(data=subcds,aes(x = Start, y = 715, xend = End, yend = 715),size = 1, color="grey90", alpha=1) +
  theme_bw()+
  coord_cartesian(xlim = c(32740000,32755000))+
  xlim(32740000,32755000) +
  scale_y_continuous(limits=c(0,800),labels=c("0","250","500","750"), breaks=c(0,250,500,750),expand=c(0.01,0)) +
  theme(legend.position="none")

png(filename="Figure7.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
dd
dev.off()
```