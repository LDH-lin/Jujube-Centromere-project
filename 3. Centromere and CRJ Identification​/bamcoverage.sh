#!/usr/bin/bash


chip1=$1
chip2=$2
input1=$3
input2=$4
genome=$5

name=$(basename "${genome}" .Hifi.final.fas)
mkdir ${name}
cd ${name}
ln -s $genome  ./${name}.Hifi.final.fas
bowtie2-build  --threads 20   ${name}.Hifi.final.fas  ${name}.Hifi.final.fas
bowtie2 -p 20 -x  ${name}.Hifi.final.fas  -1 $chip1  -2 $chip2   -S ${name}.chip.sam
samtools view -bS  -@ 20 ${name}.chip.sam  >  ${name}.chip.bam
samtools sort -@ 20 ${name}.chip.bam >  ${name}.chip.sort.bam
samtools index -@ 20   ${name}.chip.sort.bam
rm ${name}.chip.sam ${name}.chip.bam
bowtie2 -p 20 -x  ${name}.Hifi.final.fas  -1 $input1 -2 $input2 -S ${name}.input.sam
samtools view -bS  -@ 20 ${name}.input.sam > ${name}.input.bam
samtools sort -@ 20 ${name}.input.bam > ${name}.input.sort.bam
rm ${name}.input.bam ${name}.input.sam
samtools index -@ 20 ${name}.input.sort.bam
alignmentSieve -p 10 --bam  ${name}.input.sort.bam  --outFile  ${name}.input.fill.bam     --
ignoreDuplicates --minMappingQuality 20 --samFlagExclude 260
alignmentSieve -p 10 --bam  ${name}.chip.sort.bam  --outFile  ${name}.chip.fill.bam    --ign
oreDuplicates --minMappingQuality 20 --samFlagExclude 260


bamCompare --bamfile1 ${name}.chip.sort.bam  --bamfile2 ${name}.input.sort.bam   --binSize 1
000 --numberOfProcessors 40 --operation ratio   --outFileFormat bedgraph    -o ${name}.1k.be
dgraph
samtools faidx ${name}.Hifi.final.fas
cut -f 1,2 ${name}.Hifi.final.fas.fai  > ./size
bedtools makewindows -g size -w 5000 > genome_bins_5kb.bed
bedtools map -a genome_bins_5kb.bed   -b ${name}.bedgraph -c 4 -o mean >  ${name}.bamcoverag
e.bedgraph
g=$(awk '{sum += $2} END {print sum}' size)

macs3 callpeak    --broad-cutoff 0.01  --broad   -t ${name}.chip.sort.bam  -c ${name}.input.
sort.bam  -f BAM -g ${g} --outdir peak/ -n ${name}.sort -q 0.01
macs3 callpeak    --broad-cutoff 0.01  --broad   -t ${name}.chip.fill.bam  -c ${name}.input.
fill.bam  -f BAM -g ${g} --outdir peak/ -n ${name}.fill -q 0.01
