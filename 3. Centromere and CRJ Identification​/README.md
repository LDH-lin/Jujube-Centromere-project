## Centromere identification using bamcoverage.sh script
## Input files: ChIP_R* and Input_R* (post-fastp quality control)
bamcoverage.sh Chip_R1 Chip_R2 Input_R1 Input_R2 genome.fa 

## CRJ enrichment analysis (Alignment Score calculation)
auto_blast.py -q CRJ.fa -d genome.fa -e 1e-5 -th 30 -t blastn -m blast -o CRJ.blast  ## Align CRJ to genome using BLASTn
AlignmentScore.py CRJ.blast CRJ.kmer 100

## CRJ-ChIP-seq co-enrichment analysis
awk '{print$2"\t"$11"\t"$12}' CRJ.blast >CRJ.bed
bed_sort.py CRJ.bed CRJ.sort.bed
bedtools intersect -b CRJ.sort.bed -a CRJ.50k.bed -c > CRJ.50k.bedgraph
bedGraphToBigWig CRJ.50k.bedgraph DZ.size CRJ.50k.bw
computeMatrix reference-point -R hap1.crz.bed -S CRJ.50k.bw \
--beforeRegionStartLength 1000 --regionBodyLength 2000 --afterRegionStartLength 1000 -o test 
plotHeatmap -m test -o test_heatmap.pdf --colorMap RdBu --legendLocation best --refPointLabel "Center"