# ​​Repetitive Sequence Identification​
## Step1. TE analysis for gene annotation and centromere
```

BuildDatabase -name DZ    DZ.genome.fa  > build.log  2>&1
RepeatModeler -database DZ.genome.fa  -threads 30 -engine ncbi -LTRStruct -genomeSampleSizeMax 1000000000 > repeatmodeler.log 2>&1
RepeatMasker -lib DZ-families.fa ‐xsmall  -gff  -e ncbi -dir .  -pa 30 DZ.genome.fa  > repeatmasker.log 2>&1
```

## Step2. Identification of LTR-RTs using ltr.retriever.sh
```
ltr.retriever.sh DZ /data/DZ.genome.fa DZ.genome.fa
```
## Step3. Tandem Repeat Identification​​
```
trf DZ.genome.fa   2 7 7 80 10 50 500 -l 20 -d -m -h
trf2gff -i DZ.genome.fa.2.7.7.80.10.50.500.dat -o DZ.genome.fa.2.7.7.80.10.50.500.gff
```
## Step3. TE analysis for constructing phylogenetic trees of LTR-RT structural domains across multiple species
```

"DZ.genome.ltr.fa" ## LTR-RT sequences extracted from ltr.retriever identification results
"DZ" ## Prefix for generated files
TEsorter -db rexdb-plant -p 20 DZ.genome.ltr.fa  ## Identify LTR-RT structural domains using TEsorter
concatenate_domains.py DZ.genome.ltr.fa.rexdb-plant.cls.pep RT GAG INT RH     >   DZ.rexdb-plant.pep.full.aln ## Extract structural sequences based on identification results and perform multiple sequence alignment
iqtree2 -s DZ.rexdb-plant.pep.full.aln -m MFP -bb 1000  -bnni ## Construct phylogenetic tree
```

# Assessment of repetitive sequence proportion in centromeric regions

```
DZ.out -- RepeatMasker .out file for centromeric region
DZ.genome.fa.2.7.7.80.10.50.500.gff -- TRF gff3 file for centromeric region
DZ.cen.bed -- Centromeric region bed file
repeat.count.py speice.out  speice.trf.gff3 speice.cen.bed output.count ## Calculate the proportion of repetitive sequences in centromeric regions
```

