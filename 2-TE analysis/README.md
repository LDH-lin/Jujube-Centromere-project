# Step1. TE analysis for gene annotation and centromere
BuildDatabase -name ${name}    ${name}  > build.log  2>&1
RepeatModeler -database ${name}  -threads 30 -engine ncbi -LTRStruct -genomeSampleSizeMax 1000000000 > repeatmodeler.log 2>&1
RepeatMasker -lib ../01_db/${name}-families.fa ‐xsmall  -gff  -e ncbi -dir .  -pa 30 ${genome}  > repeatmasker.log 2>&1

# Step2. Identification of LTR-RTs using ltr.retriever.sh

# Step3. TE analysis for constructing phylogenetic trees of LTR-RT structural domains across multiple species
ltr.fa="genome.ltr.fa" ## LTR-RT sequences extracted from ltr.retriever identification results
pre="DZ" ## Prefix for generated files
TEsorter -db rexdb-plant -p 20 ${ltr.fa}  ## Identify LTR-RT structural domains using TEsorter
concatenate_domains.py ${ltr.fa}.fa.rexdb-plant.cls.pep RT GAG INT RH     >   ${pre}.rexdb-plant.pep.full.aln ## Extract structural sequences based on identification results and perform multiple sequence alignment
iqtree2 -s ${pre}.rexdb-plant.pep.full.aln -m MFP -bb 1000  -bnni ## Construct phylogenetic tree

# Step4. Assessment of repetitive sequence proportion in centromeric regions
## ${speice}.out -- RepeatMasker .out file for centromeric region
## ${speice}.trf.gff3 -- TRF gff3 file for centromeric region
## ${speice}.cen.bed -- Centromeric region bed file
repeat.count.py ${speice}.out  ${speice}.trf.gff3 ${speice}.cen.bed output.count ## Calculate the proportion of repetitive sequences in centromeric regions

