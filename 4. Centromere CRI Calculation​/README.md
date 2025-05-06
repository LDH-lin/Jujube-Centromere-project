# Extract LTR sequences and calculate CLTRI (Centromere-Linked Tandem Repeat Index)
## Step1: Extract LTR sequences for each species
```
Speice.fa.pass.list.gff3 -- LTR identification results from ltr.retriever
awk '$3=="long_terminal_repeat" {print$1"\t"$4"\t"$5} ' Speice.fa.pass.list.gff3 > Speice.ltr.bed
bedtools getfasta -fi Speice.genome.fa -bed Speice.ltr.bed -fo Speice.ltr.fa
```

## Step2: Perform self-alignment of LTR sequences for each species
```
auto_blast -q Speice.ltr.fa -d Speice.genome.fa -t blastn -m blast -o Speice.ltr.blast -e 1e-5 
CLTRI.count.py Speice.ltr.blast Speice.cen.bed 30 > Sspeice.list 
awk '{print $1"\t"$2/$3}' Speice.list > Speice.ltr.CLTRI
```

## Step3: Phylogenetic tree construction based on manually curated CLTRs
```
Speice.cltr.fa -- Consensus CLTR sequences across species
Speice.pass.list -- ltr.retriever filtered pass.list file
insertion.time.py Speice.cltr.fa Speice.cltr.time.fa ## Extract insertion time for each CLTR
mafft --auto speice.cltr.time.fa --adjustdirection --maxiterate 1000 > Speice.cltr.time.aln
iqtree2 -s Speice.rexdb-plant.pep.full.aln -m MFP -bb 1000 -bnni 
```

## Step4: Identify relationships between TRs and CLTRs in jujube
```
jujube.trf.fa -- TRF sequences from three jujube cultivars
jujube.J*.fa -- CLTR sequences from J1/J2 branches
cd-hit -c 0.8 -i jujube.J*.fa -o jujube.J*.cd.fa
auto_blast -q jujube.J*.cd.fa -d jujube.trf.fa -t blastn -m blast -e 1e-5 -o jujube.J*.blast
```
