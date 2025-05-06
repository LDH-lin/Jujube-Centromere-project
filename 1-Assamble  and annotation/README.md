# haplotype assembly, DZ means "Dongzao"
# current directory: 
```

hifiasm -o DZ_hap -t 40 --ul-cut 15000 -D10 --hom-cov 71 --h1 DZ_Hic_clean.R1.fastq.gz --h2 DZ_Hic_clean.R1.fastq.gz --ul ONT_UL_pass_reads.fasta.gz Hifi.ccs.fastq.gz >DZ_hap.log 2>DZ_hap.err.log
```

# the above command generated two preliminary haplotypes: DZ_hap.hic.hap1.p_ctg.gfa and DZ_hap.hic.hap2.p_ctg.gfa
# run purge_dups to remove filter the raw assemply

sh run_purge_dups.sh

# This will generated cleaned contigs named DZ_hap1.purged.fas and DZ_hap2.purged.fas
```
#run quarTeT to link contigs to chromosomes using our previous T2T assembly of DZ
python quartet.py AssemblyMapper -r DZ.t2t.final.fas -q DZ_hap1.purged.fas -i 95 -p DZ_hap1 -t 10 --plot --overwrite
python quartet.py AssemblyMapper -r DZ.t2t.final.fas -q DZ_hap2.purged.fas -i 95 -p DZ_hap2 -t 10 --plot --overwrite

#fill the gaps
python quartet.py GapFiller -d DZ_hap1.draftgenome.fasta -g Hifi.ccs.fastq.gz -i 70 -t 50 -p DZ_hap1.fgap --enablejoin --overwrite
python quartet.py GapFiller -d DZ_hap2.draftgenome.fasta -g Hifi.ccs.fastq.gz -i 70 -t 50 -p DZ_hap2.fgap --enablejoin --overwrite

#Through the above steps, we finshed the final assembly of the haplotype-resolved gapfree T2T genome of Dongzao.
#The two haplotypes were named DZ_hap1.genome.fas and DZ_hap2.genome.fas

# Repeat annotation of two haplotype genome
# First run repeatmodeler
BuildDatabase -name DZ_hap1 DZ_hap1.genome.fas
RepeatModeler -database DZ_hap1 -pa 30 -genomeSampleSizeMax 500000000
cat DZ_hap1-families.fas RepeatMasker.lib > DZ.lib
#run Repeatmasker
RepeatMasker -e rmblast -pa 10 -lib DZ.lib -dir DZ_hap1.Rm_results -gff -a DZ_hap1.genome.fas
RepeatMasker -e rmblast -pa 10 -lib DZ.lib -dir DZ_hap2.Rm_results -gff -a DZ_hap2.genome.fas

##Protein gene prediction and functional annotation

sh run_annotation.sh sample.lst
```
