#!/bin/sh
current=`pwd`
cpu=60
pwd_gseq=$pwd/00-assembly
pwd_rnaseq=$pwd/RNA-seq-for-annotation
repeat_path=$pwd/01_repeatmasker
script=$pwd/00-script
ulimit -n 100000
input_list=$1

cat $input_list | while read name
do 

  rm -rf $name
  mkdir $name
  cd $name
  rm -rf 01-STAR
  mkdir 01-STAR
  cd 01-STAR
  date
  echo "Run STAR..."
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  STAR --runThreadN 80 --runMode genomeGenerate --genomeSAindexNbases 13 --genomeDir $current/$name/01-STAR/index \
       --genomeFastaFiles $name.Hifi.final.fas
  STAR --runThreadN 80 --genomeDir $current/$name/01-STAR/index \
       --readFilesCommand zcat --readFilesIn $pwd_rnaseq/${name}.RNA.R1.clean.fq.gz $pwd_rnaseq/${name}.RNA.R2.clean.fq.gz \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ./$name \
       --outWigType bedGraph \
       --outSAMstrandField intronMotif
  stringtie ${name}Aligned.sortedByCoord.out.bam -o $name.gtf -l $name -p $cpu
  gffread -w $name.fas -g $name.Hifi.final.fas $name.gtf
  mv $name.fas $name.total.transcrips.fas
  cd ..

  date
  echo "Done STAR, Run PASA..."
  rm -rf 02-PASA
  mkdir 02-PASA
  cd 02-PASA
  ln -s ../01-STAR/$name.total.transcrips.fas .
  cp /data/02_work/yangm/02_software/01_package/PASApipeline/alignAssembly.config .
  sed -i "s/DATABASE=.*/DATABASE=\/tmp\/$name.db/" alignAssembly.config
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  $PASAHOME/bin/seqclean $name.total.transcrips.fas
  rm -rf /tmp/$name.db
  $PASAHOME/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $name.Hifi.final.fas \
       -t $name.total.transcrips.fas.clean -T -u $name.total.transcrips.fas \
       --ALIGNERS blat,gmap,minimap2 --CPU $cpu > $name.pasa.log 2> $name.pasa.error.log
  ln -s $name.db.pasa_assemblies.gff3 $name.pasa.gff3
  #rm -rf pblat_outdir $name.Hifi.final.fas.gmap $name.Hifi.final.fas.mm2
  cd ..

  date
  echo "Done PASA, Run GeMoMa..."
  rm -rf 03-GeMoMa
  mkdir 03-GeMoMa
  cd 03-GeMoMa
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/apple.genome.fna
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/apple.genome.gff3
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/Ath.fna
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/Ath.gff3
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/Jujube.Hifi.fna
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/Jujube.Hifi.gff3
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/P.armeniaca.fna
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/P.armeniaca.gff3
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/populus.fna
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/populus.gff3
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/P.persica.fna
  ln -s /data/02_work/yangm/03_project/01_pangenome26/03_ref_species/P.persica.gff3
  java -jar /data/02_work/yangm/02_software/01_package/GeMoMa-1.9/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=$cpu \
           AnnotationFinalizer.r=NO o=true p=false t=$name.Hifi.final.fas \
	   s=own i=A1 g=apple.genome.fna a=apple.genome.gff3 \
	   s=own i=A2 g=Ath.fna a=Ath.gff3 \
	   s=own i=A3 g=Jujube.Hifi.fna a=Jujube.Hifi.gff3 \
	   s=own i=A4 g=P.armeniaca.fna a=P.armeniaca.gff3 \
	   s=own i=A5 g=populus.fna a=populus.gff3 \
	   s=own i=A6 g=P.persica.fna a=P.persica.gff3 \
	   > $name.gemoma.log 2>$name.gemoma.error.log

  rm -rf GeMoMa_temp unfiltered*
  perl /data/02_work/yangm/02_software/01_package/EVidenceModeler/EvmUtils/misc/GeMoMa_gff_to_gff3.pl final_annotation.gff > $name.Gemoma.gff3
  cd ..

  date
  echo "Done GeMoMa, Run Augustus..."
  rm -rf 04-Augustus
  mkdir 04-Augustus
  cd 04-Augustus
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  
  num=$(grep ">" $name.Hifi.final.fas | wc -l)

  perl /data/02_work/yangm/01_script/blast-spliting-blast.pl $name.Hifi.final.fas aa $num F
  
  for i in $(seq 1 $num)
  do 
     nohup augustus --species=DZ --genemodel=complete --gff3=on $name.Hifi.final.fas.$i >$name.split.$i.gff3 2>$name.split.$i.log&
  done

  while ps -u yangm | grep -q "augustus"
  do
     sleep 60
     continue
  done

  for i in $(seq 1 $num)
  do 
    rm -rf $name.split.$i
    rm -rf $name.split.$i.log
    perl /data/02_work/yangm/02_software/01_package/EVidenceModeler/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl $name.split.$i.gff3 > $name.augustus.$i.gff3
  done
  
  cat $name.augustus.*.gff3 > $name.augustus.gff3
  perl -ane 's/model/$F[0].model/g; print;' $name.augustus.gff3 > tmp
  perl -ane 's/ID=gene/ID=$F[0].gene/g; print;' tmp > tmp1
  perl -ane 's/Parent=gene/Parent=$F[0].gene/g; print;' tmp1 > tmp2
  mv tmp2 $name.augustus.gff3
  rm tmp tmp1
  rm $name.split.*
  rm $name.augustus.*.gff3 
  cd ..


  date
  echo "Done Augustus, Run EVM..."
  rm -rf 05-EVM
  mkdir 05-EVM
  cd 05-EVM
  ln -s ../04-Augustus/$name.augustus.gff3 .
  ln -s ../03-GeMoMa/$name.Gemoma.gff3 .
  ln -s ../02-PASA/$name.pasa.gff3 .
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  echo "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
  echo "PROTEIN\tGeMoMa\t5" >> weights.txt
  echo "TRANSCRIPT\tassembler-$name.db\t9" >> weights.txt
  
  $EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome $name.Hifi.final.fas \
     --gene_predictions $name.augustus.gff3 \
     --protein_alignments $name.Gemoma.gff3 \
     --transcript_alignments $name.pasa.gff3 \
     --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out \
     > $name.partition_EVM.log 2> $name.partition_EVM.error.log

  $EVM_HOME/EvmUtils/write_EVM_commands.pl --genome $name.Hifi.final.fas \
      --weights $current/$name/05-EVM/weights.txt \
      --gene_predictions $name.augustus.gff3 --protein_alignments $name.Gemoma.gff3 \
      --transcript_alignments $name.pasa.gff3 \
      --output_file_name evm.out --partitions partitions_list.out > commands.list 2> $name.write_EVM.error.log

  $EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list > run-execute.log 2>run-execute.error.log

  $EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out >recombine.log 2>recombine.err.log

  $EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $name.assemble.fa >convert.log 2>convert.err.log

  find . -regex ".*evm.out.gff3" -exec cat {} \; > $name.EVM.gff3
  
  cp $name.EVM.gff3 $name.pre.gff3
  sed -i s/EVM%20prediction%20//g $name.pre.gff3
  sed -i s/evm.model/t/g $name.pre.gff3
  sed -i s/evm.TU/g/g $name.pre.gff3
  perl $script/format_gff3.pl $name.pre.gff3 $name.pre.gff3.1
  perl $script/fill_0_num.pl $name.pre.gff3.1 $name.pre.gff3
  rm -rf $name.pre.gff3.1
  find . -type d -name "*Chr*" | xargs rm -rf
  rm -rf convert-* recombine-* run-*
  cd ..


  date
  echo "Done EVM, Run filter..."
  rm -rf 06-filter
  mkdir 06-filter
  cd 06-filter
  ln -s ../05-EVM/$name.pre.gff3 .
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  perl /data/02_work/yangm/01_script/gff3_gene_extractor_v2.pl -i $name.pre.gff3 -g $name.Hifi.final.fas > tmp.log 2> tmp.err.log
  cp $name.pre.gff3.coding.primarySeq $name.pre.coding
  rm -rf $name.pre.gff3.*
  perl /data/02_work/yangm/01_script/nuclear2aa_1.pl $name.pre.coding $name.pre.coding.1 $name.pre.coding.2 $name.pre.coding.3
  grep ">" $name.pre.coding.2 | cut -c 2-100 | perl $script/classify_gff3.pl -i - -j $name.pre.gff3 -o $name.pre.gff3.nopseudo -n 1
  perl $script/sel_PSI_ratio.pl $repeat_path/$name/${name}.Rm_results/$name.fas.masked $name.pre.gff3.nopseudo $name.Nratio.lst
  perl $script/sel_subseq_gene.pl $name.pre.coding.2 $name.Nratio.lst $name.Nratio.pro $name.outNr.pro
  
  ##run emapper annotation for $name.Nratio.pro
  mkdir ${name}_Nratio
  cd ${name}_Nratio
  ln -s ../$name.Nratio.pro .
  rm -rf ${name}_mmseqs ${name}_tmp
  mkdir ${name}_mmseqs ${name}_tmp
  emapper.py -i $name.Nratio.pro -m mmseqs --cpu 30 --num_servers 15 --no_file_comments -o $name.Nratio --output_dir ./${name}_mmseqs --temp_dir=./${name}_tmp --usemem > $name.Nratio.emapper.log 2>$name.Nratio.emapper.err.log
  cd ${name}_mmseqs
  perl $script/remove.transposons.emapper.pl $name.Nratio.emapper.annotations $name.discard.lst $name.Nratio.emapper.annotations.left
  mv $name.Nratio.emapper.annotations.left $name.Nratio.emapper.annotations
  rm $name.Nratio.emapper.hits $name.Nratio.emapper.seed_orthologs
  cut -f 1,12,13 $name.Nratio.emapper.annotations | awk '$2!="-"' > $name.Nratio.emapper.kegg
  awk '{print $1}' $name.Nratio.emapper.annotations | sort | uniq | grep -v "query" > $name.Nratio.recover.lst
  cd ..
  rm -rf ${name}_tmp 
  cd ..

  ##run emapper annotation for $name.outNr.pro
  mkdir ${name}_outNr
  cd ${name}_outNr
  ln -s ../$name.outNr.pro .
  rm -rf ${name}_mmseqs ${name}_tmp
  mkdir ${name}_mmseqs ${name}_tmp
  emapper.py -i $name.outNr.pro -m mmseqs --cpu 30 --num_servers 15 --no_file_comments -o $name.outNr --output_dir ./${name}_mmseqs --temp_dir=./${name}_tmp --usemem > $name.outNr.emapper.log 2>$name.outNr.emapper.err.log
  cd ${name}_mmseqs
  perl $script/remove.transposons.emapper.pl $name.outNr.emapper.annotations $name.outNr.discard.lst $name.outNr.emapper.annotations.left
  mv $name.outNr.emapper.annotations.left $name.outNr.emapper.annotations
  rm $name.outNr.emapper.hits $name.outNr.emapper.seed_orthologs
  cut -f 1,12,13 $name.outNr.emapper.annotations | awk '$2!="-"' > $name.outNr.emapper.kegg
  perl $script/sel_subseq_gene.pl ../$name.outNr.pro $name.outNr.discard.lst $name.outNr.discard.fas $name.outNr.left.fas
  grep ">" $name.outNr.left.fas | awk '{print $1}' | sed "s/>//" > $name.outNr.left.lst
  cd ..
  rm -rf ${name}_tmp
  cd ..
  
  ##combine emapper annotation and extract the final gff
  mkdir ${name}_emappy
  cd ${name}_emappy
  ln -s ../${name}_Nratio/${name}_mmseqs/$name.Nratio.recover.lst .
  ln -s ../${name}_outNr/${name}_mmseqs/$name.outNr.left.lst
  cat $name.Nratio.recover.lst $name.outNr.left.lst > $name.finalgene.lst
  perl $script/classify_gff3.pl -i $name.finalgene.lst -j ../$name.pre.gff3 -o $name.final.gff3 -n 1
  ln -s ../${name}_Nratio/${name}_mmseqs/$name.Nratio.emapper.annotations
  ln -s ../${name}_outNr/${name}_mmseqs/$name.outNr.emapper.annotations
  head -n 1 $name.Nratio.emapper.annotations > $head.txt
  cat $name.Nratio.emapper.annotations $name.outNr.emapper.annotations | sort -k1,1 | grep -v "query" > temp.txt
  cat $head.txt temp.txt > $name.final.emapper.annotations
  cut -f 1,12,13 $name.final.emapper.annotations | awk '$2!="-"' > $name.final.emapper.annotations.kegg
  rm -rf $head.txt temp.txt
  cd ..
  ln -s ${name}_emappy/$name.final.gff3 .
  cd ..

  date
  echo "Done filter, running sequence extraction...."
  rm -rf 07-extract
  mkdir 07-extract
  cd 07-extract
  ln -s ../06-filter/$name.final.gff3 .
  ln -s $pwd_gseq/$name.Hifi.final.fas .
  perl /data/02_work/yangm/01_script/gff3_gene_extractor_v2.pl -i $name.final.gff3 -g $name.Hifi.final.fas -o $name.gene -s
  sed -i "s/\*$//" $name.gene.protein.primarySeq
  rm -rf *.gene *.intron *.transcript *.exon *.coding *.protein *.utr* $name.gene.transcript.primarySeq $name.gene.intron.primarySeq

:<<BLOCK  
  ##interproscan annotation
  mkdir ${name}_ipr
  cd ${name}_ipr
  ln -s ../$name.gene.protein.primarySeq .
  perl /data/02_work/yangm/01_script/blast-spliting-blast.pl $name.gene.protein.primarySeq aa 10 F
  
  for i in $(seq 1 10)
  do
    nohup /data/02_work/yangm/02_software/01_package/interproscan-5.57-90.0/interproscan.sh -cpu 8 -goterms -i $name.gene.protein.primarySeq.$i -dp > $name.ipr.log 2>$name.ipr.err.log &
  done
  
  while ps -u yangm | grep -q "interproscan"
  do
    sleep 60
    continue
  done
  
  cat *.tsv > all
  rm -rf $name.gene.protein.primarySeq*
  perl $script/fill_0_num.pl all all.1
  sort -k1,1 all.1 > $name.ipr.tsv
  rm -rf all all.1
  perl -ane '@tmp=$_=~/GO:\d+/g;foreach $tmp(@tmp){push @{$hash{$F[0]}}, $tmp;} END {foreach $key(keys %hash){$go=join " ", @{$hash{$key}}; print "$key\t$go\n";}}' $name.ipr.tsv > $name.ipr.tsv.1.go
  perl -ane '@tmp=$_=~/GO:\d+/g;foreach $tmp(@tmp){push @{$hash{$F[0]}}, $tmp;} END {foreach $key(keys %hash){foreach $go(@{$hash{$key}}) {print "$key\t$go\n";}}}' $name.ipr.tsv > $name.ipr.tsv.2.go
  cd ..

BLOCK

  cd ../../
  date
  echo "Done extract and interproscan annotation, pipeline finished"

done

date
echo "End predicting $name....."
echo "####################"
