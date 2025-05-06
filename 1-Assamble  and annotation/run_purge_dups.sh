#!/bin/sh

input=$1

cat $input | while read name
  do
     if [ ! -d "$DIRECTORY" ]; then
         mkdir $name	      
     fi

     cd $i
     ln -sf ../$name.contig.fas .
     ln -sf ../ccs.fofn .
     python /work/purge_dups-master/scripts/pd_config.py -n config.josn $name.contig.fas ccs.fofn
     python /work/purge_dups-master/scripts/run_purge_dups.py -p bash config.josn /work/purge_dups-master/bin $name
     cd ..
  done
