#!/bin/bash

if [ "$#" -lt 1 ]
then
  echo "usage: rpS3_trckr.sh <list>"
  echo "where <list> is a /n delimited list of paths to prodigal-predicted gene files
  echo "output folder is current folder"
  exit 1
fi

# path to folder with HMMs
hmm=${SCRATCH}/rpS3_abundance_1_AHGGTNBCX2/rpS3_Diamond2019.hmm

# Input variable
die_liste=$1

## ++++++++++++++ FUNCTIONS ++++++++++++++ ##
# 1. Extract S3
extr_rp(){
  
  # retrive sample name from file
  sample_name=$(ls $sample | awk -F \/ '{print$NF}')
  echo -e "... working on sample ${sample_name} ..."
  hmmsearch --tblout /dev/stdout -o /dev/null --cpu 1 --notextw ${hmm} ${sample} | grep -v "^#" > ${sample_name}.hmm_results.txt
  awk '{print$1}' ${sample_name}.hmm_results.txt >> ${sample_name}.all_hits.txt
  pullseq -i ${sample} -n ${sample_name}.all_hits.txt > ${sample_name}.all_hits.fasta
  awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ${sample_name}.all_hits.fasta > ${sample_name}.$domain.all_hits_woLB.fasta
  for gene in $(cat ${sample_name}.all_hits.txt); do grep -w -A 1 $gene ${sample_name}.all_hits_woLB.fasta | sed 1d | wc; done | awk '{print$3}' > ${sample_name}.$domain.length.tmp
  awk '{print$6}' ${sample_name}.hmm_results.txt > ${sample_name}.scores.tmp
  paste ${sample_name}.all_hits.txt ${sample_name}.scores.tmp ${sample_name}.length.tmp > ${sample_name}_hits.txt

  # filter tables by score and length
  awk '{ if ($2 >= 40) print $1 }' ${sample_name}_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 > ${sample_name}_final_rpS3.ids
  
  # retrieve final set of rpS3 genes for sample
  cat ${sample_name}_final_rpS3.ids | pullseq -i ${sample} -N > ${sample_name}_final_rpS3.fa
  # clean up
  rm ${sample_name}.*.all_hits.fasta ${sample_name}.*.all_hits_woLB.fasta ${sample_name}.*.scores.tmp ${sample_name}.*.length.tmp ${sample_name}*_hits.txt ${sample_name}_selected1 ${sample_name}_selected2 ${sample_name}_selected3 ${sample_name}_crit1 ${sample_name}_final_rpS3.ids #${sample_name}.*.all_hits.txt
  echo -e "... done with ${sample_name} ..."
}

# 2. Extract scaffolds
extr_scaff(){
  scaff_file=$(ls $sample | sed "s/\.faa/.fa/g")
  pullseq -i ${scaff_file} -n all_rpS3_scaffold.ids >> all_rpS3_scaffold.fna
}

# 3. Extract clusters
extr_cluster(){
  grep "^>" ${cluster} | awk '{print$1}' | awk -F_ 'NF{NF-=1};1' | sed "s/\ /_/g" | sed "s/>//g" | pullseq -i ../all_rpS3_scaffold.fna -N > ${cluster}_scaffolds.fna
  usearch -sortbylength ${cluster}_scaffolds.fna -fastaout ${cluster}_scaffolds.sorted.faa
  grep -m 1 ">" ${cluster}_scaffolds.sorted.faa | sed "s/>//g" | pullseq -i ${cluster}_scaffolds.sorted.faa -N >> ../final_rpS3_scaffolds.fna
}

## ++++++++++++++ END OF FUNCTIONS ++++++++++++++ ##

# run rpS3 search
for sample in $(cat $die_liste); do
  extr_rp
done

# determine the best clustering parameter and run the culstering
cat *_final_rpS3.fa > all_rpS3.fa
mkdir -p rpS3_clusters
usearch -sortbylength all_rpS3.fa -fastaout all_rpS3_sorted.faa
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' all_rpS3.fa | awk 
max=$(bioawk -c fastx '{print ">" $name ORS length($seq)}' CS2RZSR1G_contigs.fa | grep -vE '>' | sort | tail -n1)
min=$(bioawk -c fastx '{print ">" $name ORS length($seq)}' CS2RZSR1G_contigs.fa | grep -vE '>' | sort | head -n1)
var=$(echo "scale=3; $min/$max" | bc)
usearch -cluster_fast all_rpS3_sorted.faa -id 0.99 -centroids all_rps3_centroids.faa -uc all_rps3_99_clusters.uc -clusters rpS3_clusters/rpS3_ -query_cov $var -target_cov $var

# extract scaffolds
echo -e "... extracting scaffolds with rpS3 genes ..."
grep ">" all_rpS3.fa | sed "s/>//g" | awk '{print$1}' | awk -F_ 'NF{NF-=1};1' | sed "s/\ /_/g" > all_rpS3_scaffold.ids
for sample in $(cat $die_liste); do
  extr_scaff
done
rm all_rpS3_scaffold.ids

# extract scaffolds of clusters
cd rpS3_clusters
ls -1 rpS3_* > all_genes.txt
split all_genes.txt -l 6 batch.
for batch in $(ls -1 batch.*); do
  for cluster in $(cat $batch); do
    extr_cluster
  done
  wait
done
rm all_genes.txt batch.*
cd ..


echo -e "... all done! output file is final_rpS3_scaffolds.fna"


# store some intermediate files
mkdir -p intermediate_files
mv *_final_rpS3.fa intermediate_files/.
mv *hmm_results.txt intermediate_files/.
mv all_rpS3_sorted.faa intermediate_files/.
mv all_rps3_centroids.faa intermediate_files/.
mv all_rps3_99_clusters.uc intermediate_files/.
mv all_rpS3_scaffold.fna intermediate_files/.
