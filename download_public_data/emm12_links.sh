#!/bin/bash

set -euo pipefail

while read line; do ffq --ftp ${line} | jq -r '.[] | .url'; done < list_emm12_ABC.txt > list_links_emm12.txt

mkdir -p input

while read line
do
  if [ ${line: -10} == '1.fastq.gz' ]
  then
    NAME=$(basename $line _1.fastq.gz)
    curl -o input/${NAME}_S1_L001_R1_001.fastq.gz $line
  elif [ ${line: -10} == '2.fastq.gz' ]
  then
    NAME=$(basename $line _2.fastq.gz)
    curl -o input/${NAME}_S1_L001_R2_001.fastq.gz $line
  fi
done < list_links_emm12.txt
