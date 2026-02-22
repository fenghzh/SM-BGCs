#!/bin/bash
set -e

database_path=/home/u220220933040/BGC/diamond_database/MIBIG_database.dmnd

cat ./SRA_ID_10.txt | while read oneline
do

    #oneline=(${oneline//'\n'/})
    echo $oneline
    prefetch  $oneline --location NCBI
    sra_data=$(ls ./ | grep "^$oneline.*lite" )
    fastq-dump --split-3 $sra_data
    num=$(ls ./ | grep "^$oneline.*.fastq" | wc -l)
    if ((num ==1)); then
        a=$(ls ./ | grep "^$oneline.*.fastq" )
        b=(${a//fastq/1P.fastq })
        trimmomatic SE $a $b ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 -threads 40
    fi
    if ((num ==2)); then
        a=$(ls ./ | grep "^$oneline.*1.fastq" )
        b=$(ls ./ | grep "^$oneline.*2.fastq" )
        filter1P=(${a//1.fastq/1P.fastq })
        filter1U=(${a//1.fastq/1U.fastq })
        filter2P=(${a//1.fastq/2P.fastq })
        filter2U=(${a//1.fastq/2U.fastq })
        echo $a $b  $filter1P $filter1U $filter2P $filter2U
        trimmomatic PE -phred33 $a $b $filter1P $filter1U $filter2P $filter2U ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 -threads 40 
    fi
    fastq_filter=$(ls ./ | grep "^$oneline.*1P.fastq" )
    blast_result=(${fastq_filter//1P.fastq/_diamond.m8})
    diamond blastx -d $database_path -q $fastq_filter -o $blast_result -k 1 -p 40

    phlan_out=(${fastq_filter//1P.fastq/_metaphlan.txt})
    bowtie_out=(${fastq_filter//1P.fastq/.bowtie2.bz2})
    metaphlan $fastq_filter --bowtie2out $bowtie_out  --input_type fastq  --nproc 40  -o  $phlan_out  --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /home/u220220933040/BGC/metaphlan_database/database
    rm $a $b $filter1P $filter1U $filter2P $filter2U $sra_data

done
