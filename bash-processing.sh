#!/bin/bash
## example command ./processing.sh -i "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/all-RNAseqs" -c "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/all-RNAseqs-copy" -s "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/rmats-cohort-lists/all.RNAseqs.txt"

set -euo pipefail

while getopts i:c:s: flag
do 
    case "${flag}" in
        i) initdir=${OPTARG};;
        c) cpdir=${OPTARG};;
        s) samples=${OPTARG};;
    esac
done

#INIT_DIR="/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/all-RNAseqs"
#CP_DIR="/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/all-RNAseqs-copy"
#SAMPLE_IDS="/projects/rmorin/projects/gambl-repos/gambl-kcoyle/rmats-cohort-lists/all.RNAseqs.txt"

NEW_IDS_TSV=$cpdir"/sample_ids.txt"


file_pattern=".MATS.JC.txt"

#cp -r $initdir $cpdir

#cat $samples | perl -ne 'chomp; /([a-zA-Z0-9._\-]*)\.bam/g; print "$1\n"' | cut -d"." -f1 | tr '\n' '\t' > $cpdir/sample_ids.txt

events="A3SS A5SS SE RI"


for EVENT in ${events}
do
    echo $EVENT
    file=${cpdir}/${EVENT}${file_pattern}

    echo $file.IJC
    ## Take IJC data & store in new file
    cat $NEW_IDS_TSV  | perl -ne 's/\t/_IJC\t/g;print;'| perl -ne 's/\t$//g;print;' > $file.IJC
    echo "" >> $file.IJC
    WC=$(wc -l $file | awk '{print $1}')
    awk 'BEGIN{FS="\t"} NR>1 {print $13}' $file | tr ',' '\t' >> $file.IJC 
    WC_IJC=$(wc -l $file.IJC | awk '{print $1}')
    if [ "$WC" -eq "$WC_IJC" ]; then
        echo "File is OK"
    else
        echo "Original $WC is not equal to $WC_IJC in $file.IJC"
    fi

    ## Take SJC data & store in new file
    echo $file.SJC
    cat $NEW_IDS_TSV  | perl -ne 's/\t/_SJC\t/g;print;' | perl -ne 's/\t$//g;print;' > $file.SJC
    echo "" >> $file.SJC
    WC=$(wc -l $file | awk '{print $1}')
    awk 'BEGIN{FS="\t"} NR>1 {print $14}' $file | tr ',' '\t' >> $file.SJC 
    WC_SJC=$(wc -l $file.SJC | awk '{print $1}')
    if [ "$WC" -eq "$WC_SJC" ]; then
        echo "File is OK"
    else
        echo "Original $WC is not equal to $WC_SJC in $file.SJC"
    fi

    ## Take IncLevel data & store in new file
    echo $file.Inc
    cat $NEW_IDS_TSV  | perl -ne 's/\t/_Inc\t/g;print;' | perl -ne 's/\t$//g;print;' > $file.Inc
    echo "" >> $file.Inc
    WC=$(wc -l $file | awk '{print $1}')
    awk 'BEGIN{FS="\t"} NR>1 {print $21}' $file | tr ',' '\t' >> $file.Inc 
    WC_INC=$(wc -l $file.Inc | awk '{print $1}')
    if [ "$WC" -eq "$WC_INC" ]; then
        echo "File is OK"
    else
        echo "Original $WC is not equal to $WC_INC in $file.Inc"
    fi
    
    ## Take Identifiers & join & store in new file
    echo $file.ID
    echo "Event" > $file.ID
    WC=$(wc -l $file | awk '{print $1}')
    awk -v e=$EVENT 'NR>1 {print e"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11}' $file | perl -ne 's/\t$//g;print;' >> $file.ID
    WC_ID=$(wc -l $file.ID | awk '{print $1}')
    if [ "$WC" -eq "$WC_ID" ]; then
        echo "File is OK"
    else
        echo "Original $WC is not equal to $WC_ID in $file.ID"
    fi

    ## Join all files together with paste
    paste $file.ID $file.IJC $file.SJC $file.Inc |  perl -ne 's/\t$//g;print;' > $cpdir/$EVENT.final.txt
    
done
