#!/usr/bin/env bash

if test $# -ne 1
then
    echo "No bam file given."
    exit 1
fi

if ! test -f $BAMFILE
then
   echo "Input file with name $BAMFILE not found."
   exit 1
fi

BAMFILE=$1
SORTEDBAMFILE=$(basename $BAMFILE).sorted
INDEXBAMFILE=$(basename $SORTEDBAMFILE).bai
WORKDIR=sections-$(basename $BAMFILE)-$(date +%Y%m%d%H%M)

echo "> bamfile: $BAMFILE"
echo "> sorted-bamfile: $SORTEDBAMFILE"
echo "> index: $INDEXBAMFILE"
echo "> work-directory: $WORKDIR"
echo ""

pushd .
mkdir -p $WORKDIR
cd $WORKDIR
echo "> Current dir is $(pwd)"
echo "> Sorting $BAMFILE to $SORTEDBAMFILE"
samtools sort $BAMFILE $SORTEDBAMFILE
echo "> Indexing $SORTEDBAMFILE"
samtools index $SORTEDBAMFILE.bam
# Not that good...
samtools view -H $SORTEDBAMFILE.bam | grep @SQ | sed 's/@SQ\t//g' - | sed 's/SN://g' - | sed 's/LN://g' - | sed 's/\t/ /g' - | sort -n > regions.txt
echo "> Splitting sections"
while read line
do
    # echo "> Processing line $line"
    values=( $line )
    ID=${values[0]}
    MAXLEN=${values[1]}
    # echo "ID=$ID MAXLEN=$MAXLEN"
    for (( STARTP=0; STARTP<$MAXLEN; STARTP+=1000000 ))
    do
	OVERLAPSTARTP=$(($(($STARTP - 200000)) < 0 ? 0 : $(($STARTP - 200000))))
        ENDP=$(($(($STARTP + 1000000 -1)) < $MAXLEN ? $(($STARTP + 1000000 -1)) : $MAXLEN))
        echo "$ID:$OVERLAPSTARTP-$ENDP" >> sections.txt
    done
done < regions.txt
# ln -s $BAMFILE .
echo "> splitting reads"
while read line
do
    values=( `echo $line | tr ":" " "` )
    ID=${values[0]}
    SECT=${values[1]}
    START=( `echo $SECT | tr -d " "` )
    mkdir -p chr$ID
    samtools view -hb $SORTEDBAMFILE.bam $line > "chr$ID/chr$ID-$START.bam"
done < sections.txt
popd

exit 0
