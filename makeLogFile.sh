#! /bin/bash

META=$1

OUTPUT_DIR=`grep "Output Directory:" $META | awk '{print $3}'`

## Determine what log file should be called
n=1
LOG=""
while [ -n $LOG ]; do
    if [ ! -e ${OUTPUT_DIR}pipeline.$n.log ]; then
        LOG=${OUTPUT_DIR}pipeline.$n.log
        break
    else
        let n++
    fi
done

eval batches=(`awk '{ if ($2 == "gencall" || $2 == "zcall" || $2 == "qc") print $1 }' $META`)
echo $batches >> $LOG
for b in $batches; do
    echo $b >> $LOG
done

#while read line; do
#done <$META

## Initiate log file
echo "------------------------------" > $LOG
echo "Exome Chip QC and Analysis Pipeline v3.0" >> $LOG
echo "Jackie Goldstein [jigold@broadinstitute.org]" >> $LOG
echo "" >> $LOG
echo "META file used: $META" >> $LOG
echo "" >> $LOG

date=`date +"%m/%d/%Y %T"`
echo "Analysis started at $date" >> $LOG
echo "" >> $LOG
echo "Output Directory: $OUTPUT_DIR" >> $LOG
echo "------------------------------" >> $LOG
echo "" >> $LOG



# Return log file name back to main.sh
echo $LOG
