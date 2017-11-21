

for SAMPLE_NAME in `ls -1 /home/yk42w/project/yk42w/data/RNAseq/eBLSamples/|grep eBL`;
do
ls -al /home/yk42w/project/yk42w/data/RNAseq/eBLSamples/$SAMPLE_NAME;

for type in 2;
do

bsub -q long -n 4 -R rusage[mem=20000] -R "select[tmp>1000]" -W 03:00 \
-e "$SAMPLE_NAME"_"$type"_err.%J.txt -o "$SAMPLE_NAME"_"$type"_out.%J.txt \
~/codes/ViralGenomeAssembly/scripts/kallisto-exp-quantify.sh \
"$SAMPLE_NAME"_"$type"_TotalRNA_EBV stranded 4 \
~/project/yk42w/data/RNAseq/eBLSamples/$SAMPLE_NAME/"$SAMPLE_NAME"_cutAdpted_1.fastq.gz \
~/project/yk42w/data/RNAseq/eBLSamples/$SAMPLE_NAME/"$SAMPLE_NAME"_cutAdpted_2.fastq.gz \
EBV $type

sleep 1;

done
done

for SAMPLE_NAME in `ls -1 /home/yk42w/project/yk42w/data/RNAseq/eBLSamples/PolyA`;
do
  for type in 1 2;
  do
    bsub -q long -n 4 -R rusage[mem=20000] -R "select[tmp>1000]" -W 03:00 \
    -e "$SAMPLE_NAME"_"$type"_err.%J.txt \
    -o "$SAMPLE_NAME"_"$type"_out.%J.txt \
    ~/codes/ViralGenomeAssembly/scripts/kallisto-exp-quantify.sh \
    "$SAMPLE_NAME"_"$type"_PolyA-RNA_EBV stranded 4 \
    ~/project/yk42w/data/RNAseq/eBLSamples/PolyA/$SAMPLE_NAME/eBLRNAseqBatch3_"$SAMPLE_NAME"_R1.fastq.gz \
    ~/project/yk42w/data/RNAseq/eBLSamples/PolyA/$SAMPLE_NAME/eBLRNAseqBatch3_"$SAMPLE_NAME"_R2.fastq.gz EBV $type;

    sleep 1; \
  done
done

for SAMPLE_NAME in BL717_8-1-16;
do
  for type in 2;
  do
    bsub -q long -n 4 -R rusage[mem=20000] -R "select[tmp>1000]" -W 03:00 \
    -e "$SAMPLE_NAME"_"$type"_err.%J.txt -o "$SAMPLE_NAME"_"$type"_out.%J.txt \
    ~/codes/ViralGenomeAssembly/scripts/kallisto-exp-quantify.sh \
    "$SAMPLE_NAME"_"$type"_TotalRNA_EBV stranded 4 \
    /home/yk42w/project/yk42w/data/RNAseq/FL-BL-Batch-1/07OCT16_HFVF3BB/BL717_8-1-16_CTCTAGAT_S53_L007_R1_001.fastq.gz \
    /home/yk42w/project/yk42w/data/RNAseq/FL-BL-Batch-1/07OCT16_HFVF3BB/BL717_8-1-16_CTCTAGAT_S53_L007_R2_001.fastq.gz EBV $type;
    sleep 1;
  done
done
