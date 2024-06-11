# 2024.6.10

dir="/home4/ssyi/Mouse_Placenta/Public_data/ATAC_GSE116854"
##################################################### 0.raw
mkdir 0.raw
cd 0.raw
#download
cat SRA_Acc_List.txt | while read id;     
    do nohup prefetch ${id} --max-size 50G & 
done
for i in */*.sralite;do   # prefetch版本高下载的sralite，没有保存碱基质量信息
    nohup fastq-dump --split-3 --gzip $i &
done
#rename
ls -v *.fastq.gz >> Rename.txt
vi Rename.txt
cat Rename.txt | while read id;
do
    arr=(${id})
    raw_name=${arr[0]}
    new_name=${arr[1]}
    nohup mv ${raw_name} ${new_name} &
done


##################################################### 1.qc
mkdir 1.qc
cd 1.qc
ls -v $dir/0.raw/*.fastq.gz | while read id; do(nohup fastqc -t 10 $id -o $dir/1.qc &) ; done   

multiqc -f $dir/1.qc -o $dir/1.qc -n raw_multiqc_report



##################################################### 2.trim_galore
mkdir 2.trim_galore
for i in $dir/0.raw/*.R1.fastq.gz
do
    r1=$i
    r2=${i/.R1.fastq.gz/.R2.fastq.gz}
    nohup trim_galore -j 10 -q 25 --length 36 --paired --fastqc --trim-n \
    -o $dir/2.trim_galore $r1 $r2 &
done

multiqc -f $dir/2.trim_galore -o $dir/2.trim_galore -n clean_multiqc_report



##################################################### 3.bowtie2，align to mm10
mkdir 3.bowtie2
cd 3.bowtie2

bowtie2_index="/home/share/bowtie2_index/mm10"
for i in $dir/2.trim_galore/*R1_val_1.fq.gz;
do
    sample=$(basename $i)
    sample=${sample%%.*}
    r1=$i
    r2=${i/R1_val_1.fq.gz/R2_val_2.fq.gz}
    nohup bowtie2 -p 10 -x $bowtie2_index \
    --end-to-end --very-sensitive -X 2000 \
    --no-mixed --no-discordant --no-unal \
    --time -1 $r1 -2 $r2 -S $dir/3.bowtie2/$sample.sam > $dir/3.bowtie2/$sample.log 2>&1 &
done

ls -v *.log | while read id;
do
    sample=$(basename $id ".log")
    clean_pairs=$(grep "paired" $id |sed 's/^ *//' |cut -d " " -f 1)
    aligned_rate=$(grep "overall" $id |cut -d " " -f 1 |tr "\n" "\t" |sed "s/\t$/\n/")
    one=$(grep "concordantly" $id |grep "exactly" |sed 's/^ *//' |cut -d " " -f 1)
    one_rate=$(grep "concordantly" $id |grep "exactly" |sed 's/^ *//' |cut -d "(" -f 2|cut -d ")" -f 1)
    onemore=$(grep "concordantly" $id |grep ">" |sed 's/^ *//' |cut -d " " -f 1)
    onemore_rate=$(grep "concordantly" $id |grep ">" |sed 's/^ *//' |cut -d "(" -f 2|cut -d ")" -f 1)
    aligned_pairs=$(expr $one + $onemore)
    echo -e ${sample}"\t"${clean_pairs}"\t"${one}"\t"${one_rate}"\t"${onemore}"\t"${onemore_rate}"\t"${aligned_pairs}"\t"${aligned_rate} >> bt2_align_stat.txt
done
sed -i '1s/^/sample\ttotal_pairs\tproperly_paired\tratio\tmore_1_time\tratio\taligned_pairs\taligned_rate\n/' bt2_align_stat.txt



##################################################### 4.samtools
mkdir 4.samtools
cd 4.samtools
# sam2bam & sort
for i in $dir/3.bowtie2/*sam; do
    sample=${i##*/}  
    sample=${sample%.*}
    nohup samtools sort -@ 10 -O bam -o $dir/4.samtools/$sample.sorted.bam $i &
done
# index
for i in *.sorted.bam; do
    nohup samtools index -@ 10 $i &
done



##################################################### 5.sambamba，去除dup
mkdir 5.sambamba 
cd 5.sambamba
for i in $dir/4.samtools/*.sorted.bam;do
    sample=$(basename $i ".sorted.bam")
    nohup sambamba markdup -r -t 8 -p --overflow-list-size 600000 $i $dir/5.sambamba/${sample}.rmdup.bam > $dir/5.sambamba/${sample}.rmdup.log &
done
# --overflow-list-size 600000针对：sambamba markdup too many open files issue
### samtools flagstats比对结果
for i in $dir/5.sambamba/*.rmdup.bam;do
    sample=$(basename $i ".rmdup.bam")
    nohup samtools flagstat -@ 10 $i > $sample.flagstat &
done
# 统计去完dup后的read pairs
for i in *flagstat;do sample=$(basename $i ".flagstat");
    sample=$(basename $i ".flagstat")
    dedup_pairs=$(grep "read1" $i|cut -f1 -d " ")
    echo -e ${sample}"\t"${dedup_pairs} >> dedup_pairs.stat
done
sed -i '1s/^/sample\tdedup_pairs\n/' dedup_pairs.stat

### merged bam文件
mkdir ./merged_bam
samtools merge ./merged_bam/MII_ATAC.merge.bam MII_ATAC_rep1.rmdup.bam MII_ATAC_rep2.rmdup.bam
for i in *.merge.bam; do
    sample=$(basename $i ".bam")
    nohup samtools sort -@ 10 -O bam -o ./$sample.sorted.bam $i &
done
for i in *.merge.sorted.bam; do
    nohup samtools index -@ 10 $i &
done
# merge后去除dup
for i in ./*.merge.sorted.bam;do
    sample=$(basename $i ".merge.sorted.bam")
    nohup sambamba markdup -r -t 8 -p --overflow-list-size 600000 $i ./${sample}.rmdup.bam > ./${sample}.rmdup.log &
done



##################################################### 6.deeptools
mkdir 6.deeptools
cd 6.deeptools
for i in $dir/5.sambamba/*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup bamCoverage -p 10 -b $i -bs 50 --normalizeUsing RPKM -of "bigwig" -o ./$sample.bw &
done

### merged bam文件
mkdir merged_bam
cd merged_bam
for i in $dir/5.sambamba/merged_bam/*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup bamCoverage -p 10 -b $i -bs 50 --normalizeUsing RPKM -of "bigwig" -o ./$sample.bw &
done
