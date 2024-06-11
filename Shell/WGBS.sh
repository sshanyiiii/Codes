# 2024.6.10

dir="/home4/ssyi/Mouse_Placenta/Public_data/Gao_WGBS_GSE108711"
##################################################### 0.raw
mkdir 0.raw
cd 0.raw
cat SRA_Acc_List.txt | while read id;     
    do nohup prefetch ${id} --max-size 50G & 
done
for i in */*.sralite;do
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
cd 2.trim_galore
for i in $dir/0.raw/*.R1.fastq.gz
do
    r1=$i
    r2=${i/.R1.fastq.gz/.R2.fastq.gz}
    nohup trim_galore -j 10 --paired -q 20 --length 50 \
    --trim-n --illumina \
    -o $dir/2.trim_galore $r1 $r2 &
done

multiqc -f $dir/2.trim_galore -o $dir/2.trim_galore -n clean_multiqc_report



##################################################### 3.Bismark
mkdir 3.bismark
mkdir -p $dir/3.bismark/lambda

bismark_index="/home/share/bismark_index"
### 比对lambda,作为内参查看bisulfite转化效率
for i in $dir/2.trim_galore/*R1_val_1.fq.gz;
do
    sample=$(basename $i)
    sample=${sample%%.*}
    r1=$i
    r2=${i/R1_val_1.fq.gz/R2_val_2.fq.gz}
    nohup bismark -p 10 -N 1 ${bismark_index}/lambda \
    -1 $r1 -2 $r2 -B ${sample} --unmapped \
    -o $dir/3.bismark/lambda --temp_dir $dir/3.bismark/lambda \
    > $dir/3.bismark/lambda/${sample}.lambda.log 2>&1 &
done
#基本没比对上lambda基因组，只有200个C位点

### 比对mm10
mkdir -p $dir/3.bismark/mm10
for i in $dir/2.trim_galore/*R1_val_1.fq.gz;
do
    sample=$(basename $i)
    sample=${sample%%.*}
    r1=$i
    r2=${i/R1_val_1.fq.gz/R2_val_2.fq.gz}
    nohup bismark -p 10 -N 1 ${bismark_index}/mm10 \
    -1 $r1 -2 $r2 -B ${sample} --unmapped \
    -o $dir/3.bismark/mm10 --temp_dir $dir/3.bismark/mm10 \
    > $dir/3.bismark/mm10/${sample}.mm10.log 2>&1 &
done
#比对mm10比对率好低，17%



##################################################### 4.bismark_dedup，用bismark去除duplicate
mkdir 4.bismark_dedup
### 双端去dup，加-p
for i in $dir/3.bismark/mm10/*_pe.bam; 
do 
    sample=$(basename $i ".bam")
    nohup deduplicate_bismark -p --output_dir $dir/4.bismark_dedup --bam $i \
    > $dir/4.bismark_dedup/$sample\.dedup.log &
done
### 统计去除dup的比例
ls -v *_pe.deduplication_report.txt | while read i;
do
    sample=${i%_pe*}
    sample=${sample/2C_scCOOL_/}
    alignments=$(grep "alignments" ${i} | grep "analysed" | cut -d : -f 2 | awk '$1=$1')
    deduplicated_hit=$(grep "removed" ${i} | grep "%" | cut -d : -f 2|cut -d \( -f 1| awk '$1=$1')
    deduplicated_ratio=$(grep "removed" ${i} | grep "%" | cut -d \( -f 2 | cut -d \) -f 1 | awk '$1=$1')
    leftover=$(grep "leftover" ${i} | grep "%" | cut -d : -f 2|cut -d \( -f 1| awk '$1=$1')
    leftover_ratio=$(grep "leftover" ${i} | grep "%" | cut -d \( -f 2 | cut -d " " -f 1 | awk '$1=$1')
    echo -e ${sample}"\t"${clean_reads}"\t"${deduplicated_hit}"\t"${deduplicated_ratio}"\t"${leftover}"\t"${leftover_ratio} >> PE_deduplicate_stat.txt
done
sed -i '1s/^/sample\talignments\tdeduplicated_hit\tdeduplicated_ratio\tleftover\tleftover_ratio\n/' PE_deduplicate_stat.txt



###### 5.bismark_methylation_extractor，提取deduplicate_bismark去重的bam文件
mkdir 5.bismark_methylation_extractor
cd 5.bismark_methylation_extractor
# https://www.jianshu.com/p/877a5716c24a
# https://liuyujie0136.gitbook.io/sci-tech-notes/bioinformatics/bismark-swdmr
### 双端数据统计，-p
for i in $dir/4.bismark_dedup/*_pe.deduplicated.bam;
do
    sample=$(basename $i ".deduplicated.bam")
    nohup bismark_methylation_extractor -p \
    --parallel 5 --buffer_size 10G \
    --no_overlap \
    --comprehensive \
    --bedGraph \
    --counts \
    --report \
    -cytosine_report \
    --genome_folder /home/share/bismark_index/mm10 \
    -o $dir/5.bismark_methylation_extractor $i \
    > $dir/5.bismark_methylation_extractor/$sample\_met_ext.log 2>&1 &
done
