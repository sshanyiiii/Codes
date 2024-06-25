MED6_CUTTAG="/home3/ssyi/Embryo_Med6/20231108_CUT_TAG"
##################################################### 1.qc
mkdir 1.qc
ls -v $MED6_CUTTAG/0.0.raw/*.fastq.gz | while read id; do(nohup fastqc -t 10 $id -o $MED6_CUTTAG/1.qc &) ; done

multiqc -f $MED6_CUTTAG/1.qc -o $MED6_CUTTAG/1.qc -n raw_multiqc_report


##################################################### 2.trim_glore
mkdir 2.trim_galore
for i in $MED6_CUTTAG/0.0.raw/*.R1.fastq.gz
do
    r1=$i
    r2=${i/.R1.fastq.gz/.R2.fastq.gz}
    nohup trim_galore -j 10 -q 25 --length 50 --trim-n --paired --fastqc \
    -o $MED6_CUTTAG/2.trim_galore $r1 $r2 &
done

multiqc $MED6_CUTTAG/2.trim_galore -o $MED6_CUTTAG/2.trim_galore -n clean_multiqc_report


##################################################### 3.bowtie2
mkdir 3.bowtie2
bowtie2_index="/home/share/bowtie2_index/mm10"
for i in $MED6_CUTTAG/2.trim_galore/*R1_val_1.fq.gz;
do
    sample=$(basename $i)
    sample=${sample%%.*}
    r1=$i
    r2=${i/R1_val_1.fq.gz/R2_val_2.fq.gz}
    nohup bowtie2 -p 10 -x $bowtie2_index \
    --end-to-end --very-sensitive \
    -I 10 -X 700 \
    --no-mixed --no-discordant --no-unal \
    --time -1 $r1 -2 $r2 -S $MED6_CUTTAG/3.bowtie2/$sample.sam > $MED6_CUTTAG/3.bowtie2/$sample.log 2>&1 &
done
# 统计比对率
ls -v *.log | while read id;
do
    sample=$(basename $id ".log")
    rate=$(grep "overall" $id |cut -d " " -f 1 |tr "\n" "\t" |sed "s/\t$/\n/")
    echo -e ${sample}"\t"${rate} >> bt2_align_rate.txt
done
sed -i '1s/^/sample\taligned_rate\n/' bt2_align_rate.txt
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
# sam2bam & sort (thread_control_1.sh)
for i in $MED6_CUTTAG/3.bowtie2/*sam; do
    sample=${i##*/}  
    sample=${sample%.*}
    nohup samtools sort -@ 8 -O bam -o $MED6_CUTTAG/4.samtools/$sample.sorted.bam $i &
done
# index
for i in *.sorted.bam; do
    nohup samtools index -@ 10 $i &
done
# 把H3K27me3加测的和前面的merge到一起
mkdir batch_H3K27me3
nohup samtools merge ../W4_L2C_ctrl_H3K27me3_rep1.sorted.bam W4_L2C_ctrl_H3K27me3_rep1.sorted.bam W4_L2C_ctrl_H3K27me3_rep1_plus.sorted.bam &
nohup samtools merge ../W4_L2C_ctrl_H3K27me3_rep2.sorted.bam W4_L2C_ctrl_H3K27me3_rep2.sorted.bam W4_L2C_ctrl_H3K27me3_rep2_plus.sorted.bam &
nohup samtools merge ../W4_L2C_Zp3_Cre_H3K27me3_rep1.sorted.bam W4_L2C_Zp3_Cre_H3K27me3_rep1.sorted.bam W4_L2C_Zp3_Cre_H3K27me3_rep1_plus.sorted.bam &
nohup samtools merge ../W4_L2C_Zp3_Cre_H3K27me3_rep2.sorted.bam W4_L2C_Zp3_Cre_H3K27me3_rep2.sorted.bam W4_L2C_Zp3_Cre_H3K27me3_rep2_plus.sorted.bam &


##################################################### 5.sambamba，去除dup
mkdir 5.sambamba 
for i in $MED6_CUTTAG/4.samtools/*.sorted.bam;do
    sample=$(basename $i ".sorted.bam")
    nohup sambamba markdup -r -t 8 -p $i $MED6_CUTTAG/5.sambamba/${sample}.rmdup.bam > $MED6_CUTTAG/5.sambamba/${sample}.rmdup.log &
done
# samtools flagstats比对结果
for i in $MED6_CUTTAG/5.sambamba/*.rmdup.bam; do
    sample=$(basename $i ".rmdup.bam")
    nohup samtools flagstat -@ 5 $i > $sample.flagstat &
done
for i in *.flagstat;do
	sample=$(basename $i ".flagstat")
	dedup_pairs=$(grep "read1" $i | cut -d ' ' -f1)
	echo -e ${sample}"\t"${dedup_pairs} >> flagstat.txt
done
sed -i '1s/^/sample\tdedup_pairs\n/' flagstat.txt
# 建index
for i in *.rmdup.bam; do
    nohup samtools index -@ 10 $i &
done
### 样本相关性较好，把rep merge到一起
mkdir merged_bam
nohup samtools merge ./merged_bam/W4_L2C_ctrl_H3K4me3_rep12.rmdup.bam W4_L2C_ctrl_H3K4me3_rep1.rmdup.bam W4_L2C_ctrl_H3K4me3_rep2.rmdup.bam &
nohup samtools merge ./merged_bam/W4_L2C_Zp3_Cre_H3K4me3_rep12.rmdup.bam W4_L2C_Zp3_Cre_H3K4me3_rep1.rmdup.bam W4_L2C_Zp3_Cre_H3K4me3_rep2.rmdup.bam &
nohup samtools merge ./merged_bam/W4_L2C_ctrl_H3K4me1_rep12.rmdup.bam W4_L2C_ctrl_H3K4me1_rep1.rmdup.bam W4_L2C_ctrl_H3K4me1_rep2.rmdup.bam &
nohup samtools merge ./merged_bam/W4_L2C_Zp3_Cre_H3K4me1_rep12.rmdup.bam W4_L2C_Zp3_Cre_H3K4me1_rep1.rmdup.bam W4_L2C_Zp3_Cre_H3K4me1_rep2.rmdup.bam &
nohup samtools merge ./merged_bam/W4_L2C_ctrl_H3K27ac_rep12.rmdup.bam W4_L2C_ctrl_H3K27ac_rep1.rmdup.bam W4_L2C_ctrl_H3K27ac_rep2.rmdup.bam &
nohup samtools merge ./merged_bam/W4_L2C_Zp3_Cre_H3K27ac_rep12.rmdup.bam W4_L2C_Zp3_Cre_H3K27ac_rep1.rmdup.bam W4_L2C_Zp3_Cre_H3K27ac_rep2.rmdup.bam &


##################################################### 6.deeptools
mkdir 6.deeptools
for i in $MED6_CUTTAG/5.sambamba/*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup bamCoverage  -p 10 -b $i -bs 50 --normalizeUsing RPKM -of "bigwig" -o ./$sample.bw &
done
### 相关性/PCA
ls -v *.bw | while read file;
do
    sample=${file%.bw*}
    sample=${sample/W4_L2C_/}
    sample=${sample/Zp3_/}
    echo -e ${file}"\t"${sample} >> bw_files.txt
done
nohup multiBigwigSummary bins \
    -p 10 \
    --bwfiles $(cat bw_files.txt | cut -f 1 | tr '\n' ' ') \
    --labels $(cat bw_files.txt | cut -f 2 | tr '\n' ' ') \
    --outFileName W4_L2C_HM_signal_10kbin. \
    --outRawCounts W4_L2C_HM_signal_10kbin.tab &
nohup plotCorrelation \
    -in W4_L2C_HM_signal_10kbin.npz \
    --corMethod spearman \
    --skipZeros \
    --plotTitle "Spearman Correlation" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_SpearmanCorr_HM_10kbin.pdf \
    --removeOutliers &
nohup plotCorrelation \
    -in W4_L2C_HM_signal_10kbin.npz  \
    --corMethod pearson \
    --skipZeros \
    --plotTitle "Pearson Correlation" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    --removeOutliers -o heatmap_PearsonCorr_HM_10kbin.pdf &
nohup plotPCA \
    -in W4_L2C_HM_signal_10kbin.npz \
    --transpose \     # 记得加这个
    --plotHeight 15 \
    --plotWidth 10\
    --plotTitle "PCA on 10kb bin" \
    -o PCA_HM_10kbin.pdf &   
# H3K4me3, H3K4me1在5most TSS/promoter上的信号分布
nohup computeMatrix reference-point \
    -p 10 \
    --referencePoint TSS \
    -b 2000 -a 2000 \
    -bs 50 \
    -R /home1/ssyi/annotate/mm/NCBI_RefSeq/5most_mm10_gene.bed \
    -S W4_L2C_ctrl_H3K4me3_rep1.bw W4_L2C_ctrl_H3K4me3_rep2.bw W4_L2C_Zp3_Cre_H3K4me3_rep1.bw W4_L2C_Zp3_Cre_H3K4me3_rep2.bw \
    --samplesLabel ctrl_rep1 ctrl_rep2 Zp3_cre_rep1 Zp3_cre_rep2 \
    --skipZeros \
    -o matrix_H3K4me3_TSS.gz \
    --outFileSortedRegions regions_H3K4me3_TSS.bed &
nohup plotProfile -m matrix_H3K4me3_TSS.gz \
    -out H3K4me3_TSS.pdf \
    --perGroup \
    --yAxisLabel "H3K4me3 signal" \
    --plotTitle "H3K4me3 on pomoter (+-2kb)" \
    --plotHeight 10 \
    --plotWidth 15 \
    # --colors "#fdbb84" "#2b8cbe" \
    --colors "#b2182b" "#ef8a62" "#2166ac" "#67a9cf" \
    --plotFileFormat pdf --dpi 720 &
# scale-region
nohup computeMatrix scale-regions \
    -p 10 \
    -b 3000 -a 3000 \
    -bs 50 \
    -R /home1/ssyi/annotate/mm/NCBI_RefSeq/5most_mm10_gene.bed \
    --regionBodyLength 5000 \
    -S W4_L2C_ctrl_H3K4me3_rep1.bw W4_L2C_ctrl_H3K4me3_rep2.bw W4_L2C_Zp3_Cre_H3K4me3_rep1.bw W4_L2C_Zp3_Cre_H3K4me3_rep2.bw \
    --samplesLabel ctrl_rep1 ctrl_rep2 Zp3_cre_rep1 Zp3_cre_rep2 \
    --skipZeros \
    -o matrix_H3K4me3_scale_region.gz \
    --outFileSortedRegions regions_H3K4me3_scale_region.bed &
nohup plotProfile -m matrix_H3K4me3_scale_region.gz \
    -out H3K4me3_scale_region.pdf \
    --perGroup \
    --startLabel "TSS" \
    --endLabel "TES" \
    --yAxisLabel "H3K4me3 signal" \
    --plotTitle "H3K4me3 on genebody" \
    --plotHeight 10 \
    --plotWidth 15 \
    # --colors "#fdbb84" "#2b8cbe" \
    --colors "#b2182b" "#ef8a62" "#2166ac" "#67a9cf" \
    --plotFileFormat pdf --dpi 720 &


##################################################### 7.masc2
mkdir 7.macs2
for i in $MED6_CUTTAG/5.sambamba/*H3K4me3*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup macs2 callpeak -t $i -f BAMPE -g mm --keep-dup all \
    -n ${sample} --outdir $MED6_CUTTAG/7.macs2 > $MED6_CUTTAG/7.macs2/$sample.log &
done
for i in $MED6_CUTTAG/5.sambamba/*H3K4me1*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup macs2 callpeak -t $i -f BAMPE -g mm --broad --keep-dup all \
    -n ${sample} --outdir $MED6_CUTTAG/7.macs2 > $MED6_CUTTAG/7.macs2/$sample.log &
done
for i in $MED6_CUTTAG/5.sambamba/*H3K27me3*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup macs2 callpeak -t $i -f BAMPE -g mm --broad --keep-dup all \
    -n ${sample} --outdir $MED6_CUTTAG/7.macs2 > $MED6_CUTTAG/7.macs2/$sample.log &
done
for i in $MED6_CUTTAG/5.sambamba/*H3K27ac*.rmdup.bam;
do
    sample=$(basename $i ".rmdup.bam")
    nohup macs2 callpeak -t $i -f BAMPE -g mm --keep-dup all \
    -n ${sample} --outdir $MED6_CUTTAG/7.macs2 > $MED6_CUTTAG/7.macs2/$sample.log &
done
### rep之间的peak取交集
mkdir intersect_peak
bedtools intersect -a W4_L2C_ctrl_H3K27ac_rep1_peaks.narrowPeak -b W4_L2C_ctrl_H3K27ac_rep2_peaks.narrowPeak > intersect_peak/W4_L2C_ctrl_H3K27ac_overlap_rep12_peaks.bed
bedtools intersect -a W4_L2C_Zp3_Cre_H3K27ac_rep1_peaks.narrowPeak -b W4_L2C_Zp3_Cre_H3K27ac_rep2_peaks.narrowPeak > intersect_peak/W4_L2C_Zp3_Cre_H3K27ac_overlap_rep12_peaks.bed
bedtools intersect -a W4_L2C_ctrl_H3K4me1_rep1_peaks.broadPeak -b W4_L2C_ctrl_H3K4me1_rep2_peaks.broadPeak > intersect_peak/W4_L2C_ctrl_H3K4me1_overlap_rep12_peaks.bed
bedtools intersect -a W4_L2C_Zp3_Cre_H3K4me1_rep1_peaks.broadPeak -b W4_L2C_Zp3_Cre_H3K4me1_rep2_peaks.broadPeak > intersect_peak/W4_L2C_Zp3_Cre_H3K4me1_overlap_rep12_peaks.bed

### 合并merged bam ctrl和Zp3_Cre组的peak
mkdir merged_peak
cd merged_peak
cut -f1,2,3 ../W4_L2C_ctrl_H3K4me1_rep12_peaks.broadPeak > W4_L2C_ctrl_H3K4me1_rep12_peaks.bed
cut -f1,2,3 ../W4_L2C_Zp3_Cre_H3K4me1_rep12_peaks.broadPeak > W4_L2C_Zp3_Cre_H3K4me1_rep12_peaks.bed
cut -f1,2,3 ../W4_L2C_ctrl_H3K4me3_rep12_peaks.narrowPeak > W4_L2C_ctrl_H3K4me3_rep12_peaks.bed
cut -f1,2,3 ../W4_L2C_Zp3_Cre_H3K4me3_rep12_peaks.narrowPeak > W4_L2C_Zp3_Cre_H3K4me3_rep12_peaks.bed
# peak cat到同一个文件
cat W4_L2C_ctrl_H3K4me1_rep12_peaks.bed W4_L2C_Zp3_Cre_H3K4me1_rep12_peaks.bed > cat_ctrl_Zp3_Cre_H3K4me1_rep12_peaks.bed
cat W4_L2C_ctrl_H3K4me3_rep12_peaks.bed W4_L2C_Zp3_Cre_H3K4me3_rep12_peaks.bed > cat_ctrl_Zp3_Cre_H3K4me3_rep12_peaks.bed
# peak排序
sort -k1,1 -k2,2n cat_ctrl_Zp3_Cre_H3K4me1_rep12_peaks.bed > sort_cat_ctrl_Zp3_Cre_H3K4me1_rep12_peaks.bed
sort -k1,1 -k2,2n cat_ctrl_Zp3_Cre_H3K4me3_rep12_peaks.bed > sort_cat_ctrl_Zp3_Cre_H3K4me3_rep12_peaks.bed
# 合并peak坐标
bedtools merge -i sort_cat_ctrl_Zp3_Cre_H3K4me1_rep12_peaks.bed > merge_cat_ctrl_Zp3_Cre_H3K4me1_rep12_peaks.bed
bedtools merge -i sort_cat_ctrl_Zp3_Cre_H3K4me3_rep12_peaks.bed > merge_cat_ctrl_Zp3_Cre_H3K4me3_rep12_peaks.bed
# 计算所有rep在合并peak上的信号值
nohup computeMatrix scale-regions \
    -p 10 \
    -b 5000 -a 5000 \
    --regionBodyLength 5000 \
    -bs 50 \
    -R merge_cat_ctrl_Zp3_Cre_H3K4me3_rep12_peaks.bed \
    -S $MED6_CUTTAG/6.deeptools/W4_L2C_ctrl_H3K4me3_rep1.bw $MED6_CUTTAG/6.deeptools/W4_L2C_ctrl_H3K4me3_rep2.bw $MED6_CUTTAG/6.deeptools/W4_L2C_Zp3_Cre_H3K4me3_rep1.bw $MED6_CUTTAG/6.deeptools/W4_L2C_Zp3_Cre_H3K4me3_rep2.bw \
    --samplesLabel ctrl_rep1 ctrl_rep2 Zp3_cre_rep1 Zp3_cre_rep2 \
    --skipZeros \
    -o matrix_merged_peak_H3K4me3_scale_region.gz \
    --outFileSortedRegions merged_peak_H3K4me3_scale_region.bed &
nohup plotProfile -m matrix_merged_peak_H3K4me3_scale_region.gz \
    -out merged_peak_H3K4me3_scale_region.pdf \
    --startLabel "Start" \
    --endLabel "End" \
    --yAxisLabel "H3K4me3 signal" \
    --plotTitle "H3K4me3 on merged peak" \
    --plotHeight 10 \
    --plotWidth 15 \
    # --colors "#fdbb84" "#2b8cbe" \
    --colors "#b2182b" "#ef8a62" "#2166ac" "#67a9cf" \
    --perGroup \
    --plotFileFormat pdf --dpi 720 &


##################################################### 8.peak_anno
mkdir 8.peak_anno
for i in $MED6_CUTTAG/7.macs2/merged_bam/*.homer.*Peak;
do
    sample=${i##*/}
    sample=${sample%_*}
    nohup annotatePeaks.pl $i mm10 -cpu 5 1>$sample.ann 2>$sample.stat &
done

for i in $MED6_CUTTAG/7.macs2/intersect_peak/*bed; 
do 
    sample=$(basename $i ".bed")
    nohup annotatePeaks.pl $i mm10 -cpu 5 1>$sample.ann 2>$sample.stat &
done


