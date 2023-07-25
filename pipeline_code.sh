conda activate pipeline_macvs

cd /media/rin/FGFR3fusion
for sample in TL-20-FEDEFE_T_CSQ1 TL-21-8NDXZMXY_T_CSQ1 TL-22-HXTS2TE8_T_CSQ1 TL-21-76F2KNCB_T_CSQ1 TL-21-Y6JNT5E9_T_CSQ1 TL-22-AFGW9GPP_T_DSQ1 TL-22-TDNPZF2T_T_CSQ1 TL-20-BB6020_T_DSQ1 TL-22-BG7DFVNW_T_CSQ1 TL-22-4ZEBBBBQ_T_CSQ1 TL-22-IZF22MT5_T_CSQ1 TL-22-WB5D2ECV_T_DSQ1 TL-22-WYHF3858_T_CSQ1 TL-22-VNMGQW33_T_CSQ1 TL-22-KCPE65I2_T_CSQ1 TL-22-Z24N63NZ_T_CSQ1 TL-23-263ZB4MS_T_CSQ1 TL-23-G7DCC66X_T_CSQ1 TL-23-I4MHR3WM_T_CSQ1 TL-22-CB792G9R_T_CSQ1 TL-22-4ZU93YNJ_T_CSQ1 TL-20-BB6020_N_DSQ1
do
    if [ ! -e /Volumes/UTUCproject/TO_DRIVE/${sample}/${sample}.dedup.numAligned.txt ]
    then
        cd ${sample}

        trim_galore --fastqc --paired -j 10 ${sample}_1.fastq ${sample}_2.fastq

        bwa mem -R '@RG\tID:1\tSM:'${sample}'\tPL:illumina\tLB:lib1\tPU:unit1' -t 15 /Volumes/UTUCproject/TO_DRIVE/hg19.fa ${sample}_1_val_1.fq ${sample}_2_val_2.fq | samtools view -bS -F 0x900 - > ${sample}.aligned.bam
        samtools view -c -F 0x4 ${sample}.aligned.bam > ${sample}.numAligned.txt

        rm -R ${sample}_1_val_1.fq ${sample}_2_val_2.fq

        mkdir tmp
        samtools sort -@ 10 -m 5G -T tmp/aln.sorted -o ${sample}.aligned.bam ${sample}.aligned.bam
        samtools index ${sample}.aligned.bam

        picard MarkDuplicates --INPUT ${sample}.aligned.bam --OUTPUT ${sample}.dedup.aligned.bam --METRICS_FILE ${sample}.metrics.txt --ASSUME_SORTED True --REMOVE_DUPLICATES True
        samtools sort -@ 10 -m 3G -T tzmp/aln.sorted -o ${sample}.dedup.aligned.bam ${sample}.dedup.aligned.bam
        samtools index ${sample}.dedup.aligned.bam

        samtools view -c -F 0x4 ${sample}.dedup.aligned.bam > ${sample}.dedup.numAligned.txt
        cd ..
    fi
done

cd /media/rin/FGFR3fusion
for sample in TL-20-FEDEFE_T_CSQ1 TL-21-8NDXZMXY_T_CSQ1 TL-22-HXTS2TE8_T_CSQ1 TL-21-76F2KNCB_T_CSQ1 TL-21-Y6JNT5E9_T_CSQ1 TL-22-AFGW9GPP_T_DSQ1 TL-22-TDNPZF2T_T_CSQ1 TL-20-BB6020_T_DSQ1 TL-22-BG7DFVNW_T_CSQ1 TL-22-4ZEBBBBQ_T_CSQ1 TL-22-IZF22MT5_T_CSQ1 TL-22-WB5D2ECV_T_DSQ1 TL-22-WYHF3858_T_CSQ1 TL-22-VNMGQW33_T_CSQ1 TL-22-KCPE65I2_T_CSQ1 TL-22-Z24N63NZ_T_CSQ1 TL-23-263ZB4MS_T_CSQ1 TL-23-G7DCC66X_T_CSQ1 TL-23-I4MHR3WM_T_CSQ1 TL-22-CB792G9R_T_CSQ1 TL-22-4ZU93YNJ_T_CSQ1 TL-20-BB6020_N_DSQ1
do
	cd ${sample}
	/home/rin/bam-readcount/build/bin/bam-readcount -f /media/rin/FGFR3fusion/hg19.fa  ${sample}.dedup.aligned.bam  chr17:7577547-7577547 -w 0 > ${sample}_TP53_VAF.txt
	/home/rin/bam-readcount/build/bin/bam-readcount -f /media/rin/FGFR3fusion/hg19.fa  ${sample}.dedup.aligned.bam  chr17:41226461-41226461 -w 0 > ${sample}_BRCA1_VAF.txt
	/home/rin/bam-readcount/build/bin/bam-readcount -f /media/rin/FGFR3fusion/hg19.fa  ${sample}.dedup.aligned.bam  chr7:140481402-140481402 -w 0 > ${sample}_BRAF_VAF.txt
	/home/rin/bam-readcount/build/bin/bam-readcount -f /media/rin/FGFR3fusion/hg19.fa  ${sample}.dedup.aligned.bam  chr10:123279605-123279605 -w 0 > ${sample}_FGFR2_VAF.txt
	/home/rin/bam-readcount/build/bin/bam-readcount -f /media/rin/FGFR3fusion/hg19.fa  ${sample}.dedup.aligned.bam  chr1:115256529-115256529 -w 0 > ${sample}_NRAS_VAF.txt
	/home/rin/bam-readcount/build/bin/bam-readcount -f /media/rin/FGFR3fusion/hg19.fa  ${sample}.dedup.aligned.bam  chr4:1808911-1808911 -w 0 > ${sample}_fusion_denominator_VAF.txt
	cd ..
done

#RACC3 coordinates chr4:1741279-1742970 


#running gatk now
#run only once
conda activate gatk_env
cd /media/rin/FGFR3fusion
gatk CreateSequenceDictionary -R hg19.fa

conda activate gatk_env
cd /media/rin/FGFR3fusion
for sample in TL-20-FEDEFE_T_CSQ1 TL-21-8NDXZMXY_T_CSQ1 TL-22-HXTS2TE8_T_CSQ1 TL-21-76F2KNCB_T_CSQ1 TL-21-Y6JNT5E9_T_CSQ1 TL-22-AFGW9GPP_T_DSQ1 TL-22-TDNPZF2T_T_CSQ1 TL-20-BB6020_T_DSQ1 TL-22-BG7DFVNW_T_CSQ1 TL-22-4ZEBBBBQ_T_CSQ1 TL-22-IZF22MT5_T_CSQ1 TL-22-WB5D2ECV_T_DSQ1 TL-22-WYHF3858_T_CSQ1 TL-22-VNMGQW33_T_CSQ1 TL-22-KCPE65I2_T_CSQ1 TL-22-Z24N63NZ_T_CSQ1 TL-23-263ZB4MS_T_CSQ1 TL-23-G7DCC66X_T_CSQ1 TL-23-I4MHR3WM_T_CSQ1 TL-22-CB792G9R_T_CSQ1 TL-22-4ZU93YNJ_T_CSQ1
do
	cd ${sample}
	gatk Mutect2 -R /media/rin/FGFR3fusion/hg19.fa -I ${sample}.dedup.aligned.bam -I /media/rin/FGFR3fusion/TL-20-BB6020_N_DSQ1/TL-20-BB6020_N_DSQ1.dedup.aligned.bam -normal TL-20-BB6020_N_DSQ1 -O ${sample}.GATK.somatic_unfiltered.dedup.vcf --native-pair-hmm-threads 30
	gatk FilterMutectCalls -R /media/rin/FGFR3fusion/hg19.fa -V ${sample}.GATK.somatic_unfiltered.dedup.vcf -O ${sample}.GATK.somatic_filtered.dedup.vcf
	cd ..	
done



