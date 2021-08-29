#!/bin/bash
gunzip -k ../vcf_phase3_hg38_v2/ALL.$1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
mv ../vcf_phase3_hg38_v2/ALL.$1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf ./in.vcf
sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' in.vcf > out.vcf
rm in.vcf
../vcftools-0.1.16/vcftools \
--gzvcf out.vcf \
--keep EUR_indv_phase3.txt \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--mac 1 \
--IMPUTE \
--phased \
--max-missing 1 \
--out chr21_EUR_panel
#--bed ../wgEncodeDukeMapabilityUniqueness20bp/wgEncodeDukeMapabilityUniqueness20bp.bed \
#--maf 0.1 \

rm out.vcf
