#!/bin/bash
bcftools view \
../vcf_phase3_hg38_v2/ALL.$1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
--samples-file samples_per_panel/$2_panel.txt \
--force-samples \
--exclude-types indels,mnps,ref,bnd,other \
--min-alleles 2 \
--max-alleles 2 \
--min-ac 1:minor \
--phased \
--exclude 'AN!=2*N_SAMPLES' \
--output-file $1_$2_panel.bcf \
--output-type u
# --regions-file ../wgEncodeDukeMapabilityUniqueness20bp/wgEncodeDukeMapabilityUniqueness20bp.bed \
# --min-af 0.1:minor \

bcftools convert \
$1_$2_panel.bcf \
--haplegendsample $1_$2_panel

rm $1_$2_panel.bcf
mkdir ./$2_panel.hg38.BCFtools
mv $1_$2_panel.legend.gz ./$2_panel.hg38.BCFtools
mv $1_$2_panel.hap.gz ./$2_panel.hg38.BCFtools
mv $1_$2_panel.samples ./$2_panel.hg38.BCFtools
