


process combine_multiple_samples {
    cpus 2
    memory '64 GB'
    input:
    tuple val(tool), path(all_vcfs)
    output:
      tuple path("${out}"), path("${out}.tbi")
    publishDir params.out_dir
    script:
    out = "${tool}_1.vcf.gz"
    """
bcftools merge -m none *bcf | bcftools annotate -x INFO/SR,FORMAT/PR,FORMAT/PL,FORMAT/RCR,FORMAT/RCL,FORMAT/DR,FORMAT/GL,FORMAT/DV,FORMAT/DV,FORMAT/RC,FORMAT/RR,FORMAT/RV,FORMAT/RDCN,FORMAT/SR,FORMAT/SQ,FORMAT/RP,FORMAT/QR,FORMAT/AO,FORMAT/RS,FORMAT/DHFC,FORMAT/QA,FORMAT/AS,FORMAT/ASC,FORMAT/DHSP,FORMAT/RO,FORMAT/DHFFC,FORMAT/AP,FORMAT/DP,FORMAT/AB,FORMAT/DHBFC,FORMAT/DHSP,FORMAT/SHQ,FORMAT/FT -Ob -o merge.bcf
bcftools filter -i "ALT='<DEL>'"  -Oz -o dels.vcf.gz merge.bcf
bcftools index -t dels.vcf.gz
truvari collapse --pctsim 0 -P 0.75 -i dels.vcf.gz -o d_mrg.vcf -c d_col.vcf -f /dataB/aux/38/Homo_sapiens_assembly38.fasta >& del.log
bcftools sort -m 10G d_mrg.vcf -Ob -o d_mrg.bcf
bcftools index d_mrg.bcf
bcftools filter -e "ALT='<DEL>'"  -Oz -o others.vcf.gz merge.bcf
bcftools index -t others.vcf.gz
truvari collapse -p 0.9 -P 0.8 -i others.vcf.gz -o oth_mrg.vcf -c oth_col.vcf -f /dataB/aux/38/Homo_sapiens_assembly38.fasta >& oth.log
bcftools sort -m 60G  oth_mrg.vcf -Ob -o oth_mrg.bcf
bcftools index oth_mrg.bcf
bcftools concat -a d_mrg.bcf oth_mrg.bcf -Ob -o combined.bcf
bcftools index  combined.bcf
bcftools sort -m 60G  combined.bcf -Oz -o ${out}
bcftools index -t ${out}
    """
}

