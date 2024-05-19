
exclude_list=file(params.exclude_list)


sex_map=file(params.sex_map)


process filter_pass {
   cpus 4
   input:
      path(vcfgz)
   output:
      tuple path(out), path(tbi)
   storeDir params.truvari_out
   script:
    out = "res"+(vcfgz[0].simpleName.substring(0))+"_1.vcf.gz"
    tbi = "${out}.tbi"
   """
     bcftools view -f PASS *.vcf.gz -Oz -o $out
     bcftools index -t  --threads 4 $out
    """
}


process cleanUp {
    input:
       tuple val(tool), path(vcf)
    output:
       tuple val(tool), path("auto/*.vcf"), emit: auto
       tuple val(tool), path("XX/*.vcf"),  emit:XX
       tuple val(tool), path("XY/*.vcf"),  emit:XY
       tuple val(tool), path("all/*.vcf"),  emit:all           
    publishDir "${params.clean}/$tool"
    tag "cleanUp/${tool}"
    script:
      """
        hostname
        mkdir all auto XX XY _auto _all _XX _XY
        vcfpassfilter.py $sex_map
        for d in auto XX XY all; do
           cd _\$d
           for f in *vcf; do
              bcftools view -T ^${exclude_list} \$f | alt2type.py >  ../\$d/\$f
           done
           cd ..
        done

     """
}


process cleanUpStitch {
    input:
       tuple val(tool), path(vcf)
    output:
       tuple val(tool), path("auto_s/*.vcf"), emit: auto
       tuple val(tool), path("XX_s/*.vcf"),  emit:XX
       tuple val(tool), path("XY_s/*.vcf"),  emit:XY       
    storeDir "${params.clean}/$tool"
    tag "cleanUp/${tool}"
    script:
      """
        mkdir auto XX XY _auto _XX _YY
        vcfpassfilter.py $sex_map
        for d in auto XX YY all; do
           mkdir -p \$d
           cd \$d
           for f in *vcf; do
              stitch.py -a $f ../\${d}_s/$f
           done
           cd ..
        done
     """
}
