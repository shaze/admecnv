#!/usr/bin/env nextflow

nextflow.enable.dsl=2

 
Channel.fromPath(params.delly_path).set { delly_vcfs }
Channel.fromPath(params.gridss_path).set { gridss_vcfs }
Channel.fromPath(params.smoove_path).set { smoove_vcfs }
Channel.fromPath(params.manta_path).set { manta_vcfs }

validators = Channel.fromPath(params.cnvpytor_path)
dcnvpaths=Channel.fromPath(params.dcnvpaths)





yml = file("kn.yml")


refi = file("${params.ref}.fai")
ref  = file (params.ref)

chroms = file("chroms")

d = delly_vcfs.map { it -> ["delly", it] }.groupTuple()
g = gridss_vcfs.map { it -> ["gridss", it] }.groupTuple()
s = smoove_vcfs.map { it -> ["smoove", it] }.groupTuple()
m = manta_vcfs.map { it -> ["manta", it] }.groupTuple()


all_unfiltered  = d.mix(s).mix(g).mix(m)





params.out = "final"

support_vals = 2..4

include { cleanUp } from "./modules/cleanup.nf"
include { doSurvivor as doSurvivorSXX  } from "./modules/mod_survivor.nf" addParams(sub_tag: "-XY", enabled: false)
include { doSurvivor as doSurvivorSXY  } from "./modules/mod_survivor.nf" addParams(sub_tag: "-XX", enabled: false)
include { doSurvivor as doSurvivorA  } from "./modules/mod_survivor.nf" addParams(sub_tag: "-auto", enabled: false)



d = delly_vcfs.map { it -> ["delly", it] }.groupTuple() 
g = gridss_vcfs.map { it -> ["gridss", it] }.groupTuple()
s = smoove_vcfs.map { it -> ["smoove", it] }.groupTuple()
m = manta_vcfs.map { it -> ["manta", it] }.groupTuple()




chroms = ((1..22).collect { String.format("chr%d",it) })+["chrX","chrY"]


process annotate_short {
    cpus 2
    memory '80.GB'
    input:
       tuple val(support), path(short_sv), path(si)
       each chrom
    output:
       path("${annot}.tsv")
    script:
       abbrev="abbrev.vcf"
       annot =chrom+"_annot"
       """
         bcftools view -r $chrom $short_sv |
              cut -f 1-9 | grep -v '##'  | sd "SUPP_VEC=[01]+" ""   > $abbrev
         AnnotSV -SVinputFile $abbrev -outputDir ./ -outputFile $annot
       """
}



process filter_long {
    memory '16.GB'
    input:
       tuple path(long_sv), path(li)
    output:
       tuple path(filter_long), path(fi)
    script:
       filter_long=long_sv.simpleName+"_filtered.vcf.gz"
       fi = filter_long+".tbi"
       a = 10000
       b = params.max_long_sv_len
       """
       bcftools filter -i \
        "((INFO/SVLEN>$a)&&(INFO/SVLEN<$b)) || ((INFO/SVLEN<-$a)&&(INFO/SVLEN>-$b))"  \
       -Oz -o $filter_long  $long_sv
       tabix $filter_long
    """
}
    
process annotate_long {
    memory '24.GB'
    input:
       tuple path(long_sv), path(li)
    output:
      path("${lonf_annot}.tsv")
    publishDir "${params.out_dir}/", mode: 'copy'
    script:
       long_annot=long_sv.simpleName+"_long_annot"
       """
         AnnotSV -SVinputFile $long_sv   -outputDir ./ -outputFile $long_annot
       """
}


process visualise_short  {
    input:
       path(annotated)
    output:
    tuple path("consensus.tsv"), path("*html"), path("*xlsm")
    publishDir "${params.out_dir}/", mode: 'copy'
    """
         head -n1 chr1_annot.tsv > consensus.tsv
         tail -q -n+2  chr[1-9]_annot.tsv chr1[0-9]_annot.tsv chr2[0-2]_annot.tsv chr[XY]_annot.tsv  >> consensus.tsv
         knotAnnotSV2XL.pl --configFile $yml --annotSVfile consensus.tsv  --outDir ./ --genomeBuild hg38
         knotAnnotSV.pl --configFile $yml --annotSVfile consensus.tsv  --outDir ./ --genomeBuild hg38
    """
}

process visualise_long  {
    input:
       path(annotated)
    output:
      tuple path("*html"), path("*xlsm")
    publishDir "${params.out_dir}/", mode: 'copy'
    """
         knotAnnotSV2XL.pl --configFile $yml --annotSVfile $annotated  --outDir ./ --genomeBuild hg38
         knotAnnotSV.pl --configFile $yml --annotSVfile $annotated  --outDir ./ --genomeBuild hg38
    """
}


process indivVCFs {
   cpus 4
   input:
      tuple val(tool), path(vcfs)
   output:
     tuple val(tool), path("pass")
   clusterOptions = '--constraint=avx2'
   script:
      """
      mkdir pass
      for x in *vcf; do
         bcftools reheader --threads 4 -f $refi \$x | grep -F -v  "./1" |\
                 bcftools view -f PASS  -Oz -o pass/\${x}.gz -
         bcftools index -t --threads 4 pass/\${x}.gz
      done
      """
}


process svimmer {
    memory "40.G"
    input:
      tuple val(tool), path(pass)
    output:
       tuple val("${tool[0]}"),  path("all.vcf.gz")
    script:
      out = "${tool[0]}_1.vcf"
     """
      ls pass/*vcf.gz  > vcf_list
      svimmer vcf_list \$(cat $chroms)  | bgzip > all.vcf.gz
      """
}    

process reheader {
    cpus 4
    input:
      tuple val(label), path(vcf)
    output:
      tuple val(label), path("temp.vcf.gz")
    script:
      """ 
      hostname
      bcftools reheader --threads 4 -f $refi $vcf > temp.vcf.gz
      """
}

process decolonise {
    cpus 4
    input:
      tuple val(label), path(vcfgz)
    output:
      tuple val(label), path("nocolons.vcf.gz")
    script:
    """
      decolonise.py $vcfgz nocolons.vcf
      bgzip -@4 nocolons.vcf
    """
}

process sortVCF {
    memory "20.G"
    input:
      tuple val(label), path(nocolon)
    output:
      tuple path(out), path("${out}.tbi")
    script:
      out = "${label}_1.vcf.gz"
      """
        bcftools sort  -Oz -o ${out} $nocolon
        tabix ${out}
    """
}


process combineVCFs {
    memory "20.G"
    cpus 4
    input:
      tuple val(support), path(vcfs)
    output:
      tuple val(support), path(out), path("${out}.tbi")
    publishDir params.out_dir
    script:
       out = "${params.out}_${support}.vcf.gz"    
    """
       bcftools   concat --threads 4  *-auto*.vcf.gz *-sex*.vcf.gz -Oz -o $out
       tabix $out
    """
}


process mergeSexVCFs { // combine XX and XY calls
    input:
       tuple val(support), path(vcfs)
    output:
       tuple val(support), path ("${out}"), path("${out}.tbi")
    script:
       out = vcfs[0].simpleName.replace("-XX_","").replace("-XY_","")+"-sex.vcf.gz"
      """
         bcftools query -l *-XX_*vcf.gz > xx
         bcftools query -l *-XY_*vcf.gz > xy
         bcftools merge -m none -0  -Ob *vcf.gz -o merged.bcf
         sort xx xy > elts
         bcftools view -S elts merged.bcf -Oz -o ${out}
         tabix $out
    """
}


process sortBySample {
    // Ensure output VCF columns in alphabetic order
    input:
      tuple val(support), path(vcf), path(index)
    output:
      tuple val(support), path(ovcf), path(oindex)
    script:
      ovcf   = "s-${vcf}"
      oindex = "s-${index}"
      """
         bcftools query -l $vcf | sort > inds
         bcftools view -S inds $vcf -Oz -o ${ovcf}
         tabix ${ovcf}
      """
}

process removeLongSVs {
    input:
       tuple val(support), path(vcf), path(idx)
    output:
       tuple val(support), path(shorter), path(index)
    script:
      shorter = "s-${vcf}"
      index = "${shorter}.tbi"
      """
        bcftools filter -e '(INFO/SVLEN>10000) | (INFO/SVLEN<-10000)'   $vcf -Oz -o $shorter
        tabix $shorter
      """
}





new_line = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">'

process cnv_reheader {
    input:
     path(cnv)
    output:
     path(rcnv)
    script:
      base = cnv.simpleName
      rcnv = "${base}.dcnv.vcf"
      """
      grep "#" $cnv > header
      sed  '/ID=CNSD/a\\
$new_line' header > new_head
      stitch.py $cnv stitched.vcf
      bcftools reheader -h new_head  -o $rcnv stitched.vcf
    """
}


process validate_one_sample {
   maxForks 12
   cpus 2
   input:
    tuple val(base), path(pytor), path('gridss.vcf'), path('manta.vcf'), \
	path('dcnv.vcf'), path('dsv.vcf'), path('smoove.vcf')
   output:
      path(out)
   script:
    out = "${base}.vcf.gz"
    """
      for x in dcnv.vcf manta.vcf  dsv.vcf smoove.vcf; do  #gridss not included as too short
         echo \$x
         validate_call.py -v $pytor --out v-\${x}  \$x
         bgzip v-\${x}
         tabix v-\${x}.gz
      done
      bcftools concat -a  --threads 2   v-dcnv.vcf.gz v-manta.vcf.gz v-dsv.vcf.gz v-smoove.vcf.gz  -Ob -o base.bcf
      bcftools sort base.bcf -Ov | stitch.py - - | bgzip > $out
    """
}


process gz_and_index {
    maxForks 12
    input:
      path(vcf)
    output:
      tuple  path(vcfgz), path(idx)
    script:
      res=vcf.simpleName  // file may be gzipped already
      vcfgz = "${res}.vcf.gz"
      idx   = "${vcfgz}.csi"
     """
        pwd  > p
        which bcftools >& b
        hostname > hostname
        if [[ $vcf =~ gz ]]; then
            echo "already gz"
        else
            bcftools view -Oz -o $vcfgz $vcf
        fi
        bcftools index $vcfgz
        ls > files
     """
}

    
process combine_longs {
    cpus 2
    memory '14 GB'
    input:
    path(all_vcfs)
    output:
      tuple path("${out}"), path("${out}.tbi")
    publishDir "${params.out_dir}/", mode: 'copy'
    script:
    out = "longer_svs.vcf.gz"
    """
bcftools merge -m none *vcf.gz | bcftools annotate -x INFO/SR,FORMAT/PR,FORMAT/PL,FORMAT/RCR,FORMAT/RCL,FORMAT/DR,FORMAT/GL,FORMAT/DV,FORMAT/DV,FORMAT/RC,FORMAT/RR,FORMAT/RV,FORMAT/RDCN,FORMAT/SR,FORMAT/SQ,FORMAT/RP,FORMAT/QR,FORMAT/AO,FORMAT/RS,FORMAT/DHFC,FORMAT/QA,FORMAT/AS,FORMAT/ASC,FORMAT/DHSP,FORMAT/RO,FORMAT/DHFFC,FORMAT/AP,FORMAT/DP,FORMAT/AB,FORMAT/DHBFC,FORMAT/DHSP,FORMAT/SHQ,FORMAT/FT -Ob -o merge.bcf
bcftools filter -i "ALT='<DEL>'"  -Oz -o dels.vcf.gz merge.bcf
bcftools index -t dels.vcf.gz
truvari collapse -p 0 -P 0.75 -i dels.vcf.gz -o d_mrg.vcf -c d_col.vcf -f /dataB/aux/38/Homo_sapiens_assembly38.fasta >& del.log
bcftools sort -m 10G d_mrg.vcf -Ob -o d_mrg.bcf
bcftools index d_mrg.bcf
bcftools filter -e "ALT='<DEL>'"  -Oz -o others.vcf.gz merge.bcf
bcftools index -t others.vcf.gz
truvari collapse -p 0.9 -P 0.8 -i others.vcf.gz -o oth_mrg.vcf -c oth_col.vcf -f /dataB/aux/38/Homo_sapiens_assembly38.fasta >& oth.log
bcftools sort -m 10G  oth_mrg.vcf -Ob -o oth_mrg.bcf
bcftools index oth_mrg.bcf
bcftools concat -a d_mrg.bcf oth_mrg.bcf -Ob -o combined.bcf
bcftools index  combined.bcf
bcftools sort -m 10G  combined.bcf -Oz -o ${out}
bcftools index -t ${out}
    """
}


process merge_long_short  {
    input:
       tuple path(long_sv), path(li), val(support), path(short_sv), path(si)
    output:
       tuple path(result_gz), path(result_idx)
    publishDir params.out_dir
    script:
      result = "${params.out}.vcf"
      result_gz = "${result}.gz"
      result_idx= "${result_gz}.tbi"
      """
        bcftools query -l $short_sv > short_sample_names
        bcftools view -S short_sample_names $long_sv -Ov -o long.vcf.gz
        tabix long.vcf.gz
        bcftools concat -a $short_sv long.vcf.gz  -d all -Ob -o $result_gz
        bcftools index -t $result_gz
    """
}

workflow getShorterSVs {
    take:   XX
    take:   XY
    take:   auto
    take:   support_vals
    main:
      doSurvivorSXX(XX, support_vals)
      doSurvivorSXY(XY, support_vals)
      mergedXY =doSurvivorSXX.out.\
        mix(doSurvivorSXY.out).\
        groupTuple().map { it -> [it[0], it[1]+it[2]] }
      mergeSexVCFs(mergedXY)
      doSurvivorA(auto, support_vals) | \
	sortBySample | \
	mix(mergeSexVCFs.out) |\
        groupTuple() |\
        map { it -> [it[0], it[1]+it[2]] }|\
        combineVCFs |
	removeLongSVs |
        filter { it -> it[0]==2 } |\
        set {result }
    
    emit:
       result
}



workflow longer_svs  {
    take: all_candidates
    main:
    cnv_reheader(dcnvpaths) | mix(all_candidates) | \
        map { it -> [it.simpleName , it] } | groupTuple() | \
	map { it -> [it[0], it[1][0], it[1][1], it[1][2], it[1][3], it[1][4]] } | set{ all_candidates_plus } 
      long_candidates = validators.map { it -> [it.simpleName, it] }.join(all_candidates_plus)
      validate_one_sample(long_candidates)  | \
	gz_and_index | \
	flatten |  \
	toList | \
	combine_longs | filter_long | annotate_long 
    emit:
      annotate_long.out

}

workflow {
  main:
   params.get_one_caller_support=true
    split = cleanUp(all_unfiltered)
    shorter= getShorterSVs(split.XX, split.XY, split.auto, support_vals) 
    annotate_short(shorter,chroms).toList() | visualise_short
    split.all.map {  it -> it[1] }.flatten()  |  longer_svs | visualise_long
}
