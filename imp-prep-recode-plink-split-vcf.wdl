task plink {
  File bim
  File fam
  File bed
  String basename 
  File Chromosome
  File Exclude
  File Force_Allele1
  File Position
  File Strand_Flip
  
  command <<<
  	mkdir bfiles
	mv ${bim} bfiles/
    mv ${bed} bfiles/
    mv ${fam} bfiles/
    
    
    plink --bfile bfiles/${basename} --exclude ${Exclude} --make-bed --out TEMP1
    plink --bfile TEMP1 --update-map ${Chromosome} --update-chr --make-bed --out TEMP2
    plink --bfile TEMP2 --update-map ${Position} --make-bed --out TEMP3
    plink --bfile TEMP3 --flip ${Strand_Flip} --make-bed --out TEMP4
    plink --bfile TEMP4 --a2-allele ${Force_Allele1} --make-bed --out updated
    plink --bfile updated --real-ref-alleles --make-bed --out updated
    plink --bfile updated --real-ref-alleles --recode vcf --out ${basename} 
    
  >>>

  runtime {
    disks: "local-disk 100 HDD"
    memory: "32 GB"
    docker: "jrose77/plinkdocker"
  }

  output {
    File out_vcf = "${basename}.vcf"
  }
}

task split_vcf {
  Int chr
  File in_vcf
  String basename

  command <<<
  	mkdir in_files/
    mv ${in_vcf} in_files/
    cd in_files && bgzip ${basename}.vcf && bcftools index ${basename}.vcf.gz && cd ..
    bcftools view -r ${chr} in_files/${basename}.vcf.gz| bgzip >chr${chr}.vcf.gz
    bcftools index chr${chr}.vcf.gz
  >>>

  runtime {
    disks: "local-disk 100 HDD"
    memory: "32 GB"
    docker: "vandhanak/bcftools:1.3.1"
  }

  output {
    File chr_vcf = "chr${chr}.vcf.gz"
    File chr_vcf_ind = "chr${chr}.vcf.gz.csi"
  }
}

workflow plink_workflow {
  String basename
  Int chr
  File bim
  File fam
  File bed
  File Chromosome
  File Exclude
  File Force_Allele1
  File Position
  File Strand_Flip

  call plink {
    input:
      	basename=basename,
        bim=bim,
        fam=fam,
        bed=bed,
        Chromosome=Chromosome,
        Exclude=Exclude,
        Force_Allele1=Force_Allele1,
        Position=Position,
        Strand_Flip=Strand_Flip
  }
  call split_vcf {
  	input:
    	in_vcf=plink.out_vcf,
        chr=chr,
        basename=basename
  }

  output {
	File chr_vcf=split_vcf.chr_vcf
    File chr_vcf_ind=split_vcf.chr_vcf_ind
  }
}