task merge {
  File out
  Int mem
  Int disk
  
  command <<<
  	tar -xvf ${out}  
    # Index all VCFs
    bcftools index out/chr1.dose.vcf.gz
    bcftools index out/chr2.dose.vcf.gz
    bcftools index out/chr3.dose.vcf.gz
    bcftools index out/chr4.dose.vcf.gz
    bcftools index out/chr5.dose.vcf.gz
    bcftools index out/chr6.dose.vcf.gz
    bcftools index out/chr7.dose.vcf.gz
    bcftools index out/chr8.dose.vcf.gz
    bcftools index out/chr9.dose.vcf.gz
    bcftools index out/chr10.dose.vcf.gz
    bcftools index out/chr11.dose.vcf.gz
    bcftools index out/chr12.dose.vcf.gz
    bcftools index out/chr13.dose.vcf.gz
    bcftools index out/chr14.dose.vcf.gz
    bcftools index out/chr15.dose.vcf.gz
	bcftools index out/chr16.dose.vcf.gz
	bcftools index out/chr17.dose.vcf.gz
    bcftools index out/chr18.dose.vcf.gz
    bcftools index out/chr19.dose.vcf.gz
    bcftools index out/chr20.dose.vcf.gz
    bcftools index out/chr21.dose.vcf.gz
    bcftools index out/chr22.dose.vcf.gz
    # Create a list of all VCFs present 
    find . -name '*.dose.vcf.gz' > list.tsv
	bcftools concat -f list.tsv |bgzip > topmed_out_ALL.vcf.gz 
	bcftools index topmed_out_ALL.vcf.gz

    
  >>>

  runtime {
    disks: "local-disk ${disk} HDD"
    memory: "${mem} GB"
    docker: "vandhanak/bcftools:1.3.1"
  }

  output {
	File topmed_out_ALL= "topmed_out_ALL.vcf.gz"
    File topmed_out_ALL_ind= "topmed_out_ALL.vcf.gz.csi"
  }
}

task download {
  String link
  String password
  Int mem
  Int disk
  
  command <<<
  
    ls -lah
    mkdir out && cd out
    curl -sL ${link} | bash
    for zipped_file in *.zip; do unzip -P "${password}" $zipped_file; done
    rm *.zip
    ls -lah 
    
    cd ../
    tar -cvzf out.tar.gz out

  >>>

  runtime {
    disks: "local-disk ${disk} HDD"
    memory: "${mem} GB"
    docker: "mgibio/cle:latest"
  }

  output {
	File out = "out.tar.gz"
  }
}


workflow plink_workflow {
  String link
  String password

  
  # runtime parameters 
  Int disk = 250
  Int mem = 32

  call download {
    input:
      	link=link,
        password=password,
        mem = mem,
        disk =disk
  }
  call merge {
  	input:
    	out=download.out,
        mem = mem,
        disk =disk
  
  }
  output{
    File TOPMed_out = merge.topmed_out_ALL
    File TOPMed_out_ind = merge.topmed_out_ALL_ind
    		
    }

}