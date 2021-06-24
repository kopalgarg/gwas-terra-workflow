task plink_freq {
  File bim
  File fam
  File bed
  String basename 
  
  command <<<
  	mkdir bfiles
	mv ${bim} bfiles/
    mv ${bed} bfiles/
    mv ${fam} bfiles/
    plink --bfile bfiles/${basename} --freq --out ${basename}
  >>>

  runtime {
    disks: "local-disk 100 HDD"
    memory: "12 GB"
    docker: "jrose77/plinkdocker"
  }

  output {
    File out_freq = "${basename}.frq"
  }
}


workflow plink_workflow {
  String basename
  File bim
  File fam
  File bed

  call plink_freq {
    input:
      	basename=basename,
        bim=bim,
        fam=fam,
        bed=bed
  }
  output {
	File freq_file=plink_freq.out_freq
  }
}