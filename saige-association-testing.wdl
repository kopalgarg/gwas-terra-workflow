task split_vcf {

	File inFile
    Int chr
    Int? mem
    Int cpu
    Int preempt
    Int disk
    File inFileIndex
    
    command <<<

        bcftools view -r chr${chr} ${inFile} |bgzip > chr${chr}.vcf.gz
        bcftools index chr${chr}.vcf.gz 

    >>>

    output {
        File vcf = "chr${chr}.vcf.gz"
        File vcf_ind= "chr${chr}.vcf.gz.csi"
    }

    runtime {
        disks: 			"local-disk ${disk} HDD"
        cpu: 			"${cpu}"
        memory: 		"${mem} GB"
        preemptible: 	"${preempt}"
        docker: 		"vandhanak/bcftools:1.3.1"
    }

    




}
task fit_null {

	# input data
    	# flags and fields
    String phenotype
    String trait_type
    String? covariate_list
    String covariate_line = if defined(covariate_list) then '--covarColList=' + covariate_list else ''
    String id_col 
    
		# files
    File covariates_file
    File inbed
    File inbim
    File infam
    
    # runtime parameters 
    Int cpu
    Int disk
    Int mem
    Int preempt
    


    command <<<

        basename=$(basename ${inbed} .bed)
        directory=$(dirname ${inbed})
        plink_path="$directory"/"$basename"

        step1_fitNULLGLMM.R \
            --plinkFile=$plink_path \
            --phenoFile=${covariates_file} \
            --phenoCol=${phenotype} \
            ${covariate_line} \
            --invNormalize=F \
            --sampleIDColinphenoFile=${id_col} \
            --traitType=quantitative \
            --nThreads=${cpu} \
            --outputPrefix=step1 

    >>>

    output {
        File gmmat = "step1.rda"
        File variance_ratio = "step1.varianceRatio.txt"
    }

    runtime {
        disks: 			"local-disk ${disk} HDD"
        cpu: 			"${cpu}"
        memory: 		"${mem} GB"
        preemptible: 	"${preempt}"
        docker: 		"finngen/saige:0.39.1.fg"
    }
}


task spa {
	
    # input files
	File inVCF
    File inVCF_index
    
    # input options
    Int chr
    File samples
    File varianceRatio
    File gmmat
    String genotype
    
    # runtime parameters 
    Int cpu
    Int disk
    Int mem
    Int preempt

	command <<<
    
    	step2_SPAtests.R \
        	--vcfFile=${inVCF}  \
        	--vcfFileIndex=${inVCF_index} \
        	--vcfField=${genotype} \
        	--chrom=chr${chr} \
        	--minMAF=0.01 \
        	--minMAC=1 \
        	--sampleFile=${samples} \
        	--GMMATmodelFile=${gmmat} \
        	--varianceRatioFile=${varianceRatio} \
        	--SAIGEOutputFile=${chr}.saige2.out.SAIGE.txt.gz \
        	--numLinesOutput=2 \
        	--IsOutputAFinCaseCtrl=TRUE

    
    >>>

    output {
        File sum_stat="${chr}.saige2.out.SAIGE.txt.gz"
    }


    runtime {
        disks: 			"local-disk ${disk} HDD"
        cpu: 			"${cpu}"
        memory: 		"${mem} GB"
        preemptible: 	"${preempt}"
        docker: 		"finngen/saige:0.39.1.fg"
    }

}
workflow saige {


	# input files
    
    File covariates_file
    String phenotype
    File inbed
    File inbim
    File infam    
    File samples
    File inVCF
    File inVCF_index
    
    # input options
    String? covariate_list
    Int chr
    String genotype
    String id_col
    
    # runtime parameters 
    Int? cpu = 4 
    Int? disk = 100
    Int? mem = 32
    Int? preempt = 3
    
		call split_vcf {
    	input:
        	cpu=cpu,
            mem=mem,
            preempt=preempt,
            disk=disk,
            chr=chr,
            inFile=inVCF,
            inFileIndex=inVCF_index

    }

	    call fit_null {
        input:
            covariates_file=covariates_file,
            phenotype=phenotype,
           	inbed=inbed,
            inbim=inbim,
            infam=infam,
            covariate_list=covariate_list,
            id_col=id_col,
            cpu=cpu,
            mem=mem,
            preempt=preempt,
            disk=disk
   }

   	    call spa {
        input:
        	inVCF=split_vcf.vcf,
            inVCF_index=split_vcf.vcf_ind,
            chr=chr,
            samples=samples,
            varianceRatio=fit_null.variance_ratio,
            gmmat=fit_null.gmmat,
            genotype=genotype,
            cpu=cpu,
            mem=mem,
            preempt=preempt,
            disk=disk
            
   }
   		output{
    
    		File sum_stats = spa.sum_stat
    		
    }

}