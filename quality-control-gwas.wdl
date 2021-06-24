task inc {
    Int i
    File inputVCF
    String inputVCF_name 
    Int threads 
    Float memory
    Int diskGB

    command{
      bcftools index -f ${inputVCF} > "${inputVCF_name}.csi"
      bcftools view ${inputVCF} -r "${i}" | bgzip > "chr${i}.vcf.gz"

    }

    output {
      File outputVCF= "chr${i}.vcf.gz"
      File index= "${inputVCF_name}.csi"
    }
    runtime {
        docker: "gargk/bcftools:v1.8"
        cpu: "${threads}"
        memory: "${memory} GB"
        preemptible: 3
        disks:  "local-disk " + diskGB + " HDD"

    }
  }
   task bcfstats{
      File inputVCF
	  Int threads

   
    command{
        bcftools stats -s- ${inputVCF} > "bcfstats.tsv" 
        grep "#" -v "bcfstats.tsv" | grep "PSC" > "psc.tsv"
        }
 
    runtime{
        docker: "gargk/bcftools:v1.8" 
        cpu: "${threads}"
        memory: "20 GB"
        preemptible: 3
   }
    output{
        File bcfstats = "bcfstats.tsv"
        File psc = "psc.tsv"
        }

}
    task missingness_plink{
        File inputVCF
        Int threads

        command{
            plink1.9 --vcf ${inputVCF} --const-fid 0 -missing --out data_miss 
        }
        runtime{
            docker : "biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1"
            cpu: "${threads}"
            memory: "20 GB"
        	  preemptible: 3
        }
    output{
        File data_miss = "data_miss.imiss"
        }
}
    task run_rscript {
     File statPlot
     Int threads
     File data_miss
     File psc
          
     
        command {
            mv ${statPlot} script.R
            mv ${data_miss} data_miss.imiss
            mv ${psc} psc.tsv
           
                       
            R -e "source('script.R'); a=statPlot(missingness_file='data_miss.imiss',psc_file='psc.tsv'); write_tsv(a, 'samples_to_keep.txt', col_names = FALSE)" 
         }
         runtime{
            docker : "rocker/tidyverse" 
            cpu: "${threads}"
            memory: "20 GB"
        	preemptible: 3
         
        }
        output{
        	File to_keep = "samples_to_keep.txt"	
        }

}

	task remove_qc_failed_samples{
    	File inputVCF
        File to_keep
        Int threads 
        String inputVCF_name 

    
        command{
            bcftools index ${inputVCF} > "${inputVCF_name}.csi"
            
            bcftools view ${inputVCF} -S ${to_keep} |
            bgzip > "${inputVCF_name}_modified.vcf.gz"

            }
          
        
        runtime{
            docker: "gargk/bcftools:v1.8"    
            cpu: "${threads}"
            memory: "20 GB"
            preemptible: 3
            continueOnReturnCode: true
            
        }
        output{
            File modified_vcf_01 = "${inputVCF_name}_modified.vcf.gz"
        }
        }
        

     task plink_qc_steps{
            Int threads
            File modified_vcf_01

              command{
                  plink1.9 --vcf ${modified_vcf_01}  --maf 0.01  --geno 0.001 --indep-pairwise 50 5 0.2 --const-fid --recode vcf
          }
              runtime{
                  docker : "biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1"
                  cpu: "${threads}"
                  memory: "20 GB"
                  preemptible: 3
          }
              output{
              	File outputVCF = "plink.vcf"

          }
        }
       task merge{
       	Int threads
        Array[File] vcf
        
        	command{
            	bcftools concat ${sep=' ' vcf} > "outputVCF.vcf"
            
            	}
       		runtime{
                docker: "gargk/bcftools:v1.8"    
                cpu: "${threads}"
                memory: "20 GB"
                preemptible: 3
                continueOnReturnCode: true
            
            	}
            output{
            	File mergedOutput = "outputVCF.vcf"
            
            }
            
       }	


    workflow wf {
  	 File statPlot
     #Runtime parameters 
     Int threads
     File inputVCF
    Int diskGB = ceil(size(inputVCF, "GB") + 75)

     
     Array[Int] integers = [5,6]
        scatter(i in integers) {
            call inc{
                input: 
                	inputVCF=inputVCF,
                    i=i,
                    diskGB=diskGB
                    
        }
            call bcfstats {
                input:
                    inputVCF=inc.outputVCF,
                    threads=threads
                    

             } 
            call missingness_plink {
        		input:
                	inputVCF=inc.outputVCF,
                	threads=threads
    }
        	call run_rscript{
       			 input:
        		    statPlot=statPlot,
           		    data_miss=missingness_plink.data_miss,
           		    psc=bcfstats.psc,
           		    threads=threads
    }
      	  call remove_qc_failed_samples{
      			  input:
           		   inputVCF=inc.outputVCF,
         		   to_keep=run_rscript.to_keep,
           		   threads=threads
    }
      	 call plink_qc_steps{
       			 input:
            	  threads=threads,
            	  modified_vcf_01=remove_qc_failed_samples.modified_vcf_01
    }}
    	call merge{
        		input:
                	threads=threads,
                    vcf=plink_qc_steps.outputVCF,
        }
             
      
    }
