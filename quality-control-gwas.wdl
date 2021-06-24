workflow qc_wgs{
    File inputVCF
    File statPlot
    #Runtime parameters 
    Int diskSpaceGB
    Int threads
    
    call bcftoolsSubset {
        input:
                inputVCF=inputVCF,
                threads=threads,
                diskSpaceGB=diskSpaceGB

         }  
    call missingness_plink {
        input:
                inputVCF=inputVCF,
                threads=threads,
                diskSpaceGB=diskSpaceGB
    }
    call run_rscript{
        input:
            statPlot=statPlot,
            data_miss=missingness_plink.data_miss,
            psc=bcftoolsSubset.psc,
            threads=threads,
            diskSpaceGB=diskSpaceGB
    }
   } 

   task bcftoolsSubset{
   File inputVCF
   Int threads
   Int diskSpaceGB
   
    command{
    	pwd=`pwd`
        bcftools view -i'QUAL>20 && DP>5' ${inputVCF} |
        bcftools stats -s- > "bcfstats.tsv" 
        grep "#" -v "bcfstats.tsv" | grep "PSC" > "psc.tsv"
        }
 
    runtime{
        docker: "ernfrid/bcftools-1.6:v1" 
        }
    output{
        File bcfstats = "bcfstats.tsv"
        File psc = "psc.tsv"
        }

}
    task missingness_plink{
        File inputVCF
        Int threads
        Int diskSpaceGB

        command{
            plink1.9 --vcf ${inputVCF} --const-fid 0 -missing --out data_miss 
        }
        runtime{
            docker : "biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1"
        }
    output{
        File data_miss = "data_miss.imiss"
        }
}
    task run_rscript {
     File statPlot
     Int threads
     Int diskSpaceGB
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
        }

}