#!/usr/bin/env nextflow

// processes to include: input check, filter short reads, assembly, sourmash analysis


//display a help-msgs
if (params.help) { exit 0, helpMSG() }

    if (!params.fasta && !params.fastq && !params.genomesize){
    exit 1, "input missing, use [--fastq] and [--genomesize 5m] or just [--fasta]"
    }

        if (params.fastq && !params.genomesize) {
        exit 1, "input missing, use [--fastq] and [--genomesize 5m]"
        }



//xxxxxxxxxxxxxx//
//***fastq Channel***//
//xxxxxxxxxxxxxx//


//für inputcheck: ich muss schauen das der gewünschte input ture ist... aber die anderen inputmöglichkeiten müssen auf false gesetzt werden!

if (params.fastq && params.genomesize && !params.fasta){    
    fastq_ch = Channel  .fromPath( params.fastq, checkIfExists: true)
                        .map { file -> tuple(file.baseName, file) }
                        .view()


            // split channel to get used by input_check and nanoplot
            fastq_ch.into {input_check_ch; nanoplot_ch }
}

//xxxxxxxxxxxxxx//
//***fasta Channel***//
//xxxxxxxxxxxxxx//

else if (params.fasta && !params.fastq && !params.genomesize){
fasta_ch = Channel.fromPath( params.fasta, checkIfExists: true)
                  .map { file -> tuple(file.baseName, file) }
                  .view()


            fasta_ch.into {sourmash_analysis_ch; xx}
}



    if (!params.fasta){

        //xxxxxxxxxxxxxx//
        //***valid input?***//
        //xxxxxxxxxxxxxx//



        //when inputfile is named: foo.fasta.gz     i get some weird doublenaming like foo.fasta.fasta? because ${name} = foo.fasta
        // this is also in WTP the case
        process input_check {

            input:
            tuple val(name), file(fastq) from input_check_ch

            output:
            tuple val(name), file("${name}.fastq") into filterfastq_ch
            
            script:
            """
            case "${fastq}" in
                    *.gz) 
                        zcat ${fastq} > ${name}.fastq
                        ;;

                    *.fastq)
                        ;;

                    *)
                        echo "file format not supported...(.fastq .gz is supported)"
                        exit 1
                esac
            """
        }



        //xxxxxxxxxxxxxx//
        //***nanoplot***//
        //xxxxxxxxxxxxxx//
        process nanoplot {

            publishDir 'results/nanoplot_results'

            input:
            tuple val(name), file(fastqnp) from nanoplot_ch
            output:
            tuple val(name), file("*.html")
            tuple val(name), file("*.png")
            tuple val(name), file("*.pdf")
            tuple val(name), file("${name}_read_quality.txt")
            
            script:
            """
            NanoPlot --fastq ${fastqnp} --title ${name} --color darkslategrey --N50 --plots hex --loglength -f png --store
            NanoPlot --pickle NanoPlot-data.pickle --title ${name} --color darkslategrey --N50 --plots hex --loglength -f pdf
            mv Nano*.html ${name}_read_quality_report.html
            mv NanoStats.txt ${name}_read_quality.txt
            """
        }



        //xxxxxxxxxxxxxx//
        //***filter short fastq***//
        //xxxxxxxxxxxxxx//

        // filterfastq should move into the if statement
        if (params.nofilter)
        process filtlong {

            publishDir 'results/filter_results/'
            
            input:
            tuple val(name), file(filterfastq) from filterfastq_ch

            output:
            tuple val(name), file("${name}_filtered.fastq") into assembly_ch
            
            script:
            """
            filtlong --min_length 2000 ${filterfastq} > ${name}_filtered.fastq
            """
        }



        //xxxxxxxxxxxxxx//
        //***Assembly***//
        //xxxxxxxxxxxxxx//


        assembly = file("results/assembly_results/")
        if (!assembly.exists()) { 
            process flye {

                publishDir 'results/assembly_results'
            
                input:
                tuple val(name), file(fastq) from assembly_ch

                output:
                tuple val(name), file("${name}.fasta") into sourmash_analysis_ch
                tuple val(name), file(fastq), file("${name}.gfa"), file ("${name}.gv")

                script:
                """
                flye --nano-raw ${fastq} --genome-size ${params.genomesize} --threads 4 -o Assembly_results
                mv Assembly_results/assembly.fasta ${name}.fasta
                mv Assembly_results/assembly_graph.gfa ${name}.gfa
                mv Assembly_results/assembly_graph.gv ${name}.gv
                """
            }
        }

            else{sourmash_analysis_ch = Channel.fromPath("results/assembly_results/*.fasta")
                                                            .map { file -> tuple(file.baseName, file) }
                                                            .view()
            }


    }
//else : jump direct to after sourmash after fasta input
else {
//xxxxxxxxxxxxxx//
//***sourmash***//
//xxxxxxxxxxxxxx//


    sourmash_mg_db = file ("autodownload-database/sourmash_mg_db")
    if (!sourmash_mg_db.exists()) {


        process sourmashautodownloaddatabase {

            storeDir 'autodownload-database/sourmash_mg_db'

            output:
            file("genbank-k31.lca.json.gz") into sourmash_mg_db_ch

            script:

            """
            curl -L -o genbank-k31.lca.json.gz https://osf.io/4f8n3/download
            """
        }
    }

        else { sourmash_mg_db_ch = Channel.fromPath("autodownload-database/sourmash_mg_db/genbank-k31.lca.json.gz")
        }
                                                                    
            process sourmash {
                        
                publishDir 'results/sourmash_results/'

                input:
                tuple val(name), file(assembly) from sourmash_analysis_ch
                file(databasemg) from sourmash_mg_db_ch

                output:

                tuple val(name), file("${name}_sourmashresult.csv")

                script:
                """
                sourmash compute --scaled 10000 -k 31 ${assembly} -o ${name}.sig 
                sourmash gather -k 31 ${name}.sig ${databasemg} -o ${name}_sourmashresult.csv
                """
            }

}

/*************  
* --help
*************/
def helpMSG() {

    c_blue = "\033[0;36m";
    c_light_magenta = "\033[0;95m";
    c_reset = "\033[0m";

    log.info """
    .
    ${c_light_magenta}how to use this small pipeline:${c_reset}
    nextflow run training.nf --fastq '*/*.fastq' --genomesize 5m [your genomesize]

    you can also use *.gz fastq-files


    ${c_blue}Input:${c_reset}
    ${c_blue} --fastq ${c_reset}            '*.fastq'   -> read file(s)
    ${c_blue} --fastq ${c_reset}       readfile.fastq   -> read file
    
    """.stripIndent()
}