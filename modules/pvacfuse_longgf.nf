#!/usr/bin/env nextflow

process PvacFuse_LongGF {
    cpus 1
    tag "${name}"

    input:
    tuple val(name), path(fusion_neoag_jaffal)
    tuple val(name), path(agfusions_longgf), path(problematic_transcripts_report), path(missing_exons_report), path(hla_alleles)

    output:
    tuple val(name), path("${name}_fusion_neoag_longgf")

    script:
    """
    #extracting the column alleles HLA from the input tsv file and adding HLA before the allele names to fit pvacseq format
    awk -F'\\t' '
        NR==1 {
            for (i=1; i<=NF; i++) if (\$i=="Allele") col=i
            if (!col) {
                print "ERROR: Allele column not found" > "/dev/stderr"
                exit 1
            }
            next
        }

        \$col != "" {
            allele = \$col

            # remove whitespace
            gsub(/[[:space:]]+/, "", allele)

            # remove HLA- prefix if present
            sub(/^HLA-/, "", allele)

            # keep only 2-field resolution
            if (match(allele, /^([^*]+)\\*([0-9]+:[0-9]+)/, m)) {
                gene = m[1]
                twofield = gene "*" m[2]

                # Class I → add HLA-
                if (gene ~ /^(A|B|C|E|F|G)\$/) {
                    print "HLA-" twofield
                }
                # Class II → no prefix
                else {
                    print twofield
            }
        }
    }
    ' '${hla_alleles}' | sort -u | paste -sd "," - > "${name}_HLA_alleles.txt"


    #running singularity in an example test set, where out6 is an emplty output file and pvacseq_example_data includes exmaple dataset provided by the pVACtools 
    
    HLA_ALLELES=\$(cat ${name}_HLA_alleles.txt)

    pvacfuse run \
        ${agfusions_longgf} \
        ${name} \
        \${HLA_ALLELES} \
        all \
        ${name}_fusion_neoag_longgf \
        --percentile-threshold 2 \
        --iedb-install-directory /opt/iedb

    """
    
    stub:
    """
    mkdir -p ${name}_fusion_neoag_longgf
    """
}
