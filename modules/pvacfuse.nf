#!/usr/bin/env nextflow

process PvacFuse {
    cpus 1
    tag "${name}"

    input:
    tuple val(name), path(agfusion), path(hla_alleles)

    output:
    tuple val(name), path("${name}_fusion_neoag")

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
        ${agfusion} \
        ${name} \
        \${HLA_ALLELES} \
        all \
        ${name}_fusion_neoag \
        --percentile-threshold 2 \
        --n-threads 1 \
        --iedb-install-directory /opt/iedb

    """
    
    stub:
    """
    mkdir -p ${name}_fusion_neoag
    """
}