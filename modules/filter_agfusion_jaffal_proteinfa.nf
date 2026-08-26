process FilterAgfusion_Jaffal {
    cpus 1
    tag "${name}"

    input:
    tuple val(name), path(agfusion)

    output:
    tuple val(name), path("agfusion_filtered_jaffal"), path("problematic_transcripts_report_jaffal.txt"), path("missing_exons_report_jaffal.txt")

    script:
    """
    AGFUSION_DIR="${agfusion}"
    AGFUSION="${agfusion}"/*/
    FILTERED_DIR=./agfusion_filtered_jaffal

    mkdir -p "\$FILTERED_DIR"

    > problematic_transcripts_report_jaffal.txt
    > missing_exons_report_jaffal.txt

    shopt -s nullglob
    for d in \${AGFUSION}; do
        [ -d "\$d" ] || continue

        domains=\$(ls "\$d"/*domains.csv 2>/dev/null)
        exons=\$(ls "\$d"/*exons.csv 2>/dev/null)
        protein=\$(ls "\$d"/*_protein.fa 2>/dev/null)

        rel="\${d#\$AGFUSION_DIR/}"
        rel="\${rel%/}"

        # ---- Case 1: no exons.csv at all -> exclude the WHOLE directory ----
        if [ -z "\$exons" ]; then
            echo "Missing *exons.csv in: \$d" >> missing_exons_report_jaffal.txt
            continue
        fi

        dest="\$FILTERED_DIR/\$rel"
        mkdir -p "\$dest"
        cp -r "\$d"/. "\$dest"/

        # ---- Case 2: no domains.csv -> nothing to validate against, keep as-is ----
        if [ -z "\$domains" ]; then
            continue
        fi

        # ---- Determine bad pairs (same 5'/3' coverage check as before) ----
        badpairs="\$dest/.badpairs.tmp"
        awk -F',' -v file="\$exons" -v apos="'" -v badpairs_out="\$badpairs" '
            FNR==NR {
                if(NR==1) next
                pair=\$3","\$4
                domain_pairs[pair]=1
                next
            }
            FNR==1 { next }
            {
                pair=\$3","\$4
                if(pair in domain_pairs){
                    if(\$7 ~ /5 gene/) five[pair]=1
                    if(\$7 ~ /3 gene/) three[pair]=1
                }
            }
            END{
                for(p in domain_pairs){
                    split(p,a,",")
                    if(!(five[p] && three[p])){
                        print p > badpairs_out
                        if(!five[p])
                            printf "%s : transcript pair %s,%s missing 5%s exons\\n", file, a[1], a[2], apos
                        if(!three[p])
                            printf "%s : transcript pair %s,%s missing 3%s exons\\n", file, a[1], a[2], apos
                    }
                }
            }' "\$domains" "\$exons" >> problematic_transcripts_report_jaffal.txt

        # ---- Filter protein.fa: this is what pvacfuse actually iterates over,
        #      so bad-pair FASTA records must be removed from HERE, not from
        #      domains.csv/exons.csv, to actually prevent the crash. ----
        if [ -n "\$protein" ] && [ -f "\$badpairs" ]; then
            protein_basename=\$(basename "\$protein")
            awk -v badpairs_file="\$badpairs" '
                BEGIN {
                    while ((getline line < badpairs_file) > 0) bad[line]=1
                    close(badpairs_file)
                    RS=">"; ORS=""
                }
                NF==0 { next }
                {
                    n = split(\$0, lines, "\\n")
                    header = lines[1]
                    if (match(header, /transcripts: [^,]+/)) {
                        tstr = substr(header, RSTART+13, RLENGTH-13)
                        gsub(/_/, ",", tstr)
                        if (!(tstr in bad)) print ">" \$0
                    } else {
                        print ">" \$0
                    }
                }
            ' "\$protein" > "\$dest/\$protein_basename"
        fi

        rm -f "\$badpairs"

    done
    shopt -u nullglob

    """
}

process FilterAgfusion_LongGF {
    cpus 1
    tag "${name}"

    input:
    tuple val(name), path(agfusion)

    output:
    tuple val(name), path("agfusion_filtered_longgf"), path("problematic_transcripts_report_longgf.txt"), path("missing_exons_report_longgf.txt")

    script:
    """
    AGFUSION_DIR="${agfusion}"
    AGFUSION="${agfusion}"/*/
    FILTERED_DIR=./agfusion_filtered_longgf

    mkdir -p "\$FILTERED_DIR"

    > problematic_transcripts_report_longgf.txt
    > missing_exons_report_longgf.txt

    shopt -s nullglob
    for d in \${AGFUSION}; do
        [ -d "\$d" ] || continue

        domains=\$(ls "\$d"/*domains.csv 2>/dev/null)
        exons=\$(ls "\$d"/*exons.csv 2>/dev/null)
        protein=\$(ls "\$d"/*_protein.fa 2>/dev/null)

        rel="\${d#\$AGFUSION_DIR/}"
        rel="\${rel%/}"

        # ---- Case 1: no exons.csv at all -> exclude the WHOLE directory ----
        if [ -z "\$exons" ]; then
            echo "Missing *exons.csv in: \$d" >> missing_exons_report_longgf.txt
            continue
        fi

        dest="\$FILTERED_DIR/\$rel"
        mkdir -p "\$dest"
        cp -r "\$d"/. "\$dest"/

        # ---- Case 2: no domains.csv -> nothing to validate against, keep as-is ----
        if [ -z "\$domains" ]; then
            continue
        fi

        # ---- Determine bad pairs (same 5'/3' coverage check as before) ----
        badpairs="\$dest/.badpairs.tmp"
        awk -F',' -v file="\$exons" -v apos="'" -v badpairs_out="\$badpairs" '
            FNR==NR {
                if(NR==1) next
                pair=\$3","\$4
                domain_pairs[pair]=1
                next
            }
            FNR==1 { next }
            {
                pair=\$3","\$4
                if(pair in domain_pairs){
                    if(\$7 ~ /5 gene/) five[pair]=1
                    if(\$7 ~ /3 gene/) three[pair]=1
                }
            }
            END{
                for(p in domain_pairs){
                    split(p,a,",")
                    if(!(five[p] && three[p])){
                        print p > badpairs_out
                        if(!five[p])
                            printf "%s : transcript pair %s,%s missing 5%s exons\\n", file, a[1], a[2], apos
                        if(!three[p])
                            printf "%s : transcript pair %s,%s missing 3%s exons\\n", file, a[1], a[2], apos
                    }
                }
            }' "\$domains" "\$exons" >> problematic_transcripts_report_longgf.txt

        # ---- Filter protein.fa: this is what pvacfuse actually iterates over,
        #      so bad-pair FASTA records must be removed from HERE, not from
        #      domains.csv/exons.csv, to actually prevent the crash. ----
        if [ -n "\$protein" ] && [ -f "\$badpairs" ]; then
            protein_basename=\$(basename "\$protein")
            awk -v badpairs_file="\$badpairs" '
                BEGIN {
                    while ((getline line < badpairs_file) > 0) bad[line]=1
                    close(badpairs_file)
                    RS=">"; ORS=""
                }
                NF==0 { next }
                {
                    n = split(\$0, lines, "\\n")
                    header = lines[1]
                    if (match(header, /transcripts: [^,]+/)) {
                        tstr = substr(header, RSTART+13, RLENGTH-13)
                        gsub(/_/, ",", tstr)
                        if (!(tstr in bad)) print ">" \$0
                    } else {
                        print ">" \$0
                    }
                }
            ' "\$protein" > "\$dest/\$protein_basename"
        fi

        rm -f "\$badpairs"

    done
    shopt -u nullglob

    """
}