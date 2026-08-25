process FilterAgfusion {
    cpus 1
    tag "${name}"

    input:
    tuple val(name), path(agfusion)

    output:
    tuple val(name), path("agfusion_filtered"), path("problematic_transcripts_report.txt"), path("missing_exons_report.txt")

    script:
    """
    AGFUSION_DIR="${agfusion}"
    AGFUSION="${agfusion}"/*/
    FILTERED_DIR=./agfusion_filtered

    mkdir -p "\$FILTERED_DIR"

    > problematic_transcripts_report.txt
    > missing_exons_report.txt

    shopt -s nullglob
    for d in \${AGFUSION}; do
        [ -d "\$d" ] || continue

        domains=\$(ls "\$d"/*domains.csv 2>/dev/null)
        exons=\$(ls "\$d"/*exons.csv 2>/dev/null)

        rel="\${d#\$AGFUSION_DIR/}"
        rel="\${rel%/}

        # ---- Case 1: no exons.csv at all -> exclude the WHOLE directory ----
        if [ -z "\$exons" ]; then
            echo "Missing *exons.csv in: \$d" >> missing_exons_report.txt
            continue
        fi

        dest="\$FILTERED_DIR/\$rel"
        mkdir -p "\$dest"
        cp -r "\$d"/. "\$dest"/

        # ---- Case 2: no domains.csv -> nothing to validate against, keep as-is ----
        if [ -z "\$domains" ]; then
            continue
        fi

        # ---- Case 3: both files present -> keep the directory, but strip out
        #      only the transcript-pair ROWS in domains.csv that lack full
        #      5'/3' exon coverage. exons.csv itself is left untouched. ----
        domains_basename=\$(basename "\$domains")

        awk -F',' -v file="\$exons" -v domains_out="\$dest/\$domains_basename" -v apos="'" '
            # ---- First read domains.csv: remember every row, grouped by pair ----
            FNR==NR {
                if(NR==1){ header=\$0; next }
                pair=\$3","\$4
                domain_pairs[pair]=1
                domain_lines[pair] = (pair in domain_lines) ? domain_lines[pair] ORS \$0 : \$0
                next
            }
            # ---- Now parse exons.csv ----
            FNR==1 {next}
            {
                pair=\$3","\$4
                if(pair in domain_pairs){
                    if(\$7 ~ /5 gene/) five[pair]=1
                    if(\$7 ~ /3 gene/) three[pair]=1
                }
            }
            END{
                print header > domains_out
                for(p in domain_pairs){
                    split(p,a,",")
                    if(five[p] && three[p]){
                        print domain_lines[p] >> domains_out
                    } else {
                        if(!five[p])
                            printf "%s : transcript pair %s,%s missing 5%s exons\\n", file, a[1], a[2], apos
                        if(!three[p])
                            printf "%s : transcript pair %s,%s missing 3%s exons\\n", file, a[1], a[2], apos
                    }
                }
            }' "\$domains" "\$exons" >> problematic_transcripts_report.txt

    done

    """
}
