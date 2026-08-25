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

        # ---- Case 3: both files present -> keep the directory, but strip out
        #      only the transcript-pair ROWS in domains.csv that lack full
        #      5'/3' exon coverage. exons.csv itself is left untouched. ----
        exons_basename=\$(basename "\$exons")

        awk -F',' -v file="\$exons" -v exons_out="\$dest/\$exons_basename" -v apos="'" '
            # ---- First read domains.csv: remember which pairs are "expected" ----
            FNR==NR {
                if(NR==1) next
                pair=\$3","\$4
                domain_pairs[pair]=1
                next
            }
            # ---- Now read exons.csv: remember every row, grouped by pair ----
            FNR==1 { header=\$0; next }
            {
                pair=\$3","\$4
                exon_pairs[pair]=1
                exon_lines[pair] = (pair in exon_lines) ? exon_lines[pair] ORS \$0 : \$0
                if(pair in domain_pairs){
                    if(\$7 ~ /5 gene/) five[pair]=1
                    if(\$7 ~ /3 gene/) three[pair]=1
                }
            }
            END{
                print header > exons_out
                for(p in exon_pairs){
                    keep=1
                    # Only pairs that appear in domains.csv are subject to the
                    # 5'/3' coverage check; pairs outside domains.csv pass through
                    # untouched, since they were never part of this validation.
                    if(p in domain_pairs){
                        if(!(five[p] && three[p])) keep=0
                    }
                    if(keep){
                        print exon_lines[p] >> exons_out
                    } else {
                        split(p,a,",")
                        if(!five[p])
                            printf "%s : transcript pair %s,%s missing 5%s exons\\n", file, a[1], a[2], apos
                        if(!three[p])
                            printf "%s : transcript pair %s,%s missing 3%s exons\\n", file, a[1], a[2], apos
                    }
                }
            }' "\$domains" "\$exons" >> problematic_transcripts_report_jaffal.txt
    done

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

        # ---- Case 3: both files present -> keep the directory, but strip out
        #      only the transcript-pair ROWS in domains.csv that lack full
        #      5'/3' exon coverage. exons.csv itself is left untouched. ----
        exons_basename=\$(basename "\$exons")

        awk -F',' -v file="\$exons" -v exons_out="\$dest/\$exons_basename" -v apos="'" '
            # ---- First read domains.csv: remember which pairs are "expected" ----
            FNR==NR {
                if(NR==1) next
                pair=\$3","\$4
                domain_pairs[pair]=1
                next
            }
            # ---- Now read exons.csv: remember every row, grouped by pair ----
            FNR==1 { header=\$0; next }
            {
                pair=\$3","\$4
                exon_pairs[pair]=1
                exon_lines[pair] = (pair in exon_lines) ? exon_lines[pair] ORS \$0 : \$0
                if(pair in domain_pairs){
                    if(\$7 ~ /5 gene/) five[pair]=1
                    if(\$7 ~ /3 gene/) three[pair]=1
                }
            }
            END{
                print header > exons_out
                for(p in exon_pairs){
                    keep=1
                    # Only pairs that appear in domains.csv are subject to the
                    # 5'/3' coverage check; pairs outside domains.csv pass through
                    # untouched, since they were never part of this validation.
                    if(p in domain_pairs){
                        if(!(five[p] && three[p])) keep=0
                    }
                    if(keep){
                        print exon_lines[p] >> exons_out
                    } else {
                        split(p,a,",")
                        if(!five[p])
                            printf "%s : transcript pair %s,%s missing 5%s exons\\n", file, a[1], a[2], apos
                        if(!three[p])
                            printf "%s : transcript pair %s,%s missing 3%s exons\\n", file, a[1], a[2], apos
                    }
                }
            }' "\$domains" "\$exons" >> problematic_transcripts_report_longgf.txt
    done
    """
}