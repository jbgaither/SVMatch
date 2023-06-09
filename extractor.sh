base=/igm/home/alp020/SVMatch/Gnomad/gnomad_INFO.tsv.gz

## Pulls important information from the tp vcf files
bcftools query -f "%CHROM\t%POS\t%INFO/""SVTYPE""\t%INFO/""SVLEN""\t%INFO/""MatchId""\t%INFO/""TruScore""\n" "$1"/tp-call.vcf > "$1"/tp_call_INFO.tsv
bcftools query -f "%CHROM\t%POS\t%INFO/""SVTYPE""\t%INFO/""SVLEN""\t%ID\t%INFO/""MatchId""\t%INFO/""TruScore""\n" "$1"/tp-base.vcf > "$1"/tp_base_INFO.tsv

$RPath /igm/home/alp020/SVMatch/Scripts/extractor.R "$1"

## sorts the output bed by the call chrom and position
bedtools sort -i "$1"/out.bed > "$1"/connector.bed

## Pulls the AF popMax+AF and N_HOMALT from the full gnomad vcf and sticks it into the tab file
cat "$1"/connector.bed |
while read -r CHROM_call POS_call SVTYPE_call SVLEN_call CHROM_base POS_base SVTYPE_base SVLEN_base gnomadId TruScore
  do
    echo "$CHROM_call" "$POS_call" "$POS_base" "$SVLEN_call" "$SVLEN_base" "$SVTYPE_base" "$TruScore" $(tabix "$base" "$CHROM_base":"$POS_base"-"$POS_base" | grep -w "$gnomadId" | cut -f4,5,6) "$gnomadId" | tr -s " " "\t"
done > "$1"/intermediate.tab

## Puts a header on the top of the tab file
echo "#CHROM POS_query POS_gnomad SVLEN_query SVLEN_gnomad SVTYPE TruScore AF N_HOMALT POPMAX_AF GnomadId" | tr -s " " "\t" > "$1"/out.tab
cat "$1"/intermediate.tab >> "$1"/out.tab

## Removing intermediate files
rm -r "$1"/tp_base_INFO.tsv "$1"/tp_call_INFO.tsv "$1"/out.bed "$1"/connector.bed "$1"/intermediate.tab



