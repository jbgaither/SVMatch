## This script matches structural variants with the structural variants in the gnomad SV file.
## Expects a TAB file with four columns: CHROM POS SVLEN SVTYPE.
## The CHROM column should be in the format "chr#".
## Currently only INS and DEL are supported.
## SVLEN should all be positive.
## TAB file must be sorted.
## SV's must be from build Hg38/GRCh38

# Activate conda environment
unset PYTHONPATH
. "/igm/home/alp020/miniconda3/etc/profile.d/conda.sh"
conda activate truvari



export svfile=${svfile:-null}
export base=${base:-/igm/home/alp020/Spring2023/liftGnomad/gnomad_SVSites_Hg38_barebones.vcf.gz}
export reference=${reference:-/igm/apps/genomes/Homo_sapiens/GRCh38/no_alt/GRCh38.fa}
export RPath=${RPath:-/igm/apps/R/R-3.6.2_install/bin/Rscript}
export output=${output:-./Truvari_out}
export pctsim=${pctsim:-0.7}
export pctsize=${pctsize:-0.7}
export pctovl=${pctovl:-0.7}
export refdist=${refdist:-500}
export chunksize=${chunksize:-1000}
export multimatch=${multimatch:-true}
export sizemin=${sizemin:-50}
export sizefilt=${sizefilt:-30}
export sizemax=${sizemax:-25000000000}
export giabreport=${giabreport:-false}
export debug=${debug:-false}
export prog=${prog:-false}
export unroll=${unroll:-false}
export minhaplen=${minhaplen:-50}
export typeignore=${typeignore:-false}
export duptoins=${duptoins:-false}
export uselev=${uselev:-false}
export gtcomp=${gtcomp:-false}
export bSample=${bSample:-false}
export cSample=${cSample:-false}
export passonly=${passonly:-false}
export noref=${noref:-false}
export includebed=${includebed:-false}
export extend=${extend:-0}
export help=${help:-false}
export keepfinalvcf=${keepfinalvcf:-false}

while [ $# -gt 0 ]; do

   # Here we are assigning the values of the parameters passed in to their respective variables
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        
        ## Some of the parameters are false by default but if specified they should be made true (but the user doesn't have to explicitly say they are true). They are being assigned true here.
        if [ $param = "help" ] || [ $param = "giabreport" ] || [ $param = "debug" ] || [ $param = "prog" ] || [ $param = "unroll" ] || [ $param = "typeignore" ] || [ $param = "duptoins" ] || [ $param = "uselev" ] || [ $param = "gtcomp" ] || [ $param = "bSample" ] || [ $param = "cSample" ] || [ $param = "passonly" ] || [ $param = "includebed" ] || [ $param = "multimatch" ] || [ $param = "keepfinalvcf" ]
            then
                declare $param=true
        else
            declare $param="$2"
        fi
   fi

  shift
done


if [ $help = true ]
    then
        echo Example command: bash Truvari_new.sh --svfile ./path/to/mySVfile.tab --refdist 1000 --pctsim 0.7 --sizemin 0 --Sizefilt 0 --passonly --chunksize 1000 --sizemax 250000000
elif [ $svfile = null ]
    then
        echo Must pass in a file!
else


  ############################################################################################################################################


  filename=$(echo $svfile | awk -F"/" '{print $NF}' | awk -F".tab" '{print $1}')
  location=$(echo $svfile | awk -F"$filename" '{print $1}')
  auxhead_file=${location}tempanns.tab.hdr
  
  ## Make location "./" if no path is given.
  if [ -z $location  ]
      then
          location=./
  fi
  
  ## Make output directory
  if [ $output = ./Truvari_out ]
      then
          output=${location}Truvari_out
  fi
  
  ## Zip the tab file if it isn't already.
  if [ $(echo -n $svfile | tail -c 3)  != .gz ]
        then
            bgzip -f $svfile
  fi
  
 
  ############################################################################################################################################
  ## Make the VCF file
  
  ## Reading in the tab file and creating relevant values
  zcat ${location}${filename}.tab.gz | grep -v "#" | awk -v OFS='\t' '{print $1, int( $2 ), int( $3 ), $4, int( $2+$3 ), $1"_"$4"_"NR}' > ${location}body.tab
  cat ${location}body.tab | awk -v OFS='\t' '{print $1, int( $2 ), int( $5 ), $6, $4}' > ${location}body2.tab
  
  ## Getting the FASTA sequences
  cat ${location}body2.tab |
  while read -r CHROM POS END ID SVTYPE
    do
      samtools faidx $reference $CHROM:$POS-$END
  done > ${location}temp.tab
  
  ## make the VCF body in truvari format
  $RPath /igm/home/alp020/Spring2023/SV/Truvari/REF_ALT_fixer.R ${location}

  ## Make the header for the VCF
  echo "##fileformat=VCFv4.2" > $auxhead_file
  echo '##ALT=<ID=DEL,Description="Deletion relative to the reference">' >> $auxhead_file
  echo '##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">' >> $auxhead_file
  echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">' >> $auxhead_file
  echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">' >> $auxhead_file
  echo '##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">' >> $auxhead_file
  echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">' >> $auxhead_file
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> $auxhead_file
  echo '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">' >> $auxhead_file
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">' >> $auxhead_file
  echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' >> $auxhead_file
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> $auxhead_file
  echo "#CHROM	POS ID REF ALT QUAL FILTER INFO FORMAT MYSAMP" | tr -s " " "\t" >> $auxhead_file

  ## Put the VCF header and body together
  cat $auxhead_file ${location}new_bod.tab | bgzip -f > ${location}${filename}_truvari.vcf.gz

  ## Index the VCF file
  tabix -f ${location}${filename}_truvari.vcf.gz


  ############################################################################################################################################
  ## Run truvari

  ## Remove the output directory if it exists.
  rm -r ${output}

  ## Writing the truvari command to be executed later.
  truvari='truvari bench -b $base -c ${location}${filename}_truvari.vcf.gz -f $reference -o ${output} -r $refdist -p $pctsim -B $minhaplen -P $pctsize -O $pctovl -C $chunksize -s $sizemin -S $sizefilt --sizemax $sizemax --extend $extend'

  ## Adding extra parameters if they are specified.
  if [ $giabreport = true ]; then truvari="$truvari --giabreport"; fi
  if [ $debug = true ]; then truvari="$truvari --debug"; fi
  if [ $unroll = true ]; then truvari="$truvari --unroll"; fi
  if [ $typeignore = true ]; then truvari="$truvari --typeignore"; fi
  if [ $duptoins = true ]; then truvari="$truvari --dup-to-ins"; fi
  if [ $uselev = true ]; then truvari="$truvari --use-lev"; fi
  if [ $gtcomp = true ]; then truvari="$truvari --gtcomp"; fi
  if [ $bSample = true ]; then truvari="$truvari --bSample"; fi
  if [ $cSample = true ]; then truvari="$truvari --cSample"; fi
  if [ $passonly = true ]; then truvari="$truvari --passonly"; fi
  if [ $includebed = true ]; then truvari="$truvari --includebed"; fi
  if [ $noref != false ]; then truvari="$truvari --no-ref $noref"; fi
  if [ $multimatch = true ]; then truvari="$truvari --multimatch"; fi


  ## Executing truvari
  eval $truvari


  bash /igm/home/alp020/Spring2023/SV/Truvari/extractor.sh ${output}

  ## Remove intermediate files
  if [ $keepfinalvcf = true ]
      then
        rm -r $auxhead_file ${location}body.tab ${location}body2.tab ${location}temp.tab ${location}new_bod.tab
      else
        rm -r $auxhead_file ${location}body.tab ${location}body2.tab ${location}temp.tab ${location}new_bod.tab ${location}${filename}_truvari.vcf.gz ${location}${filename}_truvari.vcf.gz.tbi
  fi
fi

