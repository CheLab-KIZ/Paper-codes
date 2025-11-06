#!/bin/bash
#2020.07.08 by ZF (Zhang, Ding et al. 2019)
#2021.12.15 modified by Wenjie Dong

#Most executables are recommended to be added into the environmental paths
#Tools pigz, BBTools, Minia, redundans, Minimap2, samtools, BESST, GapCloser and BUSCO may be used and will be automatically checked prior to formal analyses in this script
#the default starting kmer value is 21 and thus kmer values are 21, 21+20, 21+2*20....
#the statistics of assemblies generated in the asembly procedure are summerized in assembly.statistics


##show help
if [ $# -eq 0 ]; then
    echo -e "\033[33m[HELP - ${0##*/}]\nUsage1: bash $0 config_file     # Import parameters from config file.\nUsage2: bash $0 reads_file_1 reads_file_2 busco_lineage_dir     # Using default parameters\033[0m"
    exit
fi


##Setting parameters
PARAM_NAME_ARR=("_" "data_file_path" "reads_file1" "reads_file2" "dir_BUSCO_lineage" "species" "augustus_species" "read_length" "normalization_target" "threads" "memory" "mean" "stddiv" "norm" "heter" "clean")
declare -A PARAM_HASH

#set default parameters
for i in data_file_path reads_file1 reads_file2 dir_BUSCO_lineage
do
    PARAM_HASH[$i]=""
done
PARAM_HASH[species]="undefined_sp"  #species name as the prefix of resulting files
PARAM_HASH[augustus_species]="chicken"
PARAM_HASH[read_length]=100  #read length of sequencing data
PARAM_HASH[normalization_target]=10  #target normalization depth, which may be ranged from 10 to 20 or higher for low-coverage reads
PARAM_HASH[threads]=4
PARAM_HASH[memory]=4000  #memory size (MB) available for assembly
PARAM_HASH[mean]=500
PARAM_HASH[stddiv]=200
PARAM_HASH[norm]="y"  #execute the step of normalization of sequencing data
PARAM_HASH[heter]="y"  #execute the step of reduction of heterozygous contigs
PARAM_HASH[clean]="n"  #not to clean intermediate files

echo -e "\n=============== START [$(date)] ===============\n" >> log.txt

#set custom parameters
if [[ $# > 2 ]]
    then
        PARAM_HASH[reads_file1]=$1
        PARAM_HASH[reads_file2]=$2
        PARAM_HASH[dir_BUSCO_lineage]=$3
    else
        # read parameters from config file
        CONFIG_FILE=$1
        if [ -s $CONFIG_FILE ]
            then
                sed -i "s/^[ \t]*//g" $CONFIG_FILE  # remove spaces or tabs on the left of lines
                while read LINE
                do
                    if [ "${LINE::1}" != "#" ]; then
                        ARR=(${LINE//=/ })
                        PARAM_HASH[${PARAM_NAME_ARR[${ARR[0]}]}]=${ARR[1]}
                    fi
                done < $CONFIG_FILE
                echo -e "Read parameters from \"$CONFIG_FILE\" ... Done!" | tee -a log.txt
            else
                echo "<ERROR> Config file doesn't exist. Program terminated." | tee -a log.txt
                exit
        fi
fi

if [[ ${PARAM_HASH[data_file_path]} != "" ]]; then
    if [ ${PARAM_HASH[data_file_path]: -1} == "/" ]; then PARAM_HASH[data_file_path]=${PARAM_HASH[data_file_path]%.*};fi
    PARAM_HASH[reads_file1]="${PARAM_HASH[data_file_path]}/${PARAM_HASH[reads_file1]}"
    PARAM_HASH[reads_file2]="${PARAM_HASH[data_file_path]}/${PARAM_HASH[reads_file2]}"
fi

#print parameters
echo "========== Parameters ==========" | tee -a log.txt
i=1
for KEY in ${PARAM_NAME_ARR[*]:1}
do
    echo -e "$i\t$KEY = ${PARAM_HASH[$KEY]}" | tee -a log.txt
    (( i++ ))
done
echo "================================"    | tee -a log.txt


##Checking input
CHECK_OK=1
echo "Check input ..." | tee -a log.txt
i=1
for RDFILE in ${PARAM_HASH[reads_file1]} ${PARAM_HASH[reads_file2]}
do
    if [ -s $RDFILE ]
        then
            echo "Reads file $i: $RDFILE ... OK" | tee -a log.txt
        else
            echo "Reads file $i: $RDFILE ... <ERROR> File doesn't exist." | tee -a log.txt
            CHECK_OK=0
    fi
    (( i++ ))
done
if [ -s ${PARAM_HASH[dir_BUSCO_lineage]} ]
    then
        echo "BUSCO lineage directory: ${PARAM_HASH[dir_BUSCO_lineage]} ... OK" | tee -a log.txt
    else
        echo "BUSCO lineage directory: ${PARAM_HASH[dir_BUSCO_lineage]} ... <ERROR> Directory does'n exist." | tee -a log.txt
        CHECK_OK=0
fi

if [ $CHECK_OK -eq 0 ]; then
    echo "Input error. Please make sure all the input files and directories exist. Program terminated." | tee -a log.txt
    exit
fi


##Checking the package dependency
echo -e "\nChecking the package dependency..." | tee -a log.txt

#check pigz
if [ $(which pigz) ]
    then
        echo "pigz ...... OK" | tee -a log.txt
        EXE_PIGZ=$(which pigz)
        DIR_PIGZ=${EXE_PIGZ%/*}
    else
        echo "pigz ...... <ERROR> The pigz is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check BBtools
if [ $(which bbduk.sh) ]
    then
        echo "BBtools ...... OK" | tee -a log.txt
        EXE_BBTOOLS=$(which bbduk.sh)
        DIR_BBTOOLS=${EXE_BBTOOLS%/*}
    else
        echo "BBtools ...... <ERROR> Scripts of BBtools is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check Minia
if [ $(which minia) ]
    then
        echo "Minia ...... OK" | tee -a log.txt
        EXE_MINIA=$(which minia)
        DIR_MINIA=${EXE_MINIA%/*}
    else
        echo "Minia ...... <ERROR> The minia is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check Redundans
if [ "${PARAM_HASH[heter]}" == "y" ] || [ "${PARAM_HASH[heter]}" == "Y" ]
    then
        if [ $(which redundans.py) ]
            then
                echo "Redundans ...... OK" | tee -a log.txt
                EXE_REDUNDANS=$(which redundans.py)
                DIR_REDUNDANS=${EXE_REDUNDANS%/*}
            else
                echo "Redundans ...... <ERROR> The redundans.py is not found in your PATH." | tee -a log.txt
                CHECK_OK=0
        fi
    else
        echo "Reduction of heterozygous contigs will be skipped."
fi

#check Minimap2
if [ $(which minimap2) ]
    then
        echo "Minimap2 ...... OK" | tee -a log.txt
        EXE_MINIMAP2=$(which minimap2)
        DIR_MINIMAP2=${EXE_MINIMAP2%/*}
    else
        echo "Minimap2 ...... <ERROR> The minimap2 is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check Samtools
if [ $(which samtools) ]
    then
        echo "Samtools ...... OK" | tee -a log.txt
        EXE_SAMTOOLS=$(which samtools)
        DIR_SAMTOOLS=${EXE_SAMTOOLS%/*}
    else
        echo "Samtools ...... <ERROR> The samtools is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check BESST
if [ $(which runBESST) ]
    then
        echo "BESST ...... OK" | tee -a log.txt
        EXE_BESST=$(which runBESST)
        DIR_BESST=${EXE_BESST%/*}
    else
        echo "BESST ...... <ERROR> The runBESST is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check GapCloser
if [ $(which GapCloser) ]
    then
        echo "GapCloser ...... OK" | tee -a log.txt
        EXE_GAPCLOSER=$(which GapCloser)
        DIR_GAPCLOSER=${EXE_GAPCLOSER%/*}
    else
        echo "GapCloser ...... <ERROR> The GapCloser is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi

#check BUSCO
if [ $(which run_BUSCO.py) ]
    then
        echo "BUSCO ...... OK" | tee -a log.txt
        EXE_BUSCO=$(which run_BUSCO.py)
        DIR_BUSCO=${EXE_BUSCO%/*}
    else
        echo "BUSCO ...... <ERROR> The run_BUSCO.py is not found in your PATH." | tee -a log.txt
        CHECK_OK=0
fi
echo -e "\n" >> log.txt

if [ $CHECK_OK -eq 0 ]; then
    echo "Some softwares or scripts are missing. Please make sure all the executables are in your PATH. Program terminated." | tee -a log.txt
    exit
fi


##Prepare data for assembly
mkdir 0-data

#group overlapping reads into clumps and remove duplicates
echo "Clumpify and remove duplicates......" | tee -a log.txt
$DIR_BBTOOLS/clumpify.sh in1=${PARAM_HASH[reads_file1]} in2=${PARAM_HASH[reads_file2]} out1=0-data/1.clumped.fq.gz out2=0-data/2.clumped.fq.gz pigz dedupe 1>>log.txt 2>&1
echo -e "\n"  

#quality trimming
cd 0-data/
echo "Quality trimming......" | tee -a ../log.txt
$DIR_BBTOOLS/bbduk.sh in1=1.clumped.fq.gz in2=2.clumped.fq.gz out1=1.trim.fq.gz out2=2.trim.fq.gz ziplevel=5 pigz ordered qtrim=rl trimq=15 minlen=15 ecco=t maxns=5 trimpolya=10 trimpolyg=10 trimpolyc=10 1>>../log.txt 2>&1
echo -e "\n" >> ../log.txt

#normalize coverage by down-sampling reads over high-depth areas of a genome
if [ "${PARAM_HASH[norm]}" == "y" ] || [ "${PARAM_HASH[norm]}" == "Y" ]
    then
        echo "normalizing coverage of raw data......" | tee -a ../log.txt
        $DIR_BBTOOLS/bbnorm.sh in1=1.trim.fq.gz in2=2.trim.fq.gz out1=1.nor.fq.gz out2=2.nor.fq.gz target="${PARAM_HASH[normalization_target]}" min=2 histcol=2 khist=khist.txt peaks=peaks.txt 1>>../log.txt 2>&1
        echo -e "\n" >> ../log.txt
    else
        mv 1.trim.fq.gz 1.nor.fq.gz
        mv 2.trim.fq.gz 2.nor.fq.gz
fi

mv 1.nor.fq.gz 1.fq.gz && mv 2.nor.fq.gz 2.fq.gz
if [ ${PARAM_HASH[clean]} == "y" ] || [ ${PARAM_HASH[clean]} == "Y" ]; then
    rm -rf *clump* *trim* *raw* *.fastq.gz
fi
cd ..


##Multi-kmer genome assembly

#prepare the input read list
touch reads.list
echo "../0-data/1.fq.gz" >> reads.list
echo "../0-data/2.fq.gz" >> reads.list

#prepare the kmer values used for assembling
touch kmer.list
if [ ${PARAM_HASH[read_length]} -lt 120 ]
    then
        for n in 0 1 2 3; do echo "21+20*$n"|bc; done >> kmer.list
    else
        for n in 0 1 2 3 4 5; do echo "21+20*$n"|bc; done >> kmer.list
fi

#assembly with minia
#KMER values can be adjusted according to the genome size
mkdir 1-assembly && cd 1-assembly/
cp ../reads.list reads.list
for KMER in $(cat ../kmer.list)
do
    echo "Minia assembling at k=$KMER ......" | tee -a ../log.txt
    $DIR_MINIA/minia -in reads.list -kmer-size $KMER -abundance-min 2 -out k$KMER -max-memory ${PARAM_HASH[memory]}
    test -s k$KMER.contigs.fa && echo || ($DIR_MINIA/minia -in reads.list -kmer-size $KMER -abundance-min 2 -storage-type file -out k$KMER -max-memory ${PARAM_HASH[memory]})
    test -s k$KMER.contigs.fa && echo || ($DIR_MINIA/minia -in reads.list -kmer-size `expr $KMER - 2` -abundance-min 2 -out k$KMER -max-memory ${PARAM_HASH[memory]})
    test -s k$KMER.contigs.fa && echo || ($DIR_MINIA/minia -in reads.list -kmer-size `expr $KMER - 2` -abundance-min 2 -storage-type file -out k$KMER -max-memory ${PARAM_HASH[memory]})
    test -s k$KMER.contigs.fa && (echo "Finished Minia k=$KMER" | tee -a ../log.txt) || (echo "Minia assembly with k-mer value of $KMER failed......" && exit)
    echo -e "\n"  | tee -a ../log.txt
    if [ ${PARAM_HASH[clean]} == "y" ] || [ ${PARAM_HASH[clean]} == "Y" ]; then
        rm -rf dummy* *unitigs* *h5 trashme*
    fi
    cp ../reads.list reads.list
    for n in 1 2 3; do echo "./k$KMER.contigs.fa" >> reads.list; done
done
echo "Finished Multi-kmer Minia assembly" | tee -a ../log.txt
echo -e "\n"  | tee -a ../log.txt

#generate basic assembly statistics for all contig assemblies
for kmer in $(cat ../kmer.list); do $DIR_BBTOOLS/statswrapper.sh in="k$kmer".contigs.fa format=6 >> ../assembly.statistics; done

#reduction of heterozygous contigs
if [ "${PARAM_HASH[heter]}" == "y" ] || [ "${PARAM_HASH[heter]}" == "Y" ]
    then
        echo "Reduce heterozygous sequences with Redundans......" | tee -a ../log.txt
        KMER_MAX=$(sed -n '$p' ../kmer.list)
        python2 $DIR_REDUNDANS/redundans.py -v -f k$KMER_MAX.contigs.fa -o ./reduced -t ${PARAM_HASH[threads]} --log log_redundans.txt --noscaffolding --nogapclosing --identity 0.7
        cat log* >> ../log.txt
        mv reduced/scaffolds.reduced.fa k$KMER_MAX.contigs.reduced.fa && rm -rf reduced *fa.fai log* reads.list
        echo -e "\n"  | tee -a ../log.txt
        CONTIGS="k$KMER_MAX.contigs.reduced.fa"
    else
        KMER_MAX=$(sed -n '$p' ../kmer.list)
        CONTIGS="k$KMER_MAX.contigs.fa"
        echo -e "\n"  | tee -a ../log.txt
fi

##Scaffolding assembled contigs
echo "Scaffolding......" | tee -a ../log.txt
#generate the mapping file in sorted, indexed bam format
cd .. && mkdir 2-scaffolding && cd 2-scaffolding/
$DIR_MINIMAP2/minimap2 -ax sr ../1-assembly/$CONTIGS ../0-data/1.fq.gz ../0-data/2.fq.gz -t ${PARAM_HASH[threads]} | $DIR_SAMTOOLS/samtools view -bS -@ ${PARAM_HASH[threads]} | $DIR_SAMTOOLS/samtools sort -@ ${PARAM_HASH[threads]} > map.bam
$DIR_SAMTOOLS/samtools index map.bam -@ ${PARAM_HASH[threads]}

#scaffolding with BESST
$DIR_BESST/runBESST -c ../1-assembly/$CONTIGS -f map.bam -o ./ -orientation fr --iter 10000
test -s BESST_output/pass1/Scaffolds_pass1.fa && echo || (echo "Library parameters cannot be assessed by BESST, try to use alternative parameters. (mean: ${PARAM_HASH[mean]}; std_div: ${PARAM_HASH[stddiv]})" | tee -a ../log.txt) && ($DIR_BESST/runBESST -c ../1-assembly/$CONTIGS -f map.bam -o ./ -orientation fr --iter 10000 -m ${PARAM_HASH[mean]} -s ${PARAM_HASH[stddiv]})
test -s BESST_output/pass1/Scaffolds_pass1.fa && (echo -e "Scaffolding has been finished\n" | tee -a ../log.txt) || (echo "Scaffolding with BESST failed. Please try other tools for scaffolding or use the assembly in directory 1-assembly for BUSCO analyses" && exit)
echo -e "\n"  | tee -a ../log.txt
cd ..


##Gap filling
echo "Gap filling......" | tee -a log.txt
mkdir 3-gapclosing && cd 3-gapclosing/
#prepare config file for GapCloser
touch gapcloser.config
echo "[LIB]" >> gapcloser.config
echo "q1=../0-data/1.fq.gz" >> gapcloser.config
echo "q2=../0-data/2.fq.gz" >> gapcloser.config

#gap filling with GapCloser
$DIR_GAPCLOSER/GapCloser -a ../2-scaffolding/BESST_output/pass1/Scaffolds_pass1.fa -b gapcloser.config -o scaffolds.gapcloser.fa -l ${PARAM_HASH[read_length]} -t ${PARAM_HASH[threads]}
test -s scaffolds.gapcloser.fa && (echo -e "Gap filling has been finished\n" | tee -a ../log.txt) || (echo "Gap filling with GapCloser failed. Please try other tools for gap filling or use the assemblies generated in previous steps for BUSCO analyses" && exit)


#filter the scaffolds shorter than 1000 bp
$DIR_BBTOOLS/reformat.sh in=scaffolds.gapcloser.fa out=../scaffolds_${PARAM_HASH[species]}.fa minlength=1000
cd ..

#generate basic assembly statistics
$DIR_BBTOOLS/statswrapper.sh in=1-assembly/$CONTIGS,2-scaffolding/BESST_output/pass1/Scaffolds_pass1.fa,scaffolds_${PARAM_HASH[species]}.fa >> assembly.statistics

if [ ${PARAM_HASH[read_length]} -lt 120 ]
    then
        for n in 9 7 5 3; do sed -i ""$n"d" assembly.statistics; done
    else
        for n in 13 11 9 7 5 3; do sed -i ""$n"d" assembly.statistics; done
fi

#clean intermediate files
if [ ${PARAM_HASH[clean]} == "y" ] || [ ${PARAM_HASH[clean]} == "Y" ]; then
    rm -rf 2-scaffolding/map* *list
fi


##BUSCO assessment
mkdir 4-busco && cd 4-busco/
python $DIR_BUSCO/run_BUSCO.py -i ../scaffolds_${PARAM_HASH[species]}.fa -c ${PARAM_HASH[threads]} -l ${PARAM_HASH[dir_BUSCO_lineage]} -m geno -o ${PARAM_HASH[species]} -sp ${PARAM_HASH[augustus_species]}
