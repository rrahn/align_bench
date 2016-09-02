#! /bin/bash

######################################################################
# Script to bench mark parallel alignment. 
# Author: Rene Rahn <rene.rahn@fu-berlin.de>
# 
# Arguments
# $3 path to binary build with sse4
# $4 path to binary build with avx2
# $1 number of max threads
# $2 result dir
######################################################################

MAX_THREADS=$1
RES_DIR=$2
SSE4_BIN=$3
AVX2_BIN=$4
RES_FILE="$RES_DIR/res.csv"
PAR_EXEC=("native" "omp" "tbb" "native_vec" "omp_vec" "tbb_vec")
BLOCK_SIZE=(50 100 150 200)
SEQ_MIN=(1000 10000 50000 100000) #290000)
SEQ_MAX=(1100 11000 55000 110000) #310000 )
BITS=(16 32 64)
ALPHA=("dna" "aminoacid")
RUNS=1

CMD=""

# $1: cmd_args
# $2: outfile
function exec_par {
    local par_out="$RES_DIR/$2.out"
    echo "res_par=$CMD $1 -r $RUNS -o $par_out"
    res_par=$($CMD $1 -r $RUNS -o $par_out)
    diff -q $par_out $gold_out
    diff_status=$?
    if [ $diff_status -ne 0 ]; then
        res_par="-1,$res_par"
    else
        res_par="0,$res_par"
    fi
    echo $res_par >> $RES_FILE
    #Compare result of outfile.
}

# $1: cmd args
# $2: outfile
function config_bs {
    for bs in "${BLOCK_SIZE[@]}"; do
        exec_par "$1 -bs $bs" "$2-bs_$bs" 
        #echo "$1 -bs $bs $2-bs_$bs"
    done
}

# $1: cmd args
# $2: outfile
function config_threads {
    for threads in `seq 0 5 $MAX_THREADS`; do
        if [ $threads -eq 0 ]; then
            threads="1"
        fi
        config_bs "$1 -t $threads" "$2-t_$threads"
    done
}

# $1: cmd args
# $2: outfile
function run {
    # run serial execution
    CMD=$SSE4_BIN
    gold_out="$RES_DIR/$2-exec_serial.out"
    echo "res=$CMD $1 -r $RUNS -o $gold_out"
    res=$($CMD $1 -r $RUNS -o $gold_out)
    echo "0,$res" >> $RES_FILE
    
    # run parallel execution
    for exec in "${PAR_EXEC[@]}"; do
        local args="$1 -p $exec"
        local out="$2-exec_$exec"
        config_threads "$args" "$out-arch_sse4"
        if [ -f $AVX2_BIN ]; then
            if [[ $exec == *_vec ]]; then
                CMD="$AVX2_BIN"
                config_threads "$args" "$out-arch_avx2"
            fi
        fi
        CMD="$SSE4_BIN"
    done
}

# $1: bits
# $2: cmd args
# $3: outfile
function config_seq {
    for i in "${!SEQ_MIN[@]}"; do 
        local args="$2 -ml ${SEQ_MIN[$i]} -xl ${SEQ_MAX[$i]}"
        local out="$3-ml_${SEQ_MIN[$i]}-xl_${SEQ_MAX[$i]}"
        if [ "$1" = "16" ]; then
            if [ ${SEQ_MIN[$i]} -lt 16000 ]; then
                run "$args" "$out"
            fi
        else
            run "$args" "$out"
        fi 
    done
}

# $1: cmd args
# $2: outfile
function config_bits {
    for bits in "${BITS[@]}"; do
        config_seq $bits "$1 -sw $bits" "$2-sw_$bits" 
    done
}

# $1: simd instruction set
function config_alpha {
    for alpha in "${ALPHA[@]}"; do
        config_bits "-sa $alpha" "align=sa_$alpha"
    done 
}

mkdir -p "$RES_DIR"
echo "diff_out,exec,run,seqH,seqV,score_width,alphabet,score,time,block_size,threads,vector_length" > "$RES_FILE"

if [ ! -f $SSE4_BIN ]; then
    echo "No executable specified!"
    exit
fi

config_alpha

#if [ -f $SSE4_BIN ]
#then
#    CMD="$SSE4_BIN"
#    OUTFILE="SSE4"
#    config_alpha "SSE4"
#elif [ -f $AVX2_BIN ]
#then
#    CMD="$AVX2_BIN"
#    OUTFILE="AVX2"
#    config_alpha "AVX2"
#fi
#for exec in "${EXEC[@]}"
#do
#    OUTFILE="exec_$exec-"
#    if [ "$exec" == "serial" ] 
#    then
#        invoke_t 1 
#    fi
#    ARGS="-p $exec"
#    invoke_t $MAX_THREADS
#   #:
#   # do whatever on $i
#done
