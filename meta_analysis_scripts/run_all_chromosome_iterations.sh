#Submit all jobs.
for i in {1..100};
    do 
        qsub -A PAS0272 -v ITER=${i} CELL_LINE=$1 BASE_FNAME=$2 run_chromosome_iteration.sh
    done 