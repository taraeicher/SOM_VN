#Submit all jobs.
for i in {1..100};
    do 
        qsub -A PAS0272 -v ITER=${i} run_chromosome_iteration.sh
    done 