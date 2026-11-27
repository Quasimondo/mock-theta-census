#!/bin/sh
# Regenerates run_log.txt (about 5 minutes at 30 digits on 2 cores).
D=${1:-30}
( echo "##### python3 t1_kernel_free.py $D"; python3 t1_kernel_free.py $D ) > log_t1.txt 2>&1 &
( echo "##### python3 t2_kernel.py $D"; python3 t2_kernel.py $D; echo; echo "##### python3 t3_steps.py $D"; python3 t3_steps.py $D; echo; echo "##### python3 t4_series.py 300"; python3 t4_series.py 300; echo; echo "##### python3 t5_bounds.py 20"; python3 t5_bounds.py 20 ) > log_t23.txt 2>&1
wait
cat log_t1.txt > run_log.txt; echo >> run_log.txt; cat log_t23.txt >> run_log.txt; rm -f log_t1.txt log_t23.txt
