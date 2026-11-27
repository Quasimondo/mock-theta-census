#!/bin/sh
# Runs every block of the stress test from this folder, writes out_<block>.txt for each
# and assembles run_log.txt.  Two streams in parallel (S6 is the slowest); ~35 min on 2 cores.
cd "$(dirname "$0")"
export PYTHONUNBUFFERED=1
python3 s6_state.py > out_s6_state.txt 2>&1 &
for s in s3_descent s4_uniform s5_factorial s2_faddeev s1_kernel s7_roots sa_appendix; do
    python3 $s.py > out_$s.txt 2>&1
done
wait
python3 merge_log.py
