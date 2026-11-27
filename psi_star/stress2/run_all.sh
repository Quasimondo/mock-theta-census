#!/bin/sh
# Runs every block of the second round from this folder, writes out_<block>.txt for each and
# assembles run_log.txt.  Two streams in parallel; about 2 h on 2 cores (r3p.py alone ~ 40 min).
cd "$(dirname "$0")"
export PYTHONUNBUFFERED=1
( for s in r3p t14 m_z1; do python3 $s.py > out_$s.txt 2>&1; done ) &
for s in t13a t13b lemma_j t15 p3; do
    python3 $s.py > out_$s.txt 2>&1
done
wait
python3 merge_log.py
