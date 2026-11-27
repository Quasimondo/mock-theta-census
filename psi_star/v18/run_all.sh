#!/bin/sh
# complete run of the v18 checks; output is run_log.txt (a1_formula.py takes about 12 minutes, the rest about 2)
{
echo "== a1_formula.py 30 =="; python3 a1_formula.py 30
echo; echo "== b1_roots.py 45 =="; python3 b1_roots.py 45
echo; echo "== b2_dominant.py 90 =="; python3 b2_dominant.py 90
echo; echo "== c1_series.py 300 =="; python3 c1_series.py 300
echo; echo "== c2_z1.py 30 =="; python3 c2_z1.py 30
echo; echo "== d1_compare.py =="; python3 d1_compare.py
} > run_log.txt 2>&1
