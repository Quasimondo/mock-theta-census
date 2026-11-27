#!/bin/sh
# complete run of the v17 checks; output is run_log.txt
{
echo "== g1_shift.py =="; python3 g1_shift.py 2>&1 | grep -v -i warn
echo; echo "== g2_checks.py a b c =="; python3 g2_checks.py a b c 2>&1 | grep -v -i warn
echo; echo "== g2_checks.py d e =="; python3 g2_checks.py d e 2>&1 | grep -v -i warn
} > run_log.txt
