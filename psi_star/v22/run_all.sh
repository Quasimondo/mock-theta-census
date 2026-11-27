#!/bin/sh
# complete run of the v22 checks (about 15 minutes); output is run_log.txt
{
echo "== m1_contour.py (first attempt: one straight ray from i pi) =="; python3 m1_contour.py
echo; echo "== m2_contour.py (the contour of Theorem 13) =="; python3 m2_contour.py
echo; echo "== n1_zroots.py (Theorem 14) =="; python3 n1_zroots.py 2>&1 | grep -v -i warn
} > run_log.txt 2>&1
