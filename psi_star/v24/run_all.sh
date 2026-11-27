#!/bin/sh
# Addendum v24: all checks (about 5 minutes). Usage: sh run_all.sh > run_log.txt
for c in "s1_series.py 300" "s2_law.py" "s3_real.py"; do
  echo "== $c =="; python3 $c 2>&1 | grep -v -i warn; echo
done
