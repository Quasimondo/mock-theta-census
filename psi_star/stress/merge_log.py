"""Assemble run_log.txt from the per-block outputs out_<block>.txt (written by run_all.sh)."""
order = ["s1_kernel", "s2_faddeev", "s3_descent", "s4_uniform", "s5_factorial", "s6_state", "s7_roots", "sa_appendix"]
with open("run_log.txt", "w") as out:
    for k in order:
        out.write("=" * 100 + "\n%s.py\n" % k + "=" * 100 + "\n")
        try:
            out.write(open("out_%s.txt" % k).read() + "\n")
        except FileNotFoundError:
            out.write("(missing)\n\n")
