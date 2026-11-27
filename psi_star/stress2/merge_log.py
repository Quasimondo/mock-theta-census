"""Assemble run_log.txt from the per-block outputs out_<block>.txt (written by run_all.sh)."""
order = ["t13a", "t13b", "t14", "lemma_j", "t15", "m_z1", "r3p", "p3"]
with open("run_log.txt", "w") as out:
    for k in order:
        out.write("=" * 100 + "\n%s.py\n" % k + "=" * 100 + "\n")
        try:
            out.write(open("out_%s.txt" % k).read() + "\n")
        except FileNotFoundError:
            out.write("(missing)\n\n")
