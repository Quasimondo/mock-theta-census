import sys,json
for l in sys.stdin:
  try: r=json.loads(l)
  except: print(l.strip()[:200]); continue
  if 'error' in r: print(r); continue
  print(f"{r['name'][:22]:22s} {r['h']}/{r['m']}  kappa={r['kappa']:+.10f} ~{r['kappa_rat']:>7s} err={r['kappa_err']:.1e} w~{r['weight_rat']:>4s} werr={r['weight_err']:.1e} g~{r['gamma_rat']:>6s} gerr={r['gamma_err']:.0e} resid={r['resid']:.1e} dk={r['dkappa']:.1e} prec={r['maxprec']} {r['secs']}s")
