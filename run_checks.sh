#!/bin/zsh
# One command: re-run both programs, the controls and the red controls, then check every output against the expected
# values (checks/verify.py). Needs Python 3 with cypari2 (PARI/GP >= 2.15; tested with 2.15.4 and 2.17.2).
#   PY=/path/to/python zsh run_checks.sh          quick: p = 107 with both programs, all controls (about 2 minutes)
#   FULL=1 PY=/path/to/python zsh run_checks.sh   adds p = 599, 1049, 2957 with both programs (about 10 more minutes)
PY=${PY:-python3}
cd "$(dirname "$0")"
O=out/$(date +%Y%m%d_%H%M%S); mkdir -p $O
"$PY" -I -c "import cypari2" || { echo "cypari2 not importable by $PY"; exit 1; }
run() { local name=$1; shift; echo "running $name"; "$@" 2>&1 | grep -v Warning > $O/$name.log; }
run primes        "$PY" -I checks/supersingular_primes.py 9600
run saturation    "$PY" -I checks/saturation.py 12000
run sw122_pari    "$PY" -I checks/control_1483a1.py
run A_1483a1_p5   "$PY" -I programs/pm_route.py 0,1,1,2,2 5 reg "-1,0;0,1"
run A_11a1_p19    "$PY" -I programs/pm_route.py 0,-1,1,-10,-20 19
run A_389a1_p107  "$PY" -I programs/pm_route.py 0,1,1,-2,0 107 reg "-1,1;0,0"
run B_1483a1_p5   "$PY" -I programs/ref_m2ss.py 0,1,1,2,2 5 2 "-1,0;0,1"
run B_11a1_p19    "$PY" -I programs/ref_m2ss.py 0,-1,1,-10,-20 19 0 ""
run B_389a1_p107  "$PY" -I programs/ref_m2ss.py 0,1,1,-2,0 107 2 "-1,1;0,0"
run B_red_p107_scaled_by_p "$PY" -I programs/ref_m2ss.py 0,1,1,-2,0 107 2 "-1,1;0,0" 2 107
if [ "$FULL" = 1 ]; then
  for P in 599 1049 2957; do
    run A_389a1_p$P env PM_SKIP3=1 "$PY" -I programs/pm_route.py 0,1,1,-2,0 $P reg "-1,1;0,0"
    run B_389a1_p$P "$PY" -I programs/ref_m2ss.py 0,1,1,-2,0 $P 2 "-1,1;0,0" 2
  done
fi
"$PY" -I checks/verify.py $O
