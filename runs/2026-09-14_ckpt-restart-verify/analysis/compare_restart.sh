#!/bin/bash
# 연속 실행(out_cont/, run_cont.log) vs step 40 체크포인트 재개(out/, run_restart.log) 비트 비교.
# 재개 뒤에 나온 것만 비교한다: PPE_diag step>=40 행, stdout step>=50 블록, step 50/75 덤프,
# mean_field.bin(t=0.6 s 덤프), step 80 체크포인트(전체 상태).
cd "$(dirname "$0")/.." || exit 3
fail=0
chk(){ if [ "$1" = 0 ]; then echo "  같음   $2"; else echo "  다름   $2"; fail=1; fi; }

echo "[PPE_diag step>=40]"
diff <(awk -F, '$1 ~ /^[0-9]/ && $1>=40' out_cont/PPE_diag_explicit.csv) <(awk -F, '$1 ~ /^[0-9]/ && $1>=40' out/PPE_diag_explicit.csv) >/dev/null
chk $? "PPE_diag_explicit.csv ($(awk -F, '$1 ~ /^[0-9]/ && $1>=40' out/PPE_diag_explicit.csv | wc -l) 행)"

echo "[stdout step>=50: tu_max/hot-u/cold-spot/domain 줄]"
ext(){ awk '/^[0-9]+\t tu_max=/{s=$1} s>=50 && /tu_max=|\[hot-u\]|\[cold-spot\]|\[domain\]/' "$1"; }
diff <(ext run_cont.log) <(ext run_restart.log) >/dev/null
chk $? "진단 줄 $(ext run_restart.log | wc -l) 개"

echo "[바이너리]"
for f in fluid_condition_50stp.vtk fluid_condition_75stp.vtk LDM_50stp.vtk LDM_75stp.vtk fluid_75stp.vtk mean_field.bin restart.bin; do
  if [ -f "out/$f" ] && [ -f "out_cont/$f" ]; then cmp -s "out_cont/$f" "out/$f"; chk $? "$f ($(stat -c %s "out/$f") B)"
  else echo "  없음   $f"; fail=1; fi
done
echo
[ $fail = 0 ] && echo "PASS — 재개 결과가 연속 실행과 비트 단위로 같다" || echo "FAIL"
exit $fail
