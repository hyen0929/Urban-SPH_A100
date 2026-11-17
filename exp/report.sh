#!/bin/bash
# exp/report.sh — 실행 로그·진단파일에서 수치를 뽑아 RUN.md 3절에 써넣는다.
#
#   exp/report.sh <run>              RUN.md 의 AUTO:METRICS 블록을 갱신 (반복 실행 안전)
#   exp/report.sh <run> --stdout     화면에만 출력
#
# 자동으로 채우는 것은 "측정값"뿐이다. 그 아래 "### 해석" 과 "## 4. 결론" 은
# 사람이(또는 Claude 가) 쓴다 — 이 스크립트는 절대 건드리지 않는다.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

RUN="${1:-latest}"; MODE="${2:-write}"
D=$(resolve_run "$RUN") || exit 1
ID=$(basename "$D"); LOG="$D/run.log"
[ -f "$LOG" ] || die "$ID: run.log 이 없다"

STATE=$(run_state "$D")
STEP=$(last_step "$D"); SIM=$(last_simtime "$D"); WALL=$(last_walltime "$D")
NVTK=$(ls "$D/out"/*.vtk 2>/dev/null | wc -l); SIZE=$(du -sh "$D/out" 2>/dev/null | cut -f1)
BAD=$(grep -icE "nan|segmentation|cudaError|illegal|core dumped|out of memory" "$LOG" 2>/dev/null)
SPEED=$(awk -v s="${SIM:-0}" -v w="${WALL:-0}" 'BEGIN{ if(w>0) printf "%.3f", s/w; else printf "-" }')

# 진행 줄에서 시계열 요약
# 진행 줄에서 시계열 요약. 형식: "<step>\t tu_max= 8.76\ttt_max=300.00\t..."
# ('=' 뒤에 공백이 붙으므로 먼저 제거하고 key=value 로 자른다)
read -r TU_MAX TU_LAST TT_HI TT_LO UZ_MAX < <(
  grep -E '^[0-9]+[[:space:]]+tu_max' "$LOG" 2>/dev/null | awk '
    { line=$0; gsub(/=[ \t]+/,"=",line); n=split(line,f,/[ \t]+/)
      for(i=1;i<=n;i++){ split(f[i],a,"=")
        if(a[1]=="tu_max"){ v=a[2]+0; if(v>tm)tm=v; tl=v; seen=1 }
        else if(a[1]=="uz_max"){ v=a[2]+0; if(v>uz)uz=v }
        else if(a[1]=="tt_max"){ v=a[2]+0; if(!h||v>hi)  {hi=v; h=1} }
        else if(a[1]=="tt_min"){ v=a[2]+0; if(!l||v<lo)  {lo=v; l=1} } } }
    END{ if(!seen){ print "- - - - -"; exit }
         printf "%.3f %.3f %.2f %.2f %.2f", tm, tl, hi, lo, uz }' )

BLOCK=$(mktemp)
{
echo "<!-- AUTO:METRICS  exp/report.sh 가 생성 — 직접 고치지 말 것 ($(date '+%F %T')) -->"
echo
echo "| 측정 | 값 |"
echo "|---|---|"
echo "| 상태 | **$STATE** |"
echo "| 도달 step / 시뮬레이션 시간 | ${STEP:--} / ${SIM:--} s |"
echo "| 계산 시간 (wall) | ${WALL:--} s (= ${SPEED} sim-s/wall-s) |"
echo "| 속도 tu_max (최대 / 최종) | $TU_MAX / $TU_LAST m/s |"
echo "| uz_max (최대) | $UZ_MAX m/s |"
echo "| 온도 범위 tt_min–tt_max | $TT_LO – $TT_HI K |"
echo "| 출력 | $NVTK vtk, $SIZE |"
if [ "${BAD:-0}" -gt 0 ]; then
  echo "| 로그 이상 | ⚠ nan/error 문자열 **$BAD 건** |"
else
  echo "| 로그 이상 | 없음 |"
fi
echo

for f in "$D/out"/PPE_diag_*.csv; do
  [ -f "$f" ] || continue
  echo "**$(basename "$f")** (PPE 잔차 진단)"
  echo
  echo '```'
  awk -F, 'NR==1{next} {n++; bi=$3+0; rs=$4+0; dv=$6+0; sbi+=bi; srs+=rs; sdv+=dv;
                        if(bi>mbi)mbi=bi; if(rs>mrs)mrs=rs; if(dv>mdv)mdv=dv;
                        lbi=bi; lrs=rs; ldv=dv}
       END{ if(n==0){print "(데이터 없음)"; exit}
            printf "샘플 %d개\n", n
            printf "bi_max      평균 %.4g  최대 %.4g  최종 %.4g\n", sbi/n, mbi, lbi
            printf "resid_max   평균 %.4g  최대 %.4g  최종 %.4g\n", srs/n, mrs, lrs
            printf "div_un1_avg 평균 %.4g  최대 %.4g  최종 %.4g\n", sdv/n, mdv, ldv }' "$f"
  echo '```'
  echo
done

if [ -f "$D/out/IISPH_log.txt" ]; then
  echo "**IISPH_log.txt** (반복 수)"
  echo
  echo '```'
  awk '/^#/{next} {n++; it+=$3; if($3>mx)mx=$3; cv+=$5; lit=$3}
       END{ if(n==0){print "(데이터 없음)"; exit}
            printf "step %d개  iter 평균 %.1f  최대 %d  최종 %d  조기수렴률 %.1f%%\n", n, it/n, mx, lit, 100*cv/n }' \
      "$D/out/IISPH_log.txt"
  echo '```'
  echo
fi

if [ -f "$D/out/PPE_iter_log.txt" ]; then
  echo "**PPE_iter_log.txt**"
  echo
  echo '```'
  awk '/^#/{next} {n++; it+=$3; if($3>mx)mx=$3; cv+=$7}
       END{ if(n==0){print "(데이터 없음)"; exit}
            printf "step %d개  iter 평균 %.1f  최대 %d  수렴률 %.1f%%\n", n, it/n, mx, 100*cv/n }' \
      "$D/out/PPE_iter_log.txt"
  echo '```'
  echo
fi

shopt -s nullglob
figs=("$D/analysis"/*.png)
if [ ${#figs[@]} -gt 0 ]; then
  echo "**그림**"
  echo
  for f in "${figs[@]}"; do echo "- \`analysis/$(basename "$f")\`"; done
  echo
fi

if [ "${BAD:-0}" -gt 0 ]; then
  echo "**로그 이상 마지막 5줄**"
  echo
  echo '```'
  grep -inE "nan|segmentation|cudaError|illegal|core dumped|out of memory" "$LOG" | tail -5
  echo '```'
  echo
fi
echo "<!-- /AUTO:METRICS -->"
} > "$BLOCK"

if [ "$MODE" = "--stdout" ]; then cat "$BLOCK"; rm -f "$BLOCK"; exit 0; fi

python3 - "$D/RUN.md" "$BLOCK" <<'PYEOF'
import re, sys
md, blk = sys.argv[1], sys.argv[2]
new = open(blk, encoding="utf-8").read().rstrip("\n")
t = open(md, encoding="utf-8").read()
pat = re.compile(r"<!-- AUTO:METRICS.*?<!-- /AUTO:METRICS -->", re.S)
if not pat.search(t):
    sys.exit("RUN.md 에 AUTO:METRICS 마커가 없다 — 수동으로 붙여넣어라")
open(md, "w", encoding="utf-8").write(pat.sub(lambda _: new, t, count=1))
PYEOF
rm -f "$BLOCK"

cat <<EOF
  $ID 측정값을 $ID/RUN.md 3절에 기록했다.  [$STATE]

  남은 것 (사람이 쓴다):
    RUN.md "### 해석"   — 이 수치가 1절 목적에 답이 되는가
    RUN.md "## 4. 결론" — 한 줄 결론
    runs/INDEX.md       — 그 한 줄을 결론 칸에 옮겨 적기
EOF
