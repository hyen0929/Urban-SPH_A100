#!/usr/bin/env bash
# 2026-10-14 디스크 정리 — 승인된 4개 그룹(377.8 GB)만 삭제한다.
# 각 .txt 가 삭제 대상 목록이며, 목록에 없는 파일은 건드리지 않는다.
# 보존: mean_field.bin 8개, 완주 Run A 의 restart.bin 6개,
#       src29b/src50/src50w 의 LDM_*.vtk, RUN.md·env/·analysis/·*.log 전부.
set -u
cd "$(dirname "$0")/../.." || exit 1          # → /home/kdh/Urban-SPH
M=runs/CLEANUP_2026-10-14

GROUPS=(
  g1a_step0_snapshots      # 52.2 GB  평균장 런의 step-0 스냅샷 (평균장에서 복원 가능)
  g1b_restart_prev         # 61.3 GB  완주 런의 restart.bin.prev (한 세대 전 사본)
  g2_dead_runs_out         # 96.9 GB  폐기·중단 런 5개의 out/
  g3_superseded_out        # 27.4 GB  대체된 런 5개의 out/
  g3b_input_int            #  5.3 GB  Leipzig_surface_int.txt (미사용 내부표시판)
  g4a_0917_LDM             # 104.7 GB 09-17 민감도 6런의 LDM_*.vtk
  g4b_pconv_out            # 11.3 GB  pconv 5런의 out/
  g4c_vorticity            # 24.0 GB  vorticity_*.vtk (모니터링용)
)

for g in "${GROUPS[@]}"; do
  [ -f "$M/$g.txt" ] || { echo "목록 없음: $M/$g.txt — 중단"; exit 1; }
done

echo "삭제 전: $(df -h /home | awk 'NR==2{print "사용 "$3", 여유 "$4" ("$5")"}')"
printf '%s\n' "$(date '+%F %T') 시작" >> "$M/_log.txt"

for g in "${GROUPS[@]}"; do
  n=$(wc -l < "$M/$g.txt")
  before=$(df --output=avail -k /home | tail -1)
  xargs -d '\n' -a "$M/$g.txt" rm -f 2>> "$M/_errors.log"
  after=$(df --output=avail -k /home | tail -1)
  printf '  %-22s %6d개  회수 %.1f GB\n' "$g" "$n" "$(echo "($after-$before)/1048576" | bc -l)"
  printf '%s\n' "$(date '+%F %T') $g $n개" >> "$M/_log.txt"
done

# 비워진 out/ 디렉토리만 제거 (상위 런 디렉토리는 RUN.md 가 있어 남는다)
find runs -mindepth 2 -maxdepth 2 -type d -empty \( -name out -o -name out_cont \) -print -delete

echo "삭제 후: $(df -h /home | awk 'NR==2{print "사용 "$3", 여유 "$4" ("$5")"}')"
echo
echo "다음: .git 의 도달 불가 객체 259.9 GB 회수 (히스토리 보존)"
echo "  git -C $(pwd) gc --prune=now"
