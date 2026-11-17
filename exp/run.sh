#!/bin/bash
# exp/run.sh — 실험 하나를 생성·빌드·실행하고 기록을 남긴다.
#
#   exp/run.sh <슬러그> -m "목적" [옵션...]
#
# 하는 일:
#   1. runs/<날짜>_<슬러그>/ 생성 (out/, analysis/, env/)
#   2. env/ 에 재현 정보 전부 저장 (git sha, 소스 diff, Parameters.cuh, solv, 빌드/실행 명령, 바이너리)
#   3. -D 플래그로 전용 바이너리 빌드 — Parameters.cuh 는 건드리지 않는다
#   4. 백그라운드 실행, 로그는 run.log
#   5. RUN.md 생성(목적 기록) + runs/INDEX.md 에 한 줄 추가
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

usage() {
  sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
  cat <<'EOF'

옵션:
  -m, --purpose TEXT   (필수) 이 실험을 왜 돌리는지. RUN.md 1절에 기록된다.
  -i, --input CASE     INPUT_CASE_NAME (기본: Parameters.cuh 현재값)
  -D KEY=VAL           컴파일 플래그. 반복 가능. 문자열은 -D NAME='"str"'
  -g, --gpu N          CUDA_VISIBLE_DEVICES (기본 0)
  -R, --restart PATH   체크포인트에서 재개 (SOPHIA_RESTART). 같은 입자수/dt 여야 한다
  -n, --ngpu N         SOPHIA_gpu 실행 인자 = GPU 개수 (기본 1)
  -t, --tag TAG        INDEX.md 분류 태그 (예: ppe, sgs, ldm)
      --fg             포그라운드 실행 (기본: nohup 백그라운드)
      --dry-run        디렉토리·기록·빌드까지만 하고 실행은 안 함
  -h, --help

예:
  exp/run.sh ppe-consistent -m "IISPH 일치연산자 잔차가 explicit보다 작은지 확인" \
      -i CaseCS1 -t ppe -g 1 -D IMPLICIT_PPE=1 -D CONSISTENT_PPE=1

  # 발달장에서 분기해 플래그 하나만 바꿔 A/B
  exp/run.sh wall-on -m "벽 압력 BC 효과" -i Leipzig_wallbc -g 4 -D IBM_PPE_WALL=1 \
      -R runs/2026-09-15_leipzig-dt01-sponge-uy/out/restart.bin
EOF
}

SLUG=""; PURPOSE=""; INPUT_CASE=""; GPU=0; NGPU=1; TAG=""; FG=0; DRY=0; RESTART=""
DFLAGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    -m|--purpose) PURPOSE="$2"; shift 2 ;;
    -i|--input)   INPUT_CASE="$2"; shift 2 ;;
    -D)           DFLAGS+=("-D$2"); shift 2 ;;
    -D*)          DFLAGS+=("$1"); shift ;;
    -g|--gpu)     GPU="$2"; shift 2 ;;
    -R|--restart) RESTART="$2"; shift 2 ;;
    -n|--ngpu)    NGPU="$2"; shift 2 ;;
    -t|--tag)     TAG="$2"; shift 2 ;;
    --fg)         FG=1; shift ;;
    --dry-run)    DRY=1; shift ;;
    -h|--help)    usage; exit 0 ;;
    -*)           die "모르는 옵션: $1  (exp/run.sh -h)" ;;
    *)            [ -z "$SLUG" ] && SLUG="$1" || die "인자가 너무 많다: $1"; shift ;;
  esac
done

[ -n "$SLUG" ]    || { usage; exit 1; }
[ -n "$PURPOSE" ] || die "-m \"목적\" 은 필수다. 목적 없는 실험은 나중에 해석할 수 없다."
[[ "$SLUG" =~ ^[A-Za-z0-9._-]+$ ]] || die "슬러그는 영숫자/._- 만: $SLUG"

cd "$REPO_ROOT"

# ── run 디렉토리 (같은 날 같은 슬러그면 _2, _3 …) ──
BASE="$(date +%Y-%m-%d)_$SLUG"; ID="$BASE"; n=2
while [ -e "$RUNS_DIR/$ID" ]; do ID="${BASE}_$n"; n=$((n+1)); done
RUN_DIR="$RUNS_DIR/$ID"; RUNREL="runs/$ID"
mkdir -p "$RUN_DIR/out" "$RUN_DIR/analysis" "$RUN_DIR/env"

# ── 입력 케이스 기본값 = Parameters.cuh 현재값 ──
if [ -z "$INPUT_CASE" ]; then
  INPUT_CASE=$(grep -oP '^#define\s+INPUT_CASE_NAME\s+"\K[^"]+' Parameters.cuh | head -1)
  [ -n "$INPUT_CASE" ] || die "Parameters.cuh 에서 INPUT_CASE_NAME 을 못 읽었다. -i 로 지정하라."
fi
SOLV="input/$INPUT_CASE/${INPUT_CASE}_solv.txt"
[ -f "$SOLV" ] || die "solv 파일이 없다: $SOLV  (-i 케이스 이름 확인)"
[ -f "input/$INPUT_CASE/${INPUT_CASE}.txt" ] || die "입자 입력이 없다: input/$INPUT_CASE/${INPUT_CASE}.txt"
if [ -n "$RESTART" ]; then
  [ -f "$RESTART" ] || die "체크포인트가 없다: $RESTART"
  # 재개는 '같은 바이너리' 가 원칙이지만 헤더 검증은 구조체 크기·입자수·셀수·dt 만 본다.
  # 즉 플래그만 다른 빌드로는 재개된다 (A/B 에 쓰는 이유). 입력 케이스가 다르면 거부된다.
  info "체크포인트에서 재개: $RESTART"
fi

# ── 케이스/출력 플래그는 항상 주입 ──
ALL_FLAGS=("${DFLAGS[@]+"${DFLAGS[@]}"}"
           "-DINPUT_CASE_NAME=\"$INPUT_CASE\""
           "-DOUTPUT_CASE_NAME=\"$RUNREL/out\"")
EXP_STR=$(printf '%q ' "${ALL_FLAGS[@]}")
BINPATH="$RUNREL/env/SOPHIA_gpu"   # 리포 루트 기준 상대경로 — env/build.sh 가 이식 가능해진다

# ── 재현 정보 ──
BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "-")
SHA=$(git rev-parse --short HEAD 2>/dev/null || echo "-")
DIRTY=""; [ -n "$(git status --porcelain 2>/dev/null)" ] && DIRTY=" +uncommitted"
# 솔버 소스만 (입자 입력 덱은 수백 MB 라 제외)
git diff HEAD -- '*.cuh' '*.cu' Makefile > "$RUN_DIR/env/source.diff" 2>/dev/null || true
git status --porcelain       > "$RUN_DIR/env/git_status.txt" 2>/dev/null || true
cp Parameters.cuh "$RUN_DIR/env/Parameters.cuh"
cp "$SOLV"        "$RUN_DIR/env/solv.txt"
printf '%s\n' "${ALL_FLAGS[@]}" > "$RUN_DIR/env/flags.txt"
{ echo "branch=$BRANCH"; echo "sha=$(git rev-parse HEAD 2>/dev/null || echo -)";
  echo "dirty=${DIRTY:-no}"; echo "host=$(hostname)"; echo "created=$(date '+%F %T')";
  echo "restart=${RESTART:--}"; } > "$RUN_DIR/env/provenance.txt"
nvidia-smi --query-gpu=index,name,memory.total --format=csv > "$RUN_DIR/env/gpu.txt" 2>/dev/null || true

# ── 빌드 ──
info "빌드 중 → $RUNREL/env/SOPHIA_gpu"
{ echo '#!/bin/bash'; echo "# $ID 재현 빌드"; echo "set -e"; echo "cd \"\$(dirname \"\$0\")/../../..\"";
  make -n BIN="$BINPATH" EXP="$EXP_STR"; } > "$RUN_DIR/env/build.sh"
chmod +x "$RUN_DIR/env/build.sh"
make BIN="$BINPATH" EXP="$EXP_STR" > "$RUN_DIR/env/build.log" 2>&1 \
  || { tail -20 "$RUN_DIR/env/build.log"; die "빌드 실패 — 로그: $RUNREL/env/build.log"; }

# ── 실행 명령 ──
RUN_ENV=(env "CUDA_VISIBLE_DEVICES=$GPU")
PFX="CUDA_VISIBLE_DEVICES=$GPU"
if [ -n "$RESTART" ]; then RUN_ENV+=("SOPHIA_RESTART=$RESTART"); PFX="$PFX SOPHIA_RESTART=$RESTART"; fi
CMD="$PFX $RUNREL/env/SOPHIA_gpu $NGPU"
{ echo '#!/bin/bash'; echo "# $ID 재현 실행 (리포 루트에서)"; echo "cd \"\$(dirname \"\$0\")/../../..\"";
  echo "$PFX exec $RUNREL/env/SOPHIA_gpu $NGPU"; } > "$RUN_DIR/env/launch.sh"
chmod +x "$RUN_DIR/env/launch.sh"

# ── RUN.md ──
FLAGS_MD=$( ([ ${#DFLAGS[@]} -eq 0 ] && echo "(기본값)" || printf '`%s` ' "${DFLAGS[@]}") | sed 's/ $//')
python3 - "$REPO_ROOT/exp/template/RUN.md" "$RUN_DIR/RUN.md" <<PYEOF
import sys
src, dst = sys.argv[1], sys.argv[2]
sub = {
 "ID": """$ID""", "DATE": """$(date '+%F %T')""", "BRANCH": """$BRANCH""",
 "SHA": """$SHA""", "DIRTY": """$DIRTY""", "GPU": """$GPU""",
 "PURPOSE": """$PURPOSE""", "INPUT_CASE": """$INPUT_CASE""",
 "OUT_REL": """$RUNREL/out""", "FLAGS": """$FLAGS_MD""",
 "CMD": """$CMD""", "RUNREL": """$RUNREL""",
}
t = open(src, encoding="utf-8").read()
for k, v in sub.items():
    t = t.replace("{{%s}}" % k, v)
open(dst, "w", encoding="utf-8").write(t)
PYEOF

# ── INDEX.md ──
[ -f "$INDEX_MD" ] || cat > "$INDEX_MD" <<'IDXEOF'
# 실험 대장

`exp/run.sh` 가 자동으로 한 줄씩 추가한다. "결론" 칸은 실험이 끝난 뒤
각 `RUN.md` 4절에서 한 줄로 옮겨 적는다. 상태는 `exp/status.sh` 로 갱신.

| 날짜 | 실험 | 태그 | 목적 | 결론 |
|---|---|---|---|---|
IDXEOF
PURPOSE_1L=$(echo "$PURPOSE" | tr '\n' ' ' | sed 's/|/\\|/g' ); PURPOSE_1L=${PURPOSE_1L:0:90}   # 문자 단위 잘림(UTF-8 안전)
echo "| $(date +%F) | [$ID]($ID/RUN.md) | ${TAG:--} | $PURPOSE_1L | _(진행 중)_ |" >> "$INDEX_MD"

if [ "$DRY" = 1 ]; then
  info "--dry-run: 실행 생략. 준비 완료 → $RUNREL"
  exit 0
fi

# ── 실행 ──
if [ "$FG" = 1 ]; then
  info "포그라운드 실행: $CMD"
  "${RUN_ENV[@]}" "$BINPATH" "$NGPU" 2>&1 | tee "$RUN_DIR/run.log"
else
  "${RUN_ENV[@]}" setsid nohup "$BINPATH" "$NGPU" > "$RUN_DIR/run.log" 2>&1 &
  echo $! > "$RUN_DIR/env/pid"
  sleep 2
  kill -0 "$(cat "$RUN_DIR/env/pid")" 2>/dev/null || { tail -20 "$RUN_DIR/run.log"; die "시작 직후 죽었다 — $RUNREL/run.log"; }
fi

cat <<EOF

  실험 시작:  $ID
  목적:       $PURPOSE
  기록:       $RUNREL/RUN.md
  출력:       $RUNREL/out/     로그: $RUNREL/run.log
  PID:        $(cat "$RUN_DIR/env/pid" 2>/dev/null || echo -)   GPU: $GPU

  진행 확인:  exp/status.sh
  감시:       exp/watch.sh $ID
  결과 기록:  exp/report.sh $ID
EOF
