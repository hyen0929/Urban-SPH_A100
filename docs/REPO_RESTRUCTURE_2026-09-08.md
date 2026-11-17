# 리포 재구성 및 실험 워크플로우 도입 — 2026-09-08

실험을 여러 개 돌리기 전에 (1) 디렉토리가 지저분해지지 않고 (2) 모든 실험이
"왜 돌렸나 / 무엇이 나왔나"를 자기 디렉토리 안에 들고 있도록 구조를 바꾼 기록.

사용법·규약은 [EXPERIMENT_WORKFLOW.md](EXPERIMENT_WORKFLOW.md)에 있다. 이 문서는
**무엇이 어떻게 바뀌었고, 무엇이 남았는지**만 적는다.

---

## 1. 무엇이 문제였나

| 문제 | 상태 (2026-09-08 이전) |
|---|---|
| 실험 설정이 컴파일 타임 | 케이스·스킴 스위치가 전부 `Parameters.cuh` 의 `#define` → 실험마다 소스를 고치고 `make`. **어떤 소스 상태가 어떤 출력을 만들었는지 기록이 남지 않음** |
| 루트 오염 | 출력 디렉토리 25개가 리포 루트에 쌓임 (2.3 TB) |
| `.gitignore` 수동 관리 | 실험 디렉토리 이름을 60줄에 걸쳐 하나씩 등록 중 |
| 실험 기록 부재 | 목적/결과를 남긴 건 `CaseAS_compare/REPORT.md` 하나뿐 |
| 분석 스크립트 산개 | `python/` 에 `analyze_ibm_detail2.py` 류 일회용 20여 개 |

## 2. 핵심 변화 — 실험 설정이 소스에서 빠져나왔다

`Parameters.cuh` 의 실험용 `#define` **25개를 `#ifndef` 가드로 감쌌다. 값은 전부 그대로**
(`#define` 목록·값이 변경 전후 완전 동일함을 diff 로 검증, 기본 빌드 30 초 정상).

```c
#ifndef IMPLICIT_PPE
#define IMPLICIT_PPE 0     // 기본값은 예전과 동일
#endif
```

이제 `-D` 로 덮어쓸 수 있으므로 **실험할 때 `Parameters.cuh` 를 건드리지 않는다.**
작업 트리가 더러워지지 않고, GPU 여러 장에 동시 실험을 돌려도 바이너리가 충돌하지 않는다.

가드 대상: 케이스/출력(`INPUT_CASE_NAME`, `OUTPUT_CASE_NAME`), PPE(`IMPLICIT_PPE`,
`PPE_TOLERANCE`, `PPE_DIAG*`), IISPH(`CONSISTENT_PPE`, `IISPH_*`), 스킴(`RIEMANN_EULERT`),
형상(`DEACTIVATE_IBM`, `IBM_length`), MOST/ABL(`MOST_boundary`, `MOST_z0*`, `us_ini`, `ABL_h`),
온도/스칼라(`T_ini`, `T_con`, `Scalar_Diffusion_coeff`), 채널(`CHANNEL_*`).
전체 표는 [EXPERIMENT_WORKFLOW.md §3](EXPERIMENT_WORKFLOW.md).

`Makefile` 은 `BIN` / `EXP` / `ARCH` 를 받도록 바꿨다. 인자 없이 `make` 하면 예전과 동일.

## 3. 툴킷 `exp/`

```bash
exp/run.sh ppe-consistent -m "IISPH 일치연산자 잔차가 explicit보다 작은지" \
    -i CaseCS1 -t ppe -g 1 -D IMPLICIT_PPE=1 -D CONSISTENT_PPE=1

exp/status.sh          # 전체 상태 한 화면 (RUNNING/DONE/DIED/ARCHIVED, step, sim/wall time)
exp/watch.sh <run>     # 블로우업·NaN·정체·프로세스 사망 때만 종료 (watch_implicit.sh 를 일반화)
exp/report.sh <run>    # 로그·PPE 진단 CSV → RUN.md 3절 측정값 자동 기록
```

- `-m` (목적) **필수**. 없으면 거부한다 — 두 달 뒤 `out/` 만 보고 목적을 복원하는 건 불가능하다.
- `exp/report.sh` 는 `<!-- AUTO:METRICS -->` 블록만 갈아끼운다. 몇 번을 다시 돌려도
  사람이 쓴 "해석"·"결론"은 보존된다. 실행 중에도 돌릴 수 있다.

## 4. 디렉토리 규약

```
runs/<날짜>_<슬러그>/
  RUN.md        목적 · 설정 · 결과 · 결론          ← git 추적
  out/          VTK, PPE_diag*.csv, IISPH_log …    ← git 제외 (재생성 가능)
  analysis/     이 실험 전용 분석 스크립트·그림     ← git 추적
  env/          재현 정보 (정확한 nvcc 명령, Parameters 스냅샷, source.diff,
                solv 사본, 바이너리, git SHA)      ← git 추적 (바이너리 제외)
runs/INDEX.md   실험 대장 — 한 줄 = 한 실험
archive/        2026-09-08 이전 출력 전부           ← git 제외
```

`.gitignore` 가 60줄 디렉토리 목록에서 **`runs/*/out/` 한 줄**로 줄었다.
실험을 추가해도 `.gitignore` 를 다시 고칠 일이 없다.

재현은 `bash runs/<id>/env/build.sh && bash runs/<id>/env/launch.sh`.

## 5. 이사 내역

| 대상 | 이동처 |
|---|---|
| 출력 디렉토리 25개 (AIJ, CaseA*, Kori*, Leipzig*, MOST_Validation, Ph.D. Defence …) | `archive/` |
| `bench_a100.log`, `nohup.out`, `watch_implicit.sh` | `archive/` (감시는 `exp/watch.sh` 가 대체) |
| docx · pptx · `MOST_profiles.png` | `docs/attachments/` |
| `CaseAS_compare/` | `runs/2026-05-19_sgs-ddf-vs-ssm/` — `REPORT.md` → `RUN.md` 로 **이름만** 통일, 내용 그대로. 그림·스크립트는 `analysis/` 로 |
| `python/__pycache__` | 삭제 |

루트에 남은 것: 솔버 소스(`*.cuh`, `*.cu`), `Makefile`, `CLAUDE.md`, `README.md`,
그리고 디렉토리 6개(`archive`, `docs`, `exp`, `input`, `python`, `runs`).

## 6. 검증 범위

GPU 5장이 전부 다른 작업으로 100 % 사용 중이라 **실제 시뮬레이션 실행은 하지 않았다.**

- 검증함: 기본 `make` 정상(30 s), `-D` 오버라이드가 입력·출력 경로까지 반영,
  run 디렉토리·기록 생성, `env/build.sh` 재현 명령, `.gitignore` 추적 범위,
  `report.sh` 파서(실제 로그 포맷의 합성 로그 + PPE 진단 CSV로), `status.sh`
- 미검증: 실제 GPU 실행 경로(`nohup` 기동, PID 추적, `watch.sh` 실전 동작)

첫 실험을 돌릴 때 이 부분만 확인하면 된다.

## 7. 남은 일 / 추천

1. **용량이 진짜 문제다.** 실험 종료 루틴에 `python/cleanup_vtk.py --thin N` 을 넣자
   (기본 dry-run, `--apply` 필요). 결론이 `RUN.md` 에 남으면 VTK 는 지워도 된다는 것이
   이 구조의 요점이다. `archive/` 2 TB도 같은 방식으로 솎을 수 있다.
2. `python/` 의 `analyze_*.py` 20여 개는 대부분 특정 실험 일회용이다. 앞으로 일회용은
   해당 실험 `analysis/` 에, 공용 파서(`vtk_fields.py` 등)만 `python/` 에 남긴다.
   기존 것 분류는 아직 안 했다.
3. explicit vs implicit 처럼 짝으로 도는 실험이 잦으면 `exp/compare.sh <runA> <runB>`
   (두 RUN.md 측정값을 나란히 놓는 비교 리포트) 를 추가할 만하다.

## 8. 미해결 — 출력 디렉토리 4개 소실

작업 도중 루트 출력 디렉토리 4개가 파일시스템에서 사라진 것을 발견했다.

- `Leipzig_surface` (110 GB), `Transport_steady` (55 GB), `Leipzig_small_point` (6.8 GB), `Transport_advection`
- 17:48 `ls` 와 17:50 `du` 시점에는 넷 다 루트에 존재. 18:02 `archive/` 이사 루프에서
  **이 4개만** "moved" 로그가 찍히지 않았고(나머지 20개는 전부 찍힘), 이후 루트·`archive/`·
  `/home` 전체 어디에도 없다. `df` 사용량 11 T → 10 T.
- `input/` 쪽 동명 디렉토리(입력 덱)는 멀쩡하다. 사라진 것은 출력물뿐이다.
- `/home` 은 단일 파일시스템이라 `mv` 는 `rename(2)` 한 번이고, 실패하면 원본이 남지
  사라질 수 없다. cron·실행 중인 `rm`·정리 프로세스도 없었다. 이사 시점에 이미
  존재하지 않았던 것으로 보이나, 확정하지 못했다.
- 휴지통·스냅샷이 없어 복구 경로는 없다.

**확인 필요**: 17:50–18:02 사이 다른 터미널에서 지운 것인지. 아니라면 원인을 더 파야 한다.
