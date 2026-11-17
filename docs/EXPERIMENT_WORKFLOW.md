# 실험 워크플로우

목표 두 가지: **디렉토리가 지저분해지지 않을 것**, 그리고 **모든 실험이 "왜 돌렸나 /
무엇이 나왔나"를 자기 디렉토리 안에 들고 있을 것**.

---

## 1. 한 눈에

```bash
# 실험 시작 (빌드 + 실행 + 기록 생성까지 한 번에)
exp/run.sh ppe-consistent -m "IISPH 일치연산자 잔차가 explicit보다 작은지 확인" \
    -i CaseCS1 -t ppe -g 1 -D IMPLICIT_PPE=1 -D CONSISTENT_PPE=1

exp/status.sh                  # 전체 실험 상태 한 화면
exp/watch.sh ppe-consistent    # 이상 생길 때만 종료하는 감시 (블로우업·NaN·정체·프로세스 사망)
exp/report.sh ppe-consistent   # 로그·진단파일에서 수치 뽑아 RUN.md 3절에 기록
```

## 2. 디렉토리 규약

```
runs/
  INDEX.md                     실험 대장 — 한 줄 = 한 실험 (자동 추가, 결론은 손으로)
  2026-09-08_ppe-consistent/
    RUN.md                     ★ 목적 / 설정 / 결과 / 결론      ← git 추적
    out/                       VTK, PPE_diag*.csv, IISPH_log …   ← git 제외 (재생성 가능)
    analysis/                  이 실험 전용 분석 스크립트·그림   ← git 추적
    env/                       재현 정보 일체                     ← git 추적
      build.sh                 이 실험을 만든 정확한 nvcc 명령 (그대로 재실행 가능)
      launch.sh                실행 명령
      flags.txt                -D 플래그 목록
      Parameters.cuh           빌드 시점 스냅샷
      solv.txt                 런타임 설정 사본
      source.diff              빌드 시점 솔버 소스 diff (커밋 안 된 변경 포함)
      provenance.txt           브랜치·SHA·호스트·시각
      SOPHIA_gpu               그 실험을 실제로 돌린 바이너리 (git 제외)
  ...
archive/                       2026-09-08 이전 출력 디렉토리 전부 (git 제외)
```

**`.gitignore` 는 더 이상 손댈 필요가 없다.** `runs/*/out/` 한 줄이 모든 VTK를 걸러낸다.
실험을 새로 만들 때마다 디렉토리 이름을 추가하던 예전 방식은 끝났다.

## 3. 왜 이렇게 되는가 — 컴파일 타임 설정 문제

이 코드는 케이스도 스킴 스위치도 전부 `#define` 이라, 예전에는 실험할 때마다
`Parameters.cuh` 를 고치고 `make` 했다. 그래서 **어떤 소스 상태가 어떤 출력을 만들었는지
기록이 남지 않았다**.

이제 실험용 `#define` 은 `#ifndef` 가드로 감싸져 있어(값은 그대로), `-D` 로 덮어쓸 수 있다:

```c
#ifndef IMPLICIT_PPE
#define IMPLICIT_PPE 0     // 기본값은 예전과 동일
#endif
```

`exp/run.sh` 는 이걸 이용해 **`Parameters.cuh` 를 전혀 건드리지 않고** 실험 전용
바이너리를 `runs/<id>/env/` 에 빌드한다. 결과적으로

- 작업 트리가 실험 때문에 더러워지지 않는다
- 실험 여러 개를 서로 다른 GPU 에서 동시에 돌려도 바이너리가 충돌하지 않는다
- 실험의 설정 = `env/flags.txt` 한 장. 기본값과 다른 것만 거기 적혀 있다

가드가 걸린 것(= `-D` 로 바꿀 수 있는 것):

| 분류 | 매크로 |
|---|---|
| 케이스/출력 | `INPUT_CASE_NAME`, `OUTPUT_CASE_NAME` |
| PPE | `IMPLICIT_PPE`, `PPE_TOLERANCE`, `PPE_DIAG`, `PPE_DIAG_FREQ` |
| IISPH | `CONSISTENT_PPE`, `IISPH_OMEGA`, `IISPH_ETA`, `IISPH_STAG_TOL`, `IISPH_CONV_FREQ`, `IISPH_KGC` |
| 스킴 | `RIEMANN_EULERT` |
| 형상 | `DEACTIVATE_IBM`, `IBM_length` |
| MOST/ABL | `MOST_boundary`, `MOST_z0`, `MOST_z0h`, `us_ini`, `ABL_h` |
| 온도/스칼라 | `T_ini`, `T_con`, `Scalar_Diffusion_coeff` |
| 채널 검증 | `CHANNEL_FX`, `CHANNEL_U0` |
| (이전부터) 평균화·LDM | `AVG_*`, `QS_*`, `LDM_*`, `SPONGE_*` |

가드에 없는 값을 바꾸려면 `Parameters.cuh` 를 직접 고치면 된다 — 그 경우에도
`env/Parameters.cuh` 와 `env/source.diff` 에 스냅샷이 남으므로 기록은 유지된다.

## 4. RUN.md — 실험 하나의 전부

`exp/run.sh` 가 생성하고, 절반은 자동으로 채워진다.

| 절 | 누가 쓰나 | 내용 |
|---|---|---|
| 1. 목적 | `-m` 옵션 (필수) | 왜 돌리는가, 무엇을 보면 성공/실패인가 |
| 2. 설정 | 자동 | 입력 케이스, 플래그, 실행 명령, 재현 방법 |
| 3. 결과 | `exp/report.sh` 가 측정값, 해석은 사람 | 도달 step, 속도·온도 범위, PPE 잔차, 그림 목록 |
| 4. 결론 | 사람 | 한 줄 결론 → `runs/INDEX.md` 로 복사 |

`exp/report.sh` 는 `<!-- AUTO:METRICS -->` 블록만 갈아끼운다. 몇 번을 다시 돌려도
사람이 쓴 "해석"과 "결론"은 그대로다. 실행 중에도 돌릴 수 있다(중간 상태가 찍힌다).

## 5. 규칙 몇 가지

- **목적 없는 실험은 만들 수 없다.** `-m` 이 없으면 `exp/run.sh` 가 거부한다.
  두 달 뒤에 `out/` 만 보고 왜 돌렸는지 복원하는 일은 불가능하다.
- **일회성 분석 스크립트는 그 실험의 `analysis/` 에 둔다.** `python/` 은 여러 실험이
  공유하는 것(VTK 파서 등)만. 예전처럼 `analyze_ibm_detail2.py` 가 루트 근처에
  쌓이지 않게 하는 게 목적이다.
- **용량**: 출력은 금방 수백 GB가 된다. 끝난 실험은 `python/cleanup_vtk.py --thin N`
  으로 솎아낸다(기본 dry-run, `--apply` 필요). 결론이 `RUN.md` 에 적혀 있으면
  VTK 는 지워도 된다 — 그게 이 워크플로우의 요점이다.
- **실험은 항상 리포 루트에서 실행한다.** 출력 경로가 `./runs/<id>/out/` 상대경로다.

## 6. 실험이 끝나면

```bash
exp/report.sh <run>        # 3절 측정값 채우기
# RUN.md 의 "### 해석" 과 "## 4. 결론" 을 쓴다
# 그 한 줄 결론을 runs/INDEX.md 의 결론 칸에 옮긴다
git add runs/<run> && git commit -m "실험: <한 줄 결론>"
```

VTK 를 뺀 기록 전체가 커밋되므로, 나중에 `runs/INDEX.md` 만 훑으면
지금까지 무엇을 왜 시도했고 무엇이 나왔는지 한 번에 보인다.
