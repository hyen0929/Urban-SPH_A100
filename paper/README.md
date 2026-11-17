# 논문 작성용 자료 패키지

Leipzig 실도시 건물해상 LES + 라그랑지 분산 — 2026-09-08 ~ 2026-10-07 결과 일체.

```
RESULTS.md        ★ 주 문서. 모든 결과와 해석. 여기부터 읽는다
runs_index.md       수행한 런 53건 목록과 한 줄 결론
data/               RESULTS.md 의 모든 수치가 재현되는 원자료
figs/               그림 (PDF 벡터)
scripts/            분석 스크립트
```

## 데이터

| 디렉터리 | 내용 | 주요 파일 |
|---|---|---|
| `data/99_summary/` | **전체 상관 마스터 표** | `master_correlations.csv` (24행), `make_summary.py` |
| `data/08_gaussian/` | **운영 가우시안 모델 대비** | `gaussian_vs_les.csv`, `gaussian_sensitivity.csv` (36 설정) |
| `data/00_setup/` | 소스 좌표·λp·형태 지표·높이장 | `sources_n{13,29,50}.csv`, `lambda_p_by_R.csv`, `heightfield.npz` |
| `data/01_wind/` | 풍향 4방향 | `wind4_n13_transient_metrics.csv`, `wind4_n29_meanflow_metrics.csv`, `wind4_tke_vs_lambda.csv` |
| `data/02_wall/` | 벽면/지면 침적 | `wall4_n50_meanflow_metrics.csv` |
| `data/03_ensemble/` | 난류 실현 변동 | `ensemble3_spread.json` |
| `data/05_numerics/` | 입자수·dt 수렴성 | `particle_convergence.csv`, `dt_convergence.csv` |
| `data/06_stability/` | 안정 성층 | `profile_compare.npz`, `ped_bands_{neutral,stable}.csv` |
| `data/07_other/` | 반증된 가설·방출고·유동 예측자 | `street_orientation.csv`, `release_height.csv`, `flow_predictors.csv` |

`master_correlations.csv` 열: `experiment, n_sources, flow, wall_dep, metric,
rho_0, rho_90, rho_180, rho_270, mean, spread, ci_width_mean`
(`spread` = 네 풍향의 최대−최소, `ci_width_mean` = 부트스트랩 95 % 신뢰구간 폭 평균)

## 그림

| 디렉터리 | 내용 |
|---|---|
| `figs/main/` | 본문 후보 — λp 산점도, 반경 스캔, 노출 지도, 성층 연직 프로파일·평면도 |
| `figs/stability/` | 안정 성층 전체 (f14~f21, f25, f26, f30, f34~f40) |
| `figs/supplementary/` | 0° 기준 분석 36장 (f01~f37) |

## 재현

`data/99_summary/make_summary.py` 가 런 출력(`../runs/`)에서 마스터 표를 다시 만든다.
개별 분석은 `scripts/` 의 스크립트가 각각 담당한다. 런 출력 자체(VTK, 수백 GB)는
저장소에 없으므로 `runs/<런이름>/out/` 을 참조한다.

## 핵심 수치 한눈에

논지가 네 층으로 내려간다.

| 층 | 핵심 수치 |
|---|---|
| **1. 운영모델로 안 된다** | 소스 간 CV — LES **0.408** 대 가우시안 **0.015 ~ 0.102** (36 설정 전부). LES 소스 간 최대/최소 **7.3배** |
| **2. 기하가 순위를 설명한다** | 벽/지면 면적비 → 분배비 ρs **+0.76** (λp 와 기하의 상관 +0.95) |
| **3. 나머지는 풍향** | 기하 예측 대비 소스별 **3 ~ 5배** 편차, 풍향쌍 간 잔차 상관 **+0.11** |
| **4. 잡음이 아니다** | 난류 실현 변동 ρs 폭 **0.03 ~ 0.12** (풍향 변동의 1/3 ~ 1/19) |

형태 지표별 상관(참고):

| | ρs(λp75) | 풍향 폭 | CI 폭 |
|---|---|---|---|
| 벽/지면 분배비 | +0.74 | 0.05 | 0.262 |
| 벽면 침적 | +0.72 | 0.10 | 0.281 |
| 지면 침적 | −0.39 | 0.30 | 0.455 |
| 총 침적 | −0.26 | 0.29 | 0.502 |
