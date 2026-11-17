# dt=0.01 본 런 전체 시계열 (t=0 ~ 300 s)

세 런을 이어 붙인 것이다. 체크포인트에서 분기했기 때문에 상태는 끊김 없이 이어지고,
바뀐 것은 천장 스펀지 설정뿐이다.

| 구간 | step | 시각 [s] | 런 디렉토리 | 천장 스펀지 |
|---|---|---|---|---|
| 1 | 0 ~ 11000 | 0 ~ 110.005 | `2026-09-14_leipzig-dt01-full` | 없음 |
| 2 | 11001 ~ 12000 | 110.015 ~ 120.007 | `2026-09-15_leipzig-dt01-sponge` | uz만 (γ_max 0.2) |
| 3 | 12001 ~ 30000 | 120.017 ~ 299.981 | `2026-09-15_leipzig-dt01-sponge-uy` | uz+uy (γ_max 0.2) |

**구간 1·2 런의 상한 step 이후 출력은 버린 계산이다** (구간 1 은 step 13950, 구간 2 는 13120 까지 돌다가 중단).
여기 링크에는 유효 구간만 들어 있다.

- `slab/` 601장 — `fluid_condition_*stp.vtk` (y∈[0,7.5) 슬래브, 0.5 s 간격)
- `vort/` 30장 — `vorticity_*stp.vtk` (t=155~300 s, 5 s 간격)
- `full/` 3장 — `fluid_*stp.vtk` 전체장 43,046,441 입자 (t=200, 250, 300 s)
- `LDM/` 600장 — `LDM_*stp.vtk` (t=150.25~300 s, 0.25 s 간격). 방출은 t=150 s 부터.
- `PPE_diag_timeline.csv` 3001행 — 세 런의 PPE_diag 를 step 순으로 이어 붙임 + `run` 열
- `domain_timeline.csv` 598행 — run.log 의 50 step 진단(tu_max, uz_max, ftotal_max, Σk_turb, Σ½|u|²) + `run` 열

전부 상대경로 심볼릭 링크다(디스크 추가 사용 없음). 원본을 지우면 끊어진다.
mean-flow 비교용 Run B 는 `../2026-09-15_leipzig-dt01-sponge-uy-meanflow/out/` (LDM 1201장, t=0~300 s, step 번호가 Run A 와 일치).
다시 만들려면: 이 디렉토리를 지우고 `python3 build_timeline.py`.
