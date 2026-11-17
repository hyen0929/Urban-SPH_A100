# docs — 기술 노트 모음

SOPHIA_gpu 개발 과정에서 작성한 분석·정식화·조사 문서.
프로젝트 지침(`CLAUDE.md`)과 리포 소개(`README.md`)는 루트에 그대로 둔다.

## 실험 워크플로우

| 문서 | 내용 |
|------|------|
| [EXPERIMENT_WORKFLOW.md](EXPERIMENT_WORKFLOW.md) | `exp/` 툴킷과 `runs/` 규약 — 실험 생성·감시·기록 방법 |
| [REPO_RESTRUCTURE_2026-09-08.md](REPO_RESTRUCTURE_2026-09-08.md) | 위 워크플로우 도입 기록 — 무엇이 어디로 옮겨졌고 무엇이 남았는지 |

## 수치 모델 / 정식화

| 문서 | 내용 |
|------|------|
| [KGC_keps_formulation.md](KGC_keps_formulation.md) | KGC 보정행렬과 k-ε 모델의 SPH 이산화 수식 |
| [K_EPSILON_MODEL.md](K_EPSILON_MODEL.md) | 표준 k-ε 난류 모델 구현 — 활성화 경로, 모델 상수 |
| [Temperature_advection.md](Temperature_advection.md) | 온도 advection 폭주 분석과 Vila Riemann-SPH 처방 (§11이 인계 요약) |

## PPE 솔버 조사

| 문서 | 내용 |
|------|------|
| [PPE_explicit_vs_implicit.md](PPE_explicit_vs_implicit.md) | explicit / implicit PPE 비교 설계와 근거 |
| [IISPH_investigation_summary.md](IISPH_investigation_summary.md) | implicit 폭주의 원인 = 연산자 불일치, 일치연산자 도입 경위 |

## 성능

| 문서 | 내용 |
|------|------|
| [BENCHMARK_baseline.md](BENCHMARK_baseline.md) | 최적화 전 성능 기준선과 재현 프로토콜 |
| [THESIS_GPU_NOTES.md](THESIS_GPU_NOTES.md) | GPU 가속 구현 기술 노트 (학위논문 2.3장용) |

## 기타

| 문서 | 내용 |
|------|------|
| [CLAUDE_1.md](CLAUDE_1.md) | 일반 LLM 코딩 가이드라인. 프로젝트 지침이 아니며 자동 로드되지 않는다 |
| [attachments/](attachments/) | 논문·발표 자료 (docx, pptx, png) |
