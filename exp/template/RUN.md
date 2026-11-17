# {{ID}}

> {{DATE}} · 브랜치 `{{BRANCH}}` @ `{{SHA}}`{{DIRTY}} · GPU {{GPU}}

## 1. 목적

{{PURPOSE}}

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `{{INPUT_CASE}}` |
| 출력 | `{{OUT_REL}}` |
| 빌드 플래그 | {{FLAGS}} |
| 실행 | `{{CMD}}` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash {{RUNREL}}/env/build.sh && bash {{RUNREL}}/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh {{ID}}` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
