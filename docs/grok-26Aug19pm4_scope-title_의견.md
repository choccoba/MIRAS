# grok-26Aug19pm4_scope-title_의견.md

**Lane:** Grok 4.5 Multi-Agent Group (Secondary Adversarial / Biomedical Deep-Support)  
**To:** Human PI / Opus (Digital PI)  
**Date:** 2026-08-19  
**Type:** Secondary Review — Recommendation on Decision 1 & 2  
**Project:** MIRAS  
**Label:** Secondary Review

---

## 대상 결정 사항

**Decision 1 — Scope**  
원래 승인된 full scope를 완주하거나, 연구 정체성을 좁히고 누락 데이터셋 각각을 명시적으로 폐기/분리하는 versioned amendment를 발행하거나. 조용한 이행은 선택지가 아닙니다.

**Decision 2 — 연구 정체성**  
meta-analysis / NL-calibrated cross-dataset synthesis / integrative multi-dataset analysis / cross-modal triangulation

---

## 권고 요약

| 결정 | 추천안 | 한 줄 이유 |
|------|--------|-----------|
| **Decision 2** | **NL-calibrated cross-dataset synthesis** | Wang 2026과의 차별점이 가장 명확하고, 과학적으로 가장 정직하다 |
| **Decision 1** | **Versioned amendment로 범위 축소** | Full scope 완주는 현실적으로 어렵고, 정체성을 좁히는 것과 잘 맞는다 |

---

## Decision 2 — 연구 정체성

**추천: `NL-calibrated cross-dataset synthesis`**

| 옵션 | 평가 | 이유 |
|------|------|------|
| **meta-analysis** | 비추천 | 현재 보유한 dataset들이 동일한 estimand를 독립적으로 추정한 결과가 아님. SOL이 지적한 대로 “meta-analysis”라는 이름을 쓸 자격이 아직 없다. |
| **NL-calibrated cross-dataset synthesis** | **강력 추천** | 가장 정직한 표현. Wang 2026이 하지 못한 것(외부 healthy-liver reference를 다른 플랫폼임에도 불구하고 통합)이 바로 이 프레임워크의 핵심이다. |
| integrative multi-dataset analysis | 가능하나 약함 | 너무 일반적이고 차별점이 드러나지 않음. |
| cross-modal triangulation | 가능 | miRNA + mRNA를 강조할 때 보조적으로 쓸 수 있으나, 주 정체성으로는 약하다. |

**핵심 이유**  
Wang 2026 논문이 이미 다음을 보고했다:
- IT와 IC의 hepatic miRNA profile이 매우 유사하다
- IA/ENH의 변화는 상당 부분 liver-resident cell composition에 의해 설명된다
- peripheral viral load가 hepatic miRNA에 잘 반영되지 않는다

따라서 **발견 자체**는 더 이상 새로운 것이 아니다.  
MIRAS의 진짜 기여는 **“NL-calibrated framework”**를 통해 여러 데이터셋을 같은 기준으로 놓고 해석한 점이다.  
이 정체성을 선택하면 논문이 더 좁아지지만, 훨씬 방어 가능해지고 novelty가 명확해진다.

---

## Decision 1 — Scope

**추천: Versioned scientific amendment로 범위를 명시적으로 축소**

이유:

1. **실현 가능성**  
   원래 26Jun MasterPlan에 있던 많은 dataset(GSE76903, GSE77509, GSE83148, GSE84044 등)을 지금 다시 전부 분석하려면 상당한 시간과 노력이 든다. 2027년 상반기 발표 목표와 충돌할 가능성이 높다.

2. **과학적 정직성**  
   “조용히 3개 dataset만 하고 full meta-analysis라고 부르는 것”은 이미 불가능하다고 Opus가 판정했다.  
   범위를 줄이려면 **명시적 amendment**가 필요하다.

3. **정체성과의 정합성**  
   Decision 2에서 “NL-calibrated cross-dataset synthesis”를 선택한다면,  
   “현재 분석 가능한, phenotype annotation이 명확하고 NL calibration이 가능한 dataset”으로 범위를 한정하는 것이 자연스럽다.

**실무적 제안**

- 현재 이미 분석한 dataset (GSE298398, GSE162149, GSE65359 등)을 **Primary synthesis set**으로 명시
- 나머지 MasterPlan dataset은  
  - 향후 확장 가능성으로 남기거나  
  - 명시적으로 “out of current scope”로 폐기  
- 이 결정을 `MIRAS_Scope_Amendment_v1_26Aug19.md` 같은 versioned 문서로 남긴다

---

## 종합 의견

이렇게 결정하면:
- 논문 제목과 North Star 문장을 일관되게 잡을 수 있고
- Claim ceiling을 현실적으로 설정할 수 있으며
- 2027년 상반기 발표 일정을 지킬 가능성이 높아집니다.

**최종 권고**  
Decision 2를 “NL-calibrated cross-dataset synthesis”로, Decision 1을 “versioned amendment로 범위 축소”로 결정하는 것이 현재로서 가장 균형 잡힌 선택입니다.

---

**End of Secondary Review**  
Grok 4.5 Multi-Agent Group  
2026-08-19
