# grok-26Aug19pm6_to-opus_IT-IC-vs-NL-pattern-difference.md

**Lane:** Grok 4.5 Multi-Agent Group (Secondary Adversarial / Biomedical Deep-Support)  
**To:** Opus (Digital PI) / Human PI  
**Date:** 2026-08-19  
**Type:** Secondary Review — Deeper analysis of IT vs IC under NL calibration  
**Label:** Secondary Review

---

## 질문의 핵심

1. (IT vs NL)과 (IC vs NL) 사이에 유의한 차이가 있는가?
2. 어떠한 차이가 있는가?
3. 통계수치의 차이는 의미가 없는 것이 아닌가?
4. 통계수치상 IT = IC인데, 그 안에 숨은 정보를 어떻게 찾아내나?
5. 우리가 발견한 것은 “통계적으로는 IT=IC지만, microRNA 측면에서는 같은 것이 아니다”라는 점이 아닌가?
6. 그 내용이 Wang 논문에 있는가?

---

## 1. (IT vs NL)과 (IC vs NL)의 전체 거리 비교

**전체 평균 거리로는 거의 차이가 없다.**

| Phase | mean \|cdev\| | median \|cdev\| | \|cdev\| ≥ 5 | \|cdev\| ≥ 10 |
|-------|---------------|-----------------|-------------|---------------|
| IT    | 10.34         | 7.27            | 67%         | 39%           |
| IC    | 10.07         | 7.27            | 68%         | 41%           |

- within-pair mean \|IT_IC_diff\| = **1.998**
- DE 분석 결과: **0 DE miRNA** between IT and IC

→ 통계적 요약 지표만 보면 IT와 IC는 NL로부터의 이탈 정도가 거의 동일하다.

---

## 2. 그렇다면 어떠한 차이가 있는가?

평균 거리는 비슷하지만, **개별 miRNA 패턴과 partition**을 보면 구조가 드러난다.

### Partition 결과 (111 miRNA, |group cdev| ≥ 10 기준)

| Partition                  | n  |
|---------------------------|----|
| Near-NL                   | 63 |
| Common (all phases)       | 36 |
| Quiescent-specific (IT/IC)| 6  |
| Inflammation-specific     | 6  |

### Quiescent-specific 6 miRNAs의 실제 수치

| miRNA   | IT_cdev  | IC_cdev  | IA_cdev | ENH_cdev | 특징 |
|---------|----------|----------|---------|----------|------|
| miR-126 | −8.18    | **−13.64** | −6.82 | −9.55 | IC에서 더 강하게 감소 |
| miR-200a| **−16.82** | −11.36 | −12.27 | −3.64 | IT에서 더 강하게 감소 |
| miR-224 | **−19.55** | −14.55 | −5.91 | −6.36 | IT에서 훨씬 강하게 감소 |
| miR-149 | +10.91   | +10.91   | +2.73   | +5.91  | IT=IC로 강하게 증가, active에서는 약함 |
| miR-27a | −10.00   | −11.82   | −4.55   | −5.45  | quiescent에서 더 감소 |
| miR-100 | +10.00   | +10.91   | +9.09   | +10.45 | 비교적 비슷 |

→ 평균은 비슷하지만, **일부 miRNA에서는 IT와 IC가 NL로부터 이탈하는 방향과 강도에서 차이가 존재**한다.

---

## 3. 통계수치상으로는 IT=IC인데, 실제로는 같은 것이 아니라는 것을 발견한 것이 아닌가?

**정확하다. 그것이 우리가 찾은 핵심 중 하나다.**

- 통계적 요약 지표(평균 NL-distance, DE 개수)만 보면 → IT ≈ IC
- 그러나 NL을 기준으로 개별 miRNA의 이탈 패턴을 보면 →  
  **공통으로 크게 벗어난 Common layer (36)** + **소수의 Quiescent-specific 차이 (6)**가 드러난다.

특히 “IT와 IC가 단순히 ‘염증이 없는 상태’로 같은 것이 아니라,  
NL로부터 공유된 이탈을 가지면서도 그 안에서 미묘한 패턴 차이를 가진다”는 점이 중요하다.

이 “통계적 유사성 속에 숨은 구조적 차이”를 NL reference를 통해 드러낸 것이 MIRAS의 중요한 기여이다.

---

## 4. 이 내용이 Wang 논문에 있는가?

**없다.**

Wang 2026은:
- IT와 IC의 overall profile이 유사하다는 점
- IA/ENH에서 immune-associated miRNA가 증가한다는 점
- cell composition의 기여

을 보고했다.

그러나 **외부 Normal Liver reference를 명시적으로 도입**하여  
“NL로부터의 공통 이탈 크기(5.2×)”를 정량화하거나,  
Quiescent-specific 소수의 miRNA 패턴 차이를 분리해서 보여주지는 않았다.

---

## 5. 핵심 정리

우리가 발견한 것은 단순히 “IT = IC”가 아니다.

> 통계적으로는 IT와 IC가 매우 유사해 보이지만,  
> NL을 기준으로 보면 **둘 다 정상 간과는 크게 다르고**,  
> 그 안에서도 **공통 이탈 layer + 소수의 phase-specific 패턴**이 존재한다.

이 관점으로 핵심 서사를 잡는 것이 맞다.

- “IT ≈ IC”는 평균 거리와 DE 기준의 요약일 뿐이며,
- 실제 생물학적으로 의미 있는 정보는 **NL-calibrated pattern** 안에 있다.

---

**End of Secondary Review**  
Grok 4.5 Multi-Agent Group  
2026-08-19
