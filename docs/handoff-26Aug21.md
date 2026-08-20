# MIRAS Handoff — 26Aug21

**Status:** PAUSED  
**Reason:** KDDW autumn abstract deadline **2026-08-24**. Work shifts to CS1 / CS2 for submission first.  
**Authority:** Human PI decision, 26Aug21.  
**Resume condition:** After KDDW abstract submission (and CS1/CS2 completion as decided by PI).

---

## 1. What is frozen in place

| Object | Location / identity | State |
|---|---|---|
| **Working draft** | `08_Manuscripts/draft_v3/MIRAS_M2-R1_NL-reference_piece_v0.3_26Aug20pm6.md` | **Baseline.** 53,163 B; sha256 `7192767d961afd97dcc292d4e765d2f10ba0937c647353671fe042a0c9869739` |
| **v0.1** | `..._v0.1_26Aug20pm2.md` | Retained untouched (first review pin) |
| **v0.2** | Tombstoned | Bytes destroyed by in-place edit; see `..._v0.2_26Aug20pm3_TOMBSTONE.md` |
| **Frame** | `99_Logs/log_opus-miras/opus-26Aug19pm22_FRAME_M2-R1_NL-reference.md` | Frozen structure; GSE110217 exclusion deleted |
| **July Cowork corpus** | `08_Manuscripts/draft_v2/`, `99_Logs/log_cowork-miras/` | Read-only history / comparison arm |

**Immutability rule (active):** A version under review does not change. Corrections go to the next version file. Exception: frame correction rule (false primary-record statements deleted in place and logged).

---

## 2. Scientific position at pause

### Settled

- **NL-reference is the study centre.** External normal-liver coordinate for CHB phase miRNA landscapes.
- **v1 historical consensus** = GSE162149 + GSE31383 only (153 matched core).
- **GSE21279** = living-tissue validation (not primary).
- **GTEx V10** = post-mortem acquisition-stratum comparison (not merged; no GSE accession; dbGaP phs000424.v10).
- **GSE110217** = hemangioma-resection uninvolved parenchyma; **included as candidate**; construction role open; platform unresolved (V3 vs V16.0); preprocessing PENDING.
- **Hemangioma-op NL** is an admitted procurement class (GSE21279 precedent).
- **GSE132763 `HL_adult`** = study-design-defined adult normal liver; construction role unmarked. `HL_young` = developmental sensitivity only.
- **Tier operator is scale-confounded** (Spearman 0.2166 with |pctile diff|; local misorder documented). Historical labels quarantined.
- **Two percentile formulas:** construction `100*(1-rank/n)` vs application `100*(1-(rank-1)/(n-1))`. All N-2 headlines sit on the **application** scale.
- **N-0 / N-1 / N-2 exposure classes** defined; N-1 algebra holds, historical values not portable across a changed universe.
- **core-44** (July variable name `CORE42` held 44 entries) ≠ Tier1_Golden 42.
- **Partition rule:** `|quiescent_mean|≥10` / `|active_mean|≥10` → 63 / 36 / 6 / 6.
- **Ratio 5.2× / 5.0×** = IT and IC mean|cdev| over mean|IT_IC_diff|, not IT-IC vs IA-ENH.
- **Study-identity negative:** classical study-level effect-size meta-analysis of the IT–IC contrast is unavailable (one phase-stratified liver miRNA cohort). Positive label still open.
- **GSE298398** central case study; public release **2026-12-31** → full paper target **2027 H1** (Zenodo + KDDW).

### Not settled (hard gates)

| ID | Gate | Blocks |
|---|---|---|
| **G1** | Stage-aligned Rule A / Rule B comparison chain | operator choice |
| **G2** | Corrected tier operator on the 153 core | any forward tier language |
| **G3** | Independent raw→consensus regeneration **or** PI acceptance standard that levels 1–2 suffice for v1 | freeze-grade authority |

**Standing enforcement (never “closed”):**

- **E1** — every replay uses application percentile formula + `quiescent_mean`/`active_mean` partition at 10
- **E2** — historical raw-rank tier labels quarantined downstream

**Remaining work (not gated the same way):**

- **W1** Fisher OR 7.08 — cell-origin annotation + recompute, or explicit hold
- **W2** GSE110217 platform from Feature Extraction headers
- **W3** Two corrected reference variants vs historical NL-v1 at full R1.7 endpoints

**Rule at pause:** No N-2 number leaves the M2/R1 piece as a result until G1, G2, G3 are closed.

---

## 3. Orchestration state

| Lane | Role in MIRAS |
|---|---|
| Opus 5.0 | Digital PI — frame, run cards, assembly |
| GPT-5.6 SOL | Chief adversarial / numerical reconstruction |
| Cursor (Grok 4.5 ExtraFast) | Execution of run cards |
| Grok 4.5 Multi-Agent | Secondary adversarial, consistency, biomedical support |
| Human PI | Bridge to real world; final decisions |

**Open Cursor work at pause:**

- `opus-26Aug20pm7_to-cursor_M2-R1-v0.4-assembly.md` — base card
- `opus-26Aug21am0_to-cursor_M2-R1-v0.4-addendum.md` — G1/G2/G3 vs E1/E2 vs W1–W3 split; Fisher dual dependency; R1.3 level-3 freeze note

**Do not start v0.4 assembly mid-pause unless PI explicitly reopens MIRAS.** Leave cards as the resume entry point.

---

## 4. Key log trail (recent)

| Log | Content |
|---|---|
| `grok-26Aug20pm3_to-opus_M2-R1-v0.1-adversarial-review.md` | First internal-consistency pass |
| `grok-26Aug20pm5` + `pm5b` | v0.2 review + supplement (formula lock, contamination inventory) |
| `grok-26Aug20pm11_to-opus_M2-R1-v0.3-adversarial-review.md` | v0.3 closure checklist |
| `grok-26Aug21am0_to-opus_M2-R1-v0.3-review-self-audit.md` | Self-audit: “deepest” claim withdrawn; gate wording defect |
| SOL / Cursor parallel reviews | Arithmetic re-derivation from artifacts |

---

## 5. Resume checklist (when MIRAS reopens)

1. Re-pin v0.3 file size + sha256; confirm unchanged.
2. Read `opus-26Aug21am0` addendum; treat G1–G3 / E1–E2 / W1–W3 as the work board.
3. Execute **G1** first (stage-aligned Rule A/B) — blocks operator choice.
4. Then **G2** (corrected tier).
5. Then **G3** (regeneration or PI acceptance standard).
6. Only then: N-2 side-by-side at R1.7 endpoints (W3).
7. Do not re-litigate: GSE110217 inclusion, hemangioma procurement class, GTEx non-merge, N-0 DE=0 portability, frame deletion of false exclusion.

---

## 6. Parallel decision — CS1 / CS2 first

**PI intent (26Aug21):** Complete **CS1** and **CS2** projects first and submit (KDDW abstract deadline 24 Aug 2026).

MIRAS full manuscript remains on the **2027 H1** track (GSE298398 embargo to 2026-12-31). No conflict with pausing MIRAS now.

*This handoff does not define CS1/CS2 scope — those live outside this MIRAS freeze package. Point CS work at their own project folders when active.*

---

## 7. One-line summary for next session

> MIRAS paused at **M2/R1 v0.3 baseline**; hard gates **G1–G3** open; no N-2 re-export until they close; GSE298398 paper holds for 2027 H1; **CS1/CS2 take priority through KDDW abstract 24 Aug 2026**.

---

**End of handoff**  
Grok 4.5 Multi-Agent Group · 2026-08-21  
For: Human PI / Opus / SOL / Cursor on resume
