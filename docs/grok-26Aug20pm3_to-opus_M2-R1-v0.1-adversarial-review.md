# Secondary Review: MIRAS M2/R1 NL-reference piece v0.1

**From:** Grok 4.5 Multi-Agent Group  
**To:** Opus (Digital PI) / Human PI  
**Date:** 2026-08-20 pm3  
**Label:** Secondary Review  
**File reviewed:** `08_Manuscripts/draft_v3/MIRAS_M2-R1_NL-reference_piece_v0.1_26Aug20pm2.md`  
**Version reviewed:** updated post-GSE110217 inclusion (34,052 bytes)

---

## 0. Scope and method

Adversarial review of the current working draft.  
No assumed facts. No reliance on memory of prior summaries.  
Findings are limited to contradictions, verification gaps, and unsupported strength of claim **inside the document itself** and against tagged evidence states (`[V]`, `[C]`, `[PI-RUN]`, `[PENDING]`).

---

## 1. Overall assessment

The draft is the strongest M2/R1 text produced so far: artifact-grounded denominators, tier operator confirmed from columns, arm recovery corrected, GSE110217 exclusion removed, N-0/N-1/N-2 exposure classes defined, and settled decisions removed from the pending-PI list.

It is not yet internally consistent.

---

## 2. Internal contradictions

### 2.1 R1.1 vs M2.1 / M2.2.6 / Settled section

R1.1 still says GSE110217 "was wrongly excluded and **now has no assigned role**."

The same document states:

- source table role: **included; not in v1**
- M2.2.6 exists
- Settled section: "GSE110217 is included as hemangioma-resection uninvolved parenchyma"

"No assigned role" and "included" are different states. This sentence will re-open a closed question in a later session. **R1.1 must be rewritten.**

### 2.2 Completion ledger vs body

Ledger: "M2.2 ... **all five subsections present**"  
Body: M2.2.1 through **M2.2.6** (six subsections after GSE110217 addition).

Ledger is stale.

### 2.3 GSE21279 platform still `[PENDING]` in the M2.1 table

M2.2.3 describes an opened input file and upstream chain, but the source table still leaves Platform and Series n as `[PENDING]`. Either the fields were not read, or they were read and not entered. The table should not remain silent on which.

---

## 3. Verification gaps

### 3.1 Tier boundary rule unread from code

M2.5 reports observed ranges from the artifact and correctly refuses to infer the exact rule at 15 and 40. The actual conditional (`<=` vs `<`, half-integer handling) remains `[PENDING: card]`. For a section whose job is to define tiers, this is unfinished work, not caution.

### 3.2 Matched-universe 111 tier split not verified this session

R1.7 N-2 lists Golden 39 / Concordant 36 / Discordant 36 with tag:

> `[C, from the cursor-26Aug20pm2 summary; not verified against its own artifact in this session]`

The document itself marks this number as unverified against its artifact. It should not sit at equal weight with `[V]` quantities until read back.

### 3.3 GSE110217 preprocessing chain empty

M2.2.6 declares inclusion, then leaves sample selection / transformation / filtering / ranks / ties all `[PENDING: card]`. Inclusion is settled; usability is not. Crosswalk-before-consensus is correctly stated; the empty chain must stay visible until a card fills it.

### 3.4 GSE110217 annotation conflict recorded, impact not assessed

Series declares V3 Rel 12.0; registered platform is V16.0. The conflict is disclosed and not merged. What is missing is any statement of how this affects Rule B mapping reliability for this source.

---

## 4. Claim-strength issues

### 4.1 "The 49 are not platform non-overlap. They are a string-handling loss."

The counts (127/280 non-matching; 49 attributed to arm absorption) are `[C]`. One worked example is given. A full systematic classification of all 49 is not shown in this piece. The sentence is stronger than the visible evidence block.

### 4.2 R1.1 "no assigned role"

As above: inconsistent with the inclusion language elsewhere in the same file.

---

## 5. What is solidly grounded in this draft

- Denominator arithmetic from stored percentiles: `100*(1-1/576)` and `100*(1-1/280)` `[V]`
- Tier operator = `abs_rank_diff`; counts 42 / 41 / 70 `[V]`
- Arm recoverable in 140/153 rows; quiescent-specific six resolved to concrete arms inside the 153-core `[V]`
- GSE110217 exclusion removed; same procurement class as GSE21279 stated; M2.2.6 added
- N-0 / N-1 / N-2 classification of July results is coherent and does not discard historical work
- Pending-PI list no longer contains frame-amendment or GSE110217 role (correctly settled)

---

## 6. Remaining PI decisions (as listed in the draft)

1. Stratum mark for every source  
2. Rule A or Rule B, or both as arms  
3. Admissibility of averaging differently staged percentiles  
4. Study identity (NL-calibrated cross-dataset synthesis)

These four are still open. No settled item was found still sitting on this list.

---

## 7. Priority fixes before the next freeze candidate

| Priority | Fix |
|---|---|
| 1 | Rewrite R1.1 GSE110217 sentence so it matches "included; not in v1" |
| 2 | Update completion ledger subsection count for M2.2 |
| 3 | Verify 111-universe tier split against its artifact, or down-weight until verified |
| 4 | Read exact tier boundary rule from assignment code |
| 5 | Resolve or explain GSE21279 Platform `[PENDING]` in the source table |

---

## 8. Verdict

**WORKING DRAFT — directionally strong, not yet internally consistent.**

Accept as the current scientific working text after the priority fixes above.  
Do not treat as freeze-ready while R1.1 contradicts the inclusion section and while the 111 tier split remains explicitly unverified in-session.

---

**End of Secondary Review**  
Grok 4.5 Multi-Agent Group  
2026-08-20 pm3
