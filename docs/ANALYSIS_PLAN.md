# Statistical Analysis Plan

*Written and committed before any outcome was computed (see git history). Changes after outcomes were seen are
listed at the bottom with reasons.*

## Question
Among U.S. adults with a musculoskeletal (MSK) pain condition, is receiving physical therapy in one year associated
with different health care use, opioid use, pain interference and physical function in the following year,
compared with similar adults who did not receive physical therapy?

This is the evidence question behind digital MSK care: Hinge Health and peers argue that earlier, more accessible
PT-style care reduces downstream cost and opioid use. Public data can't evaluate a specific program, but it can
show what the observational evidence looks like for PT itself and how fragile it is.

## Data
AHRQ Medical Expenditure Panel Survey (MEPS), **Panel 26 Longitudinal File (HC-244, 2021-2022)**, linked by
DUPERSID to the 2021 Medical Conditions (HC-231), 2021 office-based (HC-229G) and outpatient (HC-229F) visit files,
and the 2021/2022 Prescribed Medicines files (HC-229A, HC-239A). Survey design: `LONGWT`, `VARSTR`, `VARPSU`.

## Cohort
- Panel 26 persons with `LONGWT > 0` (in scope both years), age 18-84 at end of 2021.
- At least one 2021 condition with ICD-10 (3-character, as released in MEPS) in **M15-M19** (osteoarthritis),
  **M25** (other joint disorders incl. joint pain), **M43-M54** (spine and back), **M60-M79** (soft tissue).
  Rheumatoid arthritis, gout and osteoporosis codes are not MSK *pain* for this purpose and don't qualify alone.

## Exposure
`PT = 1` if the person had **at least one 2021 office-based or outpatient visit with a physical/occupational
therapist** (`MEDPTYPE_M18 = 10`). MEPS public files don't reliably link these visits to the MSK condition, so this
is "PT use among people with MSK pain," not "PT for that condition." Secondary: dose (≥ 6 PT visits vs none).

## Outcomes (2022 unless noted)
Primary
1. Total health care expenditures (`TOTEXPY2`).
2. Any opioid prescription fill (Multum class 60 or 191 in any therapeutic-class field).

Secondary
3. Number of opioid fills. 4. Any ER visit (`ERTOTY2 > 0`). 5. Any inpatient stay (`IPDISY2 > 0`).
6. VR-12 physical component summary, round 4 (`VPCS4`). 7. Pain interference with normal work, round 4
(`ADPAIN4`, 1 = not at all … 5 = extremely).

**Negative-control outcome:** 2022 dental expenditures (`DVTEXPY2`). PT has no plausible effect on dental spending.
A clear "effect" there indicates residual confounding (for example, healthier or better-insured people being more
likely to use both PT and dental care), and every other estimate should be read with that in mind.

## Covariates (2021 / baseline)
Age, sex, race/ethnicity, poverty category, insurance coverage, region, years of education, employment (round 2),
2021 total expenditures **excluding PT spending**, 2021 opioid fill (any), number of distinct non-MSK CCSR condition
categories, specific MSK category flags (OA, spine, joint, soft tissue), diabetes, coronary heart disease, current
smoking (round 2), and round-2 VR-12 PCS, MCS, pain interference and self-rated general health.
Missing survey items (negative MEPS codes) are set to NaN, then median-imputed with a missing indicator.

## Methods
1. **Descriptives:** survey-weighted prevalence of MSK pain and PT use, with Taylor-linearized SEs (strata/PSU).
2. **Propensity score:** logistic regression of PT on the covariates.
3. **Primary estimator: overlap weights** × `LONGWT`. Overlap weights target the population with clinical
   equipoise and are bounded by construction (no extreme weights). The estimand is the ATO.
4. **Sensitivity estimators:** IPW for the ATE (PS trimmed to [0.05, 0.95]) × `LONGWT`, and augmented IPW
   (doubly robust, linear/logistic outcome models).
5. **Balance:** weighted standardized mean differences. Target |SMD| < 0.10 for every covariate.
6. **Inference:** 95% CIs from 500 bootstrap replicates, resampling PSUs within strata and re-estimating the PS
   in each replicate.
7. **Unmeasured confounding:** E-values for the primary outcomes.
8. **Decision rule for the write-up:** an association is described as "evidence of a difference" only if its 95%
   CI excludes 0 **and** the negative control does not show a comparable imbalance. Otherwise it is reported as
   inconclusive.

## Planned secondary deliverable (growth / experimentation)
Using the same data, estimate how much variance a CUPED-style baseline adjustment would remove from the main
outcomes (year-1 → year-2 correlation), and the sample size an RCT of an MSK program would need to detect a given
effect on spend, opioid use and PCS, with and without adjustment.

## Changes after outcomes were seen
*(none at the time of writing; any change will be added here with the reason)*
