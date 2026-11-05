# CRISPRseek exported-caller trace (lead, 31 Oct 2026 ~16:20 PDT; static source reading)

Source: github.com/bioc/CRISPRseek (DESCRIPTION Version 1.53.1, devel), files fetched with gh api and kept here.
R/getOfftargetWithBulge.R last changed ef901b7b (2025-02-25); R/offTargetAnalysis.R last changed 373f6f53 (2025-03-10).

1. getOfftargetWithBulge is exported (NAMESPACE line 16). Signature (getOfftargetWithBulge.R:85-97):
   DNA_bulge = 2, RNA_bulge = 2, cas_offinder_version = c("2.4.1", "3.0.0b3"); line 100 match.arg ->
   default "2.4.1" (documented at line 32: "Defaults to \"2.4.1\""). Line 126 prepCasOFFinderBulge(version)
   -> for 2.4.1, getOfftargetWithBulgeUtils.R:154 loads the bundled inst/extdata/Cas-OFFinder/cas-offinder-bulge,
   whose RNA loops are identical to upstream v1.2 beta (A5 usage_comparison.json): emitted RNA deletion
   positions 1..n-2, omitting the first and the terminal spacer base.
   => The exported function's DEFAULT path runs the 2016 wrapper with its omission.
2. offTargetAnalysis is exported (NAMESPACE line 20). findOffTargetsWithBulge defaults to FALSE (line 576);
   when TRUE, method "CasOFFinder_v3.0.0b3" is mapped to "3.0.0b3" (lines 624-628) and passed to
   getOfftargetWithBulge (line 752) -> native Cas-OFFinder 3.0.0b3, which omits the terminal RNA position
   (paper_v7 tool audit; A2_casoffinder_v3_tags/FINDING.md).
   => Every bulge route in CRISPRseek inherits a positional omission; the default of the dedicated
   bulge function is the 2016 wrapper.
Evidence state: H2 (static trace of released source; not executed in R).
