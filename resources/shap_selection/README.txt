SHAP gene selection for Figure 6g

The author recalled selecting the red positive-contribution bars from the example-cell waterfall in Figure 6d. The archived values reproduce 14 such genes among the 19 individually displayed genes. The other 81 features are aggregated in that plot. Ranking uses the absolute SHAP value for the first explained validation cell and class 3 (DMSO-IR).

positive_contributors_14.csv gives the reconstructed selection, exact SHAP values and original absolute-value ranks. example_cell_SHAP_values_100.csv preserves all 100 values and display/selection indicators, allowing the display cutoff and sign rule to be checked. provenance.json records the source and scope.

These are local positive contributions to this cell's DMSO-IR prediction, not a global ranking over cells or evidence that increasing a gene's expression improves radioprotection. The final nine-gene panel is contained in this set; its overlap with HSC module 3 is Ptma and Uqcrb. No model was retrained, no SHAP values were recomputed, and no candidate identities were changed. The missing 84-gene measured-versus-virtual table was not reconstructed.
