# Archived treatment-classification code

These source files were copied from the saved analysis project, preserving their original relative directory layout and computations. `2-model codes` contains 22 configurations; 20 were displayed in current Fig. 6b. Inclusion, training/validation use and saved-result availability are documented in `../../resources/classifier_training_protocols.csv`. File numbering is not a performance ranking.

For new execution, prepare data with `../../scripts/prepare_classifier_inputs.py`, which uses the exact archived cell split and refuses missing genes. Run the chosen model script only after supplying the original expression input and installing its dependencies. Full matrices, trained classifiers and SHAP arrays are not published here.

Historical source titles sometimes call the fixed panel “ANOVA 100 genes.” The manuscript's panel is response-associated: it was derived using regression-based F scores against transcriptional-sensitivity targets. The distributed ordered `key_gene.csv` is the verified saved panel. Those historical plot labels do not specify a new differential-expression selection method.

The nine-gene branch here is the **final** `7-best_model_9genes` analysis, not the older candidate panel stored in a differently named local directory. Both the 100- and nine-gene classifiers used the common 2,301 test cells.

`4-SHAP/SHAP.py` preserves Kernel SHAP with a 50-center background and first 100 validation cells. It expects the historical list-of-arrays multiclass return format. Its global-importance summary does not define the unavailable intermediate 14-gene selection rule. Neural-network scripts preserve the recorded best-epoch handling; the protocol table documents where checkpoint preservation was not established. These scripts are analysis records, not a newly optimized benchmark.
