Archived ATAC-seq analysis resources

This directory contains the recovered source tables and parameter evidence for
the bulk ATAC-seq motif analysis in Figure 2e,f. The resources were checked
against the retained analysis scripts and saved files. Alignment, peak calling,
differential testing and motif enrichment were not rerun during this check.

Contents
1. differential_peaks/diffpeak_control_VS_IR.csv: the saved Control-versus-IR
   annotated differential-peak table (27,898 rows). Its Fold sign is Control
   relative to IR: negative values indicate higher accessibility in IR.
2. motif_inputs/: the saved increased/decreased region tables and the complete
   background peak-region file. The two files ending in .bed retain all columns
   of the original annotated tables; their first three columns give chromosome,
   start and end. No coordinate conversion was applied during packaging.
3. homer/IR_up/ and homer/IR_down/: original knownResults.txt files and portable
   transcriptions of the corresponding saved motifFindingParameters.txt.
4. library_processing/: portable MACS2 log-header excerpts, per-library peak and
   duplication summaries, and Bowtie2 versions read from the eight BAM headers.
5. source_manifest.tsv: original-source SHA256 values, relative source names,
   packaged paths and hashes, and the transformation applied where applicable.
   Names in the source column are relative to the original xiayou archive.
   BAM text-header hashes are recorded separately in bam_software_summary.tsv;
   they do not represent hashes of the complete BAM files.
6. validation_summary.json: the completed checks and remaining design gap.

Confirmed source correspondence
The study has four pooled biological ATAC samples, with two separately prepared
technical libraries per condition. Original library aliases are N1/N2 (Ctrl),
IR1/IR2 (IR), D1/D2 (DMSO), and DIR1/DIR2 (DMSO-IR). Technical libraries are not
independent biological replicates. The biological design is reported by the
study authors; library identities are also documented by the saved commands.

All eight BAM headers record Bowtie2 2.3.3.1 and an mm10 reference. All eight
corresponding peak-table headers record MACS2 2.2.7.1, paired-end BAMPE input,
the mouse effective genome size, and a peak-calling q-value cutoff of 0.05.
The eight peak files in the export directory are byte-identical to their
corresponding processing-directory files. Picard MarkDuplicates records confirm
duplicate removal, but the exact Picard version is not established by these
records. Versions must not be inferred from a different processing branch.

Filtering the saved differential table at FDR <0.05 and Fold < -1 or Fold >1
recovers 2,285 IR-increased and 7,791 IR-decreased regions, respectively. Their
coordinate multisets match the saved HOMER inputs exactly. The common background
file contains 135,927 regions. Saved HOMER commands specify mm10, 200-bp regions,
10 threads, and this background file. The known-motif plots retain Benjamini
q <0.05 and deduplicate uppercase motif-name prefixes before parentheses.
Enrichment ratios are target motif percentage divided by background percentage.

The HOMER Log P-value column is on the natural-logarithm scale: for example, the
first increased-region p53 entry has P=1e-329 and Log P-value=-758, while the first
decreased-region GATA3 entry has P=1e-1971 and Log P-value=-4539. Figure 2e,f should
therefore label the unchanged plotted values as -ln(P), not -log10(P).

Statistical design still requiring verification
The original differential-testing sample sheet, peak-by-library count matrix,
normalization, design/contrast, model/backend settings, and software session
record have not been recovered. The differential table contains the saved
concentration, Fold, P-value and FDR fields, but these columns alone do not verify
the reported DiffBind version or the differential-model implementation.
The available BAM-merging commands belong to deepTools correlation/signal
summaries; they do not establish whether technical libraries were combined or
modeled separately for differential testing. The original analyst must supply
that information. The reported FDR values must not be interpreted as estimates
of reproducibility across independent animals or biological pools.

Portability and file integrity
Data tables and region files are copied without modification. Absolute local
paths were removed from command transcriptions and metadata excerpts; the
analytical options were preserved. Original full-source hashes and packaged
hashes are both recorded for those transformed text resources. Large BAMs and
full eight-library MACS2 peak tables are not included here. The source manifest
is a provenance record for this subdirectory; package-wide checksums are managed
by the enclosing reproducibility package.
