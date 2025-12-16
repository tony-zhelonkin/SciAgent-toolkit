# Data Processing

**Module:** `data_processing.md`
**Purpose:** Gene filtering, normalization, and differential expression workflow

---

## 1. Critical Workflow Order

**MUST follow this exact sequence:**

1. **Load raw counts** (featureCounts or similar)
2. **Annotate genes** (Ensembl → Symbol mapping) **BEFORE filtering**
3. **Filter lowly expressed genes** using `filterByExpr()` (NOT manual thresholds)
4. **TMM normalization** with `calcNormFactors()`
5. **Differential expression** with limma-voom
6. **GSEA** on gene symbols

### Common Mistakes to Avoid

- **Never** filter before gene annotation
- **Never** use manual filtering thresholds (e.g., `rowSums(counts >= 10) >= 3`)
- **Never** annotate genes after DE analysis
- **Never** use `db_species` parameter with msigdbr (version compatibility issues)

---

## 2. Gene Annotation

### 2.1 When: IMMEDIATELY After Loading, BEFORE Filtering

```r
# Load count matrix first
dge <- DGEList(counts = count_matrix, samples = sample_df)

# CRITICAL: Annotate BEFORE filtering
ensembl_ids <- gsub("\\..*$", "", rownames(dge))  # Remove version suffix

gene_symbols <- mapIds(org.Mm.eg.db,  # or org.Hs.eg.db for human
                      keys = ensembl_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Store annotation in $genes slot
dge$genes <- data.frame(
  ensembl_id = rownames(dge),
  ensembl_id_base = ensembl_ids,
  symbol = gene_symbols,
  stringsAsFactors = FALSE
)

# Update rownames to symbols (fallback to Ensembl for unmapped)
new_rownames <- ifelse(!is.na(gene_symbols), gene_symbols, rownames(dge))

# Handle duplicates by appending Ensembl ID
dup_symbols <- duplicated(new_rownames) | duplicated(new_rownames, fromLast = TRUE)
new_rownames[dup_symbols] <- paste0(new_rownames[dup_symbols], "_", ensembl_ids[dup_symbols])

rownames(dge) <- new_rownames

# NOW proceed with filtering
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
```

### 2.2 Why This Order Matters

1. **Statistical validity**: DE analysis on Ensembl IDs (unique, stable)
2. **Interpretability**: Gene symbols for downstream analysis
3. **GSEA compatibility**: MSigDB uses gene symbols
4. **Reproducibility**: Ensembl IDs stored in $genes for traceability

---

## 3. Filtering Strategy

### 3.1 Use filterByExpr() - NOT Manual Thresholds

```r
# CORRECT: Automatic, design-aware filtering
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# WRONG: Manual, arbitrary thresholds - DON'T DO THIS
MIN_COUNT <- 10
MIN_SAMPLES <- 3
keep <- rowSums(dge$counts >= MIN_COUNT) >= MIN_SAMPLES
```

### 3.2 Why filterByExpr() is Better

| Aspect | filterByExpr() | Manual Threshold |
|--------|---------------|------------------|
| Design-aware | Yes - considers group sizes | No |
| Library-size aware | Yes - adjusts for depth | No |
| Optimized | Based on empirical research | Arbitrary |
| Recommended | Official edgeR/limma docs | No |
| Validated | Reference implementations | Varies |

### 3.3 filterByExpr() Parameters

```r
# Default behavior (recommended)
keep <- filterByExpr(dge, design = design)

# Custom parameters (rarely needed)
keep <- filterByExpr(dge, design = design,
                     min.count = 10,      # Minimum count threshold
                     min.total.count = 15, # Minimum total across samples
                     large.n = 10,        # What counts as "large" group
                     min.prop = 0.7)      # Proportion requirement for large groups
```

---

## 4. Normalization

### 4.1 TMM Normalization (Standard)

```r
# After filtering
dge <- calcNormFactors(dge, method = "TMM")

# Verify normalization factors
summary(dge$samples$norm.factors)
# Should be close to 1.0 (typically 0.9-1.1)
```

### 4.2 When to Use Alternatives

| Method | When to Use |
|--------|-------------|
| TMM | Default - most RNA-seq experiments |
| RLE | High proportion of DE genes expected |
| upperquartile | Technical replicates |
| none | Already normalized data |

---

## 5. Differential Expression Workflow

### 5.1 Design Matrix Setup

```r
# No-intercept design (simpler contrasts)
design <- model.matrix(~ 0 + Group, data = dge$samples)
colnames(design) <- levels(dge$samples$Group)

# Verify design
head(design)
```

### 5.2 voom Transformation

```r
# Option 1: Standard voom (recommended for most cases)
v <- voom(dge, design, plot = TRUE)

# Option 2: With quality weights (heteroscedastic samples)
v <- voomWithQualityWeights(dge, design, plot = TRUE)

# Save voom plot for QC
pdf("03_results/plots/QC/voom_mean_variance.pdf")
v <- voom(dge, design, plot = TRUE)
dev.off()
```

### 5.3 Model Fitting

```r
# Fit linear model
fit <- lmFit(v, design)

# Define contrasts
contrast_matrix <- makeContrasts(
  Treatment_vs_Control = Treatment - Control,
  levels = design
)

# Apply contrasts
fit <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes moderation
fit <- eBayes(fit, robust = TRUE)  # robust = TRUE for outlier resistance
```

### 5.4 Extracting Results

```r
# Get DE results for a contrast
results <- topTable(fit, coef = "Treatment_vs_Control",
                    number = Inf,           # All genes
                    sort.by = "P",          # Sort by p-value
                    adjust.method = "BH")   # FDR correction

# Add gene symbols if not in rownames
results$gene <- rownames(results)

# Summary statistics
n_sig <- sum(results$adj.P.Val < 0.05)
n_up <- sum(results$adj.P.Val < 0.05 & results$logFC > 0)
n_down <- sum(results$adj.P.Val < 0.05 & results$logFC < 0)

message("DE genes: ", n_sig, " (", n_up, " up, ", n_down, " down)")
```

---

## 6. Quality Control Checks

### 6.1 Positive Control Gene

```r
# Check expected gene (e.g., knockout target)
target_gene <- "Il2ra"
target_idx <- grep(paste0("^", target_gene, "$"), rownames(results), ignore.case = TRUE)

if (length(target_idx) > 0) {
  message(target_gene, ": logFC = ", round(results$logFC[target_idx[1]], 3),
          ", FDR = ", signif(results$adj.P.Val[target_idx[1]], 3))
} else {
  warning(target_gene, " not found in results!")
}
```

### 6.2 P-value Distribution

```r
# Should show uniform distribution with spike at 0
hist(results$P.Value, breaks = 50,
     main = "P-value Distribution",
     xlab = "Raw P-value")
```

### 6.3 MA Plot

```r
# Check for systematic bias
plotMA(fit, coef = "Treatment_vs_Control")
abline(h = 0, col = "red", lty = 2)
```

---

## 7. Common Pitfalls & Solutions

| Problem | Cause | Solution |
|---------|-------|----------|
| Gene not found in results | Annotation after filtering | Annotate BEFORE filtering |
| Too few genes after filtering | Manual thresholds too strict | Use `filterByExpr()` |
| Duplicate gene symbols in GSEA | Multiple Ensembl → same symbol | Append Ensembl suffix to duplicates |
| msigdbr error: "unused arguments" | `db_species` parameter | Remove `db_species`, use `species` only |
| No DE genes found | Biological or technical issue | Check voom plot, sample clustering |

---

## 8. Reference Implementation

### 8.1 Complete Minimal Pipeline

```r
# Load and annotate
dge <- DGEList(counts = counts, samples = metadata)
dge <- annotate_genes(dge, species = "mouse")  # Custom function

# Filter and normalize
design <- model.matrix(~ 0 + Group, data = dge$samples)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# DE analysis
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, makeContrasts(Treatment - Control, levels = design))
fit <- eBayes(fit, robust = TRUE)

# Extract results
results <- topTable(fit, number = Inf)
```
