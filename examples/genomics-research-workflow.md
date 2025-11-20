# Genomics Research Workflow

## Overview

This workflow demonstrates how to analyze genes and proteins using the SciAgent Toolkit, including protein information retrieval, interaction analysis, pathway identification, and literature review for mutations and functional studies.

## Prerequisites

- ToolUniverse MCP server installed
- PubMed plugin installed (for Claude Code)
- Basic understanding of molecular biology and genomics
- Internet connection for API calls

## Use Case

You want to investigate a specific gene (e.g., BRCA1) to understand its protein product, identify protein interactions, determine involved pathways, review literature on mutations, and analyze functional implications for cancer research.

## Steps

### Step 1: Get Protein Information from UniProt

**Query:**
```
"Get comprehensive protein information for BRCA1 from UniProt"
```

**Expected Output:**
- Protein sequence
- Function description
- Subcellular location
- Post-translational modifications
- Known domains and motifs
- Gene ontology terms

**Notes:**
- UniProt provides curated protein data
- Swiss-Prot entries are manually reviewed (highest quality)
- TrEMBL entries are automated annotations

**Alternative queries:**
```
"Find all isoforms of TP53 in UniProt"
"Get protein information for EGFR including mutations"
"Search UniProt for proteins in the DNA repair pathway"
```

### Step 2: Identify Protein-Protein Interactions

**Query:**
```
"Find protein-protein interactions for BRCA1 using interaction databases"
```

**Expected Output:**
- Direct binding partners
- Interaction confidence scores
- Experimental evidence type
- Publications supporting interactions
- Functional context of interactions

**Notes:**
- Focus on high-confidence interactions
- Check experimental methods (Y2H, co-IP, etc.)
- Consider tissue/cell-type specificity

**Related queries:**
```
"What proteins does BRCA1 interact with in DNA repair?"
"Find interactions between BRCA1 and cell cycle proteins"
"Which kinases phosphorylate BRCA1?"
```

### Step 3: Identify Biological Pathways

**Query:**
```
"Identify all biological pathways involving BRCA1"
```

**Expected Output:**
- Pathway names and IDs
- Pathway diagrams
- Other proteins in the same pathways
- Cellular processes affected
- Disease associations

**Common pathways to explore:**
- DNA damage response
- Homologous recombination
- Cell cycle checkpoints
- Transcriptional regulation

**Variation queries:**
```
"Show me the DNA repair pathways involving BRCA1"
"What pathways are disrupted by BRCA1 mutations?"
"Find pathway crosstalk between BRCA1 and p53"
```

### Step 4: Search for Mutation Data

**Query:**
```
"Search for BRCA1 mutations in cancer databases and their clinical significance"
```

**Expected Output:**
- Known pathogenic mutations
- Variants of uncertain significance (VUS)
- Mutation frequencies in different cancers
- Clinical interpretations
- Functional impact predictions

**Important considerations:**
- Distinguish pathogenic from benign variants
- Check population frequency data
- Review functional studies of variants

**Related searches:**
```
"Find hotspot mutations in BRCA1 for breast cancer"
"What are the most common BRCA1 mutations in ovarian cancer?"
"Search for BRCA1 variants of uncertain significance with functional data"
```

### Step 5: Review Scientific Literature

**Query (PubMed):**
```
"Search PubMed for recent studies on BRCA1 mutations and their functional consequences in breast cancer (last 3 years)"
```

**Query (Europe PMC via ToolUniverse):**
```
"Find review articles about BRCA1 role in homologous recombination"
```

**Expected Output:**
- Recent research articles
- Review papers and meta-analyses
- Functional characterization studies
- Clinical outcome studies
- Therapeutic targeting papers

**Focus areas:**
- Mechanism of action studies
- Mutation functional studies
- Therapeutic implications
- Clinical outcomes

### Step 6: Integrate and Analyze Findings

**Query:**
```
"Based on the protein structure, interactions, pathways, mutations, and literature, what are the key functional roles of BRCA1 and how do pathogenic mutations disrupt these functions?"
```

**Expected Output:**
- Integrated functional summary
- Mutation impact analysis
- Therapeutic implications
- Research gaps identified

**Use Sequential Thinking:**
The Sequential Thinking MCP server helps with complex multi-step analysis:
```
"Using sequential thinking, analyze how BRCA1 mutations lead to cancer susceptibility based on all the data we've gathered"
```

## Complete Example Session

Here's a complete workflow for investigating a cancer-related gene:

```
1. "Get detailed protein information for TP53 from UniProt including all known domains"

2. "Find all proteins that directly interact with TP53"

3. "Identify the key pathways involving TP53 and its interacting partners"

4. "Search for TP53 mutations in colorectal cancer with their functional impacts"

5. "Find PubMed articles about TP53 mutations affecting DNA binding and their therapeutic implications"

6. "Synthesize this information: How do mutations in the DNA-binding domain affect TP53 function and what are the implications for therapy?"
```

## Variations

### For Studying Kinase Signaling

```
1. "Get information about EGFR kinase from UniProt"
2. "Find all substrates of EGFR"
3. "Identify EGFR signaling pathways"
4. "Search for activating mutations in EGFR"
5. "Find papers on EGFR inhibitor resistance mechanisms"
```

### For Metabolic Enzymes

```
1. "Get enzyme information for ALDH2 from UniProt"
2. "Find metabolic pathways involving ALDH2"
3. "Search for genetic variants affecting enzyme activity"
4. "Review literature on ALDH2 deficiency"
```

### For Transcription Factors

```
1. "Get MYC transcription factor information"
2. "Find MYC target genes"
3. "Identify MYC-regulated pathways"
4. "Search for MYC amplifications in cancer"
5. "Find papers on targeting MYC in therapy"
```

### For Novel Gene Investigation

```
1. "Search UniProt for [novel gene] and get all available annotations"
2. "Find evolutionarily conserved domains"
3. "Identify potential interacting proteins"
4. "Search for any disease associations"
5. "Review all published literature"
```

## Troubleshooting

### Issue: Gene name ambiguity

**Solution:** Use official gene symbols
```
"Search using the official HGNC symbol: TP53 (not p53)"
"Get UniProt entry for BRCA1 (human, UniProt ID: P38398)"
```

### Issue: Too many interactions

**Solution:** Filter by confidence and relevance
```
"Find high-confidence protein interactions for BRCA1 in DNA repair"
"Show only experimentally validated interactions for TP53"
```

### Issue: Limited mutation data

**Reason:** Gene may be less studied

**Solution:** Broaden search
```
"Search for any variants in [gene] across all databases"
"Find literature about [gene] function even without mutation data"
```

### Issue: Pathway information overwhelming

**Solution:** Focus on specific processes
```
"Show only DNA repair pathways involving BRCA1"
"Find cancer-related pathways for TP53"
```

## Best Practices

1. **Start with protein structure and function**: Understanding the protein helps interpret mutations
2. **Validate interactions**: Prioritize experimentally validated data
3. **Consider tissue context**: Some interactions and functions are tissue-specific
4. **Check mutation interpretation**: Use clinical databases (ClinVar) for variant significance
5. **Review recent literature**: New functional studies constantly emerge
6. **Cross-reference databases**: Confirm findings across multiple sources
7. **Document gene symbols**: Use official nomenclature (HGNC)

## Integration with Experimental Data

### Before Experiments

```
"Based on BRCA1 protein structure and interactions, which domains should I target for mutagenesis studies?"

"What are the most functionally important regions of TP53 based on mutation and structural data?"
```

### After Experiments

```
"I found that mutation X in BRCA1 disrupts binding to RAD51. Search literature for similar mutations and their effects"

"Compare my novel TP53 variant to known pathogenic mutations in the same domain"
```

## Further Reading

- [UniProt Database](https://www.uniprot.org/)
- [STRING (Protein Interactions)](https://string-db.org/)
- [KEGG Pathways](https://www.genome.jp/kegg/pathway.html)
- [ClinVar (Clinical Variants)](https://www.ncbi.nlm.nih.gov/clinvar/)
- [HGNC Gene Nomenclature](https://www.genenames.org/)

## Related Workflows

- [Drug Discovery Workflow](drug-discovery-workflow.md) - For therapeutic targeting
- [Literature Review Workflow](literature-review-workflow.md) - For comprehensive reviews

## Example Use Cases

### Cancer Research
Investigating oncogenes, tumor suppressors, and therapeutic targets

### Rare Disease Research
Understanding genetic variants and their functional impacts

### Drug Target Validation
Evaluating proteins as therapeutic targets

### Pathway Analysis
Mapping disease-associated pathway disruptions

## Notes

This workflow provides computational analysis to inform experimental research. Findings should be validated experimentally and interpreted in proper biological context. Clinical interpretations require consultation with medical genetics professionals.
