# Drug Discovery Workflow

## Overview

This workflow demonstrates how to use the SciAgent Toolkit to evaluate potential drug candidates by searching molecular databases, checking safety profiles, finding clinical trials, and reviewing relevant literature.

## Prerequisites

- ToolUniverse MCP server installed
- PubMed plugin installed (for Claude Code)
- Internet connection for API calls
- Basic understanding of drug discovery process

## Use Case

You want to find compounds similar to a known drug (e.g., aspirin), evaluate their safety profiles, check if they're in clinical trials, and review the scientific literature to identify promising candidates for further research.

## Steps

### Step 1: Find Similar Molecules in ChEMBL

**Query:**
```
"Search ChEMBL for molecules similar to aspirin"
```

**Expected Output:**
- List of structurally similar molecules
- Similarity scores
- ChEMBL IDs
- Basic molecular properties

**Notes:**
- ChEMBL contains millions of bioactive molecules
- Similarity is based on molecular structure
- You can adjust similarity thresholds if needed

**Alternative queries:**
```
"Find molecules similar to ibuprofen"
"Search for NSAIDs in ChEMBL"
"Find aspirin analogs with improved properties"
```

### Step 2: Check FDA Approval Status

**Query:**
```
"For the top 5 similar molecules from ChEMBL, check their FDA approval status"
```

**Expected Output:**
- Approval status (approved/investigational/not approved)
- Indications
- Approval dates
- Marketing information

**Notes:**
- FDA-approved drugs have extensive safety data
- Investigational drugs may be in clinical trials
- Not all compounds will have FDA records

### Step 3: Find Clinical Trials

**Query:**
```
"Search ClinicalTrials.gov for ongoing trials involving these compounds for pain management"
```

**Expected Output:**
- Active clinical trials
- Trial phases (Phase I/II/III)
- Study objectives
- Enrollment status
- Location and contact information

**Notes:**
- Focus on Phase II/III trials for more advanced candidates
- Check trial status (recruiting, active, completed)
- Look for similar indications to your research interest

**Variations:**
```
"Find Phase III trials for [compound name]"
"Search for clinical trials of NSAIDs in cancer pain"
"Find completed trials with published results for [compound]"
```

### Step 4: Review Safety Information

**Query:**
```
"Get FDA safety data and adverse event reports for [top candidate compound]"
```

**Expected Output:**
- Known adverse events
- Safety alerts
- Black box warnings
- Drug interactions
- Contraindications

**Important:**
- Pay special attention to serious adverse events
- Look for patterns in adverse event reports
- Compare safety profiles across similar compounds

### Step 5: Search Scientific Literature

**Query (Claude Code with PubMed):**
```
"Search PubMed for recent studies (last 5 years) about COX-2 inhibitors and their safety profiles"
```

**Query (Using ToolUniverse):**
```
"Search Europe PMC for publications about [compound name] mechanism of action"
```

**Expected Output:**
- Recent research articles
- Meta-analyses and reviews
- Safety studies
- Mechanism of action papers
- Clinical outcomes data

**Notes:**
- Look for recent meta-analyses for comprehensive overviews
- Check for papers on mechanism of action
- Review safety and efficacy studies

### Step 6: Synthesize Findings

**Query:**
```
"Based on the molecular similarity, FDA data, clinical trials, and literature, which of these compounds appears most promising for pain management with the best safety profile?"
```

**Expected Output:**
- Comparative analysis of candidates
- Risk-benefit assessment
- Recommendations for further investigation
- Identified research gaps

**Use Sequential Thinking:**
The Sequential Thinking MCP server is particularly useful here for structured analysis of complex data.

## Complete Example Session

Here's a complete example workflow:

```
1. "Search ChEMBL for molecules structurally similar to aspirin with similarity >0.7"

2. "For the top 10 results, get their FDA approval status and approved indications"

3. "Which of these compounds are currently in Phase II or III clinical trials?"

4. "Get detailed adverse event data from FDA for the top 3 candidates"

5. "Search PubMed for safety and efficacy meta-analyses for these candidates"

6. "Compare the safety profiles of these three candidates and recommend the most promising one for further investigation in rheumatoid arthritis treatment"
```

## Variations

### For Anti-Cancer Drugs

```
1. "Find molecules similar to cisplatin in ChEMBL"
2. "Check for clinical trials in ovarian cancer"
3. "Search for resistance mechanisms in PubMed"
4. "Evaluate combination therapy trials"
```

### For Antibiotics

```
1. "Search ChEMBL for beta-lactam derivatives"
2. "Check FDA approval status and resistance patterns"
3. "Find clinical trials for resistant infections"
4. "Review literature on novel resistance mechanisms"
```

### For Rare Diseases

```
1. "Search DrugBank for compounds targeting [specific pathway]"
2. "Check orphan drug designations"
3. "Find any clinical trials (may be limited)"
4. "Search case reports and small studies in PubMed"
```

## Troubleshooting

### Issue: Too many results from ChEMBL

**Solution:** Add more specific filters
```
"Search ChEMBL for molecules similar to aspirin with similarity >0.8 and molecular weight <400"
```

### Issue: No FDA data found

**Reason:** Compound may not be approved or in FDA databases

**Solution:** Try other databases
```
"Search DrugBank for information about [compound name]"
"Check Europe PMC for regulatory information"
```

### Issue: API rate limiting

**Solution:**
- Wait a few minutes between requests
- Reduce number of compounds per query
- Focus on most promising candidates first

### Issue: Too much literature

**Solution:** Use more specific queries
```
"Search PubMed for meta-analyses about [compound] safety published in last 2 years"
```

## Best Practices

1. **Start broad, then narrow**: Begin with similarity searches, then focus on top candidates
2. **Cross-reference data**: Verify findings across multiple databases
3. **Check publication dates**: Prioritize recent safety and efficacy data
4. **Document your process**: Keep track of search queries and results
5. **Consider mechanism**: Understanding mechanism of action helps predict safety and efficacy

## Further Reading

- [ChEMBL Database](https://www.ebi.ac.uk/chembl/)
- [FDA Drug Database](https://www.fda.gov/drugs)
- [ClinicalTrials.gov](https://clinicaltrials.gov/)
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/)

## Related Workflows

- [Literature Review Workflow](literature-review-workflow.md) - For deeper literature analysis
- [Genomics Research Workflow](genomics-research-workflow.md) - For target identification

## Notes

This workflow is for research purposes only. Drug discovery and development requires extensive laboratory validation, clinical trials, and regulatory approval. This toolkit helps with initial screening and literature review but does not replace experimental validation.
