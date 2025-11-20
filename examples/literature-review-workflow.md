# Literature Review Workflow

## Overview

This workflow demonstrates how to conduct comprehensive literature reviews using the SciAgent Toolkit, including searching multiple databases, accessing full-text articles, extracting key findings, and identifying research gaps.

## Prerequisites

- PubMed plugin installed (for Claude Code)
- ToolUniverse MCP server installed
- Internet connection for API calls
- Basic understanding of scientific literature search

## Use Case

You need to conduct a comprehensive literature review on a specific topic (e.g., CRISPR gene editing for rare diseases), find relevant articles across multiple databases, access full-text when available, extract key findings, and identify research gaps for future work.

## Steps

### Step 1: Initial Broad Search in PubMed

**Query (Claude Code with PubMed):**
```
"Search PubMed for recent papers (last 5 years) about CRISPR gene editing for rare genetic diseases"
```

**Expected Output:**
- List of relevant articles with PMIDs
- Titles and authors
- Publication dates
- Abstracts
- Citation counts

**Notes:**
- Start broad to understand the field
- Note highly cited papers
- Identify key research groups
- Look for review articles

**Search refinement:**
```
"Narrow the search to CRISPR base editing for sickle cell disease"
"Find only clinical trials using CRISPR for rare diseases"
"Search for systematic reviews and meta-analyses on CRISPR safety"
```

### Step 2: Search Alternative Databases

**Query (Using ToolUniverse):**
```
"Search Europe PMC for CRISPR gene editing articles with full-text available"
```

**Expected Output:**
- Articles with PMC full-text access
- Open access papers
- Preprints
- Additional metadata

**Also search:**
```
"Search Semantic Scholar for CRISPR papers citing highly influential works"
"Find OpenAlex publications about CRISPR with high citation counts"
"Search CrossRef for recent CRISPR review articles"
```

**Why multiple databases:**
- Different coverage and indexing
- Access to different full-text sources
- Capture preprints and grey literature
- Find citation relationships

### Step 3: Access Full-Text Articles

**Query:**
```
"For the top 10 most relevant articles, get full-text from PMC where available"
```

**Expected Output:**
- Full-text articles (when available)
- Methods sections
- Results and figures
- Discussion and conclusions
- Supplementary materials references

**Notes:**
- Not all articles have PMC full-text
- Check for author manuscripts
- Note DOIs for institutional access

**Follow-up queries:**
```
"Extract the methods used for CRISPR delivery from these full-text articles"
"Summarize the key findings from these papers"
"What cell types were used in these studies?"
```

### Step 4: Extract Key Information

**Query:**
```
"From these articles, extract:
1. Main CRISPR techniques used
2. Target diseases studied
3. Delivery methods
4. Efficacy outcomes
5. Safety concerns reported
6. Clinical trial status if mentioned"
```

**Expected Output:**
- Structured summary of key information
- Comparison across studies
- Identification of common themes
- Noted contradictions or controversies

**Use Sequential Thinking:**
```
"Using sequential thinking, analyze the evolution of CRISPR techniques for rare diseases based on these papers chronologically"
```

### Step 5: Find Related Articles and Citation Networks

**Query:**
```
"For the 3 most important papers identified, find:
1. Articles that cite them
2. Articles they cite
3. Related articles by topic"
```

**Expected Output:**
- Citation network
- Key influencing papers
- Recent work building on findings
- Seminal papers in the field

**Follow-up:**
```
"Which papers are most frequently cited across this set?"
"Find the earliest paper describing CRISPR for therapeutic use"
"Identify review articles citing multiple papers from our set"
```

### Step 6: Identify Research Gaps

**Query:**
```
"Based on all the literature reviewed, identify:
1. Research gaps in CRISPR for rare diseases
2. Understudied diseases
3. Technical challenges mentioned repeatedly
4. Missing safety or efficacy data
5. Suggested future directions from discussion sections"
```

**Expected Output:**
- Summary of identified gaps
- Opportunities for new research
- Technical challenges to address
- Clinical translation barriers

### Step 7: Create Literature Summary

**Query:**
```
"Create a structured summary of this literature review including:
- Overview of current state of CRISPR for rare diseases
- Main techniques and approaches
- Clinical progress and trials
- Safety and efficacy data
- Research gaps and future directions
- Key papers to cite in each section"
```

**Expected Output:**
- Comprehensive structured summary
- Key papers organized by topic
- Evidence-based conclusions
- Clear identification of knowledge gaps

## Complete Example Session

Here's a complete literature review workflow:

```
1. "Search PubMed for papers about CAR-T cell therapy for B-cell lymphoma published in last 3 years"

2. "Narrow to Phase II and III clinical trials only"

3. "Search Europe PMC for full-text versions of the top 10 papers"

4. "Extract efficacy endpoints (complete response, progression-free survival) from these trials"

5. "Search for papers about adverse events and toxicities in CAR-T therapy"

6. "Find systematic reviews and meta-analyses combining these trial results"

7. "Identify which patient populations and disease subtypes have limited data"

8. "Create a structured summary with: trial designs, efficacy outcomes, safety profiles, and research gaps"
```

## Variations

### For Methodology Reviews

```
1. "Search for papers describing single-cell RNA-seq protocols"
2. "Find method comparison papers"
3. "Extract key protocol steps and variations"
4. "Identify most commonly used approaches"
5. "Note reported challenges and solutions"
```

### For Disease-Specific Reviews

```
1. "Search for Alzheimer's disease biomarker papers (last 5 years)"
2. "Find validation studies for each biomarker"
3. "Extract sensitivity and specificity data"
4. "Compare biomarker performance across studies"
5. "Identify gaps in biomarker validation"
```

### For Technique Evolution Reviews

```
1. "Search for CRISPR papers chronologically from 2012-2024"
2. "Identify major technical innovations each year"
3. "Track evolution from basic research to clinical trials"
4. "Note key breakthrough papers"
5. "Map the development of CRISPR variants (base editing, prime editing, etc.)"
```

### For Rapid Topic Assessment

```
1. "Find the 5 most recent review articles on [topic]"
2. "Get full-text for these reviews"
3. "Extract main conclusions and identified gaps"
4. "Identify the most cited primary research papers"
5. "Quick summary of current state and future directions"
```

## Advanced Techniques

### Boolean Search Strategies

```
"Search PubMed for: (CRISPR OR 'gene editing') AND ('rare disease' OR 'orphan disease') AND ('clinical trial' OR 'human study')"
```

### Date-Limited Searches

```
"Find papers about immunotherapy published in 2023-2024 only"
"Search for COVID-19 vaccine papers from early 2020"
```

### Specific Study Types

```
"Find only meta-analyses and systematic reviews about [topic]"
"Search for randomized controlled trials of [intervention]"
"Find case-control studies investigating [association]"
```

### Author and Institution Searches

```
"Find all papers by [researcher name] about CRISPR"
"Search for papers from [institution] on cancer immunotherapy"
```

### Citation Tracking

```
"Find all papers citing [landmark paper PMID]"
"Which recent papers cite the original CRISPR-Cas9 paper?"
```

## Troubleshooting

### Issue: Too many results (thousands of papers)

**Solution:** Add more specific filters
```
"Narrow to clinical trials only"
"Limit to last 2 years"
"Focus on specific disease subtypes"
"Search for systematic reviews only"
```

### Issue: Too few results

**Solution:** Broaden search terms
```
"Use more general terms (gene therapy instead of AAV gene therapy)"
"Expand date range"
"Include related terms and synonyms"
"Search multiple databases"
```

### Issue: Can't access full-text

**Solutions:**
- Check PMC for open access versions
- Look for author manuscripts
- Check institutional access
- Contact authors directly
- Use legal access methods only

### Issue: Inconsistent terminology

**Solution:** Use MeSH terms and synonyms
```
"Search using MeSH terms for standardized vocabulary"
"Include common synonyms: CRISPR, Cas9, gene editing"
"Check both US and UK spelling (tumor/tumour)"
```

### Issue: Information overload

**Solution:** Systematic filtering
```
"Start with review articles and meta-analyses"
"Focus on high-impact journals for initial assessment"
"Use citation count as one filter"
"Create hierarchy: reviews → key studies → supporting studies"
```

## Best Practices

1. **Start broad, then narrow**: Begin with general searches, refine based on results
2. **Use multiple databases**: Different databases have different coverage
3. **Document search strategy**: Record exact queries used
4. **Track search dates**: Note when searches were performed
5. **Save PMIDs/DOIs**: Keep permanent identifiers for reproducibility
6. **Check for retractions**: Verify papers haven't been retracted
7. **Read critically**: Don't assume all published findings are correct
8. **Note conflicts of interest**: Check author disclosures
9. **Compare methods**: Different methods may explain different results
10. **Update regularly**: Science evolves quickly

## Organizing Your Review

### Create a Literature Matrix

Query to help organize:
```
"Create a table comparing these 10 papers with columns: Authors, Year, Methods, Sample Size, Main Findings, Limitations"
```

### Topic Categorization

```
"Categorize these papers into themes: basic mechanisms, therapeutic applications, safety studies, clinical outcomes"
```

### Timeline Creation

```
"Organize these papers chronologically and identify how the field has evolved"
```

## Tools for Literature Management

While using SciAgent Toolkit:

1. **Export citations**: Save PMIDs and DOIs
2. **Use reference managers**: Import into Zotero, Mendeley, or EndNote
3. **Create summaries**: Save AI-generated summaries with citations
4. **Track sources**: Document which database each paper came from

## Integration with Writing

### For Introduction Sections

```
"Summarize the background on [topic] suitable for a paper introduction, with key citations"
```

### For Methods Comparison

```
"Compare the methods used across these studies and suggest which approach is most suitable for [your specific use case]"
```

### For Discussion

```
"How do my findings compare to the literature? Identify supporting and contradicting studies"
```

## Further Reading

- [PubMed Search Tips](https://pubmed.ncbi.nlm.nih.gov/help/)
- [Europe PMC](https://europepmc.org/)
- [Semantic Scholar](https://www.semanticscholar.org/)
- [PRISMA Guidelines](http://www.prisma-statement.org/) for systematic reviews

## Related Workflows

- [Drug Discovery Workflow](drug-discovery-workflow.md) - For therapeutic literature
- [Genomics Research Workflow](genomics-research-workflow.md) - For gene/protein literature

## Tips for Systematic Reviews

If conducting a formal systematic review:

1. **Pre-register protocol**: Use PROSPERO
2. **Multiple reviewers**: Plan for independent screening
3. **Document everything**: Complete search strings, dates, results
4. **Use PRISMA checklist**: Follow reporting guidelines
5. **Quality assessment**: Use appropriate tools (Cochrane, GRADE, etc.)
6. **Meta-analysis**: Consider statistical synthesis if appropriate

## Common Literature Review Pitfalls

1. **Confirmation bias**: Actively search for contradicting evidence
2. **Publication bias**: Note that negative results are under-published
3. **Recency bias**: Don't ignore older seminal papers
4. **Language bias**: Consider non-English publications if relevant
5. **Citation bias**: Highly cited doesn't always mean correct
6. **Selective reporting**: Check for outcome switching in trials

## Notes

Literature review is an iterative process. Don't expect to find everything in one search. Return to databases periodically, refine search terms based on what you learn, and update your review as new papers are published.

For systematic reviews intended for publication, follow field-specific guidelines (PRISMA, MOOSE, etc.) and consider pre-registration.
