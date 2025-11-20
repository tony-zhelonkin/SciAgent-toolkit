# MCP Infrastructure Architecture

## System Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Interfaces                          │
├──────────────────────────┬──────────────────────────────────────┤
│                          │                                       │
│    Claude Code CLI       │         Codex CLI                     │
│    (Terminal/IDE)        │         (Terminal)                    │
│                          │                                       │
│    ~/.local/bin/claude   │    npm global or Homebrew             │
│                          │                                       │
└──────────────────────────┴──────────────────────────────────────┘
              │                              │
              │ MCP Protocol                 │ MCP Protocol
              │ (JSON-RPC over stdio)        │ (JSON-RPC over stdio)
              │                              │
              ├──────────────┬───────────────┤
              │              │               │
              ▼              ▼               ▼
┌─────────────────────────────────────────────────────────────────┐
│                       MCP Servers Layer                          │
├──────────────┬────────────────┬────────────────┬────────────────┤
│              │                │                │                │
│   Serena     │  Sequential    │  ToolUniverse  │    PubMed      │
│              │   Thinking     │                │                │
├──────────────┼────────────────┼────────────────┼────────────────┤
│ Command:     │ Command:       │ Command:       │ Installation:  │
│ uvx serena   │ npx sequential │ uv run tool-   │ Plugin system  │
│              │                │ universe-smcp  │                │
├──────────────┼────────────────┼────────────────┼────────────────┤
│ Functions:   │ Functions:     │ Functions:     │ Functions:     │
│ • Symbol     │ • Reasoning    │ • ChEMBL       │ • PubMed       │
│   search     │ • Problem      │ • UniProt      │   search       │
│ • Code edit  │   solving      │ • ClinicalT.   │ • PMC access   │
│ • Refactor   │ • Step-by-step │ • Europe PMC   │ • Metadata     │
│              │   analysis     │ • DrugBank     │ • Citations    │
│              │                │ • FDA data     │ • Related      │
│              │                │ • 600+ tools   │   articles     │
└──────────────┴────────────────┴────────────────┴────────────────┘
                                       │
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────┐
│                   External Data Sources                          │
├──────────────┬────────────────┬────────────────┬────────────────┤
│ NCBI/PubMed  │ Europe PMC     │ ChEMBL         │ UniProt        │
│ OpenAlex     │ Semantic       │ ClinicalTrials │ FDA Databases  │
│ CrossRef     │ Scholar        │ DrugBank       │ And more...    │
└──────────────┴────────────────┴────────────────┴────────────────┘
```

## Component Architecture

### Installation Scripts Hierarchy

```
setup_mcp_infrastructure.sh (MAIN ORCHESTRATOR)
├── install_claude.sh
│   └── Downloads and installs Claude Code
│       └── Configures PATH
│           └── Runs verification
│
├── install_codex.sh
│   └── Installs via npm or Homebrew
│       └── Provides auth instructions
│           └── Runs verification
│
├── mcp_servers/setup_serena.sh
│   └── Installs uv/uvx if needed
│       └── Tests Serena MCP server
│           └── Returns success
│
├── mcp_servers/setup_sequential_thinking.sh
│   └── Installs Node.js/npm if needed
│       └── Tests Sequential Thinking server
│           └── Returns success
│
├── mcp_servers/setup_tooluniverse.sh
│   └── Installs uv package manager
│       └── Creates Python environment
│           └── Installs ToolUniverse
│               └── Generates configs for Claude & Codex
│                   └── Creates test script
│
├── mcp_servers/setup_pubmed.sh
│   └── Provides installation instructions
│       └── Explains plugin marketplace
│           └── Lists features
│
└── config/merge_mcp_configs.sh
    └── Detects installed servers
        └── Creates unified .mcp.json
            └── Validates JSON syntax
```

## Data Flow

### 1. User Query Flow (Claude Code)

```
User Input
   │
   ▼
Claude Code CLI
   │
   ├─► Serena MCP ──────────► Code Analysis ─────► Results
   │                                                   │
   ├─► Sequential MCP ───────► Reasoning ─────────────┤
   │                                                   │
   ├─► ToolUniverse MCP ────► ChEMBL/UniProt ────────┤
   │                          API Calls               │
   │                                                   │
   └─► PubMed Plugin ───────► NCBI E-utils ──────────┤
                                                       │
                                                       ▼
                                           Synthesized Response
                                                       │
                                                       ▼
                                                   User Output
```

### 2. ToolUniverse Request Flow

```
Claude/Codex Request
   │
   ▼
ToolUniverse MCP Server
   │
   ├─► Tool Selection (from 600+ tools)
   │
   ├─► API Request Builder
   │
   ├─► External API Call
   │   ├─► ChEMBL
   │   ├─► UniProt
   │   ├─► Europe PMC
   │   ├─► ClinicalTrials.gov
   │   └─► Others...
   │
   ├─► Response Processing
   │
   ├─► Summarization (if Azure OpenAI configured)
   │
   └─► Return to Claude/Codex
```

### 3. PubMed Query Flow

```
User Query
   │
   ▼
PubMed Plugin (Claude Code)
   │
   ├─► Query Parsing
   │
   ├─► NCBI E-utilities API
   │   ├─► esearch (find PMIDs)
   │   ├─► esummary (get metadata)
   │   ├─► efetch (get full text from PMC)
   │   └─► elink (find related articles)
   │
   ├─► Response Formatting
   │
   └─► Return Results with:
       ├─► Abstracts
       ├─► Full text (if PMC available)
       ├─► Citations
       └─► Links
```

## Configuration Architecture

### Claude Code Configuration (.mcp.json)

```json
{
  "mcpServers": {
    "serena": {
      "type": "stdio",
      "command": "uvx",
      "args": ["--from", "git+...", "serena", "start-mcp-server", ...]
    },
    "sequential-thinking": {
      "type": "stdio",
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-sequential-thinking"]
    },
    "tooluniverse": {
      "type": "stdio",
      "command": "uv",
      "args": ["--directory", "./tooluniverse-env", "run", "tooluniverse-smcp-stdio"]
    }
  }
}
```

### Codex CLI Configuration (~/.codex/config.toml)

```toml
[mcp_servers.serena]
command = "uvx"
args = ["--from", "git+...", "serena", "start-mcp-server", ...]

[mcp_servers.sequential-thinking]
command = "npx"
args = ["-y", "@modelcontextprotocol/server-sequential-thinking"]

[mcp_servers.tooluniverse]
command = "uv"
args = ["--directory", "/path/to/tooluniverse-env", "run", "tooluniverse-smcp-stdio"]
startup_timeout_sec = 60
```

## Dependency Graph

```
System Dependencies
├── Python 3.10+ ──────────┬─► uv package manager ──► ToolUniverse
│                          │
│                          └─► uvx ──────────────────► Serena
│
├── Node.js 18+ ──────────► npm/npx ────────────────► Sequential Thinking
│                          │
│                          └─────────────────────────► Codex CLI (optional)
│
└── curl/bash ────────────► Native installer ───────► Claude Code
```

## Scientific Research Capabilities by Layer

### Layer 1: Code Intelligence (Serena)
```
Code Understanding
├── Symbol-level analysis
├── Cross-reference tracking
├── Intelligent refactoring
└── Context-aware editing
```

### Layer 2: Reasoning (Sequential Thinking)
```
Problem Solving
├── Step-by-step analysis
├── Decision trees
├── Complex reasoning chains
└── Multi-step technical decisions
```

### Layer 3: Scientific Tools (ToolUniverse)
```
Research Domains
├── Drug Discovery
│   ├── ChEMBL (molecular search)
│   ├── DrugBank (drug info)
│   └── FDA (approvals, safety)
│
├── Genomics & Proteomics
│   ├── UniProt (protein data)
│   ├── Gene databases
│   └── Pathway analysis
│
├── Literature Research
│   ├── Europe PMC
│   ├── Semantic Scholar
│   ├── OpenAlex
│   └── CrossRef
│
└── Clinical Research
    ├── ClinicalTrials.gov
    ├── FDA adverse events
    └── Drug interactions
```

### Layer 4: Biomedical Literature (PubMed)
```
Literature Access
├── PubMed (36M+ citations)
├── PMC (8M+ full-text)
├── NCBI databases
├── Citation networks
└── Related articles
```

## Integration Patterns

### Pattern 1: Literature-Driven Discovery

```
PubMed Search
   │
   ├─► Find relevant papers
   │
   ▼
Sequential Thinking
   │
   ├─► Analyze methodologies
   │
   ▼
ToolUniverse
   │
   ├─► Validate with clinical trials
   ├─► Check drug databases
   │
   ▼
Synthesis
   │
   └─► Evidence-based conclusions
```

### Pattern 2: Drug Development Pipeline

```
ChEMBL (ToolUniverse)
   │
   ├─► Find similar molecules
   │
   ▼
FDA Database (ToolUniverse)
   │
   ├─► Check safety profiles
   │
   ▼
ClinicalTrials.gov (ToolUniverse)
   │
   ├─► Find ongoing trials
   │
   ▼
PubMed
   │
   ├─► Literature evidence
   │
   ▼
Sequential Thinking
   │
   └─► Synthesize findings
```

### Pattern 3: Genomics Analysis

```
UniProt (ToolUniverse)
   │
   ├─► Protein information
   │
   ▼
Interaction Databases (ToolUniverse)
   │
   ├─► Protein-protein interactions
   │
   ▼
PubMed
   │
   ├─► Recent literature
   │
   ▼
Sequential Thinking
   │
   └─► Functional analysis
```

## Scalability Considerations

### Context Window Management

```
Problem: 600+ tools can overflow context
Solution: Tool filtering
   │
   ├─► --include-tools: Specific tools only
   ├─► --exclude-tool-types: Remove categories
   └─► Multiple instances: Specialized servers
```

### Performance Optimization

```
1. Lazy Loading
   └─► MCP servers start on first use

2. Parallel Queries
   └─► Multiple MCP servers run concurrently

3. Caching
   └─► ToolUniverse caches API responses

4. Summarization
   └─► Azure OpenAI condenses long outputs
```

## Security Model

```
Authentication Layers
├── Claude Code: Anthropic account
├── Codex CLI: ChatGPT or OpenAI API key
├── PubMed: No auth (public NCBI)
├── ToolUniverse: API-specific keys
└── Azure OpenAI: Optional (summarization)

Data Flow Security
├── Local execution: Scripts run on your machine
├── API calls: Direct to external services
├── No data storage: Transient processing
└── API keys: Environment variables only
```

## Extensibility Points

### Adding New MCP Servers

```
1. Create setup script in mcp_servers/
   └─► setup_newserver.sh

2. Add to orchestrator
   └─► setup_mcp_infrastructure.sh

3. Update config merger
   └─► config/merge_mcp_configs.sh

4. Document in README
   └─► config/README.md
```

### Customizing ToolUniverse

```
1. Tool filtering
   └─► --include-tools or --exclude-tool-types

2. Multiple instances
   └─► Different tool sets per instance

3. Custom hooks
   └─► Add processing pipelines

4. Environment-specific configs
   └─► Development vs. production settings
```

---

**This architecture provides a flexible, scalable foundation for AI-powered scientific research.**
