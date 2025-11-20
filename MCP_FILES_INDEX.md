# MCP Infrastructure - Complete File Index

## Created Files Summary

**Total new files created: 16**
- 9 executable scripts (.sh)
- 5 documentation files (.md)
- 2 configuration files

---

## 📂 Project Root Files

### Installation Scripts (3 files)

| File | Purpose | Status |
|------|---------|--------|
| `setup_mcp_infrastructure.sh` | **MAIN ORCHESTRATOR** - Installs everything | ✅ Executable |
| `install_claude.sh` | Installs Claude Code | ✅ Executable |
| `install_codex.sh` | Installs Codex CLI | ✅ Executable |

### Documentation (5 files)

| File | Purpose | Read First? |
|------|---------|-------------|
| `QUICK_START.md` | 5-minute quick start guide | ⭐ START HERE |
| `MCP_SETUP_README.md` | Complete user guide | ⭐⭐ MAIN GUIDE |
| `INSTALLATION_SUMMARY.md` | What was installed and why | 📊 Reference |
| `ARCHITECTURE.md` | System architecture diagrams | 🏗️ Deep dive |
| `MCP_FILES_INDEX.md` | This file - complete file listing | 📋 Index |

---

## 📂 mcp_servers/ Directory (4 files)

All files are executable scripts:

| File | Purpose | MCP Server |
|------|---------|------------|
| `setup_serena.sh` | Installs Serena MCP | Code intelligence |
| `setup_sequential_thinking.sh` | Installs Sequential Thinking MCP | Reasoning |
| `setup_tooluniverse.sh` | Installs ToolUniverse MCP | 600+ scientific tools |
| `setup_pubmed.sh` | PubMed installation instructions | Biomedical literature |

---

## 📂 config/ Directory (2 files)

| File | Purpose | Type |
|------|---------|------|
| `merge_mcp_configs.sh` | Creates unified .mcp.json | ✅ Executable |
| `README.md` | Advanced configuration guide | 📖 Documentation |

---

## 📂 Auto-Generated Files

These files are created automatically when you run the setup:

| File | Created By | Purpose |
|------|------------|---------|
| `.mcp.json` | `merge_mcp_configs.sh` | Claude Code MCP configuration |
| `tooluniverse-env/` | `setup_tooluniverse.sh` | ToolUniverse Python environment |
| `test_tooluniverse.sh` | `setup_tooluniverse.sh` | Quick ToolUniverse test script |
| `~/.codex/config.toml` | `setup_tooluniverse.sh` | Codex CLI MCP configuration |

---

## Reading Order Recommendation

### For Quick Start (5 minutes)
1. `QUICK_START.md` - Get running immediately

### For Full Understanding (30 minutes)
1. `QUICK_START.md` - Overview
2. `MCP_SETUP_README.md` - Complete guide
3. `INSTALLATION_SUMMARY.md` - What you have
4. `config/README.md` - Advanced configuration

### For Deep Dive (1 hour)
1. All of the above, plus:
2. `ARCHITECTURE.md` - System architecture
3. Individual setup scripts in `mcp_servers/`

---

## File Relationships

```
QUICK_START.md
    │
    ├─► Points to: MCP_SETUP_README.md (detailed usage)
    │
    └─► Points to: setup_mcp_infrastructure.sh (installation)
         │
         ├─► Calls: install_claude.sh
         ├─► Calls: install_codex.sh
         ├─► Calls: mcp_servers/setup_serena.sh
         ├─► Calls: mcp_servers/setup_sequential_thinking.sh
         ├─► Calls: mcp_servers/setup_tooluniverse.sh
         ├─► Calls: mcp_servers/setup_pubmed.sh
         │
         └─► Calls: config/merge_mcp_configs.sh
              │
              └─► Creates: .mcp.json
```

---

## Usage by Persona

### I'm a Researcher - I Just Want to Use It
**Start here:**
1. `QUICK_START.md`
2. Run `./setup_mcp_infrastructure.sh`
3. Try the example queries in `QUICK_START.md`

### I'm a Power User - I Want to Customize
**Start here:**
1. `MCP_SETUP_README.md`
2. `config/README.md`
3. Customize `.mcp.json` or `~/.codex/config.toml`

### I'm a Developer - I Want to Understand the System
**Start here:**
1. `ARCHITECTURE.md`
2. `INSTALLATION_SUMMARY.md`
3. Read individual setup scripts in `mcp_servers/`

### I'm Troubleshooting
**Start here:**
1. `MCP_SETUP_README.md` → Troubleshooting section
2. `config/README.md` → Troubleshooting section
3. Run `claude doctor`

---

## File Size Reference

| File | Approximate Size | Complexity |
|------|-----------------|------------|
| `setup_mcp_infrastructure.sh` | ~12 KB | High |
| `install_claude.sh` | ~3 KB | Low |
| `install_codex.sh` | ~3 KB | Low |
| `setup_serena.sh` | ~3 KB | Low |
| `setup_sequential_thinking.sh` | ~3 KB | Low |
| `setup_tooluniverse.sh` | ~8 KB | Medium |
| `setup_pubmed.sh` | ~2 KB | Low |
| `merge_mcp_configs.sh` | ~3 KB | Low |
| `QUICK_START.md` | ~8 KB | Easy read |
| `MCP_SETUP_README.md` | ~18 KB | Comprehensive |
| `ARCHITECTURE.md` | ~15 KB | Technical |
| `config/README.md` | ~12 KB | Advanced |
| `INSTALLATION_SUMMARY.md` | ~10 KB | Reference |

**Total documentation: ~63 KB**
**Total scripts: ~37 KB**

---

## Key Features by File

### setup_mcp_infrastructure.sh
- ✅ Modular installation
- ✅ Command-line flags (--skip-claude, --skip-codex, etc.)
- ✅ Verification checks
- ✅ Colored output
- ✅ Idempotent (safe to run multiple times)

### install_claude.sh
- ✅ Native installer
- ✅ PATH configuration
- ✅ Diagnostics

### install_codex.sh
- ✅ npm or Homebrew
- ✅ Platform detection
- ✅ Authentication instructions

### setup_tooluniverse.sh
- ✅ uv installation
- ✅ Python environment creation
- ✅ Dual config (Claude + Codex)
- ✅ Azure OpenAI support
- ✅ Test script generation

### MCP_SETUP_README.md
- ✅ Quick start
- ✅ Installation steps
- ✅ Usage examples
- ✅ Scientific workflows
- ✅ Troubleshooting
- ✅ Uninstallation guide

### config/README.md
- ✅ Tool filtering
- ✅ Multiple instances
- ✅ Advanced configuration
- ✅ Performance optimization
- ✅ Security model

### ARCHITECTURE.md
- ✅ System diagrams
- ✅ Data flow
- ✅ Integration patterns
- ✅ Scalability
- ✅ Extensibility

---

## Scripts Dependency Tree

```
setup_mcp_infrastructure.sh (MAIN)
├── Requires: bash, curl, python3
├── Calls:
│   ├── install_claude.sh
│   │   └── Installs: Claude Code
│   │
│   ├── install_codex.sh
│   │   └── Installs: Codex CLI (via npm or brew)
│   │
│   ├── mcp_servers/setup_serena.sh
│   │   └── Installs: uv/uvx, Serena
│   │
│   ├── mcp_servers/setup_sequential_thinking.sh
│   │   └── Installs: Node.js/npm, Sequential Thinking
│   │
│   ├── mcp_servers/setup_tooluniverse.sh
│   │   └── Installs: uv, ToolUniverse
│   │       └── Creates: test_tooluniverse.sh
│   │
│   ├── mcp_servers/setup_pubmed.sh
│   │   └── Provides: Installation instructions
│   │
│   └── config/merge_mcp_configs.sh
│       └── Creates: .mcp.json
│
└── Validates: All installations
```

---

## Quick Reference Commands

### View Main Documentation
```bash
cat QUICK_START.md              # Quick start
cat MCP_SETUP_README.md         # Full guide
cat config/README.md            # Configuration
cat ARCHITECTURE.md             # Architecture
```

### View Script Help
```bash
./setup_mcp_infrastructure.sh --help
```

### Check What's Installed
```bash
cat INSTALLATION_SUMMARY.md
```

### Run Installation
```bash
./setup_mcp_infrastructure.sh
```

---

## Comparison with Original Setup

### Original: setup_mcp_claude.sh
- ❌ Single monolithic script (347 lines)
- ❌ Only Claude Code support
- ❌ No ToolUniverse
- ❌ No PubMed
- ❌ No Codex CLI
- ❌ Limited documentation

### New: Modular Infrastructure
- ✅ 9 modular scripts (easier to maintain)
- ✅ Claude Code + Codex CLI support
- ✅ ToolUniverse (600+ tools)
- ✅ PubMed integration
- ✅ Comprehensive documentation (5 guides)
- ✅ Advanced configuration options
- ✅ Scientific workflow examples

---

## Maintenance Guide

### To Update a Component
```bash
# Update individual component
./mcp_servers/setup_tooluniverse.sh

# Or re-run full setup
./setup_mcp_infrastructure.sh
```

### To Add a New MCP Server
1. Create `mcp_servers/setup_newserver.sh`
2. Add to `setup_mcp_infrastructure.sh`
3. Update `config/merge_mcp_configs.sh`
4. Document in `config/README.md`

### To Customize Configuration
1. Edit `.mcp.json` (for Claude Code)
2. Edit `~/.codex/config.toml` (for Codex CLI)
3. See `config/README.md` for examples

---

## Success Criteria

After running the setup, you should have:

✅ Claude Code installed and working
✅ (Optional) Codex CLI installed and working
✅ Serena MCP server configured
✅ Sequential Thinking MCP server configured
✅ ToolUniverse MCP server configured with 600+ tools
✅ Instructions for PubMed plugin installation
✅ `.mcp.json` created and valid
✅ All verification checks passed

Verify with:
```bash
claude --version
claude
/mcp  # Should show all servers
```

---

## Getting Help

| Issue | Solution |
|-------|----------|
| Installation fails | Check individual script logs |
| MCP servers not loading | Run `claude doctor` |
| Configuration errors | Validate with `python3 -m json.tool .mcp.json` |
| ToolUniverse issues | Run `./test_tooluniverse.sh` |
| General questions | Read `MCP_SETUP_README.md` |
| Advanced config | Read `config/README.md` |

---

## Summary

You now have a **complete, modular, well-documented MCP infrastructure** for building AI scientist agents with:

- 🎯 **2 AI interfaces**: Claude Code and Codex CLI
- 🔧 **4 MCP servers**: Serena, Sequential Thinking, ToolUniverse, PubMed
- 🔬 **600+ scientific tools**: Drug discovery, genomics, literature research
- 📚 **5 documentation guides**: Quick start to deep architecture
- ⚙️ **9 modular scripts**: Easy maintenance and customization
- ✅ **Production-ready**: Tested, validated, idempotent

**Start with: `QUICK_START.md`**

---

**Last updated:** 2025-11-19
**Version:** 1.0
**Total files created:** 16
