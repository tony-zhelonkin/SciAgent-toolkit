# MCP Addon Templates

Addon templates define standalone MCP servers that persist across profile switches. Unlike profile-managed servers (which get overwritten on every `switch-mcp-profile.sh` call), addons are layered on top of whatever profile is active.

## How It Works

```
Profile template  -->  switch-mcp-profile.sh  -->  .mcp.json
                              |
                              +-- merge active addons
                              |
                    .claude/active-addons.json  +  templates/mcp-addons/*.addon.json
```

1. Profile switching writes the base `.mcp.json` as usual
2. After writing, the switcher checks `.claude/active-addons.json` for enabled addons
3. Enabled addons' `mcpServers` entries are merged into `.mcp.json`
4. Environment variables are substituted in the merged result

## Addon Template Format

Each addon is a `.addon.json` file with two top-level keys:

```json
{
  "_meta": {
    "name": "my-addon",
    "description": "What this addon does",
    "requires_setup": false,
    "setup_script": "setup_my_addon.sh",
    "context_tokens_estimate": "~5k",
    "runtime_dependency": "Description of what must be running",
    "env_vars": ["MY_VAR"]
  },
  "mcpServers": {
    "my-addon": {
      "type": "stdio",
      "command": "...",
      "args": ["..."],
      "env": {
        "MY_VAR": "${MY_VAR}"
      }
    }
  }
}
```

### `_meta` (addon metadata, stripped during merge)

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Addon identifier (must match filename without `.addon.json`) |
| `description` | Yes | Human-readable description |
| `requires_setup` | No | If `true`, `manage-addon.sh` warns when setup hasn't been run |
| `setup_script` | No | Script in `scripts/mcp_servers/` that installs dependencies |
| `context_tokens_estimate` | No | Approximate context window cost |
| `runtime_dependency` | No | What must be running for this addon to work |
| `env_vars` | No | Environment variables this addon needs |

### `mcpServers` (server definitions, merged into .mcp.json)

Same format as profile templates. Each key becomes a server entry in `.mcp.json`.

## Naming Conventions

- File: `<name>.addon.json`
- The `_meta.name` must match the filename stem
- Server names in `mcpServers` must not collide with profile server names
  - Safe: `jupyter`, `my-custom-tool`
  - Collision risk: `pal`, `context7`, `sequential-thinking` (used by profiles)
- If an addon defines a server name that exists in the active profile, the addon **overwrites** it

## Managing Addons

```bash
# List available and active addons
./scripts/manage-addon.sh list

# Enable an addon (merges into .mcp.json)
./scripts/manage-addon.sh enable jupyter --project-dir /path/to/project

# Disable an addon (removes from .mcp.json)
./scripts/manage-addon.sh disable jupyter --project-dir /path/to/project

# Check addon status and dependency health
./scripts/manage-addon.sh status jupyter --project-dir /path/to/project
```

## Creating a New Addon

1. Create `templates/mcp-addons/<name>.addon.json`
2. Fill in `_meta` and `mcpServers`
3. If setup is needed, create `scripts/mcp_servers/setup_<name>.sh`
4. Add any required env vars to `templates/.env.template`
5. Test: `./scripts/manage-addon.sh enable <name> --project-dir .`

## Available Addons

| Addon | Status | Description |
|-------|--------|-------------|
| `jupyter` | Placeholder | Jupyter MCP Server for notebook interaction |
