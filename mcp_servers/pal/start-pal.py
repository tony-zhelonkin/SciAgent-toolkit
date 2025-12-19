#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path

def main():
    # 1. Determine Project Root
    # We assume the current working directory is the project root when this script is called by MCP.
    project_root = Path.cwd()
    
    # 2. Load .env variables
    # We look for .env in the project root.
    env_files = [project_root / ".env", project_root / ".devcontainer" / ".env"]
    
    # Current environment
    env = os.environ.copy()
    
    for env_file in env_files:
        if env_file.exists():
            try:
                with open(env_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'): continue
                        if '=' in line:
                            k, v = line.split('=', 1)
                            # Remove surrounding quotes
                            v = v.strip()
                            if len(v) >= 2 and ((v[0] == '"' and v[-1] == '"') or (v[0] == "'" and v[-1] == "'")):
                                v = v[1:-1]
                            # Only set if not already in env (let shell override)? 
                            # actually for secrets, the .env usually holds the truth if shell is empty.
                            if k not in env:
                                env[k] = v
            except Exception as e:
                # Silently ignore errors reading .env, just proceed
                pass

    # 3. Path to the server executable in the venv
    # This script is located at .../mcp_servers/pal/start-pal.py
    # The venv is at .../mcp_servers/pal/venv
    script_dir = Path(__file__).parent.absolute()
    
    # We want to run the python interpreter from the venv, calling the module
    venv_python = script_dir / "venv" / "bin" / "python"
    
    if not venv_python.exists():
        sys.stderr.write(f"Error: Python interpreter not found at {venv_python}\n")
        sys.exit(1)

    # 4. Construct command
    # python -m pal_mcp_server [args...]
    cmd = [str(venv_python), "-m", "pal_mcp_server"] + sys.argv[1:]
    
    # 5. Execute
    try:
        # Replace current process
        os.execve(str(venv_python), cmd, env)
    except OSError as e:
        sys.stderr.write(f"Error executing PAL server: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
