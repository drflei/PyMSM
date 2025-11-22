#!/bin/bash
# Sync the MAPS submodule to the latest commit on its configured branch.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

git submodule update --init --recursive --remote MAPS
