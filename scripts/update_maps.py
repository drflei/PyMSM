#!/usr/bin/env python3
"""
scripts/update_maps.py

Helper to clone or update the MAPS subrepository (MSM_maps).
Run from the repository root (or via python -m scripts.update_maps) to ensure the MAPS directory is populated/updated.
"""
import os
import subprocess
import sys

REPO_URL = "https://github.com/drflei/MSM_maps.git"

# MAPS directory relative to repository root: repo_root/MAPS
here = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(here, os.pardir))
MAPS_DIR = os.path.join(REPO_ROOT, "MAPS")


def run(cmd, check=True):
    print('> ' + ' '.join(cmd))
    return subprocess.run(cmd, check=check)


def ensure_maps():
    maps_path = os.path.abspath(MAPS_DIR)
    if not os.path.exists(maps_path):
        print("Cloning MSM_maps into", maps_path)
        run(["git", "clone", REPO_URL, maps_path])
        # ensure the tracked branch exists locally
        try:
            run(["git", "-C", maps_path, "checkout", "master"])
        except subprocess.CalledProcessError:
            # continue if branch not present locally
            pass
    else:
        print("Pulling latest changes into", maps_path)
        # fetch and pull the configured remote/master
        run(["git", "-C", maps_path, "fetch", "--all"])
        run(["git", "-C", maps_path, "pull", "origin", "master"])


if __name__ == "__main__":
    try:
        ensure_maps()
    except subprocess.CalledProcessError as e:
        print("Error updating MAPS:", e, file=sys.stderr)
        sys.exit(2)
