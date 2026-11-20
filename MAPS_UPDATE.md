# MAPS submodule: tracking latest maps

This repository contains a MAPS submodule that points to the MSM_maps repository. The submodule has been configured to track the upstream branch `master` so you can update it to the latest maps with Git or with the provided helper script.

Options to get the latest maps:

- Using Git (recommended for developers / CI):

  From the repository root:

  ```bash
  git submodule sync --recursive
  git submodule update --init --remote MAPS
  # If you want to persist the new submodule commit in the superproject:
  git add MAPS
  git commit -m "Update MAPS submodule to latest"
  ```

- Using the helper script (convenient for users who don't want to run git submodule commands):

  From the repository root:

  ```bash
  python3 scripts/update_maps.py
  ```

Notes:
- The submodule is configured with `branch = master`. If the upstream branch name changes, update `.gitmodules` accordingly.
- The helper script calls `git clone` / `git pull` and therefore requires `git` to be installed and network access.
