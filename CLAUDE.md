# WheatRLA repository: rules for Claude

The atlas lives in two places that must stay identical:
- GitHub: the `web/` folder of this repo, published by GitHub Pages at
  https://neerajvasistha-netize.github.io/WheatNHR-Atlas/ (deploys automatically on every push to `web/`).
- Claude: the private artifact "WheatRLA Database", slug `32c869ac-708e-4b31-90c1-b5907dd3ddc9`.

## Source of truth: this repo's `web/` folder
Every change to the atlas is made to the files in `web/` first, then sent to both places:
1. Edit or add files in `web/` (index.html, data JSON such as seed_*.json).
2. Commit and push to `main` (GitHub Pages redeploys by itself).
3. Republish the artifact from the same files with the claude-artifacts `publish_files` tool
   (slug above, `path` = absolute path of each changed file in `web/`, `as` = its name).
Never edit the artifact without committing the same files here.

## Daily check (artifact changed outside this workflow)
1. Call `list_files` on the artifact, save the JSON, run `python3 -I scripts/sync/compare.py <file>`.
2. For CHANGED/NEW files up to 4 MB: fetch them with `read_file` (`paths`, at most 4 MB per call),
   copy them into `web/`, commit "Sync from Claude artifact" and push.
3. A CHANGED file over 4 MB (index.html is about 8 MB) cannot be fetched this way, and live database
   edits are not files at all. Use the page's "Export for GitHub" button (bottom left, Claude only): it saves
   index.html and seed_live.json (all live database records) to the owner's Downloads; stage them from
   ~/Downloads with the device tools, copy into web/, commit and push.

Keep large data out of index.html where possible: put records in separate JSON files that the page
loads, so updates stay small and can be synced.

## Curation files
All atlas-making material (paper lists, screening decisions, download trackers, curation notes) is kept in this
repository under `curation/` (one subfolder per module, e.g. `curation/biotroph_papers/`). Update the files there
whenever a list or tracker changes; never commit PDFs of papers.
