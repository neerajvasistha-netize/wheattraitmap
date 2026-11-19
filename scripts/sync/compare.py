"""Compare the Claude artifact's file list with the website folder (web/).

Usage: python3 -I scripts/sync/compare.py artifact_files.json
artifact_files.json = the JSON returned by the claude-artifacts list_files tool.
Prints one line per file: SAME, CHANGED (artifact differs from repo), NEW (only in artifact)
or REPO-ONLY, plus the file size, so the sync step knows what to copy.
"""
import hashlib, json, os, sys
web = os.path.join(os.path.dirname(__file__), '..', '..', 'web')
art = {f['path']: f for f in json.load(open(sys.argv[1]))['files']}
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
repo = {n for n in os.listdir(web) if not n.startswith('.')}
for path, f in sorted(art.items()):
    p = os.path.join(web, path)
    state = 'NEW' if path not in repo else ('SAME' if sha(p) == f.get('sha256') else 'CHANGED')
    print(f"{state:9} {f.get('size', 0):>10} {path}")
for n in sorted(repo - set(art)):
    print(f"{'REPO-ONLY':9} {os.path.getsize(os.path.join(web, n)):>10} {n}")
