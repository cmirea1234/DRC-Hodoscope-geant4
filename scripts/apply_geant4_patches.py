"""Check/apply the three recorded patches to a private Geant4 source tree."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

BASE = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--source', type=Path, required=True)
p.add_argument('--apply', action='store_true', help='Without this flag, only inspect and dry-run')
a = p.parse_args()
source = a.source.resolve()
manifest = json.loads((BASE/'patches/manifest.json').read_text())
original = []
for name, item in manifest.items():
    digest = hashlib.sha256((source/name).read_bytes()).hexdigest()
    if digest == item['original_sha256']:
        original.append(name)
    elif digest != item['patched_sha256']:
        raise SystemExit('Unrecognized source hash; nothing patched: '+name)
# Check every patch before modifying any file.
for name in original:
    patch = BASE/'patches'/manifest[name]['patch']
    subprocess.run(['patch','-p1','--batch','--dry-run','-i',str(patch)],cwd=source,check=True)
if not a.apply:
    print('PATCH_PREFLIGHT_PASS',len(original),'files require patches; use --apply to apply')
else:
    for name in original:
        patch = BASE/'patches'/manifest[name]['patch']
        subprocess.run(['patch','-p1','--batch','-i',str(patch)],cwd=source,check=True)
    record = {}
    for name,item in manifest.items():
        assert hashlib.sha256((source/name).read_bytes()).hexdigest()==item['patched_sha256']
        record[name] = {key:item[key] for key in ('original_sha256','patched_sha256')}
    (source/'applied_patches.json').write_text(json.dumps(record,indent=2)+'\n')
    print('PATCHED_SOURCE_HASHES_PASS',len(record))
