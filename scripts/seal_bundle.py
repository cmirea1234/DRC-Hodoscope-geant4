"""Explicitly seal an audited source bundle. Does not stage or push Git files."""
import hashlib
import json
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
excluded={'runs','local_check','.git','.venv','__pycache__'}
manifest={}
for p in sorted(BASE.rglob('*')):
    rel=p.relative_to(BASE)
    if not p.is_file() or p.name=='PAYLOAD_SHA256.json':continue
    if any(part in excluded or part.startswith('build') for part in rel.parts):continue
    if p.suffix in ['.pyc','.root','.log','.o','.so']:continue
    assert p.stat().st_size<5*1024*1024,('unexpected large Git payload',rel)
    manifest[str(rel)]=hashlib.sha256(p.read_bytes()).hexdigest()
(BASE/'PAYLOAD_SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('SEALED',len(manifest),'files')
