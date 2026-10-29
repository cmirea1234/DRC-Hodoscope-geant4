"""Module pose input: millimetres; active extrinsic rotations X then Y then Z."""
import math
ORDER=('T1','T2','H1','H2')
def validate_layout(layout):
    if layout.get('schema_version')!=2 or set(layout.get('modules',{}))!=set(ORDER):
        raise ValueError('Expected schema_version=2 and exactly T1,T2,H1,H2')
    def vector(x):
        return isinstance(x,list) and len(x)==3 and all(type(v) in (float,int) and math.isfinite(v) for v in x)
    half=layout.get('world_half_size_mm')
    if not vector(half) or any(v<=0 or v>10000 for v in half):raise ValueError('World half sizes must be positive and <=10000 mm')
    for name in ORDER:
        p=layout['modules'][name]
        if set(p)!= {'position_mm','rotation_deg'} or not all(vector(p[k]) for k in p):raise ValueError('Expected position_mm and rotation_deg vectors: '+name)
        if any(abs(a)>360 for a in p['rotation_deg']):raise ValueError('Rotation outside +/-360 degrees')
        if any(abs(a)>=b for a,b in zip(p['position_mm'],half)):raise ValueError('Module pivot outside world: '+name)
    return layout

def pose_text(layout):
    validate_layout(layout)
    lines=['KNU_LAYOUT_V1',' '.join(map(str,layout['world_half_size_mm']))]
    for i,name in enumerate(ORDER):
        p=layout['modules'][name];lines.append(' '.join(map(str,[i]+p['position_mm']+p['rotation_deg'])))
    return '\n'.join(lines)+'\n'

def rotation_matrix(angles):
    x,y,z=map(math.radians,angles);cx,sx=math.cos(x),math.sin(x);cy,sy=math.cos(y),math.sin(y);cz,sz=math.cos(z),math.sin(z)
    return [[cz*cy,cz*sy*sx-sz*cx,cz*sy*cx+sz*sx],[sz*cy,sz*sy*sx+cz*cx,sz*sy*cx-cz*sx],[-sy,cy*sx,cy*cx]]
