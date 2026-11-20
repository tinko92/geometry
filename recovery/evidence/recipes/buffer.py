"""Select entire BG mapped features; freeze a deterministic local metric chart."""
import hashlib,json,math,pathlib,sys
D=pathlib.Path(sys.argv[2]); D.mkdir(exist_ok=False)
repo=pathlib.Path(sys.argv[1]).resolve()
(D/'fixtures').mkdir(exist_ok=True)
inputs=[('road','example/data/roads.wkt',[5,2,3,10,144,958]),
        ('coast','test/algorithms/buffer/data/gr_ll.wkt',[2,7,62,4])]
def sha(b):return hashlib.sha256(b).hexdigest()
out=[]
for kind,src,rows in inputs:
    data=(repo/src).read_bytes(); lines=data.decode().splitlines()
    for row in rows:
        original=lines[row-1]; geom='LINESTRING' if kind=='road' else 'POLYGON'
        prefix,suffix=('LINESTRING(',')') if kind=='road' else ('POLYGON((', '))')
        assert original.startswith(prefix) and original.endswith(suffix)
        pts=[tuple(map(float,p.split())) for p in original[len(prefix):-len(suffix)].split(',')]
        if kind=='coast':assert pts[0]==pts[-1]
        lon,lat=zip(*pts)
        cx,cy=(min(lon)+max(lon))/2,(min(lat)+max(lat))/2
        R=6371008.8;coslat=math.cos(math.radians(cy))
        xy=[(R*coslat*math.radians(x-cx),R*math.radians(y-cy)) for x,y in pts]
        wkt=prefix+','.join(f'{x:.17g} {y:.17g}' for x,y in xy)+suffix+'\n'
        name=f'{kind}-{row}'; (D/'fixtures'/(name+'.wkt')).write_text(wkt)
        out.append(dict(feature=name,kind=kind,source=src,line=row,source_sha256=sha(data),
            wkt_sha256=sha(wkt.encode()),vertices=len(pts),origin_lon=cx,origin_lat=cy,
            bbox_lonlat=[min(lon),min(lat),max(lon),max(lat)],
            x_scale_variation=max(abs(coslat/math.cos(math.radians(y))-1) for y in lat),
            length_or_perimeter_m=sum(math.dist(a,b) for a,b in zip(xy,xy[1:])),
            radii_m=[25] if kind=='road' else [250,-250],
            selection='axis-aligned mapped segment; explicit positive control' if name=='road-5' else 'whole mapped feature, size stratum'))
(D/'fixtures.json').write_text(json.dumps(out,indent=2)+'\n')
