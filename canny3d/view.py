from enthought.mayavi import mlab
from enthought.mayavi.tools.pipeline import line_source
from enthought.tvtk.api import tvtk

import numpy as npy, canny, pdb

blue  = (0,0,1)
green = (0,1,0)
red   = (1,0,0)

def clf():
    mlab.clf()

 
def volume(vol):
    if type(vol) != npy.ndarray:
        vol = canny.sample(vol)
    m = canny.indices(vol)
    mlab.pipeline.volume(mlab.pipeline.scalar_field(m[0],m[1],m[2],vol))

    
def gradient(vol,i=None,scale=1):
    if type(vol)==npy.ndarray:
        if i == None:
            mlab.quiver3d(vol[...,0], vol[...,1], vol[...,2],scale_factor=scale)##, scale_factor=-0.2, color=(1, 1, 1))
        else:
            X = npy.array([vol[p,0] for p in i[0]])
            Y = npy.array([vol[p,1] for p in i[1]])
            Z = npy.array([vol[p,2] for p in i[2]])
            mlab.quiver3d(i[0], i[1], i[2],X,Y,Z, scale_factor=scale, color=(1, 1, 1))
    else:
        ii = canny.indices(vol)
        dx = vol.dx(ii)
        dy = vol.dy(ii)
        dz = vol.dz(ii)
        mlab.quiver3d(dx,dy,dz)

def points(points,color=(1,0,0),mode='point',scale_factor=0.1,**args):
    points = npy.array(points)
    assert 3 in points.shape and len(points.shape) == 2
    if points.shape[0] != 3:
        points = points.T
    mlab.points3d(points[0],points[1],points[2],#npy.linspace(0,1,len(points[0])),
                  scale_factor=scale_factor,mode='point',**args)

def trajectories(trajs):

    points = []
    lines  = []
    l      = 0
    scalars= []
    for tra in trajs:
        points += list(tra)
        scs = npy.arange(1,len(tra)+1)/(1.*len(tra))
        assert len(scs) == len(tra)
        scalars += list(scs)
        lines.append(range(l,l + len(tra)))
        l += len(tra)
        
    points  = npy.array(points)
    scalars = npy.array(scalars)
    
    mesh = tvtk.PolyData(points=points, lines=lines)
    mesh.point_data.scalars = scalars
    mesh.point_data.scalars.name = 'scalars'
    mlab.pipeline.surface(mesh)

def delaunay(points,opacity=0.25,alpha=None, **args):
    points = npy.array(points)
    assert 3 in points.shape and len(points.shape) == 2
    if points.shape[0] != 3:
        points = points.T
    
    s = npy.ones(points.shape[-1])*npy.inf
    src = mlab.pipeline.scalar_scatter(points[0],points[1],points[2])
    field = mlab.pipeline.delaunay3d(src,**args)
    pdb.set_trace()
    mlab.pipeline.surface(field,opacity=opacity)

def vectors(points,vol):
    points = npy.array(points)
    vectors = []
    for p in points:
        v = vol.d(tuple(p))
        v = v / v.norm()
        vectors += [p,p+v]
    vectors = npy.array(vectors).T
    lines   = npy.r_[:len(points)*2].reshape(len(points),2)
    scalars = npy.ones(2*len(points)) *npy.inf

    pdb.set_trace()
                 
    mesh = tvtk.PolyData(points=vectors.T, lines=lines)
    mesh.point_data.scalars = scalars
    mesh.point_data.scalars.name = 'scalars'
    mlab.pipeline.surface(mesh)
