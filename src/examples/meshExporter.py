#Python
import time
from functools import partial
import os

#SciPy
import h5py
import vigra
import numpy
import vtk

#chatty

class DumbWith(object):
   def __enter__(self):
       pass
   def __exit__(self ,type, value, traceback):
       pass

class DumbLogger(object):
    def __init__(self):
        pass
    def logTimed(self, *args, **kwargs):
        return DumbWith()
    def beginProgressLoop(self, *args, **kwargs):
        pass
    def endProgressLoop(self):
        pass
    def advanceProgressLoop(self):
        pass
    def log(*args, **kwargs):
        pass

class Color:
    green = None

logger = DumbLogger()

#this project
from numpy2vtk import toVtkImageData
from videoutils import exportNEURO
from video_h5utils import rH5data, wH5data, H5Path

def mkdir_p(path):
    if os.path.exists(path):
        return
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

#===---------------------------------------------------------------------------------------------------------------===

def addBboxesToSegmentation(segFile):
    print "add bounding boxes to segmentation %s" % segFile
    fname = segFile.filename
    if not H5Path(fname+"/bbox_min").exists():
        seg = rH5data(segFile, verbose=False).astype(numpy.uint32)
        dummyData = numpy.zeros(seg.shape, dtype=numpy.float32)
        fe = vigra.analysis.extractRegionFeatures(dummyData, seg, features=["Coord<Minimum >", "Coord<Maximum >"])
        f = h5py.File(fname, 'a')
        f.create_dataset("bbox_min", data = fe["Coord<Minimum >"])
        f.create_dataset("bbox_max", data = fe["Coord<Maximum >"])
        f.close()

def addBincountToSegmentation(segFile):
    print "add bincount to segmentation %s" % segFile
    fname = segFile.filename
    if not H5Path(fname+"/bincount").exists():
        seg = rH5data(segFile, verbose=False).astype(numpy.uint32)
        print "computing bincount..."
        bc = numpy.bincount(seg.flat)
        f = h5py.File(fname, 'a')
        f.create_dataset("bincount", data=bc)
        f.close()

class BlenderMeshExporter(object):
    def __init__(self, segFile, meshPrefix):
        assert isinstance(segFile, H5Path)
        self._seg = rH5data(segFile).astype(numpy.uint32)
        m = rH5data(segFile.filename, "bbox_min")
        M = rH5data(segFile.filename, "bbox_max")
        self._bbox = (m,M)
        self._bc = rH5data(segFile.filename, "bincount")
        self._meshPrefix = meshPrefix
    
    def _createMesh(self, label, slicing):
        with logger.logTimed("numpy.where(seg == %d)" % label):
            segBlock = self._seg[slicing]
            w = numpy.where(segBlock == label)
        with logger.logTimed("writing vol"):
            vol = numpy.zeros((segBlock.shape[2], segBlock.shape[1], segBlock.shape[0]), numpy.float32)
            vol[w[2], w[1], w[0]] = 1 #FIXME: coordinate transform

        def progressCallback(what, caller, eventId):
            print "%s: %f%%" % (what, 100.0*caller.GetProgress())
            
        img = toVtkImageData(vol)
        m = vtk.vtkMarchingCubes()
        m.SetInput(img)
        m.ComputeNormalsOn()
        m.SetValue(0, 1.0) 
        m.AddObserver(vtk.vtkCommand.ProgressEvent, partial(progressCallback, "marching cubes"))
        t1 = time.time()
        m.Update()
        t2 = time.time()
        print "marching cubes took %f sec." % (t2-t1) 

        scalarsOff    = vtk.vtkMaskFields()
        scalarsOff.SetInput(m.GetOutput())
        scalarsOff.CopyAttributeOff(vtk.vtkMaskFields.POINT_DATA,
                                    vtk.vtkDataSetAttributes.SCALARS)
        scalarsOff.CopyAttributeOff(vtk.vtkMaskFields.POINT_DATA,
                                    vtk.vtkDataSetAttributes.NORMALS)
        scalarsOff.CopyAttributeOff(vtk.vtkMaskFields.CELL_DATA,
                                    vtk.vtkDataSetAttributes.SCALARS)
        scalarsOff.Update()
       
        '''
        decimate = vtk.vtkQuadricDecimation()
        decimate.AddObserver(vtk.vtkCommand.ProgressEvent, partial(progressCallback, "decimate"))
        smooth   = vtk.vtkSmoothPolyDataFilter()
        smooth.AddObserver(vtk.vtkCommand.ProgressEvent, partial(progressCallback, "smooth"))
        smooth.SetInput(scalarsOff.GetOutput());
        smooth.SetNumberOfIterations(10)
        smooth.BoundarySmoothingOn();
        smooth.SetFeatureAngle(120)
        smooth.SetEdgeAngle(160);
        smooth.SetRelaxationFactor(0.25);
        decimate.SetInput(smooth.GetOutput());
        decimate.Update();
        smooth.RemoveAllInputs();
        decimate.RemoveAllInputs();
        mesh = vtk.vtkPolyData()
        mesh.DeepCopy(decimate.GetOutput())
        '''
        mesh = vtk.vtkPolyData()
        mesh.DeepCopy(scalarsOff.GetOutput())
        return mesh

    def exportMeshes(self):
        meshPrefix = self._meshPrefix 
        mkdir_p(meshPrefix)

        seg = self._seg 
        #with logger.logTimed("transpose seg"):
        #    #seg =numpy.transpose(seg, range(seg.ndim)[::-1])
        #    seg = seg.swapaxes(0,2)
        
        bboxMin = self._bbox[0]
        bboxMax = self._bbox[1]
     
        sel = numpy.where(self._bc > 0)[0]

        for i, label in enumerate(sel):
            if label == 0:
                continue 
            print "xxxxxxxxxxxxxxxxxxxxxxxxxxxx",i,"of", len(sel)
            fname =  meshPrefix+"/%06d.neuro" % label
            if not os.path.exists(fname):
                logger.log(Color.green, "will write %s" % fname)
                slicing = [slice(max(0, bboxMin[label,i]-5), min(seg.shape[i]-1, bboxMax[label,i]+5)) for i in range(3) ]
                assert numpy.all( [x.start < x.stop for x in slicing])
                assert numpy.all( [x.start >=0 for x in slicing])
                print "yyyyyyyyyyyyyyyy", slicing, seg
                assert numpy.all( [x.stop < y for x,y in zip(slicing, seg.shape)])
                mesh = self._createMesh(label, slicing)
                x = [float(t.start) for t in slicing]
                exportNEURO(mesh, fname, x[0], x[1], x[2], 1.0,0.0,0.0)
                logger.log(Color.green, "have written %s" % fname)

#===---------------------------------------------------------------------------------------------------------------===

if __name__ == "__main__":
    import socket
    if socket.gethostname() == "carsten":
        segFile = H5Path("seg.h5/seg")
    else:
        segFile = H5Path("seg.h5/volume/data")
    addBboxesToSegmentation(segFile)
    addBincountToSegmentation(segFile)

    me = BlenderMeshExporter(segFile, "/tmp/robert")
    me.exportMeshes()

