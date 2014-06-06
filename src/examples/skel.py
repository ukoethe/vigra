from volumina.view3d.volumeRendering import RenderingManager
import vtk
import h5py
import numpy
from functools import partial
import time
from videoutils._videoutils import exportOBJ
from segmentation.h5utils import H5Path, rH5data, wH5data
from segmentation import mkdir_p
import vigra
import os
from chatty import Logger, Color
logger = Logger.getInstance()

from numpy2vtk import toVtkImageData

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

class MeshExporter(object):
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
            fname =  meshPrefix+"/%06d.obj" % label
            if not os.path.exists(fname):
                logger.log(Color.green, "will write %s" % fname)
                slicing = [slice(max(0, bboxMin[label,i]-5), min(seg.shape[i]-1, bboxMax[label,i]+5)) for i in range(3) ]
                assert numpy.all( [x.start < x.stop for x in slicing])
                assert numpy.all( [x.start >=0 for x in slicing])
                assert numpy.all( [x.stop < y for x,y in zip(slicing, seg.shape)])
                mesh = self._createMesh(label, slicing)
                x = [float(t.start) for t in slicing]
                exportOBJ(mesh, fname, x[0], x[1], x[2])
                logger.log(Color.green, "have written %s" % fname)

if __name__ == "__main__":
    #fname = "seg3d_small.h5"
    #fname = "cylinder/cylinder.h5"
    fname = "vectorialDist.h5"

    addBboxesToSegmentation(H5Path(fname+"/seg"))
    addBincountToSegmentation(H5Path(fname+"/seg"))

    #me = MeshExporter(H5Path("seg3d_small.h5/seg"), "meshes")
    #me.exportMeshes()
    
    #f = h5py.File("seg3d_small.h5")

    f = h5py.File("vectorialDist.h5")
    seg = f["seg"].value[:,:,:].astype(numpy.uint32)
    seg_bincount = f["bincount"].value
    print "seg.shape=",seg.shape
  
    '''
    argsort = numpy.argsort(seg_bincount)
    for i in range(len(argsort)-1, len(argsort)-1-10, -1):
        print argsort[i], seg_bincount[argsort[i]]
    print argsort[-1]
    print seg_bincount[argsort[:-5]]
    showLabels = [109209]
    '''
    
    showLabels = []

    ew = f["ew"].value
    ev = f["ev"].value
    magnitude = f["magnitude"].value
    vectorialDist = f["dist"].value
    f.close()
    print "magnitude.shape=",magnitude.shape
    print "ew.shape=",ew.shape
    print "ev.shape=",ev.shape
    
    #ev = 100*ew*ev

    x = ew[1,:,:,:]
    x = (x-x.min())/(x.max()-x.min())
    #x = numpy.logical_and(x < 4*1E-1, x > 1E-1)

    x = x > 4*1E-1

    #y = (seg*x).astype(numpy.uint8)

    y = x
    
    print "y.shape=",y.shape
    print "will generate #vectors = ", numpy.sum(y > 0)

    # With almost everything else ready, its time to initialize the
    # renderer and window, as well as creating a method for exiting
    # the application
    renderer = vtk.vtkRenderer()
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
    renderInteractor = vtk.vtkRenderWindowInteractor()
    renderInteractor.SetRenderWindow(renderWin)
    renderer.SetBackground(1, 1, 1) # white background
    renderWin.SetSize(400, 400)

    # A simple function to be called when the user decides to quit the
    # application.
    def exitCheck(obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)

    # Tell the application to use the function as an exit check.
    renderWin.AddObserver("AbortCheckEvent", exitCheck)

    if False:
        # create the rendering manager
        mgr = RenderingManager(renderer)
        mgr.setup(y.shape[::-1])

        for i in range(254):
            mgr.addObject()

        mgr.volume[:,:,:] = y
        mgr.update()
    else:
        w = numpy.where(y > 0)
        y_copy = y
        
        numVecs = len(w[0])
        
        numpy.random.seed(0)
        ctable = (255*numpy.random.random((256,3))).astype(numpy.uint8)
        
        #positions
        pts = vtk.vtkPoints()
        pts.SetNumberOfPoints(numVecs)
        #directions
        vecArr = vtk.vtkDoubleArray()
        vecArr.SetNumberOfComponents(3)
        vecArr.SetNumberOfTuples(numVecs)
        #colors
        colors = vtk.vtkUnsignedCharArray()
        colors.SetName("colors")
        colors.SetNumberOfComponents(3)
        
        for i in range(numVecs):
            x,y,z = w[0][i], w[1][i], w[2][i]
            label = seg[x,y,z]
            
            #if label not in showLabels:
            #    continue
            
            pts.InsertPoint(i, x,y,z)
            
            #print ew[:,x,y,z]
           
            #largest eigenvalue print ew[:,x,y,z]
            #vecArr.InsertTuple3(i, ev[0,0,x,y,z], ev[0,1,x,y,z], ev[0,2,x,y,z])
            #vecArr.InsertTuple3(i, ev[0,0,x,y,z], ev[1,0,x,y,z], ev[2,0,x,y,z])
            
            #vecArr.InsertTuple3(i, ev[1,0,x,y,z], ev[1,1,x,y,z], ev[1,2,x,y,z])
            #vecArr.InsertTuple3(i, ev[0,1,x,y,z], ev[1,1,x,y,z], ev[2,1,x,y,z])
            
            #smallest eigenvalue
            vecArr.InsertTuple3(i, ev[2,0,x,y,z], ev[2,1,x,y,z], ev[2,2,x,y,z])
            #vecArr.InsertTuple3(i, ev[0,2,x,y,z], ev[1,2,x,y,z], ev[2,2,x,y,z])
           
            #v = vectorialDist[x,y,z]
            #vecArr.InsertTuple3(i, v[0], v[1], v[2]) 
            
            colors.InsertNextTupleValue(ctable[label % 256,:])
        
        uGrid = vtk.vtkUnstructuredGrid()
        uGrid.SetPoints(pts)
        uGrid.GetPointData().SetScalars(colors)
        uGrid.GetPointData().SetVectors(vecArr)
        
        arrow = vtk.vtkArrowSource()
        #arrow = vtk.vtkCubeSource()
        #arrow = vtk.vtkSphereSource()
        #arrow = vtk.vtkLineSource()
        
        arrowPoly = arrow.GetOutput()
        
        glyph = vtk.vtkGlyph3D()
        glyph.SetInput(uGrid)
        glyph.SetColorModeToColorByScalar()
        glyph.SetSource(arrowPoly)
        glyph.SetScaleModeToScaleByVector()
        
        gMapper = vtk.vtkPolyDataMapper()
        gMapper.SetInput(glyph.GetOutput())
        gMapper.ScalarVisibilityOn()
        gMapper.SetScalarRange(uGrid.GetScalarRange())
        
        gactor = vtk.vtkActor()
        gactor.SetMapper(gMapper)
        gactor.GetProperty().SetColor(1.0, 0.0, 0.0)
        
        renderer.AddActor(gactor)
        
        # set the current point as the camera's focal point
        cam = renderer.GetActiveCamera()
        print y.shape
        pt = [ttt/2.0 for ttt in y_copy.shape]
        print pt
        cam.SetFocalPoint( *pt ) #rotate around this point
        cam.Modified();
     
    #
    # show object
    #
    
    for showLabel in showLabels:
        reader = vtk.vtkOBJReader()
        reader.SetFileName("meshes/%06d.obj" % showLabel)
        reader.Update()
    
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(reader.GetOutputPort())
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1.0,0.0,0.0)
        actor.GetProperty().SetOpacity(0.1)
        
        renderer.AddActor(actor)
    
    #
    #
    #

    renderInteractor.Initialize()

    # Because nothing will be rendered without any input, we order the
    # first render manually before control is handed over to the
    # main-loop.
    renderWin.Render()

    renderInteractor.Start()
