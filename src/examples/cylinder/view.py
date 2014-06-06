import h5py, numpy
import volumina
from PyQt4.QtGui import QApplication

app = QApplication([])
v = volumina.api.Viewer()
f = h5py.File("cylinder.h5")
d = f["seg"].value

print d.min(), d.max()

#d = (d-d.min())/float(d.max()-d.min())
d = (d == 0).astype(numpy.uint8)

f = h5py.File("20_cylinder.h5", 'w')
f.create_dataset("seg", data=d.astype(numpy.uint32))
f.close()

v.addGrayscaleLayer(d)
v.showMaximized()
app.exec_()
