import math, os, numpy, PyQt4

import PyQt4.QtCore as qcore
import PyQt4.QtGui  as qt
from PyQt4.QtCore import SIGNAL

import vigra
import vigra.ufunc

try:
    from VigraQt import OverlayViewer, ImageCursor
except:
    vigra._fallbackModule('VigraQt', 
    '''    It can be obtained at 
    http://kogs-www.informatik.uni-hamburg.de/~meine/software/vigraqt/.''')
    from VigraQt import OverlayViewer, ImageCursor

import quickdialog
import weakref

class ImageViewer(OverlayViewer):
    
    activeViewers = weakref.WeakValueDictionary()

    def __init__(self, image, normalize = True, title = None, parent = None):
        OverlayViewer.__init__(self, parent)
        self.setImage(image, normalize)
        self._savedExpression = "x"
        self._lastSaveType = 2

        if title is not None:
            self.setWindowTitle(title)
        elif hasattr(image, "name"):
            self.setWindowTitle(image.name)
        else:
            for k in xrange(1,10000):
                if not ImageViewer.activeViewers.has_key(k):
                    break
            ImageViewer.activeViewers[k] = self
            self.setWindowTitle("Image %d" % k)

        self.imageCursor = ImageCursor(self)
        self.imageCursor.setVisible(False)
        self.imageCursor.setPosition(qcore.QPoint(self.image.width/2, self.image.height/2))
        self.addOverlay(self.imageCursor)

        self.zoomInAction = qt.QAction("Zoom in", self)
        self.zoomInAction.setShortcut("+")
        self.connect(self.zoomInAction, SIGNAL("triggered()"), self.zoomInPopup)

        self.zoomOutAction = qt.QAction("Zoom out", self)
        self.zoomOutAction.setShortcut("-")
        self.connect(self.zoomOutAction, SIGNAL("triggered()"), self.zoomOutPopup)

        self.saveAction = qt.QAction("Save image...", self)
        self.saveAction.setShortcut("S")
        self.connect(self.saveAction, SIGNAL("triggered()"), self.writeImage)

        self.expressionAction = qt.QAction("Apply expression...", self)
        self.expressionAction.setShortcut("E")
        self.connect(self.expressionAction, SIGNAL("triggered()"), self.applyExpression)

        self.cursorAction = qt.QAction("Line cursor", self)
        self.cursorAction.setShortcut("L")
        self.cursorAction.setCheckable(True)
        self.cursorAction.setChecked(False)
        self.connect(self.cursorAction, SIGNAL("triggered()"), self._toggleImageCursor)

        self.popup = qt.QMenu(self)
        self.popup.addAction(self.zoomInAction)
        self.popup.addAction(self.zoomOutAction)
        self.popup.addAction(self.saveAction)
        self.popup.addAction(self.expressionAction)
        self.popup.addAction(self.cursorAction)
    
    def setImage(self, image, normalize = True):
        self.image = image
        self._normalized = normalize
        OverlayViewer.setImage(self, image.qimage(normalize))
    
    def showImageCursor(self, yesOrNo = True):
        if yesOrNo != self.cursorAction.isChecked():
            self.cursorAction.trigger()

    def _toggleImageCursor(self):
        self.imageCursor.activateTool(self.cursorAction.isChecked())
        self.imageCursor.setVisible(self.cursorAction.isChecked())

    def _slideAfterZoom(self, shift):
        if self.zoomLevel() > 0:
            shift *= 1 + self.zoomLevel()
        elif self.zoomLevel() < 0:
            shift /= 1 - self.zoomLevel()
        self.slideBy(shift)

    def zoomInPopup(self):
        beforePos = self.imageCoordinate(self.mousepos)
        self.zoomUp()
        afterPos = self.imageCoordinate(self.mousepos)
        self._slideAfterZoom(afterPos - beforePos)

    def zoomOutPopup(self):
        beforePos = self.imageCoordinate(self.mousepos)
        self.zoomDown()
        afterPos = self.imageCoordinate(self.mousepos)
        self._slideAfterZoom(afterPos - beforePos)

#    def overlayPopup(self):
#        self.overlayMenu.clear()
#        index = 0
#        hideable = False
#        showable = False
#        for o in self.viewer.overlays:
#            if hasattr(o, "name"):
#                overlayName = o.name
#            else:
#                overlayName = self.viewer._defaultOverlayName(o)
#            text = "[%d] %s" % (index, overlayName)
#
#            color = None
#            if hasattr(o, "color") and isinstance(o.color, qt.QColor):
#                color = o.color
#                pmHeight = 5
#            elif hasattr(o, "fillColor") and isinstance(o.fillColor, qt.QColor):
#                color = o.fillColor
#                pmHeight = 16
#
#            if color:
#                colorPM = qt.QPixmap(16, pmHeight)
#                colorPM.fill(color)
#                iconSet = qt.QIconSet()
#                iconSet.setPixmap(colorPM, qt.QIconSet.Automatic,
#                                  qt.QIconSet.Normal, qt.QIconSet.On)
#                id = self.overlayMenu.insertItem(iconSet,
#                       text, self.viewer.toggleOverlayVisibility)
#            else:
#                id = self.overlayMenu.insertItem(
#                    text, self.viewer.toggleOverlayVisibility)
#
#            self.overlayMenu.setItemChecked(id, o.visible)
#            self.overlayMenu.setItemParameter(id, index)
#            if o.visible:
#                hideable = True
#            else:
#                showable = True
#            index += 1
#        id = self.overlayMenu.insertItem("&Hide all", self.viewer.showOverlays)
#        self.overlayMenu.setItemParameter(id, False)
#        self.overlayMenu.setItemEnabled(id, hideable)
#        id = self.overlayMenu.insertItem("&Show all", self.viewer.showOverlays)
#        self.overlayMenu.setItemParameter(id, True)
#        self.overlayMenu.setItemEnabled(id, showable)
        
    def applyExpression(self, expr = None, normalized = None):
        if expr is not None:
            self._savedExpression = expr
        else:
            d = quickdialog.QuickDialog(self,"Enter Expression")
            d.expression = quickdialog.OptionalStringInput(d, "Execute 'lambda x: ")
            d.expression.setText(self._savedExpression)
            d.expression.setFocus()
            d.addSpacing(10)
            d.norm = quickdialog.CheckBox(d, "Normalize intensity to range 0...255")
            d.norm.setChecked(self._normalized)
            if d.exec_() == 0:
                return
            self._savedExpression = d.expression.text()
            self._normalized = True if d.norm.selection() else False

        if normalized is not None:
            self._normalized = normalized

        try:
            image, normalized = self.getDisplayedImage()
        except Exception, e:
            qt.QMessageBox.critical(self, "Error Applying Expression", str(e))
            return

        OverlayViewer.setImage(self, image.qimage(normalized))

    def getDisplayedImage(self):
        """Returns the displayed image and the normalize flag
        (BYTE or NBYTE) as tuple/pair.

        Note that the returned image is the original image if no
        expression is applied, i.e. you should not change the returned
        object.  If active, the expression is applied via
        eval() on every call of getDisplayedImage()."""

        if not self._savedExpression or self._savedExpression == "x":
            self._savedExpression = "x"
            image = self.image
        else:
            for f in vigra.ufunc.__all__:
                exec 'from vigra.ufunc import %s' % f
            for f in dir(vigra.colors):
                if not f.startswith('__'):
                    exec 'from vigra.colors import %s' % f
            x = self.image
            image = eval(self._savedExpression)

        return image, self._normalized

    def writeImage(self):
        d = quickdialog.QuickDialog(self,"Write Image")

        imageFileExtensions = '*.' + ' *.'.join(vigra.impex.listExtensions().split(' '))
        d.filedialog = quickdialog.OutputFile(
            d, "Output filename:", "Image Files ("+imageFileExtensions+")")
        d.filedialog.setFocus()

        d.choices = quickdialog.HDialogGroup(d)

        d.type = quickdialog.VChoice(d.choices, "Output Pixel Type")
        d.type.addButton("Byte", "UINT8")
        d.type.addButton("Normalized to byte", "NBYTE")
        d.type.addButton("Keep type", "NATIVE")
        d.type.selectButton(1 if self._normalized else 0)
        d.type.buttonBox.setEnabled(self._lastSaveType)

        d.choices.addStretch(1)

        d.which = quickdialog.VChoice(d.choices, "Save ...")
        d.which.addButton("displayed image (zoomed, overlays)", 0)
        d.which.addButton("displayed image (1:1)", 1)
        d.which.addButton("original image", 2)
        d.connect(d.which.buttonBox, SIGNAL("clicked(int)"), \
                  d.type.buttonBox.setEnabled)
        d.which.selectButton(self._lastSaveType)

        while True:
            if d.exec_() == 0:
                return

            filename = d.filedialog.text()
            pixelType = d.type.selection()

            self._lastSaveType = d.which.selection()
            if d.which.selection():
                if d.which.selection() == 2:
                    image = self.image
                else:
                    image = self.getDisplay()[0]
                try:
                    image.write(filename, pixelType)
                except RuntimeError, e:
                    qt.QMessageBox.critical(self, "Error", str(e))
                else:
                    return
            else:
                formats = {"png": "PNG", \
                           "bmp": "BMP", \
                           "xbm": "XBM", \
                           "xpm": "XPM", \
                           "pnm": "PPM", \
                           "ppm": "PPM", \
                           "png": "PNG", \
                           "jpg": "JPEG", \
                           "jpeg": "JPEG", \
                           "tif": "TIF"}

                _, ext = os.path.splitext(filename)
                if not formats.has_key(ext[1:]):
                    f = " ".join(formats.keys())
                    qt.QMessageBox.critical(self, "Error", \
                                   "Displayed image with overlays can only be stored as\n"+f)
                else:
                    pixmap = self.getContentsPixmap()
                    pixmap.save(filename, formats[ext[1:]])
                    return

    def contextMenuEvent(self, e):
        "handles pop-up menu"
#        self.popup.setItemEnabled(
#            self.overlayMenuIndex, len(self.viewer.overlays) > 0)
        self.mousepos = e.pos()
        self.popup.exec_(e.globalPos())

    def keyPressEvent(self, e):
        "handles keys [S], [E], and possibly [Q] (for toplevel-windows)"
        if e.key() == qcore.Qt.Key_Q and not self.parent():
            self.close()
        elif e.key() == qcore.Qt.Key_S:
            self.writeImage()
        elif e.key() == qcore.Qt.Key_E:
            self.applyExpression()
        elif e.key() == qcore.Qt.Key_L:
            self.cursorAction.trigger()
        else:
            OverlayViewer.keyPressEvent(self, e)
            
class CaptionImageViewer(qt.QFrame):
    def __init__(self, image, normalize = True, title = None, parent = None):
        qt.QFrame.__init__(self, parent)
        self.viewer =ImageViewer(image, normalize, title, parent = self) 
        self.setWindowTitle(self.viewer.windowTitle())

        self._captionCoords = 0,0
        self._xplaces = int(math.log10(image.width) + 1.0)
        self._yplaces = int(math.log10(image.height) + 1.0)
        self._valueplaces = image.bands()*5
        
        self.label = qt.QLabel(self)
        font = qt.QFont()
        font.setPointSize(10)
        font.setStyleHint(qt.QFont.TypeWriter)
        self.label.setFont(font)

        self._layout = qt.QVBoxLayout(self)
        self._layout.setSpacing(5)
        self._layout.addWidget(self.viewer, 1)
        self._layout.addWidget(self.label)

        self.connect(self.viewer, SIGNAL('mouseOver(int, int)'), self.updateCaption)
        self.connect(self.viewer.cursorAction, SIGNAL('triggered()'), self._toggleCaptionSignals)

        self.updateCaption()

    def updateCaption(self, x=None, y=None):
        x = int(round(x)) if x is not None else self._captionCoords[0]
        y = int(round(y)) if y is not None else self._captionCoords[1]
        if x < 0 or x >= self.viewer.image.width or \
           y < 0 or y >= self.viewer.image.height:
            return
            
        self._captionCoords = x,y

        label = str(x).rjust(self._xplaces) + " x " + str(y).rjust(self._yplaces) +\
                " = " + str(self.viewer.image[x,y]).ljust(self._valueplaces)
        self.label.setText(label)
        self.emit(SIGNAL('captionChanged'), self.label.text())

    def updateCaptionP(self, point):
        self.updateCaption(point.x(), point.y())

    def _toggleCaptionSignals(self):
        if self.viewer.cursorAction.isChecked():
            self.disconnect(self.viewer, SIGNAL('mouseOver(int, int)'), self.updateCaption)
            self.connect(self.viewer.imageCursor, SIGNAL('positionChanged(QPoint)'), self.updateCaptionP)
        else:
            self.connect(self.viewer, SIGNAL('mouseOver(int, int)'), self.updateCaption)
            self.disconnect(self.viewer.imageCursor, SIGNAL('positionChanged(QPoint)'), self.updateCaptionP)

    def setImage(self, image, normalize = None):
        """imageWindow.setImage(image, normalize = None)

        Replace the current image with the given one.  If normalized
        is not given (or None), the normalized state is not changed."""
        
        self.viewer.setImage(image, normalize)
        self.updateCaption()

class CursorAction(qt.QAction):
    def __init__(self, name, parent):
        qt.QAction.__init__(self, name, parent)
        self.x, self.y = -1, -1
        self.zoomLevel = 0

    def trigger(self):
        qt.QAction.trigger(self)
        for v in self.viewers:
            v.viewer.cursorAction.setChecked(self.isChecked())
            v.viewer._toggleImageCursor()
            v._toggleCaptionSignals()

    def broadcastPosition(self, pos):
        if self.x ==  pos.x() and self.y == pos.y():
            return
        self.x, self.y = pos.x(), pos.y()
        for v in self.viewers:
            v.viewer.imageCursor.setPosition(pos)

    def broadcastZoom(self, level):
        if self.zoomLevel == level:
            return
        self.zoomLevel = level
        for v in self.viewers:
            v.viewer.setZoomLevel(level)
          
class ImageWindow(qt.QFrame):
    '''Display one or more images in a grid-like layout.
    '''
    def __init__(self, parent = None):
        qt.QFrame.__init__(self, parent)
        self.cursorAction = CursorAction("Connected line cursors", self)
        self.cursorAction.setCheckable(True)
        self.cursorAction.setChecked(False)
        self.addAction(self.cursorAction)
        self.cursorAction.viewers = []
        self.layout = qt.QGridLayout(self)
        
    def setImage(self, image, x=0, y=0, normalize=True, title=None):
        '''Place the given image at the given position of this window's grid layout.

           If an image already exists at this position, it is replaced.
        '''
        if self.layout.itemAtPosition(y, x):
            self.layout.itemAtPosition(y, x).widget().setImage(image, normalize)
        else:
            viewer = CaptionImageViewer(image, normalize, title, parent = self)
            self.layout.addWidget(viewer, y, x)
            self.cursorAction.viewers.append(viewer)
            if len(self.cursorAction.viewers) == 1:
                self.setWindowTitle(viewer.windowTitle())
            if self.cursorAction.x != -1:
                viewer.viewer.imageCursor.setPosition(qcore.QPoint(self.cursorAction.x, self.cursorAction.y))
            viewer.viewer.setZoomLevel(self.cursorAction.zoomLevel)
            if self.cursorAction.isChecked():
                viewer.viewer.cursorAction.trigger()
            self.disconnect(viewer.viewer.cursorAction, SIGNAL("triggered()"), viewer.viewer._toggleImageCursor)
            self.connect(viewer.viewer.cursorAction, SIGNAL("triggered()"), self.cursorAction.trigger)
            self.connect(viewer.viewer.imageCursor, SIGNAL("positionChanged(QPoint)"), self.cursorAction.broadcastPosition)
            self.connect(viewer.viewer, SIGNAL("zoomLevelChanged(int)"), self.cursorAction.broadcastZoom)
        self.updateGeometry()
        # this call is necessary to update the sizeHint() before adjustSize() is called
        qcore.QCoreApplication.processEvents()
        self.adjustSize()

def showImage(image, normalize = True, title = None):
    if isinstance(image, str):
        image = vigra.impex.readImage(image)
    v = ImageWindow()
    v.setImage(image, normalize=normalize, title=title)
    v.show()
    return v

