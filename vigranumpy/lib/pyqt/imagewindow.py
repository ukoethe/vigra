#######################################################################
#
#         Copyright 2009-2010 by Ullrich Koethe
#
#    This file is part of the VIGRA computer vision library.
#    The VIGRA Website is
#        http://hci.iwr.uni-heidelberg.de/vigra/
#    Please direct questions, bug reports, and contributions to
#        ullrich.koethe@iwr.uni-heidelberg.de    or
#        vigra@informatik.uni-hamburg.de
#
#    Permission is hereby granted, free of charge, to any person
#    obtaining a copy of this software and associated documentation
#    files (the "Software"), to deal in the Software without
#    restriction, including without limitation the rights to use,
#    copy, modify, merge, publish, distribute, sublicense, and/or
#    sell copies of the Software, and to permit persons to whom the
#    Software is furnished to do so, subject to the following
#    conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the
#    Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#    OTHER DEALINGS IN THE SOFTWARE.
#
#######################################################################

import math, os, numpy, PyQt4

import PyQt4.QtCore as qcore
import PyQt4.QtGui  as qt
from PyQt4.QtCore import SIGNAL

import vigra
import vigra.ufunc

try:
    from VigraQt import OverlayViewer, ImageCursor
except Exception, e:
    vigra._fallbackModule('VigraQt',
    '''
    %s

    If VigraQt is missing on your system, you can download it from
    http://kogs-www.informatik.uni-hamburg.de/~meine/software/vigraqt/.''' % str(e))
    from VigraQt import OverlayViewer, ImageCursor

import quickdialog
import weakref
import viewer2svg

class Crosshair(ImageCursor):
    def __init__(self, *args):
        ImageCursor.__init__(self, *args)
        self.visible = False
        self.position = qcore.QPoint(-1, -1)
    def setVisible(self, what=True):
        self.visible = what
        if what:
            ImageCursor.setPosition(self, self.position)
        else:
            ImageCursor.setPosition(self, qcore.QPoint(-1, -1))
    def setPosition(self, pos):
        self.position = pos
        if self.visible:
            ImageCursor.setPosition(self, self.position)

class ImageViewer(OverlayViewer):

    activeViewers = weakref.WeakValueDictionary()

    def __init__(self, image, normalize=True, title=None, parent=None):
        OverlayViewer.__init__(self, parent)
        self.setImage(image, normalize)
        self._savedExpression = "x"
        self._lastSaveType = 2
        self.overlays = []

        if title is not None:
            self.setWindowTitle(title)
        elif hasattr(image, "name"):
            self.setWindowTitle(image.name)
        else:
            for k in xrange(1, 10000):
                if not ImageViewer.activeViewers.has_key(k):
                    break
            ImageViewer.activeViewers[k] = self
            self.setWindowTitle("Image %d" % k)

        #self.imageCursor = ImageCursor(self) # doesn't work anymore - setVisible() is gone
        self.imageCursor = Crosshair(self)
        self.imageCursor.setVisible(False)
        self.imageCursor.setPosition(qcore.QPoint(self.image.width // 2, self.image.height // 2))
        OverlayViewer.addOverlay(self, self.imageCursor)

        self.zoomInAction = qt.QAction("Zoom in", self)
        self.zoomInAction.setShortcut("+")
        self.connect(self.zoomInAction, SIGNAL("triggered()"), self.zoomInPopup)

        self.zoomOutAction = qt.QAction("Zoom out", self)
        self.zoomOutAction.setShortcut("-")
        self.connect(self.zoomOutAction, SIGNAL("triggered()"), self.zoomOutPopup)

        self.saveAction = qt.QAction("Save image...", self)
        self.saveAction.setShortcut("S")
        self.connect(self.saveAction, SIGNAL("triggered()"), self.writeImage)

        self.svgAction = qt.QAction("Save as SVG...", self)
        self.svgAction.setShortcut("V")
        self.connect(self.svgAction, SIGNAL("triggered()"), self.writeSVG)

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
        self.popup.addAction(self.svgAction)
        self.popup.addAction(self.expressionAction)
        self.popup.addAction(self.cursorAction)

        self.overlayMenu = self.popup.addMenu("Overlays")
        self.connect(self.overlayMenu, SIGNAL("aboutToShow()"), self.overlayPopup)

    def setImage(self, image, normalize=True):
        if not hasattr(image, "qimage"):
            image = image.view(vigra.Image)
        self.image = image
        self._normalized = normalize
        OverlayViewer.setImage(self, image.qimage(normalize))

    def showImageCursor(self, yesOrNo=True):
        if yesOrNo != self.cursorAction.isChecked():
            self.cursorAction.trigger()

    def _toggleImageCursor(self):
        self.imageCursor.activateTool(self.cursorAction.isChecked())
        self.imageCursor.setVisible(self.cursorAction.isChecked())

    def addOverlay(self, overlay):
        if not hasattr(overlay, "draw"):
            raise TypeError("addOverlay: " + str(overlay) +
              "is no valid overlay with 'draw' method!")
        if overlay.parent() is None:
            overlay.setParent(self)
        overlay.visible = True
        if not hasattr(overlay, "name") or not overlay.name:
            overlay.name = self._defaultOverlayName(overlay)
        self.overlays.append(overlay)
        OverlayViewer.addOverlay(self, overlay)
        self.update()
        return len(self.overlays) - 1

    def removeOverlay(self, overlay):
        if type(overlay) == int:
            try:
                OverlayViewer.removeOverlay(self, self.overlays[overlay])
                self.overlays.pop(overlay)
                self.update()
            except IndexError, e:
                print "No such overlay."
        else:
            try:
                self.overlays.remove(overlay)
                OverlayViewer.removeOverlay(self, overlay)
                self.update()
            except ValueError, e:
                print "No such overlay."

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

    def _defaultOverlayName(self, o):
        name = str(o.__class__)
        if name[:8] == "<class '":
            name = name[8:-2]
        try:
            name = name[name.rindex(".") + 1:]
        except ValueError:
            pass
        return name

    def overlayPopup(self):
        self.overlayMenu.clear()
        index = 0
        hideable = False
        showable = False
        for o in self.overlays:
            overlayName = o.name
            text = "[%d] %s" % (index, overlayName)

            color = None
            if hasattr(o, "color") and isinstance(o.color, qt.QColor):
                color = o.color
                pmHeight = 5
            elif hasattr(o, "fillColor") and isinstance(o.fillColor, qt.QColor):
                color = o.fillColor
                pmHeight = 16

            if color:
                colorPM = qt.QPixmap(16, pmHeight)
                colorPM.fill(color)
                icon = qt.QIcon(colorPM)
                id = qt.QAction(icon, text, self)
            else:
                id = qt.QAction(text, self)

            self.overlayMenu.addAction(id)
            id.setCheckable(True)
            self.connect(id, SIGNAL('triggered()'), self.toggleOverlayVisibilityWithParam(o))
            id.setChecked(o.isVisible())
            if o.isVisible():
                hideable = True
            else:
                showable = True
            index += 1
        id = qt.QAction("&Hide all", self)
        self.overlayMenu.addAction(id)
        self.connect(id, SIGNAL('triggered()'), self.toggleOverlayVisibilityWithParam(False))
        id.setEnabled(hideable)
        id = qt.QAction("&Show all", self)
        self.overlayMenu.addAction(id)
        self.connect(id, SIGNAL('triggered()'), self.toggleOverlayVisibilityWithParam(True))
        id.setEnabled(showable)

    def toggleOverlayVisibilityWithParam(self, o):
        return lambda: self.toggleOverlayVisibility(o)

    def toggleOverlayVisibility(self, o=None):
        '''Toggle or set visibility of given overlay and update view.
           The parameter can be a boolean - which sets the visibility of
           all overlays accordingly - an overlay object or the index
           of the overlay to be hidden/re-shown. If it is omitted, all
           overlays will be toggled.
        '''
        if o is None:
            for k in self.overlays:
                k.setVisible(not k.isVisible())
        elif type(o) is bool:
            for k in self.overlays:
                k.setVisible(o)
        else:
            if type(o) is int:
                o = self.overlays[o]
            o.setVisible(not o.isVisible())
        self.update()

    def applyExpression(self, expr=None, normalized=None):
        if expr is not None:
            self._savedExpression = expr
        else:
            d = quickdialog.QuickDialog(self, "Enter Expression")
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
        d = quickdialog.QuickDialog(self, "Write Image")

        imageFileExtensions = '*.' + ' *.'.join(vigra.impex.listExtensions().split(' '))
        d.filedialog = quickdialog.OutputFile(
            d, "Output filename:", "Image Files (" + imageFileExtensions + ")")
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
                    image.writeImage(filename, pixelType)
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
                                   "Displayed image with overlays can only be stored as\n" + f)
                else:
                    pixmap = self.getContentsPixmap()
                    pixmap.save(filename, formats[ext[1:]])
                    return

    def writeSVG(self):
        d = quickdialog.QuickDialog(self, "Write Viewer Contents to SVG")

        d.filedialog = quickdialog.OutputFile(
            d, "Output filename:", "SVG Files (*.svg)")
        d.filedialog.setFocus()

        d.choices = quickdialog.HDialogGroup(d)

        d.which = quickdialog.VChoice(d.choices, "Save ...")
        d.which.addButton("all overlays", 0)
        d.which.addButton("only displayed overlays", 1)

        d.which.selectButton(self._lastSaveType)

        while True:
            if d.exec_() == 0:
                return

            self._lastSaveType = d.which.selection()
            allOVs = (d.which.selection() == 0)

            filename = d.filedialog.text()
            basename, ext = os.path.splitext(filename)

            try:
                if ext == ".SVG" or ext == ".svg":
                    viewer2svg.viewer2svg(self, basename, not allOVs)
                else:
                    viewer2svg.viewer2svg(self, filename, not allOVs)
            except RuntimeError, e:
                qt.QMessageBox.critical(self, "Error", str(e))
            return

    def contextMenuEvent(self, e):
        "handles pop-up menu"
        self.overlayMenu.setEnabled(len(self.overlays) > 0)
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
        elif e.key() == qcore.Qt.Key_Right or e.key() == qcore.Qt.Key_Left or \
          e.key() == qcore.Qt.Key_Up or e.key() == qcore.Qt.Key_Down:
            OverlayViewer.keyPressEvent(self, e)
        elif e.key() == qcore.Qt.Key_Plus or e.key() == qcore.Qt.Key_Greater:
            OverlayViewer.zoomUp(self)
        elif e.key() == qcore.Qt.Key_Minus or e.key() == qcore.Qt.Key_Less:
            OverlayViewer.zoomDown(self)
        else:
            self.emit(qcore.SIGNAL("keyPressed"), (e.key()))
            e.ignore()

    def keyReleaseEvent(self, e):
        self.emit(qcore.SIGNAL("keyReleased"), (e.key()))
        e.ignore()

    def mousePressEvent(self, e):
        imagePos = OverlayViewer.imageCoordinateF(self, qcore.QPoint(e.x(), e.y()))
        self.emit(qcore.SIGNAL("mousePressed"), (imagePos.x(), imagePos.y(), e.button()))
        OverlayViewer.mousePressEvent(self, e)
        e.ignore()

class CaptionImageViewer(qt.QFrame):
    def __init__(self, image, normalize=True, title=None, parent=None):
        qt.QFrame.__init__(self, parent)
        self.viewer = ImageViewer(image, normalize, title, parent=self)
        self.setWindowTitle(self.viewer.windowTitle())

        self._captionCoords = 0, 0
        self._xplaces = int(math.log10(self.viewer.image.width) + 1.0)
        self._yplaces = int(math.log10(self.viewer.image.height) + 1.0)
        self._valueplaces = self.viewer.image.channels * 5

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

        self._captionCoords = x, y

        label = str(x).rjust(self._xplaces) + " x " + str(y).rjust(self._yplaces) +\
                " = " + str(self.viewer.image[x, y]).ljust(self._valueplaces)
        self.label.setText(label)
        self.emit(SIGNAL('captionChanged'), self.label.text())

    def updateCaptionP(self, point):
        self.updateCaption(point.x(), point.y())

    def _toggleCaptionSignals(self):
        if self.viewer.cursorAction.isChecked():
            self.disconnect(self.viewer,
              SIGNAL('mouseOver(int, int)'), self.updateCaption)
            self.connect(self.viewer.imageCursor,
              SIGNAL('positionChanged(QPoint)'), self.updateCaptionP)
        else:
            self.connect(self.viewer,
              SIGNAL('mouseOver(int, int)'), self.updateCaption)
            self.disconnect(self.viewer.imageCursor,
              SIGNAL('positionChanged(QPoint)'), self.updateCaptionP)

    def setImage(self, image, normalize=None):
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
        if self.x == pos.x() and self.y == pos.y():
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
    def __init__(self, parent=None):
        qt.QFrame.__init__(self, parent)
        self.cursorAction = CursorAction("Connected line cursors", self)
        self.cursorAction.setCheckable(True)
        self.cursorAction.setChecked(False)
        self.addAction(self.cursorAction)
        self.cursorAction.viewers = []
        self.layout = qt.QGridLayout(self)

    def setImage(self, image, x=0, y=0, normalize=True, title=None):
        """Place the given image at the given position of this window's grid layout.

           If an image already exists at this position, it is replaced.
        """
        if self.layout.itemAtPosition(y, x):
            self.layout.itemAtPosition(y, x).widget().setImage(image, normalize)
        else:
            CIviewer = CaptionImageViewer(image, normalize, title, parent=self)
            self.layout.addWidget(CIviewer, y, x)
            self.cursorAction.viewers.append(CIviewer)
            if len(self.cursorAction.viewers) == 1:
                self.setWindowTitle(CIviewer.windowTitle())
            if self.cursorAction.x != -1:
                CIviewer.viewer.imageCursor.setPosition(
                  qcore.QPoint(self.cursorAction.x, self.cursorAction.y))
            CIviewer.viewer.setZoomLevel(self.cursorAction.zoomLevel)
            if self.cursorAction.isChecked():
                CIviewer.viewer.cursorAction.trigger()
            self.disconnect(CIviewer.viewer.cursorAction, SIGNAL("triggered()"),
                            CIviewer.viewer._toggleImageCursor)
            self.connect(CIviewer.viewer.cursorAction, SIGNAL("triggered()"),
                         self.cursorAction.trigger)
            self.connect(CIviewer.viewer.imageCursor, SIGNAL("positionChanged(QPoint)"),
                         self.cursorAction.broadcastPosition)
            self.connect(CIviewer.viewer, SIGNAL("zoomLevelChanged(int)"),
                         self.cursorAction.broadcastZoom)
            self.updateGeometry()
        # this call is necessary to update the sizeHint() before adjustSize() is called
        qcore.QCoreApplication.processEvents()
        self.adjustSize()

    def viewer(self, x=0, y=0):
        if self.layout.itemAtPosition(y, x):
            return self.layout.itemAtPosition(y, x).widget().viewer
        raise ValueError("ImageWindow.viewer(): viewer at (%d, %d) is undefined." % (x, y))

def showImage(image, normalize=True, title=None):
    if isinstance(image, str):
        image = vigra.impex.readImage(image)
    v = ImageWindow()
    v.setImage(image, normalize=normalize, title=title)
    v.show()
    return v
