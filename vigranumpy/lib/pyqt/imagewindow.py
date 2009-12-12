import math, os
import PyQt4
import PyQt4.QtCore as qcore
import PyQt4.QtGui  as qt
from PyQt4.QtCore import SIGNAL
from VigraQt import OverlayViewer, ImageCursor
import vigra.vigranumpycmodule
import quickdialog


class ImageWindow(qt.QFrame):
    def __init__(self, image, normalize = True, title = None,
                 parent = None):
        qt.QFrame.__init__(self, parent)
        self._normalized = normalize
        if title is not None:
            self.setWindowTitle(title)
        elif hasattr(image, "name"):
            self.setWindowTitle(image.name)
        else:
            self.setWindowTitle("Image Viewer")
        self.image = image
        self._savedExpression = "x"
        self._xplaces = int(math.log10(self.image.width) + 1.0)
        self._yplaces = int(math.log10(self.image.height) + 1.0)
        self._valueplaces = image.bands()*5
        self._lastSaveType = 2

        self.viewer = OverlayViewer(self)
        self.viewer.setImage(image.qimage(normalize))
        self.viewer.setContextMenuPolicy(qcore.Qt.NoContextMenu)
        
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

        self.updateCaption(0,0)
        
        self.zoomInAction = qt.QAction("Zoom in", self)
        self.zoomInAction.setShortcut("+")
        self.connect(self.zoomInAction, SIGNAL("triggered()"), self.zoomInPopup)

        self.zoomOutAction = qt.QAction("Zoom out", self)
        self.zoomOutAction.setShortcut("-")
        self.connect(self.zoomOutAction, SIGNAL("triggered()"), self.zoomOutPopup)

        self.saveAction = qt.QAction("Save image...", self)
        self.saveAction.setShortcut("S")
        self.connect(self.saveAction, SIGNAL("triggered()"), self.writeImage)
        self.addAction(self.saveAction)

        self.expressionAction = qt.QAction("Apply expression...", self)
        # FIXME: This gives: QAction::eventFilter: Ambiguous shortcut overload: E
        #        when there are several viewers in an ImageGrid
        self.expressionAction.setShortcut("E")
        self.connect(self.expressionAction, SIGNAL("triggered()"), self.applyExpression)
        self.addAction(self.expressionAction)

        self.cursorAction = qt.QAction("Line cursor", self)
        self.cursorAction.setCheckable(True)
        self.cursorAction.setChecked(False)
        self.addAction(self.cursorAction)
        self.imageCursor = ImageCursor(self.viewer)
        self.viewer.addOverlay(self.imageCursor)
        self.imageCursor.setVisible(False)
        self.imageCursor.setPosition(qcore.QPoint(self.image.width/2, self.image.height/2))

        self.popup = qt.QMenu(self)
        self.popup.addAction(self.zoomInAction)
        self.popup.addAction(self.zoomOutAction)
        self.popup.addAction(self.saveAction)
        self.popup.addAction(self.expressionAction)
        self.popup.addAction(self.cursorAction)
#        self.popup.insertItem("Apply expression...\tE", self._menuApplyExpression)
#        self.overlayMenu = qt.QPopupMenu(self.popup)
#        self.overlayMenuIndex = self.popup.insertItem("&Overlays", self.overlayMenu)
#        self.connect(self.overlayMenu, SIGNAL("aboutToShow()"), self.overlayPopup)
#        #self.connect(self.overlayMenu, SIGNAL("activated(int)"), self.toggleOverlay)
#        if self.isTopLevel():
#            self.popup.insertItem("&Close image", self.close)
#        self.viewer.connect(self.viewer, SIGNAL("mousePressed"), self.openPopup)
#        self.connect(self, SIGNAL("mousePressed"), self.openPopup)

    def updateCaption(self, x, y):
        x = int(round(x))
        y = int(round(y))
        if x < 0 or x >= self.image.width or \
           y < 0 or y >= self.image.height:
            return

        label = str(x).rjust(self._xplaces) + " x " + str(y).rjust(self._yplaces) +\
                " = " + str(self.image[x,y]).ljust(self._valueplaces)
        self.label.setText(label)
        self.emit(SIGNAL('captionChanged'), self.label.text())

    def updateCaptionP(self, point):
        self.updateCaption(point.x(), point.y())
        
    def overlayPopup(self):
        self.overlayMenu.clear()
        index = 0
        hideable = False
        showable = False
        for o in self.viewer.overlays:
            if hasattr(o, "name"):
                overlayName = o.name
            else:
                overlayName = self.viewer._defaultOverlayName(o)
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
                iconSet = qt.QIconSet()
                iconSet.setPixmap(colorPM, qt.QIconSet.Automatic,
                                  qt.QIconSet.Normal, qt.QIconSet.On)
                id = self.overlayMenu.insertItem(iconSet,
                       text, self.viewer.toggleOverlayVisibility)
            else:
                id = self.overlayMenu.insertItem(
                    text, self.viewer.toggleOverlayVisibility)

            self.overlayMenu.setItemChecked(id, o.visible)
            self.overlayMenu.setItemParameter(id, index)
            if o.visible:
                hideable = True
            else:
                showable = True
            index += 1
        id = self.overlayMenu.insertItem("&Hide all", self.viewer.showOverlays)
        self.overlayMenu.setItemParameter(id, False)
        self.overlayMenu.setItemEnabled(id, hideable)
        id = self.overlayMenu.insertItem("&Show all", self.viewer.showOverlays)
        self.overlayMenu.setItemParameter(id, True)
        self.overlayMenu.setItemEnabled(id, showable)

    def _slideAfterZoom(self, shift):
        if self.viewer.zoomLevel() > 0:
            shift *= 1 + self.viewer.zoomLevel()
        elif self.viewer.zoomLevel() < 0:
            shift /= 1 - self.viewer.zoomLevel()
        self.viewer.slideBy(shift)
    
    def showImageCursor(self, yesOrNo = True):
        if yesOrNo != self.cursorAction.isChecked():
            self.cursorAction.trigger()

    def _toggleImageCursor(self):
        self.imageCursor.activateTool(self.cursorAction.isChecked())
        self.imageCursor.setVisible(self.cursorAction.isChecked())
        if self.cursorAction.isChecked():
            self.disconnect(self.viewer, SIGNAL('mouseOver(int, int)'), self.updateCaption)
            self.connect(self.imageCursor, SIGNAL("positionChanged(QPoint)"), self.updateCaptionP)
        else:
            self.disconnect(self.imageCursor, SIGNAL("positionChanged(QPoint)"), self.updateCaptionP)
            self.connect(self.viewer, SIGNAL('mouseOver(int, int)'), self.updateCaption)            

    def zoomInPopup(self):
        beforePos = self.viewer.imageCoordinate(self.mousepos)
        self.viewer.zoomUp()
        afterPos = self.viewer.imageCoordinate(self.mousepos)
        self._slideAfterZoom(afterPos - beforePos)

    def zoomOutPopup(self):
        beforePos = self.viewer.imageCoordinate(self.mousepos)
        self.viewer.zoomDown()
        afterPos = self.viewer.imageCoordinate(self.mousepos)
        self._slideAfterZoom(afterPos - beforePos)

    def _menuApplyExpression(self):
        self.applyExpression()

    def applyExpression(self, expr = None, normalized = None):
        if expr != None:
            self._savedExpression = expr
        else:
            d = quickdialog.QuickDialog(self,"Enter Expression")
            d.expression = quickdialog.OptionalStringInput(d, "Execute 'lambda x: ")
            d.expression.setText(self._savedExpression)
            d.addSpacing(10)
            d.norm = quickdialog.CheckBox(d, "Normalize intensity to range 0...255")
            d.norm.setChecked(self._normalized)
            if d.exec_() == 0:
                return
            self._savedExpression = d.expression.text()
            self._normalized = True if d.norm.selection() else False

        if normalized != None:
            self._normalized = normalized

        try:
            image, normalized = self.getDisplay()
        except Exception, e:
            qt.QMessageBox.critical(self, "Error Applying Expression", str(e))
            return

        self.viewer.setImage(image.qimage(normalized))

    def getDisplay(self):
        """Returns the displayed image and the normalize flag
        (BYTE or NBYTE) as tuple/pair.

        Note that the returned image is the original image if no
        expression is applied, i.e. you should not change the returned
        object.  If active, the expression is applied via
        transformImage() on every call of getDisplay()."""

        if not self._savedExpression or self._savedExpression == "x":
            self._savedExpression = "x"
            image = self.image
        else:
            x = self.image
            image = eval(self._savedExpression)

        return image, self._normalized

    def keyPressEvent(self, e):
        "handles keys [S], [E], and possibly [Q] (for toplevel-windows)"
        if e.key() == qcore.Qt.Key_Q and not self.parent():
            self.close()
            e.accept()
        elif e.key() == qcore.Qt.Key_S:
            self.writeImage()
            e.accept()
        elif e.key() == qcore.Qt.Key_E:
            self.applyExpression()
            e.accept()

    def contextMenuEvent(self, e):
        "handles pop-up menu"
#        self.popup.setItemEnabled(
#            self.overlayMenuIndex, len(self.viewer.overlays) > 0)
        self.mousepos = e.pos()
        self.popup.exec_(e.globalPos())
        e.accept()

    def replaceImage(self, image, normalize = None):
        """imageWindow.replaceImage(image, normalize = None)

        Replace the current image with the given one.  If normalized
        is not given (or None), the normalized state is not changed."""
        
        self.image = image
        self.applyExpression(self._savedExpression, normalize)

    def writeImage(self):
        d = quickdialog.QuickDialog(self,"Write Image")

        imageFileExtensions = '*.' + ' *.'.join(vigranumpycmodule.impexListExtensions().split(' '))
        d.filedialog = quickdialog.OutputFile(
            d, "Output filename:", "Image Files ("+imageFileExtensions+")")
        d.filedialog.filename.setFocus()

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
                           "jpg": "JPEG", \
                           "jpeg": "JPEG"}

                _, ext = os.path.splitext(filename)
                if not formats.has_key(ext[1:]):
                    f = " ".join(formats.keys())
                    qt.QMessageBox.critical(self, "Error", \
                                   "Displayed image with overlays can only be stored as\n"+f)
                else:
                    pixmap = self.viewer.getContentsPixmap()
                    # FIXME: the next fails when the extension has 4 characters
                    pixmap.save(filename, formats[filename[-3:]])
                    return

class CursorAction(qt.QAction):
    def trigger(self):
        qt.QAction.trigger(self)
        for v in self.viewers:
            v.cursorAction.setChecked(self.isChecked())
            v._toggleImageCursor()

    def broadcastPosition(self, pos):
        self.x, self.y = pos.x(), pos.y()
        for v in self.viewers:
            v.imageCursor.setPosition(pos)
        
    
class ImageGrid(qt.QFrame):
    def __init__(self, parent = None):
        qt.QFrame.__init__(self, parent)
        self.cursorAction = CursorAction("Connected line cursors", self)
        self.cursorAction.setCheckable(True)
        self.cursorAction.setChecked(False)
        self.addAction(self.cursorAction)
        self.cursorAction.viewers = []
        self.layout = qt.QGridLayout(self)
        
    def addViewer(self, viewer, x, y):
        self.layout.addWidget(viewer, y, x)
        self.cursorAction.viewers.append(viewer)
        if self.cursorAction.isChecked():
            viewer.cursorAction.setChecked(True)
            viewer._toggleImageCursor()
            viewer.imageCursor.setPosition(qcore.QPoint(self.cursorAction.x, self.cursorAction.y))
        self.disconnect(viewer.cursorAction, SIGNAL("triggered()"), viewer._toggleImageCursor)
        self.connect(viewer.cursorAction, SIGNAL("triggered()"), self.cursorAction.trigger)
        self.connect(viewer.imageCursor, SIGNAL("positionChanged(QPoint)"), self.cursorAction.broadcastPosition)
        # this call is necessary due to a Qt bug (Issue N7468) ???
        qt.QApplication.sendPostedEvents(None, qcore.QEvent.LayoutRequest)
        # FIXME: this doesn't do the desired window resize 
        # (but an explicit call of adjustSize() after returning from this function works
        self.adjustSize()

def showImage(image, normalize = True):
    if isinstance(image, str):
        image = vigranumpycmodule.readImage(image)
    v = ImageWindow(image, normalize)
    v.show()
    return v

