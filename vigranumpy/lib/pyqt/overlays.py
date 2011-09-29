import VigraQt
from PyQt4 import QtCore, QtGui

class Overlay(VigraQt.Overlay):
    def __init__(self, color=QtCore.Qt.red, fillColor=None, width=0,
                 name=None, aa=False, parent=None):
        VigraQt.Overlay.__init__(self, parent)
        self.color = color and QtGui.QColor(color)
        self.fillColor = fillColor and QtGui.QColor(fillColor)
        self.width = width
        self.name = name
        self.setAntialiased(aa)

    def _setupPainter(self, p):
        if self.color:
            p.setPen(QtGui.QPen(self.color, self.width))
        else:
            p.setPen(QtCore.Qt.NoPen)

        if self.fillColor:
            p.setBrush(self.fillColor)
        else:
            p.setBrush(QtCore.Qt.NoBrush)

class OverlayGroup(Overlay):
    def __init__(self, overlays, color=QtCore.Qt.red, fillColor=QtCore.Qt.red,
                 name=None, parent=None):
        Overlay.__init__(self, color, fillColor, parent=parent)
        self.overlays = overlays
        self.name = name
        # Coordinate and aliasing settings from individual overlays
        # are used so following setting won't get effective
        self.setCoordinateSystem(VigraQt.Overlay.UnscaledPixel)
        if parent:
            parent.addOverlay(self)
            self.setParent(parent)

    def setParent(self, parent):
        Overlay.setParent(self, parent)
        for o in self.overlays:
            o.setParent(parent)

    def draw(self, p, r):
        for o in self.overlays:
            p.save()
            if self.parent():
                if o.coordinateSystem() & VigraQt.Overlay.Scaled:
                    p.scale(self.parent().zoomFactor(), self.parent().zoomFactor())
                if o.coordinateSystem() & VigraQt.Overlay.Pixel:
                    p.translate(0.5, 0.5)
                p.setRenderHint(QtGui.QPainter.Antialiasing, o.isAntialiased())
            o.draw(p, r)
            p.restore()

class PointOverlay(Overlay):
    def __init__(self, points, color=QtCore.Qt.red, fillColor=QtCore.Qt.red,
                 radius=0.5, colors=None, name=None, aa=False, parent=None):
        Overlay.__init__(self, color, fillColor, aa=aa, parent=parent)
        self.originalPoints = points
        self.color = color and QtGui.QColor(color)
        self.colors = colors
        self.name = name
        self.radius = radius
        if parent:
            parent.addOverlay(self)

    def draw(self, p, r):
        visibleRect = QtCore.QRectF(
          VigraQt.OverlayViewer.imageCoordinateF(self.parent(), r.topLeft()),
          VigraQt.OverlayViewer.imageCoordinateF(self.parent(), r.bottomRight()))
        w = (2.0 * self.radius + 1.0) / self.parent().zoomFactor()
        if self.colors:
            for i in range(len(self.originalPoints)):
                point = QtCore.QPointF(*self.originalPoints[i])
                if visibleRect.contains(point):
                    if len(self.colors) > i:
                        p.setBrush(QtGui.QBrush(self.colors[i]))
                        p.setPen(self.colors[i])
                    else:
                        self._setupPainter(p)
                    p.drawEllipse(point, w, w)
        else:
            self._setupPainter(p)
            for i in range(len(self.originalPoints)):
                point = QtCore.QPointF(*self.originalPoints[i])
                if visibleRect.contains(point):
                    p.drawEllipse(point, w, w)


class EdgeOverlay(Overlay):
    def __init__(self, edges, color=QtCore.Qt.red, colors=None, fillColor=None,
                 width=0, name=None, aa=False, parent=None):
        Overlay.__init__(self, color, fillColor, width, name, aa, parent)
        # input should be list of Polygons, but also accept single Polygon
        self.width = width
        self.colors = colors
        if not hasattr(edges[0][0], "__getitem__"):
            self.originalEdges = [edges]
        else:
            self.originalEdges = edges
        if parent:
            parent.addOverlay(self)

    def draw(self, p, r):
        visibleRect = QtCore.QRectF(
          VigraQt.OverlayViewer.imageCoordinateF(self.parent(), r.topLeft()),
          VigraQt.OverlayViewer.imageCoordinateF(self.parent(), r.bottomRight()))
        self._setupPainter(p)
        if self.colors:
            for j, polygon in enumerate(self.originalEdges):
                qpolf = QtGui.QPolygonF(len(polygon))
                for i, (x, y) in enumerate(polygon):
                    qpolf[i] = QtCore.QPointF(x, y)
                if qpolf.boundingRect().intersects(visibleRect):
                    if len(self.colors) > j:
                        p.setBrush(QtGui.QBrush(self.colors[j]))
                        p.setPen(self.colors[j])
                    else:
                        self._setupPainter(p)
                    p.drawPolygon(qpolf) if qpolf.isClosed() else p.drawPolyline(qpolf)
        else:
            for polygon in self.originalEdges:
                qpolf = QtGui.QPolygonF(len(polygon))
                for i, (x, y) in enumerate(polygon):
                    qpolf[i] = QtCore.QPointF(x, y)
                if qpolf.boundingRect().intersects(visibleRect):
                    p.drawPolygon(qpolf) if qpolf.isClosed() else p.drawPolyline(qpolf)


class TextOverlay(Overlay):
    def __init__(self, text, pos, color=QtCore.Qt.red, pointsize=None,
                 zoomText=False, name=None, aa=False, parent=None):
        Overlay.__init__(self, color, color, aa=aa, parent=parent)
        self.text = text
        self.pos = pos
        self.pointsize = pointsize
        self.name = name
        if zoomText:
            self.setCoordinateSystem(VigraQt.Overlay.ScaledPixel)
        else:
            self.setCoordinateSystem(VigraQt.Overlay.UnscaledPixel)
        if parent:
            parent.addOverlay(self)

    def draw(self, p, r):
        visibleRect = QtCore.QRectF(
          VigraQt.OverlayViewer.imageCoordinateF(self.parent(), r.topLeft()),
          VigraQt.OverlayViewer.imageCoordinateF(self.parent(), r.bottomRight()))

        if self.coordinateSystem() & VigraQt.Overlay.Scaled:
            position = QtCore.QPointF(*self.pos)
        else:
            position = QtCore.QPointF(*map(
              lambda x: x * self.parent().zoomFactor(), self.pos + [0.5, 0.5]))

        self._setupPainter(p)
        if self.pointsize:
            defaultFont = p.font()
            defaultFont.setPointSizeF(self.pointsize)
            p.setFont(defaultFont)
        twidth = p.fontMetrics().width(self.text)
        theight = p.fontMetrics().height()
        textRect = QtCore.QRectF(position.x() - twidth / 2.0,
          position.y() - theight / 2.0, twidth, theight)
        if self.coordinateSystem() & VigraQt.Overlay.Unscaled:
            twidth /= self.parent().zoomFactor()
            theight /= self.parent().zoomFactor()
            displayRect = QtCore.QRectF(self.pos[0] - twidth / 2.0,
              self.pos[1] - theight / 2.0, twidth, theight)
        else:
            displayRect = textRect
        if displayRect.intersects(visibleRect):
            p.drawText(textRect, QtCore.Qt.AlignCenter, self.text)

class MapOverlay(Overlay):
    def __init__(self, geomap, edgeColor=QtCore.Qt.blue,
                 nodeColor=QtCore.Qt.red, name=None, aa=False, parent=None):
        Overlay.__init__(self, aa=aa, parent=parent)
        self.edges, self.nodes = [], []
        for edge in geomap.edgeIter():
            self.edges.append(edge)
        for node in geomap.nodeIter():
            self.nodes.append(node.position())
        self.eo = EdgeOverlay(self.edges, edgeColor, parent=None)
        self.po = PointOverlay(self.nodes, nodeColor, radius=0.3, parent=None)
        self.eo.setParent(parent)
        self.po.setParent(parent)
        if parent:
            parent.addOverlay(self)

    def setParent(self, parent):
        Overlay.setParent(self, parent)
        self.eo.setParent(parent)
        self.po.setParent(parent)

    def draw(self, p, r):
        self.eo.draw(p, r)
        self.po.draw(p, r)
