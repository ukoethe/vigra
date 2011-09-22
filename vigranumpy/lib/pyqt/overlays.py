import VigraQt
from PyQt4 import QtCore, QtGui

class Overlay(VigraQt.Overlay):
    def __init__(self, parent, color = QtCore.Qt.red, fillColor = None, width = 0, name = None, aa = False):
        VigraQt.Overlay.__init__(self, parent)
        self.color = color and QtGui.QColor(color)
        self.fillColor = fillColor and QtGui.QColor(fillColor)
        self.width = width
        self.name = name
	self.setAntialiasing = aa

    def _setupPainter(self, p):
        if self.color:
            p.setPen(QtGui.QPen(self.color, self.width))
        else:
            p.setPen(QtCore.Qt.NoPen)

        if self.fillColor:
            p.setBrush(self.fillColor)
        else:
            p.setBrush(QtCore.Qt.NoBrush)


class PointOverlay(Overlay):
    def __init__(self, parent, points, color = QtCore.Qt.red, fillColor = QtCore.Qt.red,
                 radius = 0.5, colors = None, name = None, aa = False, isSuboverlay = False):
        Overlay.__init__(self, parent, color, fillColor, name = name, aa = aa)
        self.originalPoints = points
        self.color = color and QtGui.QColor(color)
        self.colors = colors
        self.name = name
        self._parent = parent
        self.radius = radius
        if not isSuboverlay:
            parent.addOverlay(self)

    def draw(self, p, r):
        visibleRect = QtCore.QRectF(VigraQt.OverlayViewer.imageCoordinateF(self._parent, r.topLeft()), 
          VigraQt.OverlayViewer.imageCoordinateF(self._parent, r.bottomRight()))
        w = (2.0 * self.radius + 1.0) / self._parent.zoomFactor()
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
    def __init__(self, parent, edges, color = QtCore.Qt.red, colors = None, fillColor = None, 
                 width = 0, name = None, aa = False, isSuboverlay = False):
        Overlay.__init__(self, parent, color, fillColor, width, name, aa)
        # input should be list of Polygons, but also accept single Polygon
        self._parent = parent
        self.width = width
        self.colors = colors
        if not hasattr(edges[0][0], "__getitem__"):
            self.originalEdges = [edges]
        else:
            self.originalEdges = edges
        if not isSuboverlay:
            parent.addOverlay(self)

    def draw(self, p, r):
        visibleRect = QtCore.QRectF(VigraQt.OverlayViewer.imageCoordinateF(self._parent, r.topLeft()), 
          VigraQt.OverlayViewer.imageCoordinateF(self._parent, r.bottomRight()))
        self._setupPainter(p)
        if self.colors:
            for j, polygon in enumerate(self.originalEdges):
                qpolf = QtGui.QPolygonF(len(polygon))
                for i, (x, y) in enumerate(polygon):
                    qpolf[i] = QtCore.QPointF(x,y)
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
                    qpolf[i] = QtCore.QPointF(x,y)
                if qpolf.boundingRect().intersects(visibleRect):
                    p.drawPolygon(qpolf) if qpolf.isClosed() else p.drawPolyline(qpolf)


class TextOverlay(Overlay):
    def __init__(self, parent, textlist, name = None, aa = False, isSuboverlay = False):
        Overlay.__init__(self, parent, name = name, aa = aa)
        self.textlist = textlist # textlist = [["text", [PosX, PosY], optional color, optional pointsize], [...], ...]
        self.name = name
        self.setCoordinateSystem(VigraQt.Overlay.UnscaledPixel)
        self._parent = parent
        if not isSuboverlay:
            parent.addOverlay(self)
        
    def draw(self, p, r):
        visibleRect = QtCore.QRectF(VigraQt.OverlayViewer.imageCoordinateF(self._parent, r.topLeft()), 
          VigraQt.OverlayViewer.imageCoordinateF(self._parent, r.bottomRight()))

        for entry in self.textlist:
            text = entry[0]
            position =  QtCore.QPointF(*map(lambda x: x * self._parent.zoomFactor(), entry[1]))
            if len(entry) > 2:
                p.setPen(QtGui.QPen(entry[2] and QtGui.QColor(entry[2]), 0))
            else:
                self._setupPainter(p)
            if len(entry) > 3:
                defaultFont = p.font()
                defaultFont.setPointSizeF(entry[3])
                p.setFont(defaultFont)
            twidth, theight = p.fontMetrics().width(text), p.fontMetrics().height()
            textRect = QtCore.QRectF(position.x() - twidth / 2.0, position.y() - theight / 2.0, twidth, theight)
            displayRect = QtCore.QRectF(entry[1][0] - twidth / 2.0, entry[1][1] - theight / 2.0, twidth, theight)
            if displayRect.intersects(visibleRect):
                p.drawText(textRect, QtCore.Qt.AlignCenter, text)

class MapOverlay(Overlay):   
    def __init__(self, parent, geomap, edgeColor = QtCore.Qt.blue, 
                   nodeColor = QtCore.Qt.red, name = None, aa = False, isSuboverlay = False):
        Overlay.__init__(self, parent, name = name, aa = aa)
        self._parent = parent
        self.edges, self.nodes = [], []
        for edge in geomap.edgeIter():
    	    self.edges.append(edge)
        for node in geomap.nodeIter():
            self.nodes.append(node.position())
        self.eo = EdgeOverlay(self._parent, self.edges, edgeColor, isSuboverlay = True)
        self.po = PointOverlay(self._parent, self.nodes, nodeColor, radius = 0.3, isSuboverlay = True)
        if not isSuboverlay:
            parent.addOverlay(self)

    def draw(self, p, r):
        self.eo.draw(p, r)
        self.po.draw(p, r)

