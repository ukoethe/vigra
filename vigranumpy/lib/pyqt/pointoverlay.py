import VigraQt
from PyQt4 import QtCore, QtGui

class Overlay(VigraQt.Overlay):
    def __init__(self, parent, color = QtCore.Qt.red, fillColor = None, width = 0, name = None):
        VigraQt.Overlay.__init__(self, parent)
        self.color = color and QtGui.QColor(color)
        self.fillColor = fillColor and QtGui.QColor(fillColor)
        self.width = width
        self.name = name

    def _setupPainter(self, p):
        if self.color:
            p.setPen(QtGui.QPen(self.color, self.width))
        else:
            p.setPen(QtCore.Qt.NoPen)

        if self.fillColor:
            p.setBrush(self.fillColor)
        else:
            p.setBrush(QtCore.Qt.NoBrush)

class PointOverlay(VigraQt.Overlay):
    def __init__(self, parent, points, color = QtCore.Qt.red, radius = 0.5, 
                 colors = None, name = None):
        VigraQt.Overlay.__init__(self, parent)
        self.originalPoints = points
        self._qpointlist = QtGui.QPolygonF(len(points))
        self.color = color and QtGui.QColor(color)
        self.colors = colors
        self.name = name
        self._parent = parent
        self.radius = radius

    def _calculatePoints(self):
        if self._qpointlist.size() < len(self.originalPoints):
            self._qpointlist = QtGui.QPolygonF(len(self.originalPoints))
        d0 = 0.5
        for i, p in enumerate(self.originalPoints):
            self._qpointlist.setPoint(i, p[0] + d0, p[1] + d0)

    def setPoints(self, points):
        self.originalPoints = points
        self._qpointlist = QtGui.QPolygonF(len(points))

    def draw(self, p, r):
        p.setBrush(QtGui.QBrush(self.color))
        p.setPen(self.color)
        w = (2.0 * self.radius + 1.0) / self._parent.zoomFactor()
        for i in range(len(self.originalPoints)):
            x, y = self.originalPoints[i]
            if self.colors:
                if len(self.colors) > i:
                    p.setBrush(QtGui.QBrush(self.colors[i]))
                    p.setPen(self.colors[i])
                else:
                    p.setBrush(QtGui.QBrush(self.color))
                    p.setPen(self.color)
            p.drawEllipse(QtCore.QPointF(x,y), w, w)


