import os 
from PyQt4 import QtCore, QtGui

def viewer2svg(viewer, basepath, onlyVisible = False, moveBy = QtCore.QPointF(0.5, 0.5)):
    outvec=[]
    outvec.append('<?xml version="1.0" standalone="no"?>\n')
    outvec.append('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n')
    outvec.append('  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')

    pngFilename = basepath + "_bg.png"
    viewer.image.writeImage(pngFilename, "")
    _, bgFilename = os.path.split(pngFilename)

    outvec.append('\n<svg version="1.1" xmlns="http://www.w3.org/2000/svg"\n')
    outvec.append('  width="21cm" height="29.7cm" preserveAspectRatio="xMinYMin meet"\n')
    outvec.append('  viewBox="-1 -1 ' + str(viewer.image.width + 1) + ' ' + \
      str(viewer.image.height + 1) + '" xmlns:xlink="http://www.w3.org/1999/xlink">\n')
    outvec.append('\n<g style="fill:none">\n')

    outvec.append('\n<image xlink:href="' + bgFilename + '"\n')
    outvec.append('  x="0" y="0" width="' + str(viewer.image.width) + '" height="' + \
      str(viewer.image.height) + '" />\n')

    ovs = []
    for ov in viewer.overlays:
        if onlyVisible and not ov.isVisible():
            continue
        ovname = viewer._defaultOverlayName(ov)
        if ovname == "MapOverlay":
            ovs.append([viewer._defaultOverlayName(ov.eo), ov.eo])
            ovs.append([viewer._defaultOverlayName(ov.po), ov.po])
        else:
            ovs.append([ovname, ov])            

    for overlay in ovs:
        if overlay[0] == "EdgeOverlay":
            overlay = overlay[1]
            color = 'rgb' + str(overlay.color.getRgb()[:3]) + '; opacity:' + str(overlay.color.getRgb()[-1] / 255.0)
            if not overlay.colors:
                for i, edge in enumerate(overlay.originalEdges):
                    outvec.append(writeEdge(edge, overlay.width, color, moveBy))
            else:
                for i, edge in enumerate(overlay.originalEdges):
                    if len(overlay.colors) > i:
                        color = overlay.colors[i] if hasattr(overlay.colors[i], "getRgb") else \
                          QtGui.QColor(overlay.colors[i])
                        color = 'rgb' + str(color.getRgb()[:3]) + '; opacity:' + str(color.getRgb()[-1] / 255.0)
                    outvec.append(writeEdge(edge, overlay.width, color, moveBy))
                
        elif overlay[0] == "PointOverlay":
            overlay = overlay[1]
            color = '  style="fill:rgb' + str(overlay.color.getRgb()[:3]) + '; opacity:' + str(overlay.color.getRgb()[-1] / 255.0) + '"/>\n'
            radius = '" r="' + str(overlay.radius if overlay.radius > 0 else 0.5) + '"\n'
            pointList = []
            for point in overlay.originalPoints:
                pointList.append(QtCore.QPointF(*point) + moveBy)
            for point in pointList:
                outvec.append('<circle cx="' + str(point.x()) + '" cy="' + str(point.y()) + radius + color)
        elif overlay[0] == "TextOverlay":
            overlay = overlay[1]
            for element in overlay.textlist:
                if len(element) == 4:
                    outvec.extend(writeText(text = element[0], position = element[1], color = element[2], size = element[3]))
                elif len(element) == 3:
                    outvec.extend(writeText(text = element[0], position = element[1], color = element[2]))
                else:
                    outvec.extend(writeText(text = element[0], position = element[1]))
        else:
            print str(overlay[0]) + " not supported yet.\n"

    outvec.append('\n</g>\n')
    outvec.append('</svg>\n')

    f = open(basepath + ".svg", 'w')
    for line in outvec:
        f.write(line)
    f.close()

def writeEdge(edge, width, color, moveBy):
    qpolf = QtGui.QPolygonF(len(edge))
    for i, (x, y) in enumerate(edge):
        qpolf[i] = QtCore.QPointF(x,y) + moveBy
    result = "\n"
    if qpolf.size() == 2:
        result += '<line x1="' + str(qpolf[0].x()) + '" y1="' + str(qpolf[0].y()) + '" '
        result += 'x2="' + str(qpolf[1].x()) + '" y2="' + str(qpolf[1].y())
    elif qpolf.size() > 2:
        result += '<polyline points="' + str(qpolf[0].x()) + '" y1="' + str(qpolf[0].y())
        for pos in range(1, qpolf.size()):
            result += ' ' + str(qpolf[pos].x()) + '" y1="' + str(qpolf[pos].y())
    result += '"\n  style="stroke:' + color + '; stroke-width:' + str(width if width > 0 else 0.5) + ';\n'
    result += '  stroke-linejoin:bevel; stroke-linecap:butt;"/>\n'
    return result

def writeText(text, position, color = None, size = None):
    if not size:
        size= "6"
    if not color:
        color = 'fill:rgb(0, 255, 0); opacity:1; stroke:rgb(0, 0, 0); stroke-width:0.3;'
    else:
        color = 'fill:rgb' + str(QtGui.QColor(color).getRgb()[:3]) + '; opacity:' + \
          str(QtGui.QColor(color).getRgb()[-1] / 255.0) + '; stroke:rgb(0, 0, 0); stroke-width:0.3;'
    style = '  style="' + color + '\n    dominant-baseline: central; ' + \
      'text-anchor: middle; font-size: ' + str(size) + 'pt; font-family: sans-serif"'

    return '\n<text x="' + str(position[0]) + '" y="' + str(position[1]) + '"\n' + \
            style + '>' + text.toUtf8().data() + '</text>\n'

