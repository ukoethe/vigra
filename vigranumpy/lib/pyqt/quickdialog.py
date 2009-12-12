import PyQt4.QtGui  as qt
import PyQt4.QtCore as qcore
from PyQt4.QtCore import SIGNAL
import os

def alignLabels(*args):
    m = 0
    for dialogElement in args:
        fontMetrics = qt.QFontMetrics(dialogElement.font())
        for line in dialogElement.label.text().ascii().split('\n'):
            labelWidth = fontMetrics.width(line)
            m = max(m, labelWidth)
    for dialogElement in args:
        dialogElement.label.setFixedWidth(m+10)

class FileDialog(qt.QFrame):
    def __init__(self, parent, label, filter):
        qt.QFrame.__init__(self, parent)
        parent.addWidget(self)
        self.filter = filter
        self.label = qt.QLabel(label)
        self.filename = qt.QLineEdit()
        self.filebrowser = qt.QPushButton("Browse...")
        self.filebrowser.setFocusPolicy(qcore.Qt.NoFocus)

        self._layout = qt.QHBoxLayout()
        self._layout.setSpacing(5)
        self._layout.addWidget(self.label)
        self._layout.addWidget(self.filename, 1)
        self._layout.addWidget(self.filebrowser)
        
        self.setLayout(self._layout)
            
    def text(self):
        return str(qcore.QFile.encodeName(self.filename.text()))

class InputFile(FileDialog):
    def __init__(self, parent, label, filter):
        FileDialog.__init__(self, parent, label, filter)
        self.connect(self.filebrowser, SIGNAL("clicked()"), self.browse)
        
    def browse(self):
        fn = qt.QFileDialog.getOpenFileName( "", self.filter, self)
        if not fn.isNull():
            self.filename.setText(fn)
        
    def validate(self):
        try:
            filename = str(qcore.QFile.encodeName(self.filename.text()))
            file = open(filename)
            file.close()
            return True
        except IOError:
            qt.QMessageBox.critical(None, "Error", "File '" + filename + "' not found")
            return False

class OutputFile(FileDialog):
    def __init__(self, parent, label, filter):
        FileDialog.__init__(self, parent, label, filter)
        self.connect(self.filebrowser, SIGNAL("clicked()"), self.browse)
        
    def browse(self):
        fn = qt.QFileDialog.getSaveFileName( self, "Save File", "", self.filter)
        if not fn.isNull():
            self.filename.setText(fn)
        
    def validate(self):
        try:
            filename = str(qcore.QFile.encodeName(self.filename.text()))
            file = open(filename)
            file.close()
            return not qt.QMessageBox.warning(
                None, "Warning", "File '" + filename + "' exists",
                "Overwrite", "Cancel")
        except IOError:
            return True

class _OptionalValueInput(qt.QFrame):
    def __init__(self, parent, label):
        qt.QFrame.__init__(self, parent)
        parent.addWidget(self)
        self.label = qt.QLabel(label)
        self.variable = qt.QLineEdit()
        self.variable.setValidator(self._QValidator(self.variable))

        self._layout = qt.QHBoxLayout()
        self._layout.setSpacing(5)
        self._layout.addWidget(self.label)
        self._layout.addWidget(self.variable, 1)
        
        self.setLayout(self._layout)
            
    def setValue(self, text):
        self.variable.setText(str(self._text2Value(text)))
    
    def value(self):
        text = self.text()
        if text == "":
            return None
        return self._text2Value(text)

    def text(self):
        return str(self.variable.text())
        
    def validate(self):
        try:
            v = self.value()
            if v == None:
                return True
        except:
            qt.QMessageBox.critical(None, "Error",
                qcore.QString("Field '%1' must contain "+self._mustContain).arg(
                    self.label.text()))
            return False
        try:
            if v < self.min:
                qt.QMessageBox.critical(None, "Error",
                    qcore.QString("Field '%1' value must be >= "+str(self.min)).arg(
                        self.label.text()))
                return False
        except AttributeError:
            pass
        try:
            if v > self.max:
                qt.QMessageBox.critical(None, "Error",
                    qcore.QString("Field '%1' value must be <= "+str(self.max)).arg(
                        self.label.text()))
                return False
        except AttributeError:
            pass
        return True
            
class OptionalIntegerInput(_OptionalValueInput):
    _QValidator = qt.QIntValidator
    _text2Value = int
    _mustContain = "an integer"

class IntegerInput(OptionalIntegerInput):
    def value(self):
        return int(self.text())

class OptionalFloatInput(_OptionalValueInput):
    _QValidator = qt.QDoubleValidator
    _text2Value = float
    _mustContain = "a float"

class FloatInput(OptionalFloatInput):
    def value(self):
        return float(self.text())

class OptionalStringInput(qt.QFrame):
    def __init__(self, parent, label):
        qt.QFrame.__init__(self, parent)
        parent.addWidget(self)
        self.label = qt.QLabel(label)
        self.variable = qt.QLineEdit()

        self._layout = qt.QHBoxLayout()
        self._layout.setSpacing(5)
        self._layout.addWidget(self.label)
        self._layout.addWidget(self.variable, 1)
        
        self.setLayout(self._layout)
            
    def setText(self, text):
        self.variable.setText(text)
    
    def text(self):
        return str(self.variable.text())
    
    def unicode(self):
        return unicode(self.variable.text())

class StringInput(OptionalStringInput):
    def __init__(self, parent, label):
        OptionalStringInput.__init__(self, parent, label)
            
    def validate(self):
        if self.text() == "":
            qt.QMessageBox.critical(
                None, "Error",
                qcore.QString("Field '%1' empty").arg(self.label.text()))
            return False
        return True

OutputVariable = StringInput
InputVariable = StringInput
OptionalInputVariable = OptionalStringInput

class CheckBox(qt.QCheckBox):
    def __init__(self, parent, label):
        qt.QCheckBox.__init__(self, label, parent)
        parent.addWidget(self)

    def selection(self):
        return self.isChecked()


class Choice(qt.QFrame):
    def __init__(self, parent, label, vertical = 0):
        qt.QFrame.__init__(self, parent)
        parent.addWidget(self)
        
        self.buttonBox = qt.QGroupBox(label, self)
        if vertical:
            self.buttonBox.layout = qt.QVBoxLayout(self.buttonBox)
        else:
            self.buttonBox.layout = qt.QHBoxLayout(self.buttonBox)
        
        self.layout = qt.QHBoxLayout(self)
        self.layout.addWidget(self.buttonBox)
        self.layout.addStretch(5)
        
        self.buttons = []
        self.results = []
    
    def addButton(self, label, result):
        self.buttons.append(qt.QRadioButton(label))
        self.buttonBox.layout.addWidget(self.buttons[-1])
        self.results.append(result)
        self.buttons[0].setChecked(True)
        
    def addSpacing(self, spacing):
        self.buttonBox.addSpace(spacing)
        
    def selectButton(self, index):
        if index >= 0 and index < len(self.buttons):
            self.buttons[index].setChecked(True)
        
    def selection(self):
        for k in range(len(self.buttons)):
            if self.buttons[k].isChecked():
                return self.results[k]
        return None # should never happen

class HChoice(Choice):
    def __init__(self, parent, label):
        Choice.__init__(self, parent, label, 0)
        
class VChoice(Choice):
    def __init__(self, parent, label):
        Choice.__init__(self, parent, label, 1)
        
class DialogGroup(qt.QFrame):
    def __init__(self, parent, vertical = 0):
        qt.QFrame.__init__(self, parent)
        parent.addWidget(self)
        if vertical:
            self.layout = qt.QVBoxLayout(self)
            self.defaultAlignment = qcore.Qt.AlignLeft
        else:
            self.layout = qt.QHBoxLayout(self)
            self.defaultAlignment = qcore.Qt.AlignTop
        self.widgets = []
               
    def addWidget(self, widget, stretch = 0, alignment = None):
        if alignment is None:
            alignment = self.defaultAlignment
        self.layout.addWidget(widget, stretch, alignment)
        self.widgets.append(widget)
        
    def addSpacing(self, spacing):
        self.layout.addSpacing(spacing)

    def addStretch(self, stretch):
        self.layout.addStretch(stretch)
        
    def addLabel(self, labelString):
        label = qt.QLabel(labelString, self)
        self.addWidget(label, 0, qcore.Qt.AlignLeft)
        
    def validate(self):
        for i in self.widgets:
            try:
                if i.validate() == 0:
                    return False
            except AttributeError:
                continue
        return True

class HDialogGroup(DialogGroup):
    def __init__(self, parent):
        DialogGroup.__init__(self, parent, 0)
        
class VDialogGroup(DialogGroup):
    def __init__(self, parent):
        DialogGroup.__init__(self, parent, 1)
       
#class DialogStack(qt.QWidgetStack):
#    def __init__(self, parent, widgetMapping = None):
#        qt.QWidgetStack.__init__(self, parent)
#        parent.addWidget(self)
#        self.widgetMapping = widgetMapping
#        self.size = 0
#            
#    def raiseWidget(self, index):
#        if self.widgetMapping:
#            qt.QWidgetStack.raiseWidget(self, self.widgetMapping[index])
#        else:
#            qt.QWidgetStack.raiseWidget(self, index)
#
#    def addWidget(self, widget):
#        qt.QWidgetStack.addWidget(self, widget, self.size)
#        self.size = self.size + 1
#        
#    def validate(self):
#        try:
#            return self.visibleWidget().validate()
#        except AttributeError:
#            pass

class QuickDialog(qt.QDialog):
    def __init__(self, parent, title):
        qt.QDialog.__init__(self, parent)

        self.layout = qt.QVBoxLayout(self)
        self.layout.addStretch(5)
        self.layout.addSpacing(20)
        
        self.insertButtons()
        
        self.widgets = []
        self.setWindowTitle(title)
        self.setOrientation(qcore.Qt.Vertical)
        self.resize(500,-1)
        
    def insertButtons(self):
        self.buttons = qt.QFrame(self)
        self.buttons.OK = qt.QPushButton("OK", self.buttons)
        self.buttons.Cancel = qt.QPushButton("Cancel", self.buttons)
        self.buttons.OK.setDefault(1)
        self.connect(self.buttons.Cancel, SIGNAL("clicked()"), self.reject)
        self.connect(self.buttons.OK, SIGNAL("clicked()"), self.tryAccept)
        
        self.buttons.layout = qt.QHBoxLayout(self.buttons)
        self.buttons.layout.addStretch(5)
        self.buttons.layout.addWidget(self.buttons.OK)
        self.buttons.layout.addWidget(self.buttons.Cancel)
        self.layout.addWidget(self.buttons)
        
    def addWidget(self, widget, stretch = 0, alignment = None):
        if alignment is None:
            alignment = qcore.Qt.AlignTop
        self.layout.insertWidget(len(self.widgets), widget, stretch, alignment)
        self.widgets.append(widget)
        
    def addSpacing(self, spacing):
        self.layout.insertSpacing(len(self.widgets), spacing)
        self.widgets.append(None)

    def addStretch(self, stretch):
        self.layout.insertStretch(len(self.widgets), stretch)
        self.widgets.append(None)

    def addLabel(self, labelString):
        label = qt.QLabel(labelString, self)
        self.addWidget(label, 0, qcore.Qt.AlignLeft)
        
    def setHelp(self, *functionSeq):
        helpString = ""
        functionList = list(*functionSeq)
        while len(functionList) > 0:
            function = functionList.pop()
            if (len(functionList) == 0) and (function.__doc__):
                helpString = helpString + function.__doc__
            elif function.__doc__:
                helpString = helpString + function.__doc__ + os.linesep + \
                    "--------------------------------------------------------"+\
                    "--------------------------------" + os.linesep
        
        if not hasattr(self.buttons, "Help"):
            self.buttons.Help = qt.QPushButton("Help", self.buttons)
            self.buttons.Help.setToggleButton(1)
            self.buttons.layout.insertWidget(3, self.buttons.Help)
            self.connect(self.buttons.Help, SIGNAL("toggled(bool)"), self.showExtension)
        
        if int(qt.qVersion()[0]) < 3:
            self.help = qt.QMultiLineEdit(self)
            self.help.setText(helpString)
            if self.help.numLines() > 20:
                self.help.setFixedVisibleLines(20)
            else:
                self.help.setFixedVisibleLines(self.help.numLines()+1)

            self.help.setReadOnly(1)
            self.help.setWordWrap(qt.QMultiLineEdit.WidgetWidth)
        else:
            self.help = qt.QVBox(self)
            self.help.text = qcore.QtextEdit(self.help)
            self.help.text.setText(helpString)
            self.help.text.setReadOnly(1)
            self.help.text.setWordWrap(qt.QMultiLineEdit.WidgetWidth)
            total_height = self.help.text.heightForWidth(self.help.width())
            if  total_height > self.help.text.height():
                self.help.text.setMinimumSize(self.help.text.width(), min(300, total_height))
                
        self.setExtension(self.help)
        
    def tryAccept(self):
        for i in self.widgets:
            try:
                if i.validate() == 0:
                    return
            except AttributeError:
                continue
        self.accept()
