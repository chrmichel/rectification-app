from PyQt5.QtWidgets import *
from mywidgets import InputWindow, Done, Verify, Output
from controller import Controller
import matplotlib.pyplot as plt


class GUI(QMainWindow):

    def __init__(self, parent=None):
        
        super().__init__(parent)
        self.setWindowTitle("Rectification app")
        self.ctrl = Controller(self)
        central = QWidget()
        self.setCentralWidget(central)
        self.iw = None
        self.ow = None
        self.grid = QGridLayout()
        self.d1 = Done()
        self.d2 = Done()

        self.btn1 = QPushButton("Get input")
        self.grid.addWidget(QLabel("Get input parameters:"), 0, 0)
        self.grid.addWidget(self.btn1, 0, 1)
        self.btn2 = QPushButton("Calc")
        self.btn2.setDisabled(True)
        self.grid.addWidget(QLabel("Run calculations:"), 1, 0)
        self.grid.addWidget(self.btn2, 1, 1)
        self.btn3 = QPushButton("Results")
        self.btn3.setDisabled(True)
        self.grid.addWidget(QLabel("Show results:"), 2, 0)
        self.grid.addWidget(self.btn3, 2, 1)
        self.btn4 = QPushButton("Use .ini")
        self.grid.addWidget(self.btn4, 0, 2)

        self.reset_btn = QPushButton("Reset")
        self.grid.addWidget(self.reset_btn, 6, 0)

        central.setLayout(self.grid)
        self._connect()
        
        self.ow = Output()
    
    def get_params(self) -> dict:
        iw = InputWindow(self)
        iw.show()
        if iw.exec_():
            self.grid.addWidget(self.d1, 0, 3)
            self.btn1.setDisabled(True)
            self.btn4.setDisabled(True)
            self.btn2.setEnabled(True)
            return iw.results
    
    def use_ini(self):
        self.grid.addWidget(self.d1, 0, 3)
        self.btn1.setDisabled(True)
        self.btn4.setDisabled(True)
        self.btn2.setEnabled(True)
    
    def calc(self, params):
        self.btn2.setDisabled(True)
        verify = Verify(self, params)
        if verify.exec_():
            self.grid.addWidget(self.d2, 1, 3)
            self.btn3.setEnabled(True)
        else:
            self.reset()

    def _connect(self):
        self.btn1.clicked.connect(self.ctrl._get_params)
        self.btn2.clicked.connect(self.ctrl._calc)
        self.btn4.clicked.connect(self.ctrl._read_config)
        self.reset_btn.clicked.connect(self.ctrl._reset)
        self.btn3.clicked.connect(self.ctrl._output)
    
    def reset(self):
        self.btn2.setDisabled(True)
        self.btn3.setDisabled(True)
        self.btn1.setEnabled(True)
        self.btn4.setEnabled(True)
        self.grid.removeWidget(self.d1)
        self.grid.removeWidget(self.d2)
        self.d1.setParent(None)
        self.d2.setParent(None)
        self.ow.temp_tab.axes.clear()
        self.ow.x_tab.axes.clear()
        self.ow.y_tab.axes.clear()

    def output(self):
        nrange = range(self.ctrl.rect.n)

        tplot = self.ow.temp_tab.axes
        tplot.plot(nrange, self.ctrl.rect.T_j)
        tplot.set_ylabel("Temperature / K")
        tplot.set_xlabel("Plate number")
        # tplot.invert_yaxis()
        
        comps = self.ctrl.rect.comps

        yplot = self.ow.y_tab.axes
        ylines = {}
        for y, c in zip(self.ctrl.rect.y_jicor.T, comps):
            line, = yplot.plot(nrange, y, label=c)
            ylines[c] = line
        
        yplot.set_ylabel("Vapor composition")
        yplot.set_xlabel("Plate number")
        # yplot.invert_yaxis()
        yplot.grid(True)
        yplot.legend()

        xplot = self.ow.x_tab.axes
        for x, c in zip(self.ctrl.rect.x_jicor.T, comps):
            xplot.plot(nrange, x, label=c)
        
        xplot.set_ylabel("Liquid composition")
        xplot.set_xlabel("Plate number")
        # xplot.invert_yaxis()
        xplot.grid(True)
        xplot.legend()


        self.ow.show()

if __name__ == '__main__':
    app = QApplication([])
    main = GUI()
    main.show()
    app.exec()