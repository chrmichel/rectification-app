from PyQt5.QtWidgets import *
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

COMP_LIST = ["Methanol", "Ethanol", "n-Propanol", "Water", "Acetone", "Benzene"]
M_MAX = 6
N_MAX = 12


class Verify(QDialog):
    def __init__(self, parent=None, params: dict={"key": ["value", 2]}):
        super().__init__(parent)
        l = QVBoxLayout()
        f = QFormLayout()
        self.setLayout(l)
        l.addWidget(QLabel("Please verify the following parameters to be used:"))
        l.addLayout(f)
        for k, v in params.items():
            f.addRow(QLabel(k), QLabel(str(v)))
        bbox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        bbox.accepted.connect(self.accept)
        bbox.rejected.connect(self.reject)
        l.addWidget(bbox)


class Done(QLabel):
    def __init__(self, parent=None):
        super().__init__("Done", parent)


class InputWindow(QDialog):

    def __init__(self, parent=None, comp_list=COMP_LIST, m=3) -> None:
        super().__init__(parent)
        self.setWindowTitle("Rectification")
        self.grid = QGridLayout()
        self.results = {}
        self.complist = comp_list
        self.m = m
        run_btn = QPushButton("Run")
        run_btn.clicked.connect(self.save_params)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.close)
        self.initGeneral()
        self.initComps()
        self.initFeed()
        self.grid.addWidget(run_btn, 1, 1)
        self.setLayout(self.grid)
    
    def initGeneral(self):
        self.form = QFormLayout()
        self.gengroup = QGroupBox("General", self)
        self.gengroup.setLayout(self.form)
        self.n = 4
        self.p = 1.013 # bar
        self.nbox = QSpinBox(self.gengroup)
        self.nbox.setRange(1, N_MAX)
        self.nbox.setValue(self.n)
        self.mbox = QSpinBox(self.gengroup)
        self.mbox.setRange(2, M_MAX)
        self.mbox.setValue(self.m)
        self.p_field = QLineEdit(str(self.p), self.gengroup)
        self.nbox.valueChanged.connect(self.nChanged)
        self.mbox.valueChanged.connect(self.mChanged)
        self.boil_check = QCheckBox(self.gengroup)
        self.cond_check = QCheckBox(self.gengroup)
        self.cond_check.setChecked(True)
        self.form.addRow("Number of plates: ", self.nbox)
        self.form.addRow("Number of components: ", self.mbox)
        self.form.addRow("Column pressure [bar]: ", self.p_field)
        self.form.addRow("Total reboiler?", self.boil_check)
        self.form.addRow("Total condenser?", self.cond_check)
        self.grid.addWidget(self.gengroup, 0, 0)
    
    def initComps(self):
        self.complay = QGridLayout()
        self.compgroup = QGroupBox("Components", self)
        self.compgroup.setLayout(self.complay)
        self.comps = [QComboBox(self.compgroup) for i in range(M_MAX)]
        self.z_fields = [QLineEdit(self.compgroup) for i in range(M_MAX)]
        self.complay.addWidget(QLabel("Name:", self.compgroup), 0, 1)
        self.complay.addWidget(QLabel("Feed fraction:", self.compgroup), 0, 2)
        for i, (c, z) in enumerate(zip(self.comps, self.z_fields)):
            c.addItem("")
            c.addItems(self.complist)
            c.setCurrentIndex(i+1)
            z.setText(str(round(1./self.m, 3)))
            z.returnPressed.connect(self.adjust_z)
            self.complay.addWidget(QLabel(f"{i+1}: ", self.compgroup), i+1, 0)
            self.complay.addWidget(c, i+1, 1)
            self.complay.addWidget(z, i+1, 2)
        self.z_fields[-1].setDisabled(True)
        self.mChanged()
        self.grid.addWidget(self.compgroup, 0, 1)
    
    def initFeed(self):
        self.feed = QGroupBox("Feed", self)
        feedlay = QFormLayout(self.feed)
        self.feed.setLayout(feedlay)
        self.nfbox = QSpinBox(self.feed)
        self.nfbox.setRange(1, self.n)
        self.nfbox.setValue(3)
        self.tf = 350 # K
        self.tf_field = QLineEdit(str(self.tf), self.feed)
        self.F = 100 # kmol/h
        self.F_field = QLineEdit(str(self.F), self.feed)
        self.D = 50 # kmol/h
        self.D_field = QLineEdit(str(self.D), self.feed)
        self.nu = 1.5
        self.nu_field = QLineEdit(str(self.nu), self.feed)
        self.f = 0
        self.f_field = QLineEdit(str(self.f), self.feed)
        feedlay.addRow("Feed plate: ", self.nfbox)
        feedlay.addRow("Feed temperature [K]: ", self.tf_field)
        feedlay.addRow("Feed flow rate [kmol/h]: ", self.F_field)
        feedlay.addRow("Distillate flow rate [kmol/h]: ", self.D_field)
        feedlay.addRow("Reflux ratio: ", self.nu_field)
        feedlay.addRow("Vapor fraction: ", self.f_field)
        self.grid.addWidget(self.feed, 1, 0)

    def mChanged(self):
        
        self.m = self.mbox.value()

        for i, (c, z) in enumerate(zip(self.comps, self.z_fields)):
            if i<self.m:
                z.setText(str(round(1./self.m, 3)))
                z.setEnabled(True)
                c.setEnabled(True)
            else:
                z.setText("")
                z.setDisabled(True)
                c.setCurrentIndex(0)
                c.setDisabled(True)
        self.z_fields[self.m-1].setDisabled(True)

    def nChanged(self):
        try:
            self.n = int(self.nbox.text())
        except ValueError:
            self.n = 4
            self.nbox.setValue(self.n)
            
        self.nfbox.setRange(1, self.n)

    def adjust_z(self):
        m = int(self.mbox.text())
        z = sum([float(zf.text()) for zf in self.z_fields[:m-1]])
        if 0<=z<=1:
            self.z_fields[m-1].setText(str(round(1-z, 3)))
        else:
            QMessageBox.warning(self, "Invalid composition", "Resetting composition", QMessageBox.Ok)
            for z in self.z_fields[:m]:
                z.setText(str(round(1./m, 3)))
  
    def save_params(self):
        self.results["n"] = int(self.nbox.text())
        m = int(self.mbox.text())
        self.results["z_f"] = [float(zf.text()) for zf in self.z_fields[:m]]
        self.results["comps"] = [str(comp.currentText()) for comp in self.comps[:m]]
        self.results["t_f"] = float(self.tf_field.text())
        self.results["n_f"] = self.nfbox.value()
        self.results["f"] = float(self.F_field.text())
        self.results["d"] = float(self.D_field.text())
        self.results["nu"] = float(self.nu_field.text())
        self.results["eta"] = float(self.f_field.text())
        self.results["boil"] = self.boil_check.isChecked()
        self.results["cond"] = self.cond_check.isChecked()
        self.results["p_tot"] = float(self.p_field.text())
        
        self.accept()


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


class Output(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout()
        tabmain = QTabWidget()
        layout.addWidget(tabmain)
        self.setLayout(layout)
        self.setWindowTitle("Result visualization")

        self.temp_tab = MplCanvas()
        self.y_tab = MplCanvas()
        self.x_tab = MplCanvas()
        tabmain.addTab(self.temp_tab, "Temperature")
        tabmain.addTab(self.y_tab, "Vapor composition")
        tabmain.addTab(self.x_tab, "Liquid composition")


if __name__ == '__main__':
    class MainWindow(QMainWindow):
        def __init__(self, parent=None) -> None:
            super().__init__(parent)
            w = InputWindow(parent=self)
            if w.exec_():
                for x in w.results:
                    print(x, '\t', w.results[x])
    
    app = QApplication([])
    main = MainWindow()
    main.show()
    app.exec()