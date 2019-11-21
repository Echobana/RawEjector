#! python3
# -*- coding: utf-8 -*-

# -------------------- IMPORTANT LAUNCH INFO --------------------------
# Spyder (Anaconda v.2019-07) does not support 'QtWidgets.qApp.quit()'
# If you use Anaconda3 v.2019-07 comment out the line 95
# If you use PyCharm v2019.2.2 use raw code (line 95 is not commented)
# Data line edit is used for data file (*.*)
# Sensor info line edit is used for sensor information (*.*)
# The requirements for sensor_info.xlsx see in the Handler header
# IDE versions're given taking into account the time of writing the code
# ----------------------------------------------------------------------

# Return file result.pdf to code launch directory
# E.S.Tsyrendorzhiev E1-92 2019


from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon
import rawsHandler as hdl


class MyWindow(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowIcon(QIcon(r'F:\Загрузки\ejector_project\ejector_raw\data.ico'))

        self.dataLabel = QtWidgets.QLabel('Data: ')
        self.sensorLabel = QtWidgets.QLabel('Sensor info: ')
        self.parametersLabel = QtWidgets.QLabel('Parameters info: ')

        self.dataEdit = QtWidgets.QLineEdit(r"F:\ejector_raw_files\2")
        self.sensorEdit = QtWidgets.QLineEdit(r"F:\ejector_raw_files\sensor_status")
        self.parametersEdit = QtWidgets.QLineEdit(r"F:\Загрузки\ejector_project\ejector_raw_files\parameters")

        self.dataBrowse = QtWidgets.QPushButton('Browse', self)
        self.dataBrowse.clicked.connect(self.btnDataClicked)

        self.sensorBrowse = QtWidgets.QPushButton('Browse', self)
        self.sensorBrowse.clicked.connect(self.btnSensorClicked)

        self.parametersBrowse = QtWidgets.QPushButton('Browse', self)
        self.parametersBrowse.clicked.connect(self.btnParametersClicked)

        self.solveBrowse = QtWidgets.QPushButton('Plot', self)
        self.solveBrowse.clicked.connect(self.btnSolveClicked)

        self.throttle_btn = QtWidgets.QPushButton('Throttle')
        self.massflow_btn = QtWidgets.QPushButton('Mass Flow')
        self.massflow_btn.setCheckable(True)
        self.throttle_btn.setCheckable(True)
        # self.throttle_btn.clicked[bool].connect(self.setType)
        # self.massflow_btn.clicked[bool].connect(self.setType)

        self.vbox_main = QtWidgets.QVBoxLayout()
        self.hbox_main = QtWidgets.QHBoxLayout()

        self.vbox_mid = QtWidgets.QVBoxLayout()
        self.vbox_mid.addWidget(self.dataLabel)
        self.vbox_mid.addWidget(self.sensorLabel)
        self.vbox_mid.addWidget(self.parametersLabel)
        # self.vbox_mid.addWidget(self.checkboxLabel)

        self.vbox_left = QtWidgets.QVBoxLayout()
        self.vbox_left.addWidget(self.dataEdit)
        self.vbox_left.addWidget(self.sensorEdit)
        self.vbox_left.addWidget(self.parametersEdit)

        self.vbox_right = QtWidgets.QVBoxLayout()
        self.vbox_right.addWidget(self.dataBrowse)
        self.vbox_right.addWidget(self.sensorBrowse)
        self.vbox_right.addWidget(self.parametersBrowse)

        self.hbox_main.addLayout(self.vbox_mid)
        self.hbox_main.addLayout(self.vbox_left)
        self.hbox_main.addLayout(self.vbox_right)
        self.vbox_main.addLayout(self.hbox_main)
        self.vbox_main.addWidget(self.solveBrowse)

        self.setLayout(self.vbox_main)

    def btnDataClicked(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self)[0]
        self.dataEdit.setText(fname)

    def btnSensorClicked(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self)[0]
        self.sensorEdit.setText(fname)

    def btnParametersClicked(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self)[0]
        self.dataEdit.setText(fname)

    def setType(self, pressed):
        source = self.sender()

        if pressed:
            val = 255
        else:
            val = 0

        if source.text() == 'Throttle':
            self.col.setThrottle(val)
        elif source.text() == 'Mass Flow':
            self.col.setMassFlow(val)

    def btnSolveClicked(self):
        data = hdl.opener(self.dataEdit.text())
        sensor_info = hdl.opener(self.sensorEdit.text())
        parameters_file = self.parametersEdit.text()

        params = hdl.set_experiment_parameters(parameters_file)
        ver_data = hdl.verification(sensor_info.T, data)
        data_dict = hdl.av(ver_data)
        p_or, p_01, p_04, p_02, ej_coeff_list, comp_ratio_list, eff_list, m1_list, m2_list = hdl.solver(data_dict, params)
        trans = hdl.transient(ver_data)
        p_or_tns, p_01_tns, p_04_tns, p_02_tns, ej_coeff_list_tns, comp_ratio_list_tns, eff_list_tns, m1_list_tns, m2_list_tns = hdl.solver(
            trans, params)

        self.graph = {
            (r'Коэффициент эжекции, k', r'$p_{02}$, кПа', 'Title_1'): (
                ej_coeff_list, p_02, ej_coeff_list_tns, p_02_tns),
            (r'$p_{04}$, кПа', r'$p_{02}$, кПа', 'Title_2'): (p_04, p_02, p_04_tns, p_02_tns),
            (r'$p_{04}$, кПа', r'$\varepsilon$', 'Title_3'): (p_04, comp_ratio_list, p_04_tns, comp_ratio_list_tns),
            (r'$p_{04}$, кПа', r'$\eta$', 'Title_4'): (p_04, eff_list, p_04_tns, eff_list_tns)
        }  # Оси и легенды графиков 1-4

        # maximize = ((r'$p_{04}$, МПа', r'$\varepsilon$', 'Title_3'), (r'$p_{04}$, МПа', r'$\eta$', 'Title_4'))
        # fig0, ax0 = plt.subplots()
        # ax0.scatter(0, 1)
        # plt.show()
        hdl.mult_plot(self.graph, sensor_info, data_dict)
        hdl.pdf_saver(r'./result.pdf')

        QtWidgets.qApp.quit()


if __name__ == '__main__':
    import sys

    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.setWindowTitle("Program")
    window.setFixedSize(350, 150)
    window.show()
    sys.exit(app.exec_())
