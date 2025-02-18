# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GCview.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


import sys,os
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import QWebEngineView
from GCviewWeb import GCviewWeb_Form


class GCview_Form(QWidget):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(551, 407)
        Form.setWindowIcon(QIcon("./logo/logo.ico"))
        self.gridLayout_2 = QtWidgets.QGridLayout(Form)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.textBrowser_2 = QtWidgets.QTextBrowser(Form)
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.gridLayout.addWidget(self.textBrowser_2, 2, 0, 1, 1)
        self.label = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout.addWidget(self.pushButton_2, 3, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.pushButton = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 4, 0, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

        # action
        self.pushButton.clicked.connect(self.web_open)
        self.pushButton_2.clicked.connect(self.read_file1)

        ## default
        self.textBrowser_2.setPlaceholderText("D:/input/test.json")

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "GCview"))
        self.label.setText(_translate("Form", "GCview"))
        self.pushButton_2.setText(_translate("Form", "Choose"))
        self.label_2.setText(_translate("Form", "Input JSON file"))
        self.pushButton.setText(_translate("Form", "Run"))


    def read_file1(self):
        openfile_name = QtWidgets.QFileDialog.getOpenFileName(self, 'choose file', '')[0]
        # print(openfile_name)
        self.textBrowser_2.setText(openfile_name)

    def web_open(self):
        # 绑定按钮点击事件，发射信号
        content = self.textBrowser_2.toPlainText()
        content = content.strip()
        print(content)
        if content == '':
            w = QWidget()
            QMessageBox.critical(w, "error",
                                 "Please add correct path!")
        else:
            self.winTable = GCviewWeb_Form(content)
            self.winTable.show()



if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = GCview_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
