# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Prodigal.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


import sys,os
import time
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from Bio import SeqIO
import psutil


class Prodigal_Form(QWidget):
    def setupUi(self, Clustal):
        Clustal.setObjectName("Clustal")
        Clustal.resize(702, 436)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Clustal.sizePolicy().hasHeightForWidth())
        Clustal.setSizePolicy(sizePolicy)
        Clustal.setStyleSheet("background-image: url(D:/Documents/Desktop/bb.png)")
        self.label = QtWidgets.QLabel(Clustal)
        self.label.setGeometry(QtCore.QRect(280, 10, 121, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Clustal)
        self.label_2.setGeometry(QtCore.QRect(30, 40, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.textBrowser_2 = QtWidgets.QTextBrowser(Clustal)
        self.textBrowser_2.setGeometry(QtCore.QRect(30, 80, 221, 31))
        self.textBrowser_2.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.label_3 = QtWidgets.QLabel(Clustal)
        self.label_3.setGeometry(QtCore.QRect(30, 140, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName("label_3")
        self.textBrowser_3 = QtWidgets.QTextBrowser(Clustal)
        self.textBrowser_3.setGeometry(QtCore.QRect(30, 180, 221, 31))
        self.textBrowser_3.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser_3.setObjectName("textBrowser_3")
        self.pushButton_2 = QtWidgets.QPushButton(Clustal)
        self.pushButton_2.setGeometry(QtCore.QRect(280, 80, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_3 = QtWidgets.QPushButton(Clustal)
        self.pushButton_3.setGeometry(QtCore.QRect(280, 180, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.textBrowser = QtWidgets.QTextBrowser(Clustal)
        self.textBrowser.setGeometry(QtCore.QRect(30, 340, 221, 71))
        self.textBrowser.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser.setObjectName("textBrowser")
        self.label_4 = QtWidgets.QLabel(Clustal)
        self.label_4.setGeometry(QtCore.QRect(60, 300, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label_4.setFont(font)
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.pushButton = QtWidgets.QPushButton(Clustal)
        self.pushButton.setGeometry(QtCore.QRect(580, 370, 81, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(22)
        self.pushButton.setFont(font)
        self.pushButton.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.pushButton.setObjectName("pushButton")
        self.verticalGroupBox = QtWidgets.QGroupBox(Clustal)
        self.verticalGroupBox.setGeometry(QtCore.QRect(420, 80, 271, 111))
        self.verticalGroupBox.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.verticalGroupBox.setObjectName("verticalGroupBox")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalGroupBox)
        self.verticalLayout.setObjectName("verticalLayout")
        self.radioButton = QtWidgets.QRadioButton(self.verticalGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton.sizePolicy().hasHeightForWidth())
        self.radioButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.radioButton.setFont(font)
        self.radioButton.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.radioButton.setAutoExclusive(True)
        self.radioButton.setObjectName("radioButton")
        self.verticalLayout.addWidget(self.radioButton)
        self.radioButton_2 = QtWidgets.QRadioButton(self.verticalGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton_2.sizePolicy().hasHeightForWidth())
        self.radioButton_2.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.radioButton_2.setFont(font)
        self.radioButton_2.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.radioButton_2.setAutoExclusive(True)
        self.radioButton_2.setObjectName("radioButton_2")
        self.verticalLayout.addWidget(self.radioButton_2)
        self.verticalGroupBox_2 = QtWidgets.QGroupBox(Clustal)
        self.verticalGroupBox_2.setGeometry(QtCore.QRect(420, 249, 271, 111))
        self.verticalGroupBox_2.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.verticalGroupBox_2.setObjectName("verticalGroupBox_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalGroupBox_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.radioButton_3 = QtWidgets.QRadioButton(self.verticalGroupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton_3.sizePolicy().hasHeightForWidth())
        self.radioButton_3.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.radioButton_3.setFont(font)
        self.radioButton_3.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.radioButton_3.setAutoExclusive(True)
        self.radioButton_3.setObjectName("radioButton_3")
        self.verticalLayout_2.addWidget(self.radioButton_3)
        self.radioButton_4 = QtWidgets.QRadioButton(self.verticalGroupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton_4.sizePolicy().hasHeightForWidth())
        self.radioButton_4.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.radioButton_4.setFont(font)
        self.radioButton_4.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.radioButton_4.setAutoExclusive(True)
        self.radioButton_4.setObjectName("radioButton_4")
        self.verticalLayout_2.addWidget(self.radioButton_4)
        self.radioButton_5 = QtWidgets.QRadioButton(self.verticalGroupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton_5.sizePolicy().hasHeightForWidth())
        self.radioButton_5.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.radioButton_5.setFont(font)
        self.radioButton_5.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.radioButton_5.setAutoExclusive(True)
        self.radioButton_5.setObjectName("radioButton_5")
        self.verticalLayout_2.addWidget(self.radioButton_5)
        self.label_9 = QtWidgets.QLabel(Clustal)
        self.label_9.setGeometry(QtCore.QRect(30, 230, 161, 22))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_9.setFont(font)
        self.label_9.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_9.setObjectName("label_9")
        self.textEdit = QtWidgets.QTextEdit(Clustal)
        self.textEdit.setGeometry(QtCore.QRect(30, 260, 221, 31))
        self.textEdit.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textEdit.setObjectName("textEdit")
        self.label_7 = QtWidgets.QLabel(Clustal)
        self.label_7.setGeometry(QtCore.QRect(420, 40, 171, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_7.setFont(font)
        self.label_7.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(Clustal)
        self.label_8.setGeometry(QtCore.QRect(420, 210, 181, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_8.sizePolicy().hasHeightForWidth())
        self.label_8.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_8.setFont(font)
        self.label_8.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_8.setObjectName("label_8")

        self.retranslateUi(Clustal)
        QtCore.QMetaObject.connectSlotsByName(Clustal)

        # button action
        self.pushButton.clicked.connect(self.build_tree)
        self.pushButton_2.clicked.connect(self.read_file1)
        self.pushButton_3.clicked.connect(self.read_file2)

        self.textEdit.setPlaceholderText(" Coden: 11")
        self.radioButton.setChecked(True)
        self.radioButton_3.setChecked(True)


    def retranslateUi(self, Clustal):
        _translate = QtCore.QCoreApplication.translate
        Clustal.setWindowTitle(_translate("Clustal", "Prodigal"))
        self.label.setText(_translate("Clustal", "Prodigal"))
        self.label_2.setText(_translate("Clustal", "Input fasta file"))
        self.label_3.setText(_translate("Clustal", "Output fasta file"))
        self.pushButton_2.setText(_translate("Clustal", "Choose"))
        self.pushButton_3.setText(_translate("Clustal", "Choose"))
        self.label_4.setText(_translate("Clustal", "Status"))
        self.pushButton.setText(_translate("Clustal", "Run"))
        self.radioButton.setText(_translate("Clustal", "Single"))
        self.radioButton_2.setText(_translate("Clustal", "Meta"))
        self.radioButton_3.setText(_translate("Clustal", "GFF"))
        self.radioButton_4.setText(_translate("Clustal", "GBK"))
        self.radioButton_5.setText(_translate("Clustal", "SCO"))
        self.label_9.setText(_translate("Clustal", "Translation tables"))
        self.label_7.setText(_translate("Clustal", "Single / Meta"))
        self.label_8.setText(_translate("Clustal", "Output file format"))


    def read_file1(self):
        openfile_name = QtWidgets.QFileDialog.getOpenFileName(self, 'choose file', '')[0]
        print(openfile_name)
        self.textBrowser_2.setText(openfile_name)

    def read_file2(self):
        openfile_name = QtWidgets.QFileDialog.getSaveFileName(self, "choose file", "./")[0]
        print(openfile_name)
        self.textBrowser_3.setText(openfile_name)


    def build_tree(self):
        try:
            fasta = self.textBrowser_2.toPlainText()
            out = self.textBrowser_3.toPlainText()
            path = os.path.dirname(out)

            def check_process_running(process_name): # 检查进程是否运行
                for process in psutil.process_iter(['name']):
                    if process.info['name'] == process_name:
                        return True
                return False

            try:
                coden = self.textEdit.text()
                if isinstance(coden, int) == True:
                    coden = coden
                else:
                    QMessageBox.critical(self, "error", "Coden value error!")
            except:
                coden = 11

            if any([len(fasta), len(out)]) == False:
                QMessageBox.warning(self, "warning", "Please add correct file path!", QMessageBox.Cancel)
            else:
                self.textBrowser.setText('Running! please wait')
                QApplication.processEvents()  # 逐条打印状态

                if self.radioButton.isChecked() and self.radioButton_3.isChecked():
                    os.popen(r".\tools\prodigal\prodigal.exe -f gff -g %s -o %s -p single -i %s"
                             % (coden, out, fasta))
                elif self.radioButton.isChecked() and self.radioButton_4.isChecked():
                    os.popen(r".\tools\prodigal\prodigal.exe -f gbk -g %s -o %s -p single -i %s"
                             % (coden, out, fasta))
                elif self.radioButton.isChecked() and self.radioButton_5.isChecked():
                    os.popen(r".\tools\prodigal\prodigal.exe -f sco -g %s -o %s -p single -i %s"
                             % (coden, out, fasta))
                elif self.radioButton_2.isChecked() and self.radioButton_3.isChecked():
                    os.popen(r".\tools\prodigal\prodigal.exe -f gff -g %s -o %s -p meta -i %s"
                             % (coden, out, fasta))
                elif self.radioButton_2.isChecked() and self.radioButton_4.isChecked():
                    os.popen(r".\tools\prodigal\prodigal.exe -f gbk -g %s -o %s -p meta -i %s"
                             % (coden, out, fasta))
                elif self.radioButton_2.isChecked() and self.radioButton_5.isChecked():
                    os.popen(r".\tools\prodigal\prodigal.exe -f sco -g %s -o %s -p meta -i %s"
                             % (coden, out, fasta))
                else:
                    QMessageBox.critical(self, "error", "Please choose Substitution model!")

                process_name = 'prodigal.exe'
                time.sleep(5)
                while True: # 判断 iqtree.exe 是否运行完成
                    if check_process_running(process_name):
                        print(f"The process {process_name} is running.")
                        time.sleep(30)
                        continue
                    else:
                        print(f"The process {process_name} is not running.")
                        self.textBrowser.setText('Finished!!!')
                        break
        except:
            QMessageBox.critical(self, "error", "Check fasta file format!")

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Clustal = QtWidgets.QWidget()
    ui = Prodigal_Form()
    ui.setupUi(Clustal)
    Clustal.show()
    sys.exit(app.exec_())