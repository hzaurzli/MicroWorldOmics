# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Muscle.ui'
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


class WorkThread(QThread):
    # 自定义信号对象
    trigger = pyqtSignal(str)

    def __int__(self):
        # 初始化函数
        super(WorkThread, self).__init__()

    def run(self):
        def check_process_running(process_name):  # 检查进程是否运行
            for process in psutil.process_iter(['name']):
                if process.info['name'] == process_name:
                    return True
            return False

        try:
            path = os.path.abspath('.')
            if '\\' in path:
                path = path.strip().split('\\')
                path = '/'.join(path)

            os.popen(r"muscle -align %s -output %s"
                     % (fasta, out))

            process_name = 'muscle'
            time.sleep(3)
            while True:  # 判断 iqtree 是否运行完成
                if check_process_running(process_name):
                    print(f"The process {process_name} is running.")
                    time.sleep(10)
                    continue
                else:
                    print(f"The process {process_name} is not running.")
                    break

            self.trigger.emit('Finished!!!')

        except Exception as ex:
            self.trigger.emit('Some errors have occurred, %s!' % ex)


class Muscle_Form(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.work = WorkThread()

    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(609, 400)
        Form.setWindowIcon(QIcon("./logo/logo.ico"))
        self.gridLayout_4 = QtWidgets.QGridLayout(Form)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_4 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label_4.setFont(font)
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 0, 0, 1, 1)
        self.textBrowser = QtWidgets.QTextBrowser(Form)
        self.textBrowser.setObjectName("textBrowser")
        self.gridLayout_3.addWidget(self.textBrowser, 1, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout_3, 1, 1, 1, 1)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout, 0, 0, 1, 2)
        self.gridLayout_5 = QtWidgets.QGridLayout()
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label_3 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName("label_3")
        self.gridLayout_5.addWidget(self.label_3, 3, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.gridLayout_5.addWidget(self.label_2, 0, 0, 1, 1)
        self.textBrowser_2 = QtWidgets.QTextBrowser(Form)
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.gridLayout_5.addWidget(self.textBrowser_2, 1, 0, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout_5.addWidget(self.pushButton_2, 2, 0, 1, 1)
        self.textBrowser_3 = QtWidgets.QTextBrowser(Form)
        self.textBrowser_3.setObjectName("textBrowser_3")
        self.gridLayout_5.addWidget(self.textBrowser_3, 4, 0, 1, 1)
        self.pushButton_3 = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.gridLayout_5.addWidget(self.pushButton_3, 5, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout_5, 1, 0, 1, 1)
        self.gridLayout_6 = QtWidgets.QGridLayout()
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.pushButton = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(22)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout_6.addWidget(self.pushButton, 0, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout_6, 2, 1, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

        # button action
        self.pushButton.clicked.connect(self.calculation)
        self.pushButton_2.clicked.connect(self.read_file1)
        self.pushButton_3.clicked.connect(self.read_file2)

        ## default
        self.textBrowser_2.setPlaceholderText("D:/input/test.fa")
        self.textBrowser_3.setPlaceholderText("D:/output/test.aln")

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Muscle"))
        self.label_4.setText(_translate("Form", "Status"))
        self.label.setText(_translate("Form", "Muscle"))
        self.label_3.setText(_translate("Form", "Output fasta file"))
        self.label_2.setText(_translate("Form", "Input fasta file"))
        self.pushButton_2.setText(_translate("Form", "Choose"))
        self.pushButton_3.setText(_translate("Form", "Choose"))
        self.pushButton.setText(_translate("Form", "Run"))


    def read_file1(self):
        openfile_name = QtWidgets.QFileDialog.getOpenFileName(self, 'choose file', '')[0]
        print(openfile_name)
        self.textBrowser_2.setText(openfile_name)

    def read_file2(self):
        openfile_name = QtWidgets.QFileDialog.getSaveFileName(self, "choose file", "./")[0]
        print(openfile_name)
        self.textBrowser_3.setText(openfile_name)

    def finished(self, str):
        self.textBrowser.setText(str)

    def calculation(self):
        try:
            global fasta, out, path
            fasta = self.textBrowser_2.toPlainText()
            out = self.textBrowser_3.toPlainText()
            path = os.path.dirname(out)

            def is_fasta(filename):
                with open(filename, "r") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    return any(fasta)

            if 0 in [len(fasta), len(out)]:
                QMessageBox.warning(self, "warning", "Please add correct file path!", QMessageBox.Cancel)
            else:
                if is_fasta(fasta) == False:
                    QMessageBox.critical(self, "error", "Check fasta file format!")
                else:
                    self.textBrowser.setText('Running! please wait')
                    QApplication.processEvents()  # 逐条打印状态

                    # 启动线程, 运行 run 函数
                    self.work.start()
                    # 传送信号, 接受 run 函数执行完毕后的信号
                    self.work.trigger.connect(self.finished)

        except:
            QMessageBox.critical(self, "error", "Check fasta file format!")


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Muscle_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
