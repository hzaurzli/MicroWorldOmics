# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ShinyIGV.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import os, sys, re, math, time
from ShinyWeb import ShinyWeb_Form


# QtWidgets.QWidget 要与 ui 窗口一致 QWidget 对应 QWidget; QMainWindow 对应 QMainWindow
class winTest(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('My Browser')
        self.setStyleSheet("background-image: url(./logo/backgroundpage.png)")

    """对QDialog类重写，实现一些功能"""

    def closeEvent(self, event):
        """
        重写closeEvent方法，实现dialog窗体关闭时执行一些代码
        :param event: close()触发的事件
        :return: None
        """
        try:
            with os.popen('netstat -aon|findstr "55368"') as res:
                res = res.read().split('\n')
            result = []
            for line in res:
                temp = [i for i in line.split(' ') if i != '']
                if len(temp) > 4:
                    result.append({'pid': temp[4], 'address': temp[1], 'state': temp[3]})

            if len(result) != 0:
                for task in result:
                    pid = task['pid']
                    os.popen('taskkill -pid %s -f' % (pid))

        except:
            return None  # 设置正常退出


class WorkThread(QThread):
    # 自定义信号对象
    trigger = pyqtSignal(str)

    def __int__(self):
        # 初始化函数
        super(WorkThread, self).__init__()

    def run(self):
        path = os.path.abspath('.')
        if '\\' in path:
            path = path.strip().split('\\')
            path = '/'.join(path)

        os.popen(path + '/Shiny/R-4.3.2/bin/Rscript ' +
                 path + '/Shiny/ShinyScript/ShinyIGV/ShinyIGV.R')

        time.sleep(25)
        self.trigger.emit('ShinyApp has been started!!!')


class ShinyIGV_Form(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.work = WorkThread()

    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(633, 434)
        Form.setWindowIcon(QIcon("./logo/logo.ico"))
        Form.setStyleSheet("background-image: url(./logo/green_back.png);")
        self.gridLayout_4 = QtWidgets.QGridLayout(Form)
        self.gridLayout_4.setVerticalSpacing(0)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.gridLayout_5 = QtWidgets.QGridLayout()
        self.gridLayout_5.setVerticalSpacing(0)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout_5.addWidget(self.label, 0, 0, 1, 3)
        self.label_4 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label_4.setFont(font)
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.gridLayout_5.addWidget(self.label_4, 3, 0, 1, 3)
        self.textBrowser = QtWidgets.QTextBrowser(Form)
        self.textBrowser.setStyleSheet("background-image: url(./logo/white.png)")
        self.textBrowser.setObjectName("textBrowser")
        self.gridLayout_5.addWidget(self.textBrowser, 5, 0, 2, 3)
        self.label_3 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.gridLayout_5.addWidget(self.label_3, 7, 0, 1, 3)
        self.pushButton = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(25)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout_5.addWidget(self.pushButton, 8, 0, 1, 3)
        self.pushButton_2 = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(25)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout_5.addWidget(self.pushButton_2, 9, 0, 1, 3)
        self.gridLayout_4.addLayout(self.gridLayout_5, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

        # button action
        self.pushButton.clicked.connect(self.start)
        self.pushButton_2.clicked.connect(self.open)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "ShinyIGV"))
        self.label.setText(_translate("Form", "ShinyIGV"))
        self.label_4.setText(_translate("Form", "Status"))
        self.label_3.setText(_translate("Form", "Please close this window before starting another Shiny Apps!!!"))
        self.pushButton.setText(_translate("Form", "Start App"))
        self.pushButton_2.setText(_translate("Form", "Open Web"))


    def start(self):
        self.textBrowser.setText(
            'Starting ShinyApp!!!')
        QApplication.processEvents()  # 逐条打印状态

        # 启动线程, 运行 run 函数
        self.work.start()
        # 传送信号, 接受 run 函数执行完毕后的信号
        self.work.trigger.connect(self.finished)

    def open(self):
        self.winTable = ShinyWeb_Form(port=55368)
        self.winTable.show()


    def finished(self, str):
        self.textBrowser.setText(str)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    WT = QtWidgets.QWidget()
    WT = winTest()
    ui = ShinyIGV_Form()
    ui.setupUi(WT)
    WT.show()
    sys.exit(app.exec_())
