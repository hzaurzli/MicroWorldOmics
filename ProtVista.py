# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ProtVista.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


import sys,os,shutil
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import QWebEngineView


class winTest(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('My Browser')
        self.setStyleSheet("background-image: url(D:/Documents/Desktop/bb.png)")

    """对QDialog类重写，实现一些功能"""

    def closeEvent(self, event):
        """
        重写closeEvent方法，实现dialog窗体关闭时执行一些代码
        :param event: close()触发的事件
        :return: None
        """
        try:
            if os.path.exists(os.path.dirname(content_final) + '/external_P05067.json') == True:
                print(content_final)
                os.remove(os.path.dirname(content_final) + '/external_P05067.json')
            else:
                event.ignore()  # 设置正常退出
        except:
            return None  # 设置正常退出


class Provista_Form(QWidget):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(691, 431)
        Form.setStyleSheet("background-image: url(D:/Documents/Desktop/bb.png)")
        self.textBrowser_2 = QtWidgets.QTextBrowser(Form)
        self.textBrowser_2.setGeometry(QtCore.QRect(120, 200, 341, 31))
        self.textBrowser_2.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(280, 30, 141, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(120, 160, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.pushButton_2 = QtWidgets.QPushButton(Form)
        self.pushButton_2.setGeometry(QtCore.QRect(500, 200, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton = QtWidgets.QPushButton(Form)
        self.pushButton.setGeometry(QtCore.QRect(300, 300, 81, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

        # action
        self.pushButton.clicked.connect(self.web_open)
        self.pushButton_2.clicked.connect(self.read_file1)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "ProtVista"))
        self.label.setText(_translate("Form", "ProtVista"))
        self.label_2.setText(_translate("Form", "Input JSON file"))
        self.pushButton_2.setText(_translate("Form", "Choose"))
        self.pushButton.setText(_translate("Form", "Run"))


    def read_file1(self):
        openfile_name = QtWidgets.QFileDialog.getOpenFileName(self, 'choose file', '')[0]
        print(openfile_name)
        self.textBrowser_2.setText(openfile_name)


    def web_open(self):
        # 绑定按钮点击事件，发射信号
        content = self.textBrowser_2.toPlainText()
        content = content.strip()
        content_file = os.path.basename(content)

        bool_1 = content_file.endswith(".json")
        bool_2 = content_file.endswith(".JSON")
        bool_3 = content_file.endswith(".Json")
        print(bool_1 or bool_2 or bool_3 == True)

        if bool_1 or bool_2 or bool_3 == True:
            path_1 = os.path.abspath('.')
            path_1 = path_1.strip().split('\\')
            path = '/'.join(path_1) + "/html/protvista/"

            shutil.copy(content, path)
            global content_final
            if os.path.exists(path + 'external_P05067.json'):
                os.remove(path + 'external_P05067.json')
                os.rename(path + content_file, path + 'external_P05067.json')
                content_final = path + 'external_P05067.json'
            else:
                os.rename(path + content_file, path + 'external_P05067.json')
                content_final = path + 'external_P05067.json'

            print(content_final)
            self.winTable = MainWindow(content_final)
            self.winTable.show()
        else:
            QMessageBox.critical(self, "error", "Please check format(suffix:json or JSON)!")


class MainWindow(QMainWindow):
    def __init__(self, content):
        super().__init__()
        self.setWindowTitle('My Browser')
        self.showMaximized()

        self.path_1 = os.path.abspath('.')
        self.path_1 = self.path_1.strip().split('\\')
        self.path = '/'.join(self.path_1)

        #####放入WebEngineView
        self.webview = WebEngineView()
        self.webview.load(QUrl(self.path + "/html/protvista/ProtVista.html"))
        self.setCentralWidget(self.webview)

        #####web页面加载完毕，调用函数
        self.webview.page().loadFinished.connect(lambda: self.run_js(content))

    def run_js(self, content):
        js_string_1 = '''
            var yourDiv = document.getElementById('yourDiv');
            var ProtVista = require('ProtVista');
            var instance = new ProtVista({
            el: yourDiv,
            uniprotacc: "P05067", 
            defaultSources: true,
            customDataSource: {
                 url: "'''

        js_string_2 = '''",
                source: 'myLab',
                useExtension: true
                }
            });
        '''
        content_path = os.path.dirname(content)
        content = content_path + '/external_'

        js_string = js_string_1 + content + js_string_2
        print(js_string)
        self.webview.page().runJavaScript(js_string)


class WebEngineView(QWebEngineView):
    windowList = []

    # 重写createwindow()
    def createWindow(self, QWebEnginePage_WebWindowType):
        new_webview = WebEngineView()
        new_window = MainWindow()
        new_window.setCentralWidget(new_webview)
        # new_window.show()
        self.windowList.append(new_window)  # 注：没有这句会崩溃
        return new_webview


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    # Form = QtWidgets.QWidget()
    WT = winTest()
    ui = Provista_Form()
    ui.setupUi(WT)
    WT.show()
    sys.exit(app.exec_())