# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Phylotreejs.ui'
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



class Phylotreejs_Form(QWidget):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(551, 407)
        Form.setWindowIcon(QIcon("./logo/logo.ico"))
        self.gridLayout_2 = QtWidgets.QGridLayout(Form)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setVerticalSpacing(0)
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
        self.textBrowser_2.setPlaceholderText("D:/input/test.nwk")

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Phylotree"))
        self.label.setText(_translate("Form", "Phylotree"))
        self.pushButton_2.setText(_translate("Form", "Choose"))
        self.label_2.setText(_translate("Form", "Input Tree file (nwk)"))
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
            self.winTable = MainWindow(content)
            self.winTable.show()


class MainWindow(QMainWindow):
    def __init__(self, content):
        super().__init__()
        self.setWindowTitle('Phylotree')
        self.setWindowIcon(QIcon("./logo/logo.ico"))
        self.showMaximized()

        self.path_1 = os.path.abspath('.')
        self.path_1 = self.path_1.strip().split('\\')
        self.path = '/'.join(self.path_1)

        #####放入WebEngineView
        self.webview = WebEngineView()
        self.webview.load(QUrl(self.path + "/html/phylotree/index.html"))
        self.setCentralWidget(self.webview)

        #####web页面加载完毕，调用函数
        self.webview.page().loadFinished.connect(lambda: self.run_js(content)) # 信号函数传参

    def run_js(self, content):
        js_string_1 = '''
        $(document).ready(function() {
            var file_name = "'''

        js_string_2 = '''";
            d3.text(file_name, function(error, newick) {
              var height = 550,
                width = 290,
                tree = d3.layout.phylotree()
                  .svg(d3.select("#tree_display"))
                  .options({
                    'left-right-spacing': 'fit-to-step',
                    'top-bottom-spacing': 'fit-to-size',
                    'selectable': true,
                    'reroot': false,
                    'hide': false,
                    'show-scale': false
                  })
                  .align_tips(true)
                  .font_size(14)
                  .size([height, width])
                  .node_circle_size(0);
              tree(newick).layout();
              function my_node_style_text(node) {
                node['text-italic'] = !node['text-italic'];
                d3.layout.phylotree.trigger_refresh(tree);
              }
              function my_menu_title(node) {
                if (node['text-italic']) {
                  return "Remove Italics";
                }
                return "Custom function";
              }
              tree.get_nodes()
                .filter(d3.layout.phylotree.is_leafnode)
                .forEach(function(tree_node) {
                  d3.layout.phylotree.add_custom_menu(tree_node, // add to this node
                    my_menu_title, // display this text for the menu
                    function() { //Custom function
                      console.log(tree_node['name'])
                    },
                    d3.layout.phylotree.is_leafnode
                  );
                });
            })
          })
        '''

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
    Form = QtWidgets.QWidget()
    ui = Phylotreejs_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())