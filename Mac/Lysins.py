# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Lysins.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


import sys,os
import time
from models.Lysins.Feature import all_feature, readFasta
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from numpy import loadtxt, savetxt
import os,psutil
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import numpy as np
import joblib
from sklearn.linear_model import LogisticRegression as LR
import time,os,sys,re,math
from sklearn.ensemble import RandomForestClassifier as RF



class WorkThread(QThread):
    # 自定义信号对象
    trigger = pyqtSignal(str)

    def __int__(self):
        # 初始化函数
        super(WorkThread, self).__init__()

    def run(self):
        try:
            global out

            def model_calculation(fasta, out, path):
                def get_base_proba(test_features, feature_index):
                    base_feature = []
                    for idx, (k, v) in zip(feature_index, clf_feature_order.items()):
                        features = test_features[:, idx]
                        for j in v:
                            model = eval(j)
                            base_proba = model.predict_proba(features)[:, -1]
                            base_feature.append(base_proba)
                    return np.array(base_feature).T

                def meta_model(X, y, path):
                    meta_clf = LR(solver='liblinear', random_state=Randon_seed)
                    meta_clf.fit(X, y)
                    joblib.dump(meta_clf, path)

                def meta_pred(fastafile, path_meta, path):
                    seqs = readFasta(fastafile)
                    test_full_features, feature_index = all_feature(seqs, path)
                    base_feature = get_base_proba(test_full_features, feature_index)
                    meta_clf = joblib.load(path_meta)
                    result = meta_clf.predict_proba(base_feature)
                    return result

                clf_feature_order = {
                    "AAC": ["AAC_ERT", "AAC_ANN", "AAC_XGB", "AAC_KNN", "AAC_LR"],
                    "BPNC": ["BPNC_ERT", "BPNC_ANN", "BPNC_XGB", "BPNC_KNN", "BPNC_LR"],
                    "CTD": ["CTD_ERT", "CTD_ANN", "CTD_XGB", "CTD_KNN", "CTD_LR"],
                    "AAE": ["AAE_ERT", "AAE_ANN", "AAE_XGB", "AAE_KNN", "AAE_LR"],
                    "AAI": ["AAI_ERT", "AAI_ANN", "AAI_XGB", "AAI_KNN", "AAI_LR"],
                    "GAAC": ["GAAC_ERT", "GAAC_ANN", "GAAC_XGB", "GAAC_KNN", "GAAC_LR"],
                }

                test_result = meta_pred(fasta, path + '/meta/Meta.m', path)[:, -1]
                np.savetxt(out, test_result, fmt='%.4f', delimiter=',')

            global AAC_ERT, AAC_ANN, AAC_XGB, AAC_KNN, AAC_LR
            global BPNC_ERT, BPNC_ANN, BPNC_XGB, BPNC_KNN, BPNC_LR
            global CTD_ERT, CTD_ANN, CTD_XGB, CTD_KNN, CTD_LR
            global AAE_ERT, AAE_ANN, AAE_XGB, AAE_KNN, AAE_LR
            global AAI_ERT, AAI_ANN, AAI_XGB, AAI_KNN, AAI_LR
            global GAAC_ERT, GAAC_ANN, GAAC_XGB, GAAC_KNN, GAAC_LR

            AAC_ERT = joblib.load(path + '/base/AAC_ERT.m')
            AAC_ANN = joblib.load(path + '/base/AAC_ANN.m')
            AAC_XGB = joblib.load(path + '/base/AAC_XGB.m')
            AAC_KNN = joblib.load(path + '/base/AAC_KNN.m')
            AAC_LR = joblib.load(path + '/base/AAC_LR.m')

            BPNC_ERT = joblib.load(path + '/base/BPNC_ERT.m')
            BPNC_ANN = joblib.load(path + '/base/BPNC_ANN.m')
            BPNC_XGB = joblib.load(path + '/base/BPNC_XGB.m')
            BPNC_KNN = joblib.load(path + '/base/BPNC_KNN.m')
            BPNC_LR = joblib.load(path + '/base/BPNC_LR.m')

            CTD_ERT = joblib.load(path + '/base/CTD_ERT.m')
            CTD_ANN = joblib.load(path + '/base/CTD_ANN.m')
            CTD_XGB = joblib.load(path + '/base/CTD_XGB.m')
            CTD_KNN = joblib.load(path + '/base/CTD_KNN.m')
            CTD_LR = joblib.load(path + '/base/CTD_LR.m')

            AAE_ERT = joblib.load(path + '/base/AAE_ERT.m')
            AAE_ANN = joblib.load(path + '/base/AAE_ANN.m')
            AAE_XGB = joblib.load(path + '/base/AAE_XGB.m')
            AAE_KNN = joblib.load(path + '/base/AAE_KNN.m')
            AAE_LR = joblib.load(path + '/base/AAE_LR.m')

            AAI_ERT = joblib.load(path + '/base/AAI_ERT.m')
            AAI_ANN = joblib.load(path + '/base/AAI_ANN.m')
            AAI_XGB = joblib.load(path + '/base/AAI_XGB.m')
            AAI_KNN = joblib.load(path + '/base/AAI_KNN.m')
            AAI_LR = joblib.load(path + '/base/AAI_LR.m')

            GAAC_ERT = joblib.load(path + '/base/GAAC_ERT.m')
            GAAC_ANN = joblib.load(path + '/base/GAAC_ANN.m')
            GAAC_XGB = joblib.load(path + '/base/GAAC_XGB.m')
            GAAC_KNN = joblib.load(path + '/base/GAAC_KNN.m')
            GAAC_LR = joblib.load(path + '/base/GAAC_LR.m')

            model_calculation(fasta=fasta, out=path + '/tmpLysinActivity.txt', path=path)

            with open(fasta) as fa:
                fa_dict = {}
                for line in fa:
                    line = line.replace('\n', '')
                    if line.startswith('>'):
                        seq_name = line[1:]
                        fa_dict[seq_name] = ''
                    else:
                        fa_dict[seq_name] += line.replace('\n', '')
            fa.close()

            lis = []
            with open(path + '/tmpLysinActivity.txt') as ac:
                for i in ac:
                    i = i.replace('\n', '')
                    lis.append(i)
            ac.close()

            for i_1 in range(0, len(lis)):
                key = list(fa_dict.keys())[i_1]
                val = [fa_dict.get(key, [])] + [lis[i_1]]
                fa_dict[key] = val

            with open(out, 'w') as f:
                for key in fa_dict:
                    lines = key + '\t' + fa_dict[key][0] + '\t' + fa_dict[key][1] + '\n'
                    print(lines)
                    f.write(lines)
            f.close()

            os.remove(path + '/tmpLysinActivity.txt')

            self.trigger.emit('Finished!!!')

        except Exception as ex:
            self.trigger.emit('Some errors have occurred, %s!' % ex)

class Lysins_Form(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.work = WorkThread()

    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(683, 533)
        Form.setWindowIcon(QIcon("./logo/logo.ico"))
        self.gridLayout_6 = QtWidgets.QGridLayout(Form)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setSpacing(0)
        self.gridLayout.setObjectName("gridLayout")
        self.textBrowser_2 = QtWidgets.QTextBrowser(Form)
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.gridLayout.addWidget(self.textBrowser_2, 2, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.label_3 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 4, 0, 1, 1)
        self.textBrowser_3 = QtWidgets.QTextBrowser(Form)
        self.textBrowser_3.setObjectName("textBrowser_3")
        self.gridLayout.addWidget(self.textBrowser_3, 5, 0, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout.addWidget(self.pushButton_2, 3, 0, 1, 1)
        self.pushButton_3 = QtWidgets.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.gridLayout.addWidget(self.pushButton_3, 6, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_7.setFont(font)
        self.label_7.setAlignment(QtCore.Qt.AlignCenter)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 1, 1, 1, 1)
        self.label_4 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label_4.setFont(font)
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 7, 0, 1, 1)
        self.textBrowser = QtWidgets.QTextBrowser(Form)
        self.textBrowser.setObjectName("textBrowser")
        self.gridLayout.addWidget(self.textBrowser, 8, 0, 1, 1)
        self.tableWidget = QtWidgets.QTableWidget(Form)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        self.gridLayout.addWidget(self.tableWidget, 2, 1, 4, 1)
        self.label_6 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 6, 1, 1, 1)
        self.pushButton_4 = QtWidgets.QPushButton(Form)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_4.sizePolicy().hasHeightForWidth())
        self.pushButton_4.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.pushButton_4.setFont(font)
        self.pushButton_4.setObjectName("pushButton_4")
        self.gridLayout.addWidget(self.pushButton_4, 7, 1, 1, 1)
        self.pushButton = QtWidgets.QPushButton(Form)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton.sizePolicy().hasHeightForWidth())
        self.pushButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(22)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 8, 1, 1, 1)
        self.gridLayout_6.addLayout(self.gridLayout, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

        # button action
        self.pushButton.clicked.connect(self.calculation)
        self.pushButton_2.clicked.connect(self.read_file1)
        self.pushButton_3.clicked.connect(self.read_file2)
        self.pushButton_4.clicked.connect(self.table_read)

        ## default
        self.textBrowser_2.setPlaceholderText("D:/input/test.fa")
        self.textBrowser_3.setPlaceholderText("D:/output/result.txt")

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Lysins activity"))
        self.label_2.setText(_translate("Form", "Input fasta file"))
        self.label.setText(_translate("Form", "Lysins activity"))
        self.label_3.setText(_translate("Form", "Output result file"))
        self.pushButton_2.setText(_translate("Form", "Choose"))
        self.pushButton_3.setText(_translate("Form", "Choose"))
        self.label_7.setText(_translate("Form", "Result table"))
        self.label_4.setText(_translate("Form", "Status"))
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Form", "ID"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Form", "Seq"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Form", "Pred"))
        self.label_6.setText(_translate("Form", "If the program is finished, click \'Table\' to display the result"))
        self.pushButton_4.setText(_translate("Form", "Table"))
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
            global fasta
            global out, Randon_seed

            fasta = self.textBrowser_2.toPlainText()
            out = self.textBrowser_3.toPlainText()

            def is_fasta(filename):
                with open(filename, "r") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    return any(fasta)

            Randon_seed = 100
            print(len(out))
            if 0 in [len(fasta), len(out)]:
                QMessageBox.warning(self, "warning", "Please add correct file path!", QMessageBox.Cancel)
            else:
                if is_fasta(fasta) == False:
                    QMessageBox.critical(self, "error", "Check fasta file format!")
                else:
                    self.textBrowser.setText('Running! please wait')
                    QApplication.processEvents()  # 逐条打印状态

                    # 加载模型
                    global path

                    path_1 = '.'
                    path_1 = os.path.abspath(path_1)
                    path_1 = path_1.split('\\')
                    path_1 = '/'.join(path_1)
                    path = path_1 + '/models/Lysins/Models'
                    print(path)

                    # 启动线程, 运行 run 函数
                    self.work.start()
                    # 传送信号, 接受 run 函数执行完毕后的信号
                    self.work.trigger.connect(self.finished)

        except:
            QMessageBox.critical(self, "error", "Check fasta file format!")

    def table_read(self):
        try:
            f = open(out)
            count = 0
            for line in f:
                count = count + 1

            nrows = int(count)
            ncols = 3
            self.tableWidget.setRowCount(nrows)  # 设置行数
            self.tableWidget.setColumnCount(ncols)

            f = open(out)
            row_num = 0
            for line in f:
                print(line)
                li = line.strip().split('\t')
                col_num = 0
                for i in li:
                    item = QTableWidgetItem(i)
                    print(item)
                    self.tableWidget.setItem(row_num, col_num, item)
                    print(row_num,col_num)
                    col_num = col_num + 1
                row_num = row_num + 1
        except:
            QMessageBox.critical(self, "error", "Please run program first!!!")


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Lysins_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

