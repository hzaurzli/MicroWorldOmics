# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'PhaTYP.ui'
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
import pandas as pd
import numpy as np
import pickle as pkl
import subprocess
import shutil
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess as sub
import torch
from datasets import DatasetDict,Dataset
import argparse
import pandas as pd
import pyarrow as pa
import numpy as np
import pickle as pkl
from scipy.special import softmax
from transformers import AutoTokenizer
from transformers import DataCollatorWithPadding
from transformers import AutoModelForSequenceClassification
from transformers import BertTokenizer, LineByLineTextDataset
from transformers import BertConfig, BertForMaskedLM, DataCollatorForLanguageModeling
from transformers import TrainingArguments, Trainer
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from Bio import SeqIO
import psutil


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
            if os.path.exists(out_tmp):
                os.remove(out_tmp)
            else:
                event.ignore()  # 设置正常退出
        except:
            return None  # 设置正常退出


class WorkThread(QThread):
    # 自定义信号对象
    trigger = pyqtSignal(str)

    def __int__(self):
        # 初始化函数
        super(WorkThread, self).__init__()

    def run(self):
        def read_abc(position, output):
            with open(output, 'w') as abc:
                with open(position) as tab:
                    for line in tab:
                        line = line.replace('\n', '')
                        li = line.split('\t')
                        if li[0] != li[1]:
                            lines = li[0] + '\t' + li[1] + '\t' + li[10] + '\n'
                            abc.write(lines)
                tab.close()
            abc.close()

        def check_process_running(process_name):  # 检查进程是否运行
            for process in psutil.process_iter(['name']):
                if process.info['name'] == process_name:
                    return True
            return False

        if not os.path.isdir(out_fn):
            os.makedirs(out_fn)

        if not os.path.isdir(transformer_fn):
            os.makedirs(transformer_fn)

        rec = []
        for record in SeqIO.parse(contigs, 'fasta'):
            if len(record.seq) > int(length):
                rec.append(record)
        SeqIO.write(rec, out_fn + '/filtered_contigs.fa', 'fasta')

        try:
            p = sub.Popen(
                r'.\tools\prodigal\prodigal.exe -i ' + out_fn + '/filtered_contigs.fa -a ' + out_fn + '/test_protein.fa -f gff -p meta')
            p.wait()
        except:
            QMessageBox.critical(self, "error", "check fasta file format!")

        try:
            # create database
            p = sub.Popen(
                r'.\tools\diamond\diamond.exe makedb --threads 2 --in ./models/PhaTYP/database/database.fa -d ' + out_fn + '/database.dmnd',
                shell=False)
            p.wait()
            # running alignment
            p = sub.Popen(
                r'.\tools\diamond\diamond.exe blastp --threads 2 --sensitive -d ' + out_fn + '/database.dmnd -q ' + out_fn + '/test_protein.fa -o ' + out_fn + '/results.tab -k 1',
                shell=False)
            p.wait()

            read_abc(out_fn + '/results.tab', out_fn + '/results.abc')
        except:
            QMessageBox.critical(self, "error", "check fasta file format!")

        # Load dictonary and BLAST results D:\tools\Pycharm\pyqt\models\PhaTYP\database
        proteins_df = pd.read_csv('./models/PhaTYP/database/proteins.csv')
        proteins_df.dropna(axis=0, how='any', inplace=True)
        pc2wordsid = {pc: idx for idx, pc in enumerate(sorted(set(proteins_df['cluster'].values)))}
        protein2pc = {protein: pc for protein, pc in
                      zip(proteins_df['protein_id'].values, proteins_df['cluster'].values)}
        blast_df = pd.read_csv(out_fn + "/results.abc", sep='\t', names=['query', 'ref', 'evalue'])

        # Parse the DIAMOND results
        contig2pcs = {}
        for query, ref, evalue in zip(blast_df['query'].values, blast_df['ref'].values,
                                      blast_df['evalue'].values):
            conitg = query.rsplit('_', 1)[0]
            idx = query.rsplit('_', 1)[1]
            pc = pc2wordsid[protein2pc[ref]]
            try:
                contig2pcs[conitg].append((idx, pc, evalue))
            except:
                contig2pcs[conitg] = [(idx, pc, evalue)]

        # Sorted by position
        for contig in contig2pcs:
            contig2pcs[contig] = sorted(contig2pcs[contig], key=lambda tup: tup[0])

        # Contigs2sentence
        contig2id = {contig: idx for idx, contig in enumerate(contig2pcs.keys())}
        id2contig = {idx: contig for idx, contig in enumerate(contig2pcs.keys())}
        sentence = np.zeros((len(contig2id.keys()), 300))
        sentence_weight = np.ones((len(contig2id.keys()), 300))
        for row in range(sentence.shape[0]):
            contig = id2contig[row]
            pcs = contig2pcs[contig]
            for col in range(len(pcs)):
                try:
                    _, sentence[row][col], sentence_weight[row][col] = pcs[col]
                    sentence[row][col] += 1
                except:
                    break

        # propostion
        rec = []
        for key in blast_df['query'].values:
            name = key.rsplit('_', 1)[0]
            rec.append(name)
        counter = Counter(rec)
        mapped_num = np.array([counter[item] for item in id2contig.values()])

        rec = []
        for record in SeqIO.parse(out_fn + '/test_protein.fa', 'fasta'):
            name = record.id
            name = name.rsplit('_', 1)[0]
            rec.append(name)
        counter = Counter(rec)
        total_num = np.array([counter[item] for item in id2contig.values()])
        proportion = mapped_num / total_num

        # Store the parameters
        pkl.dump(sentence, open(transformer_fn + '/sentence.feat', 'wb'))
        pkl.dump(id2contig, open(transformer_fn + '/sentence_id2contig.dict', 'wb'))
        pkl.dump(proportion, open(transformer_fn + '/sentence_proportion.feat', 'wb'))
        pkl.dump(pc2wordsid, open(transformer_fn + '/pc2wordsid.dict', 'wb'))

        feat = pkl.load(open(transformer_fn + '/sentence.feat', 'rb'))
        pcs = pkl.load(open('./models/PhaTYP/database/pc2wordsid.dict', 'rb'))
        id2pcs = {item: key for key, item in pcs.items()}
        text = []
        label = []
        for line in feat:
            sentence = ""
            flag = 0
            for i in range(len(line) - 2):
                if line[i] - 1 == -1:
                    flag = 1
                    sentence = sentence[:-1]
                    break
                sentence = sentence + id2pcs[line[i] - 1] + ' '
            if flag == 0:
                sentence = sentence[:-1]
            text.append(sentence)
            label.append(1)

        feat_df = pd.DataFrame({'label': label, 'text': text})
        feat_df.to_csv(transformer_fn + '/bert_feat.csv', index=None)

        out_dir = os.path.dirname(out)
        if out_dir != '':
            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)

        id2contig = pkl.load(open(transformer_fn + '/sentence_id2contig.dict', 'rb'))
        bert_feat = pd.read_csv(transformer_fn + '/bert_feat.csv')

        SENTENCE_LEN = 300  # len
        NUM_TOKEN = 45583  # PC

        CONFIG_DIR = "./models/PhaTYP/config"
        OUTPUT_DIR = "finetune"

        # load the token configuration
        tokenizer = BertTokenizer.from_pretrained(CONFIG_DIR, do_basic_tokenize=False)

        def preprocess_function(examples):
            return tokenizer(examples["text"], truncation=True)

        train = pa.Table.from_pandas(bert_feat)
        test = pa.Table.from_pandas(bert_feat)
        train = Dataset(train)
        test = Dataset(test)

        data = DatasetDict({"train": train, "test": test})

        tokenized_data = data.map(preprocess_function, batched=True)
        data_collator = DataCollatorWithPadding(tokenizer=tokenizer)
        model = AutoModelForSequenceClassification.from_pretrained("./models/PhaTYP/model/", num_labels=2)

        training_args = TrainingArguments(
            output_dir='./models/PhaTYP/results',
            overwrite_output_dir=False,
            do_train=True,
            do_eval=True,
            learning_rate=2e-5,
            num_train_epochs=10,
            per_device_train_batch_size=32,
            per_device_eval_batch_size=32,
            weight_decay=0.01,
        )

        trainer = Trainer(
            model=model,
            args=training_args,
            train_dataset=tokenized_data["train"],
            eval_dataset=tokenized_data["test"],
            tokenizer=tokenizer,
            data_collator=data_collator,
        )

        with torch.no_grad():
            pred, label, metric = trainer.predict(tokenized_data["test"])

        prediction_value = []
        for item in pred:
            prediction_value.append(softmax(item))
        prediction_value = np.array(prediction_value)

        all_pred = []
        all_score = []
        for score in prediction_value:
            pred = np.argmax(score)
            if pred == 1:
                all_pred.append('temperate')
                all_score.append(score[1])
            else:
                all_pred.append('virulent')
                all_score.append(score[0])

        pred_csv = pd.DataFrame({"Contig": id2contig.values(), "Pred": all_pred, "Score": all_score})
        pred_csv.to_csv(out, index=False)

        process_name = 'diamond.exe'
        time.sleep(3)
        while True:  # 判断 iqtree.exe 是否运行完成
            if check_process_running(process_name):
                print(f"The process {process_name} is running.")
                time.sleep(10)
                continue
            else:
                print(f"The process {process_name} is not running.")
                break

        self.trigger.emit('Finished!!!' + '\n' + 'lysogen_prediction.csv is your result!!!')

class PhaTYP_Form(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.work = WorkThread()

    def setupUi(self, Clustal):
        Clustal.setObjectName("Clustal")
        Clustal.resize(702, 467)
        Clustal.setStyleSheet("background-image: url(D:/Documents/Desktop/bb.png)")
        self.label = QtWidgets.QLabel(Clustal)
        self.label.setGeometry(QtCore.QRect(250, 10, 191, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Clustal)
        self.label_2.setGeometry(QtCore.QRect(60, 50, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_2.setFont(font)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.textBrowser_2 = QtWidgets.QTextBrowser(Clustal)
        self.textBrowser_2.setGeometry(QtCore.QRect(60, 80, 181, 31))
        self.textBrowser_2.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.label_3 = QtWidgets.QLabel(Clustal)
        self.label_3.setGeometry(QtCore.QRect(60, 120, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName("label_3")
        self.textBrowser_3 = QtWidgets.QTextBrowser(Clustal)
        self.textBrowser_3.setGeometry(QtCore.QRect(60, 160, 181, 31))
        self.textBrowser_3.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser_3.setObjectName("textBrowser_3")
        self.pushButton_2 = QtWidgets.QPushButton(Clustal)
        self.pushButton_2.setGeometry(QtCore.QRect(260, 80, 61, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_3 = QtWidgets.QPushButton(Clustal)
        self.pushButton_3.setGeometry(QtCore.QRect(260, 160, 61, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.textBrowser = QtWidgets.QTextBrowser(Clustal)
        self.textBrowser.setGeometry(QtCore.QRect(60, 340, 191, 91))
        self.textBrowser.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textBrowser.setObjectName("textBrowser")
        self.label_4 = QtWidgets.QLabel(Clustal)
        self.label_4.setGeometry(QtCore.QRect(100, 290, 101, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(19)
        self.label_4.setFont(font)
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.pushButton = QtWidgets.QPushButton(Clustal)
        self.pushButton.setGeometry(QtCore.QRect(560, 400, 81, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(22)
        self.pushButton.setFont(font)
        self.pushButton.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.pushButton.setObjectName("pushButton")
        self.pushButton_4 = QtWidgets.QPushButton(Clustal)
        self.pushButton_4.setGeometry(QtCore.QRect(350, 350, 331, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.pushButton_4.setFont(font)
        self.pushButton_4.setObjectName("pushButton_4")
        self.label_6 = QtWidgets.QLabel(Clustal)
        self.label_6.setGeometry(QtCore.QRect(350, 320, 331, 20))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_6.setObjectName("label_6")
        self.tableWidget = QtWidgets.QTableWidget(Clustal)
        self.tableWidget.setGeometry(QtCore.QRect(355, 70, 321, 241))
        self.tableWidget.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        self.label_5 = QtWidgets.QLabel(Clustal)
        self.label_5.setGeometry(QtCore.QRect(60, 200, 161, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        self.label_5.setFont(font)
        self.label_5.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_5.setObjectName("label_5")
        self.textEdit = QtWidgets.QTextEdit(Clustal)
        self.textEdit.setGeometry(QtCore.QRect(60, 230, 181, 31))
        self.textEdit.setStyleSheet("background-image: url(D:/Documents/Desktop/white.png)")
        self.textEdit.setObjectName("textEdit")

        self.retranslateUi(Clustal)
        QtCore.QMetaObject.connectSlotsByName(Clustal)

        # button action
        self.pushButton.clicked.connect(self.calculation)
        self.pushButton_2.clicked.connect(self.read_file1)
        self.pushButton_3.clicked.connect(self.read_file2)
        self.pushButton_4.clicked.connect(self.table_read)

        self.textEdit.setPlaceholderText(" Contig length filter: 3000")

    def retranslateUi(self, Clustal):
        _translate = QtCore.QCoreApplication.translate
        Clustal.setWindowTitle(_translate("Clustal", "PhaTYP"))
        self.label.setText(_translate("Clustal", "PhaTYP"))
        self.label_2.setText(_translate("Clustal", "Input fasta file"))
        self.label_3.setText(_translate("Clustal", "Output fasta file"))
        self.pushButton_2.setText(_translate("Clustal", "Choose"))
        self.pushButton_3.setText(_translate("Clustal", "Choose"))
        self.label_4.setText(_translate("Clustal", "Status"))
        self.pushButton.setText(_translate("Clustal", "Run"))
        self.pushButton_4.setText(_translate("Clustal", "Table"))
        self.label_6.setText(_translate("Clustal", "If the program is finished, click \'Table\' to display the result"))
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Clustal", "Contig"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Clustal", "Pred"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Clustal", "Score"))
        self.label_5.setText(_translate("Clustal", "Contig length"))


    def read_file1(self):
        openfile_name = QtWidgets.QFileDialog.getOpenFileName(self, 'choose file', '')[0]
        print(openfile_name)
        self.textBrowser_2.setText(openfile_name)


    def read_file2(self):
        openfile_name = QtWidgets.QFileDialog.getExistingDirectory(self, "choose file", "./")
        print(openfile_name)
        self.textBrowser_3.setText(openfile_name)

    def finished(self, str):
        self.textBrowser.setText(str)

    def calculation(self):
        try:
            def is_fasta(filename):
                with open(filename, "r") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    return any(fasta)

            global out, contigs, out_fn, transformer_fn, length
            contigs = self.textBrowser_2.toPlainText()
            out_fn = self.textBrowser_3.toPlainText()
            transformer_fn = out_fn + '/phatyp'

            out = out_fn + '/lysogen_prediction.csv'

            if any([len(contigs), len(out_fn)]) == False:
                QMessageBox.warning(self, "warning", "Please add correct file path!", QMessageBox.Cancel)
            else:
                if is_fasta(contigs) == False:
                    QMessageBox.critical(self, "error", "Check fasta file format!")
                else:
                    self.textBrowser.setText(
                        'Running! please wait(5-8mins)' + '\n' + 'If no response,never close window!!!')
                    QApplication.processEvents()  # 逐条打印状态

                    try:
                        length = int(self.textEdit.toPlainText())
                    except:
                        length = 3000

                    # 启动线程, 运行 run 函数
                    self.work.start()
                    # 传送信号, 接受 run 函数执行完毕后的信号
                    self.work.trigger.connect(self.finished)

        except:
            QMessageBox.critical(self, "error", "Check fasta file format!")

    def table_read(self):
        try:
            global out_tmp
            out_p = os.path.dirname(out)
            out_tmp = out_p + '/lysogen_prediction_tmp.csv'

            with open(out_tmp, 'w') as w:
                f = open(out)
                count = 0
                for line in f:
                    if count == 0:
                        print('1')
                    else:
                        w.write(line)
                        print(line)
                    count = count + 1
            w.close()

            f = open(out_tmp)
            count = 0
            for line in f:
                count = count + 1

            nrows = int(count)
            print(nrows)
            ncols = 3
            self.tableWidget.setRowCount(nrows)  # 设置行数
            self.tableWidget.setColumnCount(ncols)

            f = open(out_tmp)
            row_num = 0
            for line in f:
                print(line)
                li = line.strip().split(',')
                col_num = 0
                for i in li:
                    item = QTableWidgetItem(i)
                    print(item)
                    self.tableWidget.setItem(row_num, col_num, item)
                    print(row_num, col_num)
                    col_num = col_num + 1
                row_num = row_num + 1

        except:
            QMessageBox.critical(self, "error", "Please run program first!!!")


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    # Clustal = QtWidgets.QWidget()
    WT = winTest()
    ui = PhaTYP_Form()
    ui.setupUi(WT)
    WT.show()
    sys.exit(app.exec_())