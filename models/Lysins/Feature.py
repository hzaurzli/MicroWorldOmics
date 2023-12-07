import re, os, sys, math
import pandas as pd
import joblib
import numpy as np
from sklearn.ensemble import RandomForestClassifier as RF
import argparse

def readFasta(file):
    if os.path.exists(file) == False:
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open(file) as f:
        records = f.read()

    if re.search('>', records) == None:
        print('The input file seems not in fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    myFasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append(sequence)

    return myFasta


letters = list('ACDEFGHIKLMNPQRSTVWY')


# AAE
def AAE_1(fastas):
    length = float(len(fastas))
    amino_acids = dict.fromkeys(letters, 0)
    encodings = []
    for AA in amino_acids:
        hits = [a.start() for a in list(re.finditer(AA, fastas))]
        p_prev = 0
        p_next = 1
        sum = 0
        while p_next < len(hits):
            distance = (hits[p_next] - hits[p_prev]) / length
            sum += distance * math.log(distance, 2)
            p_prev = p_next
            p_next += 1
        amino_acids[AA] = -sum
        encodings.append(amino_acids[AA])
    return encodings


def AAE(seq):
    encodings = []
    for fastas in seq:
        fastas_NT5 = "%s" % fastas[:5]
        fastas_CT5 = "%s" % fastas[-5:]
        encodings_full = AAE_1(fastas)
        encodings_CT5 = AAE_1(fastas_CT5)
        encodings_NT5 = AAE_1(fastas_NT5)
        encodings.append(encodings_full + encodings_NT5 + encodings_CT5)
    return encodings


# AAI
def AAI_1(fastas,path):
    encodings = []
    fileAAindex1 = open(path + '/' + R'Features/AAIindex/AAindex_1.txt')
    fileAAindex2 = open(path + '/' + R'Features/AAIindex/AAindex_2.txt')
    records1 = fileAAindex1.readlines()[1:]
    records2 = fileAAindex2.readlines()[1:]
    AAindex1 = []
    AAindex2 = []
    for i in records1:
        AAindex1.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
    for i in records2:
        AAindex2.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
    index = {}
    for i in range(len(letters)):
        index[letters[i]] = i
    fastas_len = len(fastas)
    for i in range(len(AAindex1)):
        total = 0
        for j in range(fastas_len):
            temp = AAindex1[i][index[fastas[j]]]
            total = total + float(temp)
        encodings.append(total / fastas_len)
    for i in range(len(AAindex2)):
        total = 0
        for j in range(fastas_len):
            temp = AAindex2[i][index[fastas[j]]]
            total = total + float(temp)
        encodings.append(total)
    return encodings


def AAI(seqs, path):
    encodings = []
    for fastas in seqs:
        fastas_NT5 = "%s" % fastas[:5]
        fastas_CT5 = "%s" % fastas[-5:]

        encodings_full = AAI_1(fastas, path)
        encodings_CT5 = AAI_1(fastas_CT5, path)
        encodings_NT5 = AAI_1(fastas_NT5, path)
        encodings.append(encodings_full + encodings_NT5 + encodings_CT5)
    return encodings


# BPNC
def BPNC(seqs):
    encodings = []
    for fastas in seqs:
        fastas_NT5 = "%s" % fastas[:5]
        fastas_CT5 = "%s" % fastas[-5:]
        full = fastas_NT5+fastas_CT5
        encoding = []
        for AA in full:
            for AA1 in letters:
                tag = 1 if AA == AA1 else 0
                encoding.append(tag)
        encodings.append(encoding)
    return encodings


# CTD_function
group1 = {'hydrophobicity_PRAM900101': 'RKEDQN', 'normwaalsvolume': 'GASTPDC', 'polarity': 'LIFWCMVY',
          'polarizability': 'GASDT', 'charge': 'KR', 'secondarystruct': 'EALMQKRH', 'solventaccess': 'ALFCGIVW'}
group2 = {'hydrophobicity_PRAM900101': 'GASTPHY', 'normwaalsvolume': 'NVEQIL', 'polarity': 'PATGS',
          'polarizability': 'CPNVEQIL', 'charge': 'ANCQGHILMFPSTWYV', 'secondarystruct': 'VIYCWFT',
          'solventaccess': 'RKQEND'}
group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW', 'normwaalsvolume': 'MHKFRYW', 'polarity': 'HQRKNED',
          'polarizability': 'KMHFRYW', 'charge': 'DE', 'secondarystruct': 'GNPSD', 'solventaccess': 'MSPTHY'}
groups = [group1, group2, group3]
propertys = ('hydrophobicity_PRAM900101', 'normwaalsvolume', 'polarity', 'polarizability', 'charge', 'secondarystruct',
             'solventaccess')


def Count_C(sequence1, sequence2):
    sum = 0
    for aa in sequence1:
        sum = sum + sequence2.count(aa)
    return sum


def Count_D(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >= 1 else 1 for i in cutoffNums]
    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence))
                    break
        if myCount == 0:
            code.append(0)
    return code


def CTD(seqs):
    encodings = []
    for seq in seqs:
        code = []
        code2 = []
        CTDD1 = []
        CTDD2 = []
        CTDD3 = []
        aaPair = [seq[j:j + 2] for j in range(len(seq) - 1)]
        for p in propertys:
            c1 = Count_C(group1[p], seq) / len(seq)
            c2 = Count_C(group2[p], seq) / len(seq)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]

            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code2 = code2 + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
            CTDD1 = CTDD1 + [value / float(len(seq)) for value in Count_D(group1[p], seq)]
            CTDD2 = CTDD2 + [value / float(len(seq)) for value in Count_D(group2[p], seq)]
            CTDD3 = CTDD3 + [value / float(len(seq)) for value in Count_D(group3[p], seq)]
        encodings.append(code + code2 + CTDD1 + CTDD2 + CTDD3)
    return encodings


# DPC
diPeptides = [aa1 + aa2 for aa1 in letters for aa2 in letters]


def DPC(seqs):
    encodings = []
    for seq in seqs:
        AADict = {}
        for aa in range(len(letters)):
            AADict[letters[aa]] = aa

        tmpCode = [0] * 400
        for j in range(len(seq) - 2 + 1):
            tmpCode[AADict[seq[j]] * 20 + AADict[seq[j + 1]]] = tmpCode[AADict[seq[j]] * 20 + AADict[seq[j + 1]]] + 1
        if sum(tmpCode) != 0:
            tmpDPC = [i / sum(tmpCode) for i in tmpCode]
        encodings.append(tmpDPC)
    return encodings


# GTPC
group = {'alphatic': 'GAVLMI', 'aromatic': 'FYW', 'postivecharge': 'KRH', 'negativecharge': 'DE', 'uncharge': 'STCPNQ'}
groupKey = group.keys()
tripeptide = [g1 + '.' + g2 + '.' + g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]


def GTPC(seqs):
    encodings = []
    for seq in seqs:
        index = {}
        myDict = {}
        for key in groupKey:
            for aa in group[key]:
                index[aa] = key
        for t in tripeptide:
            myDict[t] = 0
        sum = 0
        for j in range(len(seq) - 2):
            myDict[index[seq[j]] + '.' + index[seq[j + 1]] + '.' + index[seq[j + 2]]] = myDict[index[seq[j]] + '.' + index[
                seq[j + 1]] + '.' + index[seq[j + 2]]] + 1
            sum = sum + 1
        code = []
        if sum == 0:
            for t in tripeptide:
                code.append(0)
        else:
            for t in tripeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)
    return encodings


def savetsv(encodings, file='feature.tsv'):
    with open(file, 'w') as f:
        if encodings == 0:
            f.write('Descriptor calculation failed.')
        else:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i) - 1])) + '\n')
    return None


#[0-1]normalization
def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range

def all_feature(fastas,path):
    # feature_AAE = normalization(np.array(AAE(fastas)))
    # feature_AAI = normalization(np.array(AAI(fastas)))
    # feature_BPNC = normalization(np.array(BPNC(fastas)))
    # feature_CTD = normalization(np.array(CTD(fastas)))
    # feature_DPC = normalization(np.array(DPC(fastas)))
    # feature_GTPC = normalization(np.array(GTPC(fastas)))

    feature_AAE = np.array(AAE(fastas))
    feature_AAI = np.array(AAI(fastas,path))
    feature_BPNC = np.array(BPNC(fastas))
    feature_CTD = np.array(CTD(fastas))
    feature_DPC = np.array(DPC(fastas))
    feature_GTPC = np.array(GTPC(fastas))

    _ = [feature_AAE.shape[1], feature_AAI.shape[1],
         feature_BPNC.shape[1], feature_CTD.shape[1],
         feature_DPC.shape[1], feature_GTPC.shape[1]]
    feature_subset_index = []
    idx = 0
    for __ in _:
        feature_subset_index.append([fea_idx + idx for fea_idx in range(__)])
        idx += __

    encodings = np.hstack((feature_AAE, feature_AAI,
                             feature_BPNC, feature_CTD,
                             feature_DPC, feature_GTPC))
    return encodings, feature_subset_index

def input_args():
    """
    Usage:
    python Features.py -p pos.fasta -n neg.fasta -o ./data/Features/Data1.csv
    """

    parser = argparse.ArgumentParser(usage="Usage Tip;",
                                     description = "DeepLysin Feature Extraction")
    parser.add_argument("--file", "-f", required = True,
                        help = "input file(.fasta)")
    parser.add_argument("--method", "-m", default='ALL',
                        help = "select the sunsets of features(example:AAE or ALL)")
    parser.add_argument("--label", "-l", required=True, choices=[0,1], type=int,
                        help = "sample label,1 or 0")
    parser.add_argument("--model_path", "-m", required=True, help="Model path")
    parser.add_argument("--out", "-o", required=True, help="output path and filename")
    return parser.parse_args()


if __name__ == '__main__':
    args = input_args()
    Select_feature = {
        'AAE':'AAE(fastas)',
        'AAI':'AAI(fastas)',
        'BPNC':'BPNC(fastas)',
        'CTD':'CTD(fastas)',
        'DPC':'DPC(fastas)',
        'GTPC':'GTPC(fastas)',
    }
    fastas = readFasta(args.file)
    if args.method == 'ALL':
        method,feature_subset_index = all_feature(fastas, os.path.abspath(Args.model_path))
        print(feature_subset_index)
    else:
        method = Select_feature.get(args.method, None)
        method = eval(method)
    encodings = pd.DataFrame(method)
    column_1 = encodings.columns.tolist()
    column_1.insert(0,"class")
    encodings["class"] = args.label
    encodings = encodings.reindex(columns=column_1)
    encodings.to_csv(args.out,index=False,sep=',')


