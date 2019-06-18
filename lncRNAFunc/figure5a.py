#!/usr/bin/env python
# coding:utf-8
__author__ = 'yangrui'

import os
import pandas as pd
import numpy as np
from keras.models import load_model
from scipy.stats import pearsonr, spearmanr


lncRNAdb = '/public/noncode/yangrui/rna/hgdb/lncRNA.transcript.bed'


def find_lncRNAs_from_k562(rbp, way):
    infile = '{}/{}.reg.75.pad_{}.txt'.format(rbp, rbp, way)
    bedfile = '{}/{}.reg.75.{}.bed'.format(rbp, rbp, way)
    fout = open(bedfile, 'w')
    with open(infile) as f:
        head = f.readline().strip('\n').split('\t')
        cidx, sidx, eidx, stidx = head.index('chrom'), head.index('start'), head.index('end'), head.index('strand')
        lidx, ridx, idx = head.index('left_length'), head.index('right_length'), head.index('signal')
        for line in f:
            li = line.strip('\n').split('\t')
            chrom, start, end, strand = li[cidx], int(li[sidx]), int(li[eidx]), li[stidx]
            left, right, signal = int(li[lidx]), int(li[ridx]), li[idx]
            name = ':'.join([chrom, str(start), str(end), strand, signal])
            start += left
            end -= right
            fout.write('\t'.join([chrom, str(start), str(end), signal, name, strand])+'\n')
    fout.close()

    intersectBed = '{}/{}.reg.75.{}.intersect.bed'.format(rbp, rbp, way)
    cmd = 'bedtools intersect -a {bed1} -b {bed2} -wa -wb > {out}'.format(bed1=bedfile, bed2=lncRNAdb, out=intersectBed)
    print(cmd)
    os.system(cmd)
    return intersectBed


def extract_predict_info(filename):
    hdict = {}
    with open(filename) as f:
        head = f.readline().strip().split('\t')
        cidx, sidx, eidx, stidx, idx = head.index('chrom'), head.index('start'), head.index('end'), head.index('strand'), head.index('signal')
        for line in f:
            li = line.strip('\n').split('\t')
            chrom, start, end, strand, signal = li[cidx], li[sidx], li[eidx], li[stidx], li[idx]
            key = ':'.join([chrom, start, end, strand, signal])
            if key not in hdict.keys():
                hdict[key] = li
            else:
                print(key, 'replicate')
    return hdict, head


def find_lncRNAs_from_k562_main(rbps):
    # 筛选出每个rbp样本中lncRNA的peak的CLIP-seq signal及模型预测值

    for rbp in rbps:
        print(rbp)
        testfile = find_lncRNAs_from_k562(rbp, 'test')
        validfile = find_lncRNAs_from_k562(rbp, 'valid')
        fp1 = open(testfile)
        fp2 = open(validfile)
        lines1 = fp1.readlines()
        lines2 = fp2.readlines()
        fp1.close()
        fp2.close()
        tdict, head = extract_predict_info('{}/{}.reg.75.pad_{}.normal.txt'.format(rbp, rbp, 'test'))
        vdict, head = extract_predict_info('{}/{}.reg.75.pad_{}.normal.txt'.format(rbp, rbp, 'valid'))
        fout = open('{}/{}.lncRNA.reg.inpt'.format(rbp, rbp), 'w')
        fout.write('\t'.join(head + ['gene', 'gid', 'tid']) + '\n')
        for line in lines1 + lines2:
            li = line.strip('\n').split('\t')
            start, end = int(li[1]), int(li[2])
            start2, end2 = int(li[7]), int(li[8])
            if start >= start2 and end <= end2:
                key = li[4]
                gene, gid, tid = li[9], li[10], li[12]
                if key in tdict.keys():
                    value = tdict[key]
                    value[1], value[2] = str(start), str(end)
                    fout.write('\t'.join(value + [gene, gid, tid]) + '\n')
                elif key in vdict.keys():
                    value = vdict[key]
                    value[1], value[2] = str(start), str(end)
                    fout.write('\t'.join(value + [gene, gid, tid]) + '\n')
                else:
                    print(key)
        fout.close()


class Sample:

    def __init__(self, input):
        data = pd.read_csv(input, sep='\t')
        self.df = data
        self.seqs = list(data['center_seq'])
        self.strands = list(data['strand'])
        self.genes = list(data['gene'])
        self.gdis = list(data['gid'])
        self.tids = list(data['tid'])
        self.chroms = list(data['chrom'])
        self.starts = list(data['start'])
        self.ends = list(data['end'])

    def load_data(self, width):
        y = self.df['signal'].apply(lambda x: float(x))
        y = np.array(list(y))

        self.df['x_center'] = self.df['center_seq'].apply(lambda x: self._get_array(x))
        x_pad = self.df['seq'].apply(lambda row: self._get_array(self._padseq(row)))

        x_center = self.df['x_center']
        x_center = np.array(list(x_center)).reshape(-1, width, 4)
        x_pad = np.array(list(x_pad)).reshape(-1, width+150*2, 4)

        return [x_center, x_pad], y

    def _padseq(self, seq):
        ival = 375 - len(seq)
        x = ival // 2
        return 'N' * x + seq + 'N' * (ival - x)

    def _get_array(self, seq):
        base_dict = {'A': 0, 'G': 1, 'C': 2, 'T': 3}
        lst = []
        for c in seq:
            item = ['0'] * 4
            if c != 'N':
                item[base_dict[c]] = '1'
            else:
                item = ['0.25'] * 4
            lst.append(item)
        return lst


def compute_pearson(model, x_test, y_test, rbp, genes, gids, tids, chroms, starts, ends, strands):
    outfile = '{}/{}.lncRNA.reg.predict'.format(rbp, rbp)
    fout = open(outfile, 'w')
    y_pred_test = model.predict(x_test)[:, 0]
    # res_test = pearsonr(y_pred_test, y_test)
    # res_test2 = spearmanr(y_pred_test, y_test)

    for i in range(len(y_pred_test)):
        fout.write('\t'.join([str(y_pred_test[i]), str(y_test[i]), chroms[i], str(starts[i]), str(ends[i]), strands[i], genes[i], gids[i], tids[i]])+'\n')
    fout.close()

    return outfile


def searchmodel(models):
    model = ''
    epoch = 0
    for item in models:
        tmp = int(item.split('-')[0].split('.')[-1])
        if tmp > epoch:
            epoch = tmp
            model = item
    return model


def plot_file(rbp, infile):
    outfile = '{}/{}.lncRNA.scatter.plot'.format(rbp, rbp)
    hdict = {}
    with open(infile) as f:
        for line in f:
            li = line.strip().split('\t')
            chrom, start, end, strand = li[2], li[3], li[4], li[5]
            key = ''.join([chrom, start, end, strand])
            hdict[key] = [li[0], li[1]]

    y_pred, y_true = [], []
    fout = open(outfile, 'w')
    fout.write('\t'.join(['y_pred', 'y_true'])+'\n')
    for value in hdict.values():
        fout.write('\t'.join([value[0], value[1]])+'\n')
        y_pred.append(float(value[0]))
        y_true.append(float(value[1]))
    fout.close()

    pcc = pearsonr(y_pred, y_true)
    print('\t'.join([str(round(pcc[0], 3)), rbp, str(len(y_true))]))


def figure5a(rbp):
    way = 'seq_5.1'
    width = 75
    prefix = '{rbp}/{rbp}.reg'.format(rbp=rbp)
    testfile = '{}/{}.lncRNA.reg.inpt'.format(rbp, rbp)
    test = Sample(testfile)
    X_test, y_test = test.load_data(width)
    genes = test.genes
    gids = test.gdis
    tids = test.tids
    chroms = test.chroms
    starts = test.starts
    ends = test.ends
    strands = test.strands

    try:
        modeldir = 'model_{}_{}'.format(way, width)
        models = os.listdir('{prefix}/{modeldir}'.format(prefix=prefix, modeldir=modeldir))
        importModel = '{prefix}/{modeldir}/{model}'.format(prefix=prefix, modeldir=modeldir, model=searchmodel(models))
    except:
        importModel = '{}/model_{}_{}/model.hdf5'.format(prefix, way, width)

    model = load_model(importModel)

    out = compute_pearson(model, X_test, y_test, rbp, genes, gids, tids, chroms, starts, ends, strands)
    out = '{}/{}.lncRNA.reg.predict'.format(rbp, rbp)
    plot_file(rbp, out)


if __name__ == '__main__':
    fp = open('rbps.filter.txt')
    rbps = [item.strip('\n').split()[0] for item in fp.readlines()]
    fp.close()
    find_lncRNAs_from_k562_main(rbps)

    for rbp in rbps:
        figure5a(rbp)

    cmd = 'Rscript figure5a.r figure5a.txt'
    os.system(cmd)