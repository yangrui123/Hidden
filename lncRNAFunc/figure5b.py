#!/usr/bin/env python
# coding:utf-8
__author__ = 'yangrui'

import sys
import os
import pandas as pd
import numpy as np
from glob import glob
from scipy.stats import pearsonr, spearmanr
import math



class Sample:

    def __init__(self, input):
        self.df = pd.read_csv(input, sep='\t')
        self.genes = list(self.df['gid'])
        self.chroms = list(self.df['chrom'])
        self.starts = list(self.df['start'])
        self.ends = list(self.df['end'])
        self.strands = list(self.df['strand'])

    def load_data(self):
        for seq in self.df['left']:
            if len(seq) != 150:
                print(seq)
        self.df['x_center'] = self.df['center'].apply(lambda x: self._get_array(x))
        x_pad = self.df['seq'].apply(lambda row: self._get_array(row))

        x_center = self.df['x_center']
        x_center = np.array(list(x_center)).reshape(-1, 75, 4)
        x_pad = np.array(list(x_pad)).reshape(-1, 375, 4)

        return [x_center, x_pad]

    def _get_array(self, seq):
        seq = seq.upper()
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


def getfasta(bedfile, outfile, length=None):
    cmd = "bedtools getfasta -fi /public/noncode/yangrui/rna/data/hg19/ucsc.hg19.fasta -bed {bedfile} " \
          "-fo {outfile} -s -name -tab".format(bedfile=bedfile, outfile=outfile)
    # print(cmd)
    os.system(cmd)
    return outfile


def center_pad(start, end):
    lst = []
    for i in range(start, end - 10, 50):
        lst.append([i, min(end, i + 75)])
    return lst


def splitbed(bedfile, rbp):
    hdict = {}
    with open(bedfile) as f:
        for line in f:
            li = line.strip().split('\t')
            # print(li)
            chrom, start, end, strand = li[0], li[1], li[2], li[5]
            gene, gid = li[3], li[4]
            hdict[gid] = (chrom, start, end, strand, gene)

    idx = []
    for i, gid in enumerate(hdict.keys()):
        if i % 100 == 0:
            fout = open('rbps/{}/lncRNA_{}.bed'.format(rbp, i // 100 + 1), 'w')
            idx.append(i//100+1)
        item = hdict[gid]
        start, end = int(item[1]), int(item[2])
        chrom, strand = item[0], item[3]
        regions = center_pad(start, end)
        for (s, e) in regions:
            left = max(s - 150, start)
            right = min(e + 150, end)
            name = ':'.join([chrom, str(left), str(right), gid, strand, str(s - left), str(e - s), str(right - e)])
            fout.write('\t'.join([chrom, str(left), str(right), name, '.', strand]))
            fout.write('\n')
    fout.close()
    return idx


def padfile(seqfile, outfile):
    fout = open(outfile, 'w')
    fout.write('\t'.join(['gid', 'center', 'left', 'right', 'seq', 'chrom', 'start', 'end', 'strand'])+'\n')
    with open(seqfile) as f:
        for line in f:
            li = line.strip().split('\t')
            infos = li[0].split(':')
            gene = infos[3]
            chrom, start, end, strand = infos[0], int(infos[1]), int(infos[2]), infos[4]
            seq = li[1]
            left, right = int(infos[-3]), int(infos[-1])
            start += left
            end -= right
            left_seq = 'N' * (150 - left) + seq[:left]
            center_seq = seq[left: len(seq) - right]
            ival = 75 - len(center_seq)
            x = ival // 2
            center_seq = 'N' * x + center_seq + 'N' * (ival - x)
            right_seq = seq[len(seq) - right:] + 'N' * (150 - len(seq[len(seq) - right:]))
            ival = 375 - len(seq)
            x = ival // 2
            seq = 'N' * x + seq + 'N' * (ival - x)
            fout.write('\t'.join([gene, center_seq, left_seq, right_seq, seq, chrom, str(start), str(end), strand]) + '\n')
    fout.close()


def creat_input(rbp):
    bedfile = 'rbps/{}/{}.bed'.format(rbp, rbp)
    idx = splitbed(bedfile, rbp)
    for i in idx:
        print(i)
        bedfile = 'rbps/{}/lncRNA_{}.bed'.format(rbp, i)
        seqfile = getfasta(bedfile, 'rbps/{}/lncRNA_{}.seq'.format(rbp, i))
        padfile(seqfile, 'rbps/{}/lncRNA_{}.pad.seq'.format(rbp, i))


def read_lncRNA():
    lncRNAdb = '../hgdb/lncRNA.transcript.max.bed'
    hdict = {}
    with open(lncRNAdb) as f:
        for line in f:
            li = line.strip().split('\t')
            gid = li[4]
            hdict[gid] = li
    return hdict


def extract_lncRNA_for_validation(rbp):
    dbdict = read_lncRNA()
    filename = '../k562_v3/{}/{}.lncRNA.reg.predict'.format(rbp, rbp)
    fp = open(filename)
    lines = fp.readlines()
    fp.close()
    gdict = {}
    for line in lines:
        li = line.strip().split('\t')
        gene, gid, tid = li[-3], li[-2], li[-1]
        items = dbdict[gid]
        chrom, start, end, strand = items[0], items[1], items[2], items[5]
        gdict[gid] = (chrom, start, end, gene, gid, strand, tid)

    outdir = 'rbps/{}'.format(rbp)
    if not os.path.exists(outdir):
        os.system('mkdir -p {}'.format(outdir))
    outfile = '{}/{}.bed'.format(outdir, rbp)
    fout = open(outfile, 'w')
    for value in gdict.values():
        fout.write('\t'.join(value)+'\n')
    fout.close()


def model_predict(model, filename, outfile):
    fout = open(outfile, 'w')
    print('reading file ......')
    test = Sample(filename)
    genes = test.genes
    chroms = test.chroms
    starts = test.starts
    ends = test.ends
    strands = test.strands
    X_test = test.load_data()

    print('predicting .....')
    y_pred = model.predict(X_test, batch_size=500)[:, 0]

    for i in range(len(genes)):
        fout.write('\t'.join([chroms[i], str(starts[i]), str(ends[i]), genes[i], str(y_pred[i]), strands[i]])+'\n')
    fout.close()


def score(rbp):
    from keras.models import load_model
    modelpath = 'models/{}.model.hdf5'.format(rbp)
    print('loading model ......')
    model = load_model(modelpath)

    seqs = glob('rbps/{}/lncRNA_*.pad.seq'.format(rbp))
    outs = []
    for seq in seqs:
        idx = seq.split('.')[0].split('_')[-1]
        print(seq)
        outfile = 'rbps/{}/predicts_{}.reg.bed'.format(rbp, idx)
        model_predict(model, seq, outfile)
        outs.append(outfile)
    outfile = 'rbps/{}/predict.bed'.format(rbp)
    cmd = 'cat '
    for res in outs:
        cmd = cmd + res + ' '
    cmd = cmd + ' > ' + outfile
    os.system(cmd)


def extract_regions(rbp):
    filename = '../k562_v3/{}/{}.lncRNA.reg.predict'.format(rbp, rbp)
    fp = open(filename)
    lines = fp.readlines()
    fp.close()
    gdict = {}
    for line in lines:
        li = line.strip().split('\t')
        gid = li[-2]
        chrom, start, end, strand = li[2], li[3], li[4], li[5]
        signal = li[1]
        gene = li[-3]
        if gid not in gdict.keys():
            gdict[gid] = set()
        gdict[gid].add((chrom, start, end, gene, gid, strand, signal))

    outfile = 'rbps/{}/lncRNA.signal'.format(rbp)
    fout = open(outfile, 'w')
    for gid in gdict.keys():
        for item in gdict[gid]:
            fout.write('\t'.join(list(item))+'\n')
    fout.close()


def intersect(rbp):
    bed1 = 'rbps/{}/lncRNA.signal'.format(rbp)
    bed2 = 'rbps/{}/predict.bed'.format(rbp)
    outbed = 'rbps/{}/lncRNA.intersect'.format(rbp)
    cmd = 'bedtools intersect -a {bed1} -b {bed2} -wa -wb > {out}'.format(bed1=bed1, bed2=bed2, out=outbed)
    os.system(cmd)
    return outbed


def compute_pcc(intersectBed, rbp):
    hdict = {}
    with open(intersectBed) as f:
        for line in f:
            li = line.strip().split('\t')
            chrom, start, end, strand, gid = li[0], li[1], li[2], li[5], li[4]
            chrom1, start1, end1, strand1, gid1 = li[7], li[8], li[9], li[12], li[10]
            signal = float(li[6])
            pred = float(li[11])
            key = ':'.join([chrom, start, end, strand, gid])
            length = min(int(end1), int(end)) - max(int(start), int(start1))
            if key not in hdict.keys() or(hdict[key][0] < length):
                hdict[key] = [length, signal, pred]
    y_true = []
    y_pred = []
    for value in hdict.values():
        y_true.append(value[1])
        y_pred.append(value[2])
    pcc = pearsonr(y_true, y_pred)
    scc = spearmanr(y_true, y_pred)
    outfile = 'rbps/{}/pcc.txt'.format(rbp)
    print(round(pcc[0], 3), '\t', round(scc[0], 3), '\t', rbp)
    fout = open(outfile, 'w')
    fout.write('\t'.join(['pcc', 'p_value', 'scc', 's_p_value'])+'\n')
    fout.write('\t'.join([str(round(pcc[0], 3)), str(round(pcc[1], 6)), str(round(scc[0], 3)), str(round(scc[1], 6))])+'\n')


def compute_pcc_lncRNA(intersectBed, rbp):
    hdict = {}
    with open(intersectBed) as f:
        for line in f:
            li = line.strip().split('\t')
            chrom, start, end, strand, gid = li[0], li[1], li[2], li[5], li[4]
            chrom1, start1, end1, strand1, gid1 = li[7], li[8], li[9], li[12], li[10]
            signal = float(li[6])
            pred = float(li[11])
            key = ':'.join([chrom, start, end, strand, gid])
            length = min(int(end1), int(end)) - max(int(start), int(start1))
            if key not in hdict.keys() or (hdict[key][0] < length):
                hdict[key] = [length, signal, pred]

    gdict = {}
    for key in hdict.keys():
        gid = key.split(':')[-1]
        if gid not in gdict.keys():
            gdict[gid] = []
        gdict[gid].append((hdict[key][1], hdict[key][2]))

    outfile = 'rbps/{}/lncRNA_pcc.txt'.format(rbp)
    fout = open(outfile, 'w')
    fout.write('\t'.join(['lncRNA', 'pcc', 'p_value', 'scc', 's_value', 'number'])+'\n')
    # y_s = []
    # y_p = []
    for gid in gdict.keys():
        y_true = []
        y_pred = []
        num = len(gdict[gid])
        for signal, pred in gdict[gid]:
            y_true.append(signal)
            y_pred.append(pred)
        # y_s.append(max(y_true))
        # y_p.append(max(y_pred))
        pcc = pearsonr(y_true, y_pred)
        scc = spearmanr(y_true, y_pred)

        if scc[1] > 0.05 or math.isnan(scc[1]) or pcc[1] > 0.05 or math.isnan(pcc[1]):
            continue
        fout.write('\t'.join([gid, str(round(pcc[0], 3)), str(round(pcc[1], 6)),
                              str(round(scc[0], 3)), str(round(scc[1], 6)), str(num)]) + '\n')

    # pcc = pearsonr(y_s, y_p)
    # scc = spearmanr(y_s, y_p)
    # print(round(pcc[0], 3), '\t', round(scc[0], 3), '\t', rbp)


def main(rbp):
    extract_lncRNA_for_validation(rbp)
    creat_input(rbp)
    score(rbp)
    extract_regions(rbp)
    intersectBed = intersect(rbp)
    intersectBed = 'rbps/{}/lncRNA.intersect'.format(rbp)
    # compute_pcc(intersectBed, rbp)
    compute_pcc_lncRNA(intersectBed, rbp)    #max


if __name__ == '__mian__':
    fp = open('rbps.filter.txt')
    rbps = [item.strip('\n').split()[0] for item in fp.readlines()]
    fp.close()

    for rbp in rbps:
        main(rbp)
