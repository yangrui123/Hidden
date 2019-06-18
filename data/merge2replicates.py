#!/usr/bin/env python
# coding:utf-8
__author__ = 'yangrui'

import os


def combine2peakfiles(peak1, peak2, outfile):
    fout = open(outfile, 'w')
    hdict = {}
    with open(peak2) as f:
        for line in f:
            li = line.strip().split('\t')
            chrom, start, end, pvalue, signal, strand = li[0], li[1], li[2], li[3], li[4], li[5]
            key = ':'.join([chrom, start, end, strand])
            hdict[key] = float(signal)
    with open(peak1) as f:
        for line in f:
            li = line.strip().split('\t')
            chrom, start, end, pvalue, signal, strand = li[0], li[1], li[2], li[3], li[4], li[5]
            key = ':'.join([chrom, start, end, strand])
            signal = (float(signal) + hdict[key]) / 2
            fout.write('\t'.join([chrom, start, end, '', '', strand, str(signal), pvalue]))
            fout.write('\n')
    return outfile


def sortbed(bed, outfile):
    cmd = 'bedtools sort -i {bed} > {out}'.format(bed=bed, out=outfile)
    print(cmd)
    os.system(cmd)
    return outfile


def intersectbeds(bed1, bed2, outfile):
    cmd = 'bedtools intersect -a {bed1} -b {bed2} -wa -wb > {out}'.format(bed1=bed1, bed2=bed2, out=outfile)
    print(cmd)
    os.system(cmd)
    return outfile


def extractregion(bedfile, outfile):
    fout = open(outfile, 'w')
    idx = 1
    with open(bedfile) as f:
        for line in f:
            li_1 = line.strip('\n').split('\t')[:8]
            li_2 = line.strip('\n').split('\t')[10:-2]
            start1, end1 = int(li_1[1]), int(li_1[2])
            start2, end2 = int(li_2[1]), int(li_2[2])
            start = max(start1, start2)
            end = min(end1, end2)
            chrom, strand, pvalue = li_1[0], li_1[5], li_1[7]
            name = 'mergePeaks{}'.format(idx)
            idx += 1
            fout.write('\t'.join([chrom, str(start), str(end), name, pvalue, strand]))
            fout.write('\n')

    fout.close()

    return outfile


def bamcount(bam, outfile):
    cmd = 'samtools view -c -F 4 {bam} > {outfile}'.format(bam=bam, outfile=outfile)
    print(cmd)
    os.system(cmd)
    return outfile


def computesignal(bedfile, expbam, mockbam, readcount, mockcount, outfile):
    cmd = 'perl overlap_peakfi_with_bam_PE.pl {expbam} {inptbam} {peakfile} {readcount} {mockcount} {outfile}'.format(
        expbam=expbam, inptbam=mockbam, peakfile=bedfile, outfile=outfile, readcount=readcount, mockcount=mockcount)
    print(cmd)
    os.system(cmd)
    return outfile


def merge2peaks(bed1, bed2, bam1, bam2, mockbam, prefix):
    # sort two replicates
    bed1 = sortbed(bed1, prefix + '_1.sorted.bed')
    bed2 = sortbed(bed2, prefix + '_2.sorted.bed')

    # two replicates intersect
    bedfile = intersectbeds(bed1, bed2, prefix + '.intersect.bed')

    bedfile = extractregion(bedfile, prefix + '.merge.bed')

    bedfile = sortbed(bedfile, prefix + '.sort.merge.bed')

    # count reads number
    rep1_count = bamcount(bam1, prefix + '.exp.read_count1')
    rep2_count = bamcount(bam2, prefix + '.exp.read_count2')
    mock_count = bamcount(mockbam, prefix + '.mock.read_count')
    #
    peak1 = computesignal(bedfile, bam1, mockbam, rep1_count, mock_count, prefix+'.peak1')
    peak2 = computesignal(bedfile, bam2, mockbam, rep2_count, mock_count, prefix + '.peak2')

    peak = combine2peakfiles(peak1, peak2, prefix+'.peak')
    return peak


if __name__ == '__main__':
    bed1 = ''
    bed2 = ''
    bam1 = ''
    bam2 = ''
    mockbam = ''
    prefix = ''
    merge2peaks(bed1, bed2, bam1, bam2, mockbam, prefix)