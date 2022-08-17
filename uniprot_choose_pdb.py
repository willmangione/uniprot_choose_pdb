#!/util/common/python/py37/anaconda-5.3.1/bin/python
from sklearn.cluster import DBSCAN
import numpy as np
import sys
import argparse


parser = argparse.ArgumentParser(description='This script will select the lowest resolution PDB chains corresponding to a query UniProt ID.' 
                                 'If using locally or on a separate cluster, first run'
                                 '"wget http://protinfo.compbio.buffalo.edu/cando/data/v2.2+/mappings/pdb_2_uniprot.csv"'
                                 'and then use the following flag: "-l pdb_2_uniprot.csv" for all subsequent runs.')
parser.add_argument('-u', '--uniprot', help='The query UniProt ID', required=False)
parser.add_argument('-b', '--batch', help='Path to file containing UniProt IDs on each line', required=False)
parser.add_argument('-o', '--out', help='Name of file to save output of batch', required=False)
parser.add_argument('-l', '--local', help='Run with a local mapping (not on CCR)', required=False)
parser.add_argument('-d', '--detailed', help='Output more information about the selections',
                    action='store_true', required=False)
parser.add_argument('-v', '--visual', help='Output a visual of sequence coverage for each chain selected',
                    action='store_true', required=False)
parser.add_argument('-y', '--polyprotein', help='Polyprotein sequences use the lower cutoff of 30 residues',
                    action='store_true', required=False)

args = vars(parser.parse_args())

if args['batch']:
    if not args['out']:
        print('Please enter an out file name with the -o flag')
        quit()


def overlap(x, y, t, ov=False, count=False):
    min1 = x[0]
    max1 = x[1]
    min2 = y[0]
    max2 = y[1]
    if not ov:
        return abs((max(0, min(max1, max2) - max(min1, min2)) - (max1 - min1))) / t
    if ov:
        if count:
            return abs((max(0, min(max1, max2) - max(min1, min2))))
        else:
            return abs((max(0, min(max1, max2) - max(min1, min2)))) / t

def matrix(xs, ln):
    pairs = []
    for x in xs:
        ps = []
        for y in xs:
            ps.append(overlap(x, y, ln))
        pairs.append(ps)
    return pairs


def merge(x, seg, xs, segs):
    i = xs.index(x)
    y = segs[i]
    y1, y2 = y[0], y[1]
    seg1, seg2 = seg[0], seg[1]

    if y1 < seg1:
        if seg1 - y2 != 1:
            #print('cannot merge {}'.format(x))
            return None
        else:
            new = [y1, seg2]
    else:
        if y1 - seg2 != 1:
            #print('cannot merge {}'.format(x))
            return None
        else:
            new = [seg1, y2]
    segs[i] = new
    return xs, segs


#f = open('pdb_chain_uniprot_human.tsv', 'r')
if args['local']:
    f = open(args['local'], 'r')
else:
    f = open('/projects/academic/rams/wmangion/mappings/proteins/pdb_2_uniprot.csv', 'r')
d = {}
lengths = {}
resolu = {}
nmr = []
cluster_sizes = []
li = 0

res_meths = ['X-RAY DIFFRACTION', 'ELECTRON MICROSCOPY', 'NEUTRON DIFFRACTION', 'FIBER DIFFRACTION',
             'ELECTRON CRYSTALLOGRAPHY', 'POWDER DIFFRACTION']
nmr_meths = ['SOLID-STATE NMR', 'INFRARED SPECTROSCOPY', 'SOLUTION NMR', 'THEORETICAL MODEL',
             'SOLUTION SCATTERING']

def check_method(meth_list):
    for m in meth_list:
        if m in res_meths:
            return True
    return False


bad_i = 0
bad_i2 = 0
for l in f.readlines()[1:]:
    #print(li)
    #ls = l.strip().split('\t')
    ls = l.strip().split(',')
    pdb = ls[0] + ls[1]
    uni = ls[2]
    seq_length = float(ls[9])
    lengths[uni] = seq_length
    meth = ls[11]
    if not meth:
        continue
    elif meth in res_meths:
        try:
            resolu[pdb] = float(ls[12])
        except ValueError:
            bad_i += 1
            nmr.append(pdb)
    elif meth[0] == '"':
        i = 0
        curr = meth
        mult = [meth]
        while curr[-1] != '"':
            i += 1
            mult.append(ls[11+i])
            curr = ls[11+i]
        if check_method(mult):
            try:
                resolu[pdb] = float(ls[12+i])
            except:
                print('oh boy', i, ls)
    elif meth in nmr_meths:
        nmr.append(pdb)
    else:
        print('idk', meth)
        continue
    #if uni != uniprot:
    #    continue
    x = int(ls[7])
    y = int(ls[8])
    if uni in d:
        if pdb in d[uni][0]:
            result = merge(pdb, [x, y], d[uni][0], d[uni][1])
            if result == None:
                continue
            else:
                L = result[0]
                Is = result[1]
                d[uni] = [L, Is]
        else:
            d[uni][0].append(pdb)
            d[uni][1].append([x,y])
    else:
        d[uni] = [[pdb], [[x,y]]]
    li += 1


def representative(xs):
    curr = xs[0]
    n = 0
    r = 10000
    for x in xs:
        diff = x[1][1] - x[1][0]
        try:
            res = resolu[x[0]]
            if res < r:
                curr = x
                r = res
            elif res == r:
                if diff > n:
                    curr = x
                    n = diff
        except KeyError:
            if n == 0 and r == 10000:
                n = diff
                curr = x
    cluster_sizes.append((curr[0], len(xs)))
    return curr


def contained(i, xs, allx, thresh=1.0, upper_cutoff=100):
    length = xs[i]
    if length < upper_cutoff:
        return True
    cutoff = length * thresh
    for j in range(len(xs)):
        if i == j:
            continue
        else:
            if xs[j] >= cutoff:
                if length > allx[j][j]:
                    pass
                else:
                    return True
            else:
                pass
    return False


# --> SCRIPT <--


if args['uniprot']:
    uniprot = args['uniprot']

    try:
        length = lengths[uniprot]
        if length == 0:
            print('UniProt length is 0, please check.')
            quit()
    except KeyError:
        print('UniProt ID {} cannot be found in this mapping file.'.format(uniprot))
        quit()

    it3 = matrix(d[uniprot][1], length)
    if args['polyprotein']:
        m = DBSCAN(eps=0.1, min_samples=2)
    else:
        m = DBSCAN(eps=0.4, min_samples=2)
    them = m.fit(np.array(it3))


    d2 = {}

    for li in range(len(them.labels_)):
        cluster = them.labels_[li]
        pro = d[uniprot][0][li]
        rang = d[uniprot][1][li]
        if cluster not in d2:
            d2[cluster] = [(pro, rang)]
        else:
            d2[cluster].append((pro, rang))

    bigs = []
    for i in list(d2.keys()):
        bigs.append(representative(d2[i]))

    filt = []
    for big in bigs:
        f2 = []
        for b in bigs:
            f2.append(overlap(big[1], b[1], 1, ov=True))
        filt.append(f2)
    #print(bigs)

    #for l in filt:
    #    print(l)

    #print(filt)
    #print(bigs)
    #quit()

    segs = []

    if length >= 150 and not args['polyprotein']:
        for li in range(len(filt)):
            if not contained(li, filt[li], filt, thresh=0.5, upper_cutoff=100):
                p = bigs[li][0]
                seg = bigs[li][1]
                try:
                    res = resolu[bigs[li][0]]
                except KeyError:
                    if bigs[li][0] in nmr:
                        print('Warning: no X-Ray Diffraction structure available, using NMR structure')
                    else:
                        print('No resolution information present')
                segs.append(bigs[li])
                #print(p, seg, res)
    else:
        for li in range(len(filt)):
            if not contained(li, filt[li], filt, thresh=0.5, upper_cutoff=30):
                p = bigs[li][0]
                seg = bigs[li][1]
                try:
                    res = resolu[bigs[li][0]]
                except KeyError:
                    if bigs[li][0] in nmr:
                        print('Warning: no X-Ray Diffraction structure available, using NMR structure')
                    else:
                        print('No resolution information present')
                segs.append(bigs[li])
                #print(p, seg, res)

    segs = sorted(segs, key=lambda x:x[1][0])


    def coverage(L, segs):
        sorted_segs = sorted(segs, key=lambda x:x[1][0])
        if len(sorted_segs) == 1:
            return str((((segs[0][1][1] - segs[0][1][0]) + 1) / L) * 100)[0:4] + '%'
        else:
            aa_count = 0
            for si in range(len(sorted_segs)):
                seg = sorted_segs[si][1]
                if si == len(sorted_segs) - 1:
                    aa_count += (seg[1] - seg[0])
                    return str((aa_count / L) * 100)[0:4] + '%'
                next_seg = sorted_segs[si+1][1]
                ov = overlap(seg, next_seg, L, ov=True, count=True)

                if ov == 0:
                    aa_count += (seg[1] - seg[0]) + 1
                else:
                    aa_count += (seg[1] - seg[0]) - ov + 1


    def visual(u, segs):

        def nticks(x, y):
            return int((y-x) / 10)

        def start_end(x, y):
            xd = int(x/10)
            yd = int(y/10)
            return xd+1, yd+1

        length = lengths[u]
        ticks = int(length / 10)
        half = int(((ticks + 2) / 2) - 3)
        spacer = (' ' * half) + u
        print(spacer)
        st = '[' + ('-' * ticks) + ']'
        print(st)

        seg_st = list((' ' * ticks) + '  ')
        mid_st = list((' ' * ticks) + '  ')
        mids = []
        for seg in segs:
            x,y = start_end(seg[1][0], seg[1][1])
            if seg_st[x-1] == ']':
                seg_st[x-1] = '|'
            else:
                seg_st[x-1] = '['
            seg_st[y] = ']'
            for i in range(x, y):
                if seg_st[i] == ' ':
                    seg_st[i] = '-'
            mid = int(((x + y) / 2) - 2)
            mids.append((seg[0], mid))

        for m in mids:
            for i in range(len(m[0])):
                mid_st[i+m[1]] = m[0][i]

        print(''.join(seg_st))
        print(''.join(mid_st))


    if args['detailed']:
        for s in segs:
            for cl in cluster_sizes:
                if s[0] == cl[0]:
                    try:
                        res = resolu[s[0]]
                    except KeyError:
                        res = 'N/A'
                    print(s[0])
                    print('  |--> residue {} to {},'.format(s[1][0], s[1][1]),
                          'resolution = {}, chosen from a cluster of {} chains.'.format(res, cl[1]))
        print('Total sequence coverage: {} (of {} total residues).'.format(coverage(lengths[uniprot],
                                                                                    segs), int(lengths[uniprot])))
    else:
        for s in segs:
            print(s[0])

    if args['visual']:
        visual(uniprot, segs)

elif args['batch']:
    fo = open(args['out'], 'w')
    prots = []
    fb = open(args['batch'], 'r')
    for l in fb:
        prots.append(l.strip())

    for uniprot in prots:
        try:
            length = lengths[uniprot]
            if length == 0:
                fo.write('{}\n'.format(uniprot))
                continue
        except KeyError:
            fo.write('{}\n'.format(uniprot))
            continue

        it3 = matrix(d[uniprot][1], length)
        m = DBSCAN(eps=0.4, min_samples=2)
        them = m.fit(np.array(it3))

        d2 = {}

        for li in range(len(them.labels_)):
            cluster = them.labels_[li]
            pro = d[uniprot][0][li]
            rang = d[uniprot][1][li]
            if cluster not in d2:
                d2[cluster] = [(pro, rang)]
            else:
                d2[cluster].append((pro, rang))

        bigs = []
        for i in list(d2.keys()):
            bigs.append(representative(d2[i]))

        filt = []
        for big in bigs:
            f2 = []
            for b in bigs:
                f2.append(overlap(big[1], b[1], 1, ov=True))
            filt.append(f2)
        # print(bigs)

        # for l in filt:
        #    print(l)

        pdbs = []

        if length >= 150 and not args['polyprotein']:
            for li in range(len(filt)):
                if not contained(li, filt[li], filt, thresh=0.5, upper_cutoff=100):
                    p = bigs[li][0]
                    seg = bigs[li][1]
                    try:
                        res = resolu[bigs[li][0]]
                        pdbs.append(p)
                    except KeyError:
                        if bigs[li][0] in nmr:
                            pdbs.append(p+'*')
                    #    else:
                    #        print('No resolution information present')

                    # print(p, seg, res)
        else:
            for li in range(len(filt)):
                if not contained(li, filt[li], filt, thresh=0.5, upper_cutoff=30):
                    p = bigs[li][0]
                    seg = bigs[li][1]
                    #try:
                    #    res = resolu[bigs[li][0]]
                    #except KeyError:
                    #    if bigs[li][0] in nmr:
                    #        print('Warning: no X-Ray Diffraction structure available, using NMR structure')
                    #    else:
                    #        print('No resolution information present')
                    pdbs.append(p)
                    # print(p, seg, res)

        fo.write('{}\t{}\n'.format(uniprot, '\t'.join(pdbs)))

