import csv

import numpy as np
from sklearn.cluster import KMeans


def load_data(filename):
    d = []
    with open(filename, 'r') as f:
        csvreader = csv.reader(f,
                               delimiter='\t',
                               quoting=csv.QUOTE_NONE,
                               lineterminator='\n')
        csvreader.next()
        for row in csvreader:
            d.append(map(float, row[1:]))

    return np.array(d)


def parse_argv(argv):
    # use sys.argv[1:] to substitute argv
    import argparse

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('item')
    parser.add_argument('out')
    parser.add_argument('ubout')
    parser.add_argument('-i', '--maxiter', default=100, type=int)
    parser.add_argument('-n', '--numclusters', type=int, default=1)
    parser.add_argument('-m', '--numiters', type=int, default=1)
    args = parser.parse_args(argv)
    return args


if __name__ == '__main__':
    import sys
    args = parse_argv(sys.argv[1:])
    data = load_data(args.item)
    estimator = KMeans(init='k-means++',
                       max_iter=args.numiters,
                       n_clusters=args.numclusters,
                       n_init=args.numiters)
    estimator.fit(data)
    with open(args.out, 'w') as f:
        f.write('{0}\t{1}\n'.format(*(estimator.cluster_centers_.shape)))
        for cluster in estimator.cluster_centers_:
            f.write('\t'.join(map(str, cluster)))
            f.write('\n')

    # upper bound
    m = {}
    for idx, l in enumerate(estimator.labels_):
        m.setdefault(l, []).append(data[idx])

    with open(args.ubout, 'w') as f:
        f.write('{0}\t{1}\n'.format(*(estimator.cluster_centers_.shape)))
        for l in m:
            ub = np.max(m[l], axis=0)
            f.write('{0}\t'.format(sum(ub)))
            f.write('\t'.join(map(str, ub)))
            f.write('\n')
