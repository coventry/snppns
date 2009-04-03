import sys, os, cPickle, csv, tempfile
from collections import defaultdict as ddict
import SNPPosterior

numcols, subpath = sys.argv[1:3]
numcols = int(numcols)

subject_order = [l.strip() for l in open(subpath)]

callsuffix = '.pickle.results'
callfiles = sys.argv[3:]
assert all(p.endswith(callsuffix) for p in callfiles)

def colbounds(callpath):
    callpath = os.path.basename(callpath)
    boundstring = callpath.split('.')[-3]
    bounds = map(int, boundstring.split('-'))
    refseq = callpath[:-len('.'+boundstring+callsuffix)]
    return refseq, bounds

callpaths = ddict(list)
for path in callfiles:
    refseq, bounds = colbounds(path)
    callpaths[refseq].append((bounds, path))

def col_iterator(boundspaths):
    for bounds, path in sorted(boundspaths):
        start = bounds[0]
        for idx, column in enumerate(cPickle.load(open(path))):
            yield idx+start, column

def col_group_iterator(boundspaths):
    col_i = col_iterator(boundspaths)
    while 1:
        rv = []
        for colidx, r in enumerate(col_i):
            rv.append(r)
            if colidx == numcols-1:
                break
        if rv:
            yield rv[0][0], rv[-1][0], rv
        else:
            raise StopIteration

report_iterators = dict(
    (n, col_group_iterator(bp)) for n, bp in callpaths.items())

def report_dist(dist):
    '''Sort the genotypes by  decreasing probability, and report in
    that order those with probability greater than 1%'''
    sorted_dist = sorted(zip(dist, SNPPosterior.pcalls), reverse=True)
    return [' %s %2i' % (genotype, int(100*prob))
            if prob >= 0.01 else 6*' '
            for prob, genotype in sorted_dist]

for refseq, groups in report_iterators.items():
    for start, end, group in groups:
        csvoutput = ddict(list)
        for idx, column in group:
            # Prior  is  population-level  frequencies, dists  is  the
            # posterior  distributions for  the individual  samples at
            # this site.
            if column[0] is None:
                # There were no calls; no work to do
                continue
            prior, dists = column
            # First row is the population-level distribution
            rows = [[l] for l in report_dist(prior)]
            # Now do the sample-level distributions
            for dist in dists:
                if dist is None:
                    rows[0].append('None  ')
                    for row in rows[1:]:
                        row.append(6*' ')
                else:
                    dist *= prior
                    dist /= sum(dist)
                    for newrow, row in zip(report_dist(dist), rows):
                        row.append(newrow)
            csvoutput['Col Idx'].append([idx])
            for name, entry in zip(['Prior'] + subject_order,
                                   zip(*rows)):
                csvoutput[name].append(
                    filter(None, (l.strip() for l in entry)))
        outpath = tempfile.mktemp(dir='.')
        final_path = '%s.%i-%i.csv.gz' % (refseq, start, end+1)
        cmd = 'gzip -c > ' + outpath
        csvwriter = csv.writer(os.popen(cmd, 'w'), dialect='excel')
        if csvoutput:
            for name in (['Col Idx', 'Prior'] + subject_order):
                row = csvoutput[name]
                for rowidx in range(max(map(len, row))):
                    crow = []
                    crow.append(name if rowidx == 0 else '')
                    for entry in row:
                        crow.append(entry[rowidx] if len(entry) > rowidx else '')
                    csvwriter.writerow(crow)
        else:
            csvwriter.writerow(['No calls'])
        del csvwriter
        os.rename(outpath, final_path)
