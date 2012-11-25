#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import argparse
import random
import bisect
from multiworkers.multiworker import Controller, Worker


class MyWorker(Worker):
    def __init__(self, work_queue, result_queue, current_queue, global_params):
        Worker.__init__(self, work_queue, result_queue, current_queue, global_params)

    def __load_fraction_reads(self, offset, pick_proba):
        with open(self.global_params['infile'].name, 'r') as infile:
            infile.seek(offset)
            iread = Read.fromLine(infile.readline())
            cur_seqid = iread.seqid
            if random.random() < pick_proba:
                yield iread
            for read in (Read.fromLine(l) for l in infile):
                if read.seqid == cur_seqid:
                    if random.random() < pick_proba:
                        yield read
                else:
                    return
            #if random.random() < pick_proba:
                #yield read

    def __load_all_reads(self, offset):
        with open(self.global_params['infile'].name, 'r') as infile:
            infile.seek(offset)
            iread = Read.fromLine(infile.readline())
            cur_seqid = iread.seqid
            yield iread
            for read in (Read.fromLine(l) for l in infile):
                if read.seqid == cur_seqid:
                    yield read
                else:
                    return
            #yield read

    def __load_reads(self, offset):
        pick_proba = self.global_params['fraction']
        if pick_proba == 1.0:
            return self.__load_all_reads(offset)
        else:
            return self.__load_fraction_reads(offset, pick_proba)

    def __ascending_index(self, feats):
        features = feats[:]
        features.sort(key=lambda x: x.start)
        return features

    def __descending_index(self, feats):
        features = feats[:]
        features.sort(key=lambda x: x.end)
        return features

    def __features_encompassing_read(self, read, findex, rindex):
        up_index = bisect.bisect_right([e.start for e in findex], read.start)
        down_index = bisect.bisect_left([e.end for e in rindex], read.end)
        up_features = set(findex[:up_index])
        down_features = set(rindex[down_index:])
        return up_features & down_features

    def do(self, job):
        current_seqid_annot = self.global_params['annotation'].get(job['seqid'], [])
        findex = self.__ascending_index(current_seqid_annot)
        rindex = self.__descending_index(current_seqid_annot)
        for read in self.__load_reads(job['offset']):
            matching_features = self.__features_encompassing_read(read, findex, rindex)
            if self.global_params['exclude_ambiguous'] is False or \
                    (self.global_params['exclude_ambiguous'] and len(matching_features) < 2):
                for feature in matching_features:
                    feature.add(read, self.global_params['norm'])
        return {
            'seqid': job['seqid'],
            'offset': job['offset'],
            'counts': current_seqid_annot
        }


class MyController(Controller):
    def __init__(self, jobs, global_params, num_cpu=1, quiet=False,
                 worker_class=MyWorker, debug=False):
        jobs = (job for job in self.index_infile(global_params['infile']))
        Controller.__init__(self, jobs, global_params, num_cpu=num_cpu,
                            quiet=quiet, worker_class=worker_class,
                            debug=debug)

    def index_infile(self, infile):
        infile.seek(0)
        offset = 0
        line = infile.readline()
        cur_seqid = line.strip().split('\t')[0]
        yield {'seqid': cur_seqid, 'offset': offset, 'length': 0}
        while line:
            offset = infile.tell()
            line = infile.readline()
            if line:
                l = line.strip().split('\t')
                if l[0] != cur_seqid:
                    yield {'seqid': l[0], 'offset': offset, 'length': 0}
                    cur_seqid = l[0]
        infile.seek(0)

    def finish(self):
        output = ProfilesCollection()
        for result in self.results:
            for feature in result['counts']:
                output.append(feature)
        output_ids = [f.attributes['ID'] for f in output]
        # output features on non-processed seqids
        for seqid, features in self.global_params['annotation'].iteritems():
            for feature in features:
                if feature.attributes['ID'] not in output_ids:
                    output.append(feature)
        output.sort(key=lambda x: x.attributes['ID'])
        print >>self.global_params['outfile'], output


class BaseFeature:
    def __init__(self, seqid, start, end, strand):
        self.seqid = seqid
        self.start = min(start, end)
        self.end = max(start, end)
        self.strand = strand

    def __hash__(self):
        return hash((self.seqid, self.start, self.end, self.strand))

    def __eq__(self, other):
        return self.seqid == other.seqid and \
                self.start == other.start and \
                self.end == other.end and \
                self.strand == other.strand

    def __str__(self):
        return '\t'.join(map(str, [
            self.seqid,
            self.start,
            self.end,
            self.strand,
        ]))

    @classmethod
    def fromLine(cls, line):
        l = line.strip().split('\t')
        try:
            seqid, start, end, strand = l
        except ValueError:
            seqid, start, end, strand = l[:3] + ['=']
        return cls(seqid, int(start), int(end), strand)
    
    @property
    def length(self):
        return abs(self.end - self.start)


class Read(BaseFeature):
    def __init__(self, seqid, start, end, strand):
        BaseFeature.__init__(self, seqid, start, end, strand)
        self.midpos = self.start + self.length / 2


class Feature(BaseFeature):
    def __init__(self, seqid, source, method, start, end, score, strand, phase,
                 attributes):
        BaseFeature.__init__(self, seqid, start, end, strand)
        self.source = source
        self.method = method
        self.phase = phase
        self.attributes = attributes
        self.score = 0.0

    def __str__(self):
        return '\t'.join(map(str, [
            self.seqid,
            self.source,
            self.method,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.phase,
            self.formatAttributes(),
        ]))

    @classmethod
    def fromLine(cls, line):
        l = line.strip().split('\t')
        attributes = {}
        if len(l) > 8:
            a = l[8].split(';')
            for i in a:
                k, v = i.split('=')
                attributes[k] = v
        return cls(l[0], l[1], l[2], int(l[3]), int(l[4]), l[5], l[6], l[7], attributes)

    def formatAttributes(self):
        return ';'.join(['='.join([k, v]) for k, v in self.attributes.iteritems()])

    def add(self, read, norm):
        self.increment(read, norm)

    def increment(self, read, norm):
        norm_factor = self.length if norm else 1
        self.score += 1.0 / norm_factor


class ProfilesCollection(list):
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)

    def __str__(self):
        return '\n'.join([str(i) for i in self])


def parseAnnotation(filename):
    annotation = {}
    with open(filename, 'r') as infile:
        for feature in (Feature.fromLine(line) for line in infile):
            if feature.seqid in annotation:
                annotation[feature.seqid].append(feature)
            else:
                annotation[feature.seqid] = [feature]
    return annotation


def main(args):
    controller = MyController(
        jobs={},
        global_params={
            'infile': args.infile,
            'outfile': args.outfile,
            'annotation': parseAnnotation(args.annotation),
            'norm': args.norm,
            'fraction': args.fraction,
            'exclude_ambiguous': args.exclude_ambiguous
        },
        num_cpu=args.num_cpu,
        quiet=args.quiet,
        worker_class=MyWorker,
        debug=args.debug
    )
    controller.start()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--infile', dest='infile',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help="Input file containing reads in bed format (seqid, start, stop)"
    )
    parser.add_argument(
        '-o', '--outfile', dest='outfile',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output file"
    )
    parser.add_argument(
        '-a', '--annotation', dest='annotation',
        help="Annotation file in GFF3 format"
    )
    parser.add_argument(
        '-N', '--num_cpu', dest='num_cpu',
        default=1,
        type=int,
        help="Number of parallel jobs to run"
    )
    parser.add_argument(
        '-m', '--exclude_ambiguous', dest='exclude_ambiguous',
        action='store_true',
        default=False,
        help='Should reads overlapping multiple features be excluded.'
    )
    parser.add_argument(
        '-t', '--overlap_method', dest='overlap_method',
        choices=['full', 'partial'],
        default='full',
        help='Method to use for determining overlap of a read with a feature'
    )
    parser.add_argument(
        '-f', '--fraction', dest='fraction',
        type=float,
        default=1.0,
        help='Fraction of the reads to count'
    )
    parser.add_argument(
        '-n', '--norm', dest='norm',
        action='store_true',
        default=False,
        help="normalize by size"
    )
    parser.add_argument(
        '-q', '--quiet', dest='quiet',
        action='store_true',
        default=False,
        help='Be quiet'
    )
    parser.add_argument(
        '-d', '--debug', dest='debug',
        action='store_true',
        default=False,
        help='Debug mode'
    )
    main(parser.parse_args())
