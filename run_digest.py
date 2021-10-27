import os
import csv
import itertools
import argparse

import numpy as np
from tabulate import tabulate
from Bio import Restriction as biorest
from Bio import Seq as bioseq

# TODO: Add CLI
# TODO: Add ability to load any file
# TODO: Update enzyme list
# TODO: Display information about buffers


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--single', '-s',
        action='store_true',
        help='Include single enzyme digests'
    )
    parser.add_argument(
        '--double', '-d',
        action='store_true',
        help='Include double enzyme digests'
    )
    parser.add_argument(
        '--num', '-n',
        default=10,
        type=int,
        help='Number of results to display'
    )
    parser.add_argument(
        '--file', '-f',
        default='plasmid.txt',
        type=str,
        help='Name of file that contains plasmid sequence'
    )

    return parser.parse_args()


def score_bands(band_list):
    # Taking the old cost function w/o modification
    cost = 0
    if len(band_list) == 0 or len(band_list) == 1:
        return np.inf
    elif len(band_list) > 6:
        cost += 5 * len(band_list)
    elif len(band_list) > 3:
        cost += len(band_list)
    elif len(band_list) == 3:
        cost += -1

    for band in band_list:
        if band >= 10000:
            cost += 100
        elif band >= 7000:
            cost += 50
        elif band >= 6000:
            cost += 25
        elif band <= 500:
            cost += 50

    spacing = np.ediff1d(bands)
    for n, sp in enumerate(spacing):
        if sp < 150:
            cost += 50
        elif band_list[n] > 5000 and sp < 1000:
            cost += 50

    return cost

def _round(nums):
    # TODO: Dirty way of returning list with descending band size
    return [int(round(n, -2)) for n in nums][::-1]


if __name__ == '__main__':
    # Load all restriction enzymes
    args = parse_args()

    if not args.single and not args.double:
        # If both are false, just do default behavior
        args.single = True
        args.double = True

    with open('enzymes.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        enzdic = dict(reader)
        enz_names = list(enzdic.keys())

    # Make restritcion objects for all single and double enzymes
    enzymes = []
    for enz in enz_names:
        enz = enz.split('_HF')[0]
        att = getattr(biorest, enz, None)
        if att is not None: enzymes.append(att)
    singles = [biorest.RestrictionBatch([i]) for i in enzymes]
    doubles = [biorest.RestrictionBatch(i) for i in itertools.combinations(enzymes, 2)]
    enzymes = singles + doubles

    # Load the plasmid
    with open(args.file, 'r') as f_open:
        plasmid_sequence = f_open.read()
    plasmid_sequence = plasmid_sequence.strip()
    plasmid_sequence = plasmid_sequence.upper()
    plasmid_sequence = bioseq.Seq(plasmid_sequence)

    digest = {}
    for enz in enzymes:
        ana = biorest.Analysis(enz, plasmid_sequence, linear=False)

        # For double enzymes, skip if one enzyme doesn't cut
        if any([len(v) == 0 for v in ana.full().values()]):
            continue

        cuts = sorted([sl for v in ana.full().values() for sl in v])
        if len(cuts) > 0:
            # The band that crosses the origin needs to be handled separately
            ori_band = len(plasmid_sequence) - cuts[-1] + cuts[0]
            oth_bands = list(np.diff(cuts))
            bands = sorted([ori_band] + oth_bands)
        digest[str(enz)] = (score_bands(bands), bands)

    digest = {k: v for k, v in sorted(digest.items(), key=lambda x: x[1])}

    # Display results
    display = args.num
    headers = ['Enzyme', 'Pattern', 'Score']
    print(f'Plasmid length: {len(plasmid_sequence)}')
    if args.single:
        single = {k: v for k, v in digest.items() if '+' not in k}
        single_table = [tuple([k, _round(v[1]), v[0]])
                        for n, (k, v) in enumerate(single.items())
                        if n < display]
        print('\n Single enzymes')
        print(tabulate(single_table, headers))

    if args.double:
        double = {k: v for k, v in digest.items() if '+' in k}
        double_table = [tuple([k, _round(v[1]), v[0]])
                        for n, (k, v) in enumerate(double.items())
                        if n < display]
        print('\n Double enzymes')
        print(tabulate(double_table, headers))
