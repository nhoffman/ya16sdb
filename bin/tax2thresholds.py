#!/usr/bin/env python3
"""Create a full taxonomic table of thresholds based on a csv of defaults.

Child tax_ids not specified in the defaults receive the parent threshold
"""
import argparse
import pandas
import sys


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'taxonomy',
        help='must have tax_id column and rank columns')
    parser.add_argument(
        'thresholds',
        help='with required columns tax_id, low, target_rank')
    parser.add_argument(
        '-o', '--out',
        default=sys.stdout,
        metavar='FILE')
    args = parser.parse_args()

    taxonomy = pandas.read_csv(
        args.taxonomy,
        comment='#',
        dtype=str,
        na_filter=True,  # False is faster
        ).set_index('tax_id')

    default_thresholds = pandas.read_csv(
        args.thresholds,
        comment='#',
        dtype=dict(tax_id=str, low=float, target_rank=str),
        na_filter=True,
        usecols=['tax_id', 'low', 'target_rank']
        )

    tax_cols = taxonomy.columns.tolist()
    tax_cols = tax_cols[tax_cols.index('root'):]

    # out output data structure
    full_tree = pandas.DataFrame(index=taxonomy.index)

    # start with the root column and move right
    for index, rank in enumerate(tax_cols):
        defaults = default_thresholds[
            default_thresholds['target_rank'] == rank]
        defaults = defaults.set_index('tax_id')

        # iterate taxonomy most to least specificity and join with defaults
        target_thresholds = []
        for specificity in reversed(tax_cols):
            thresholds = taxonomy[[specificity]].join(
                defaults[['low']], on=specificity, how='inner')
            # append just the low column
            target_thresholds.append(thresholds['low'])

        # concat and take just the first(), most specific threshold
        target_thresholds = pandas.concat(target_thresholds)
        target_thresholds = target_thresholds.groupby(
            target_thresholds.index, sort=False).first()

        full_tree[rank] = target_thresholds

        # fill in tax holes
        index = max(index-1, 0)
        blanks = full_tree[full_tree[rank].isnull()]
        full_tree[rank] = full_tree[rank].fillna(
            blanks[tax_cols[index]])

    full_tree.to_csv(args.out)


if __name__ == '__main__':
    main()
