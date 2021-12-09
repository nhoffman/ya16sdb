#!/usr/bin/env python3
"""
Add ANI_report_prokaryotes.tsv tax check data
"""
import argparse
import pandas
import ya16sdb


def extract_wgs_master_prefix(accession):
    """
    Adapted from https://gitlab.labmed.uw.edu/molmicro/reference-packages/
    plotting/-/blob/master/plotting/subcommands/get_ref_info.py#L153

    given an ncbi identifier, use a regex to
    return the 4 or 6 alphabetical code corresponding to
    NCBI convention.
    See https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
    """

    # remove 'N[Z,R]_' etc. refseq prefixes from accessions
    # see https://www.ncbi.nlm.nih.gov/books/NBK21091/table/
    # ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
    accession = accession.split("_", 1)[-1]
    match = ya16sdb.WGS_PREFIX.search(accession)
    return match.group() if match else pandas.NA


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    p.add_argument('ani', help='ANI_report_prokaryotes.txt')
    p.add_argument('asm', help='assembly_summary_genbank.txt')
    args = p.parse_args()
    info = pandas.read_feather(args.feather)
    ani = pandas.read_csv(
        args.ani,
        dtype=str,
        na_values='na',
        sep='\t',
        usecols=[
            '# genbank-accession',
            'refseq-accession',
            'best-match-type-assembly',
            'best-match-species-name',
            'best-match-species-taxid',
            'taxonomy-check-status',
            'best-match-type-category',
            'best-match-type-ANI',
            'best-match-type-qcoverage',
            'declared-type-ANI',
            'declared-type-qcoverage'
            ])
    ani = ani.rename(
        columns={
            '# genbank-accession': 'assembly_genbank',
            'refseq-accession': 'assembly_refseq'})
    # create our two map files, one from the assembly_refseq
    ani_map = ani[~ani['assembly_refseq'].isna()]
    ani_map = ani_map[['assembly_refseq', 'assembly_genbank']]
    asm_map = pandas.read_csv(
        args.asm,
        dtype=str,
        header=1,
        sep='\t',
        usecols=['wgs_master', '# assembly_accession'])
    asm_map = asm_map.rename(
        columns={'# assembly_accession': 'assembly_genbank'})
    asm_map = asm_map[~asm_map['wgs_master'].isna()]
    asm_map['wgs_prefix'] = asm_map['wgs_master'].apply(
        lambda x: extract_wgs_master_prefix(x))
    info['wgs_prefix'] = info['accession'].apply(
        lambda x: extract_wgs_master_prefix(x))
    # create two assembly_genbank columns prioritizing
    # assembly_genbank accessions from records with a wgs_prefix
    # that can be cross-referenced with the assembly_summary_genbank.txt file
    info = info.merge(asm_map, how='left', on='wgs_prefix')
    info = info.merge(
        ani_map, how='left', on='assembly_refseq', suffixes=['', '_'])
    # if a record can not be cross-referenced with assembly_summary_genbank.txt
    # then just use the assembly_refseq accession found in the genbank file
    info['assembly_genbank'] = info['assembly_genbank'].fillna(
        info['assembly_genbank_'])
    # cleanup
    drop = ['assembly_genbank_', 'assembly_refseq', 'wgs_prefix', 'wgs_master']
    info = info.drop(drop, axis='columns')
    # Now that we have our assembly_genbank column ready we
    # can merge in the rest of the ANI columns
    info = info.merge(ani, how='left', on='assembly_genbank')
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
