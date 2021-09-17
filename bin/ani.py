#!/usr/bin/env python3
"""
Add ANI_report_prokaryotes.tsv tax check data
"""
import argparse
import pandas

ani_header = [
    'assembly_genbank',  # genbank-accession
    'assembly_refseq',  # refseq-accession
    'taxid',
    'species-taxid',
    'organism-name',
    'species-name',
    'assembly-name',
    'assembly-type-category',
    'excluded-from-refseq',
    'declared-type-assembly',
    'declared-type-organism-name',
    'declared-type-category',
    'declared-type-ANI',
    'declared-type-qcoverage',
    'declared-type-scoverage',
    'best-match-type-assembly',
    'best-match-species-taxid',
    'best-match-species-name',
    'best-match-type-category',
    'best-match-type-ANI',
    'best-match-type-qcoverage',
    'best-match-type-scoverage',
    'best-match-status',
    'comment',
    'taxonomy-check-status'
    ]

asm_header = [
    'assembly_genbank',  # assembly_accession
    'bioproject',
    'biosample',
    'wgs_master',
    'refseq_category',
    'taxid',
    'species_taxid',
    'organism_name',
    'infraspecific_name',
    'isolate',
    'version_status',
    'assembly_level',
    'release_type',
    'genome_rep',
    'seq_rel_date',
    'asm_name',
    'submitter',
    'assembly_refseq',  # gbrs_paired_asm
    'paired_asm_comp',
    'ftp_path',
    'excluded_from_refseq',
    'relation_to_type_material',
    'asm_not_live_date'
    ]


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
        header=None,
        na_values='na',
        names=ani_header,
        sep='\t',
        skiprows=[0],  # header row
        usecols=[
            'assembly_genbank',
            'assembly_refseq',
            'best-match-species-name',
            'best-match-species-taxid',
            'taxonomy-check-status']
            )
    ani_map = ani[~ani['assembly_refseq'].isna()]
    ani_map = ani_map[['assembly_genbank', 'assembly_refseq']]
    asm_map = pandas.read_csv(
        args.asm,
        dtype=str,
        header=None,
        names=asm_header,
        sep='\t',
        skiprows=[0, 1],  # comment and header rows
        usecols=['assembly_genbank', 'wgs_master'])
    asm_map = asm_map[~asm_map['wgs_master'].isna()]
    asm_map['wgs_prefix'] = asm_map['wgs_master'].apply(lambda x: x[:7])
    info['wgs_prefix'] = info['accession'].apply(lambda x: x[:7])
    info = info.merge(ani_map, how='left')
    info = info.merge(asm_map, how='left', on='wgs_prefix', suffixes=['', '_'])
    info['assembly_genbank'] = info['assembly_genbank'].fillna(
        info['assembly_genbank_'])
    drop = ['assembly_genbank_', 'assembly_refseq', 'wgs_prefix', 'wgs_master']
    info = info.drop(drop, axis='columns')
    info = info.merge(ani, how='left')
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
