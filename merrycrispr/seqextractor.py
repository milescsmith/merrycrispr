#!/usr/bin/env python3

import click
import pyfaidx
import pandas as pd

from gtfparse import read_gtf
from typing import List
from copy import deepcopy


@click.command()
@click.option("--library_type", "-l",
              help="target library type.  Accepted values are 'knockout', 'repressor/activator', 'excision', or 'Cas13'",
              default=None,
              required=False,
              type=str)
@click.option("--gtf", "-g",
              help="input GTF/GFF file",
              default=None,
              required=False,
              type=str)
@click.option("--fasta", "-f",
              help="FASTA sequence file",
              default=None,
              type=str)
@click.option("--output", "-o",
              help="output file",
              default=None,
              type = str)
@click.option("--feature_type", "-t",
              help="feature type",
              default=None)
@click.option("--gene_name", "-n",
              help="gene(s) to extract",
              default=None,
              type=str)
@click.option("--bound", help="Retrieve a given number of bases on either side of the feature "
                              "instead of the sequence corresponding to a feature",
              default=None,
              type=int)
@click.option("--show_features", help="Scan a GFF file to identify the features present",
              is_flag=True)
@click.option("--show_genes",
              help="Scan a GFF file to identify the genes present",
              default=False,
              is_flag=True)
@click.option("--show_geneids", help="Scan a GFF file to identify the geneIDs present",
              default=False,
              is_flag=True)
@click.help_option()
def main(gtf, library_type, fasta, output, feature_type, gene_name, bound, show_features, show_genes, show_geneids) -> None:
    """A utility for extracting sequences from a FASTA file for a given GFF annotation"""

    if gene_name:
        gene_name = gene_name.split()
    if library_type:
        if library_type == "knockout":
            feature_type = "exon"
            extract(gtffile=gtf,
                    fastafile=fasta,
                    feature_type=feature_type,
                    outfile=output,
                    gene_name=gene_name)
        elif library_type == "repressor/activator":
            feature_type = "exon"
            if bound is None:
                bound = 100
            extract_for_tss_adjacent(gtffile=gtf,
                    fastafile=fasta,
                    feature_type=feature_type,
                    outfile=output,
                    gene_name=gene_name,
                    boundary=bound)
        elif library_type == "excision":
            if bound is None:
                bound = 100
            extract(gtffile=gtf,
                    fastafile=fasta,
                    feature_type=feature_type,
                    outfile=output,
                    gene_name=gene_name,
                    boundary=bound)
        elif library_type == "Cas13":
            extract(gtffile=gtf,
                    fastafile=fasta,
                    feature_type="CDS",
                    outfile=output,
                    gene_name=gene_name)

    elif show_features:
        if gtf:
            display_gtf_features(gtf)
        else:
            print("Please enter the name of the GFF which you wish to scan for features")
    elif show_genes:
        if gtf:
            display_gtf_genes(gtf, feature_type)
        else:
            print("Please enter the name of the GFF which you wish to scan for features")
    elif show_geneids:
        if gtf:
            display_gtf_geneids(gtf, feature_type)
        else:
            print("Please enter the name of the GFF which you wish to scan for features")
    elif not output:
        print("Please enter an output file name")
    elif not fasta:
        print("Please enter the name of the file containing matching sequences (in FASTA format)")
    elif not gtf:
        print("Please enter the name of the file containing features (in GFF format)")
    else:
        extract(gtffile=gtf,
                fastafile=fasta,
                feature_type=feature_type,
                outfile=output,
                gene_name=gene_name.split(),
                boundary=int(bound))


def extract(gtffile: str,
            fastafile: str,
            feature_type: str,
            outfile: str,
            gene_name: List[str],
            boundary: int = 0,
            **kwargs) -> None:

    print("Parsing GTF/GFF file.")
    records = read_gtf(filepath_or_buffer=gtffile, 
                       chunksize=8192 * 1024)
    if gene_name:
        records = records[records['gene_name'].isin(gene_name)]
    if feature_type:
        records = records[records['feature'] == feature_type]
        records = records[['seqname', 'feature', 'start', 'end', 'strand', 'frame', 'gene_name', f'{feature_type}_id']].drop_duplicates()
    else:
        records = records[['seqname', 'feature', 'start', 'end', 'strand', 'frame', 'gene_name']].drop_duplicates()

    if "exon_number" in kwargs.keys():
        records = records[records['exon_number'] == kwargs['exon_number']]

    print(f"{len(records)} total records found.")

    print(f"Loading the sequences in {fastafile}."
          f"Note: if this is the first time opening this file, "
          "it may take a few moments as an index is built.")
    sequences = pyfaidx.Fasta(fastafile)
    print(f"Finished loading {fastafile}")

    if boundary > 0:
        records = split_record(records, boundary)
    final_list = [match_seq(_, sequences) for _ in records.itertuples()]

    with open(outfile, 'w') as o_file:
        for entry in final_list:
            o_file.writelines(f"> {entry.fancy_name}\n{entry.seq}\n")


def display_gtf_features(gtffile: str) -> None:
    gtf = read_gtf(gtffile)
    feature_set = set(gtf['feature'])

    print(f"{len(feature_set)} features found.  These include:")
    [print(_) for _ in feature_set]


def display_gtf_genes(gtffile: str,
                      feature_type: list = None) -> None:
    gtf = read_gtf(gtffile)

    if feature_type is not None:
        gtf = gtf[gtf.feature == feature_type]

    gene_set = set(gtf['gene_name'])

    print(f"{len(gene_set)} genes found.  These include:")
    [print(_) for _ in gene_set]


def display_gtf_geneids(gtffile: str,
                        feature_type: list = None) -> None:
    gtf = read_gtf(gtffile)

    if feature_type is not None:
        gtf = gtf[gtf.feature == feature_type]

    gene_set = set(gtf['gene_id'])

    print(f"{len(gene_set)} genes found.  These include:")
    [print(_) for _ in gene_set]


def extract_for_tss_adjacent(gtffile: str,
            fastafile: str,
            outfile: str,
            gene_name: str,
            boundary: int = 100,
            **kwargs) -> None:

    # read and parse in GTF
    print("Parsing GTF/GFF file.")
    records = read_gtf(gtffile)
    if gene_name:
        records = records[records['gene_name'].isin(gene_name)]

    # break up the genome into forward and reverse strands
    # for forward strand genes, we want the lowest coordinate for exon 1 of each gene
    # for reverse strand genes, we want the highest
    forward_starts = records[records['strand'] == "+"]. \
        groupby('gene_name'). \
        apply(lambda x: x[pd.to_numeric(x['exon_number']) == 1]). \
        reset_index(drop=True). \
        groupby('gene_name'). \
        apply(lambda y: y.nsmallest(1, "start")). \
        reset_index(drop=True)

    reverse_starts = records[records['strand'] == "-"]. \
        groupby('gene_name'). \
        apply(lambda x: x[pd.to_numeric(x['exon_number']) == 1]). \
        reset_index(drop=True). \
        groupby('gene_name'). \
        apply(lambda y: y.nlargest(1, "start")). \
        reset_index(drop=True)

    predicted_tss = pd.concat([forward_starts, reverse_starts])
    predicted_tss = predicted_tss[
        ['seqname', 'feature', 'start', 'end', 'strand', 'frame', 'gene_name']].drop_duplicates()

    print(f"{len(predicted_tss)} total records found.")

    print(f"Loading the sequences in {fastafile}."
          f"Note: if this is the first time opening this file, "
          "it may take a few moments as an index is built.")
    sequences = pyfaidx.Fasta(fastafile)
    print(f"Finished loading {fastafile}")

    # for our list of predicted start sites, extract +/- an interval surrounding the TSS
    # what we want is: 
    #  -boundary --- start --- +boundary
    # easiest way to do this is to make our end coordinate the start coordinate plus the boundary
    # value and make the new start the original start - boundary 
    predicted_tss["gene_name"] += "_TSS"
    predicted_tss['end'] = predicted_tss['start'] + boundary
    predicted_tss['start'] = predicted_tss['start'] - boundary
    final_list = [match_seq(_, sequences) for _ in predicted_tss.itertuples()]

    with open(outfile, 'w') as o_file:
        for entry in final_list:
            o_file.writelines(f"> {entry.fancy_name}\n{entry.seq}\n")


def match_seq(rec: pd.Series, 
              sequences: pyfaidx.Fasta) -> pyfaidx.Sequence:
    try:
        if rec['strand'] == "-":
            rev = True
        else:
            rev = False
        seq = pyfaidx.Sequence(name=f"{rec['gene_name']}_{rec['feature']}_"
                                    f"{rec['strand']}_{rec['start']}_{rec['end']}",
                               seq=sequences.get_seq(name=rec['seqname'],
                                                     start=rec['start'],
                                                     end=rec['end'],
                                                     rc=rev).seq)
    except ValueError:
        print(f"problem with {rec['gene_name']} {rec['start']} "
              f"{rec['end']} {rec['seqname']} {rec['strand']}")
    return seq


def split_record(rec: pd.DataFrame, 
                 padding: int,
                 rev: bool=False) -> pd.DataFrame:
    rec1 = deepcopy(rec)
    rec2 = deepcopy(rec)
    if rev == False:
        rec1["end"] = rec1["start"]
        rec1["start"] -= padding
        rec1["gene_name"] += "_upstream"
        rec2["start"] = rec2["end"]
        rec2["end"] += padding
        rec2["gene_name"] += "_downstream"
    elif rev == True:
        rec1["end"] = rec1["start"]
        rec1["start"] -= padding
        rec1["gene_name"] += "_downstream"
        rec2["start"] = rec2["end"]
        rec2["end"] += padding
        rec2["gene_name"] += "_upstream"
    new_rec = pd.concat([rec1, rec2])
    return new_rec


if __name__ == "__main__":
    main()