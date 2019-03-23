""" GFFextractor.py: a module for extracting sequences matching
a certain type from a FASTA file using information
from a GFF file and writing it out as a FASTA
"""
# !/usr/bin/env python

"""GFFextractor.

Usage:
    GFFextractor.py (-g | --gtf) <gtf> (-f | --fasta) <fasta> (-o | --output) <output> [-t | --feature_type <feature>] [-n | --gene_name <gene>] [-b | --bound <boundary>]
    GFFextractor.py --show_features <features>
    GFFextractor.py --show_genes <genes>
    GFFextractor.py --show_geneids <geneids>
    GFFextractor.py (-h, --help)

Options:
  -g, --gtf <gtf>                 Input GTF/GFF file
  -f, --fasta <fasta>             FASTA sequence file
  -o, --output <output>           Output file
  -t, --feature_type <feature>    Feature type
  -n, --gene_name <gene>          Gene(s) to extract. If a list, must names must be separated by spaces.  
  -b, --bound <boundary>          Retrieve a given number of bases on either side of
                                  the feature instead of thesequence corresponding to
                                  a feature
  --show_features <features>      Scan a GFF file to identify the features present
  --show_genes <genes>            Scan a GFF file to identify the genes present
  --show_geneids <geneids>        Scan a GFF file to identify the geneIDs present
  -h, --help                      Show this message and exit.
"""

__author__='miles smith'
__email__='mileschristiansmith@gmail.com'
import sys

import click
import progressbar
import pyfaidx
from gtfparse import read_gtf


@click.command()
@click.option("--gtf", "-g",
              help="input GTF/GFF file",
              default=None,
              required=True,
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
              default=0,
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
def main(gtf, fasta, output, feature_type, gene_name, bound, show_features, show_genes, show_geneids) -> None:
    """A utility for extracting sequences from a FASTA file for a given GFF annotation"""

    if len(sys.argv) == 1:
        click.help_option()
        # print(f"GFFextractor: A utility for extracting sequences "
        #       f"from a FASTA file for a given annotation GFF file.\n"
        #       "For more information, run 'python GFFextractor --help")
        sys.exit(1)

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
            gene_name: str,
            boundary: int) -> None:

    print("Parsing GTF/GFF file.")
    records = read_gtf(gtffile)
    if gene_name:
        records = records[records['gene_name'].isin(gene_name)]
    if feature_type:
        records = records[records['feature'] == feature_type]
        records = records[['seqname', 'feature', 'start', 'end', 'strand', 'frame', 'gene_name', f'{feature_type}_id']].drop_duplicates()
    else:
        records = records[['seqname', 'feature', 'start', 'end', 'strand', 'frame', 'gene_name']].drop_duplicates()

    print(f"{len(records)} total records found.")

    print(f"Loading the sequences in {fastafile}.  "
          f"Note: if this is the first time opening this file, "
          "it may take a few moments as an index is built.")
    sequences = pyfaidx.Fasta(fastafile)
    print(f"Finished loading {fastafile}")

    seq_count = 1
    seq_widgets = ['Matching features to sequences: ', progressbar.Counter(),
                   f'/{records.index}', progressbar.Percentage(),
                   ' ', progressbar.Bar(), progressbar.Timer(), ' ', progressbar.ETA()]
    seq_progress = progressbar.ProgressBar(widgets=seq_widgets,
                                           maxval=len(records.index)).start()

    final_list = []
    for rec in records.itertuples():
        seq_progress.update(seq_count)
        seq_count += 1

        if rec.start != rec.end: # you'd be surprised
            if rec.strand == "+":
                # for a normal, say, exonic sequence
                if boundary == 0:
                    try:
                        seq = pyfaidx.Sequence(name=f"{rec.gene_name}_{getattr(rec, f'{feature_type}_id')}_{rec.strand}_{rec.start}_{rec.end}",
                                               seq=sequences[rec.seqname][rec.start:rec.end].seq)
                        final_list.append(seq)
                    except ValueError:
                        print(f"problem with {rec.gene_name} {rec.start} {rec.end} {rec.seqname} {rec.strand}")
                # for excising sequences
                else:
                    try:
                        upstream = pyfaidx.Sequence(name=f"{rec.gene_name} {getattr(rec, f'{feature_type}_id')} upstream",
                                                    seq=sequences[rec.seqname][(rec.start - boundary):rec.start].seq)
                        final_list.append(upstream)
                        downstream = pyfaidx.Sequence(name=f"{rec.gene_name} {getattr(rec,f'{feature_type}_id')} downstream",
                                                      seq=sequences[rec.seqname][(rec.end):rec.end + boundary].seq)
                        final_list.append(downstream)
                    except ValueError:
                        print(f"problem with {rec.gene_name} {rec.start} {rec.end} {rec.seqname} {rec.strand}")
            if rec.strand == "-":
                # for a normal, say, exonic sequence
                if boundary == 0:
                    try:
                        seq = pyfaidx.Sequence(name=f"{rec.gene_name}_{getattr(rec, f'{feature_type}_id')}_{rec.strand}_{rec.start}_{rec.end}",
                                               seq=sequences[rec.seqname][rec.start:rec.end].reverse.complement.seq)
                        final_list.append(seq)
                    except ValueError:
                        print(f"problem with {rec.gene_name} {rec.start} {rec.end} {rec.seqname} {rec.strand}")
                # for excising sequences
                else:
                    try:
                        downstream = pyfaidx.Sequence(name=f"{rec.gene_name} {getattr(rec, f'{feature_type}_id')} downstream",
                                                      seq=sequences[rec.seqname][(rec.start - boundary):rec.start].reverse.complement.seq)
                        final_list.append(downstream)
                        upstream = pyfaidx.Sequence(name=f"{rec.gene_name} {getattr(rec, f'{feature_type}_id')} upstream",
                                                    seq=sequences[rec.seqname][(rec.end):rec.end + boundary].reverse.complement.seq)
                        final_list.append(upstream)
                    except ValueError:
                        print(f"problem with {rec.gene_name} {rec.start} {rec.end} {rec.seqname} {rec.strand}")

    seq_progress.finish()


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


if __name__ == "__main__":
    main()
