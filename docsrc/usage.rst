Usage
------

MerryCRISPR has three subfunctions: `new-species`_, `prep-sequences`_, and `create-library`_::

    Usage: merrycrispr [OPTIONS] COMMAND [ARGS]...

    Options:
      --help  Show this message and exit.

    Commands:
      create-library  Build a CRISPR library
      new-species     Import a new species from Ensembl.
      prep-sequences  Generate target sequences to search for spacers

new-species
~~~~~~~~~~~
::

    Usage: merrycrispr new-species [OPTIONS]

      Import a new species from Ensembl.  Optionally, build a Bowtie index for
      the files.

    Options:
      -s, --species_value TEXT      Species value to retrieve sequences for.  Must
                                    be a value for the attribute chosen, i.e.
                                    'canis_familiaris' if using the 'name'
                                    attribute
      -a, --species_attribute TEXT  Species attribute to search when retrieving
                                    sequences. Acceptable values include
                                    'taxon_id', 'accession', 'aliases',
                                    'division', 'groups', 'release', 'name',
                                    'strain', 'strain_collection', 'display_name',
                                    'assembly', and 'common_name'
      -d, --dest TEXT               Location in which to place downloaded files
                                    and, if selected, new Bowtie index.
      -b, --build_bowtie            Build a bowtie index for the retrieved species
      -w, --show_available_builds   Display a list of the genome builds available
                                    from Ensembl
      --help                        Show this message and exit.

The `new-species` subcommand is used to import sequences from Ensembl. Currently, `new-species` is only able to
retrieve the most current release ofa genome build (version 97 as this is written). To import files, you will need to
know the value for one of the attributes for that genome's release.  In plain English, it is probably best if you use
`accession` or `assembly` for the `--species_attribute` argument, look up the value on
`Ensembl <https://uswest.ensembl.org/info/about/species.html>`_, and use that for the `--species_value`
argument.  Alternatively, you can run `new-species` with just the `--show_available` argument to see possible values.
While you can use the `common_name` (i.e "dog") or species `name` (i.e. "Canis familiaris`), those are
non-unique for some values; if that is the case, `new-species` will exit and as you to make a more precise request.

The `--dest` argument is optional.  If you do not provide it, the files will be placed in a subdirectory of the current
working directory called "/mc_resources".

Because it can take a up to several hours to build the index, by default `new-species` does not build one on species
import.  Passing the `--build_bowtie` flag when you run the command will tell it to build the index.

If there are already matching files in the `dest` directory, `new-species` will check to see if the file sizes match
those on Ensembl and, if they do, will not redownload the files.  This allows you to run the program some time after
downloading to build the index.

An example of using `new-species`::

    merrycrispr new-species --species_value homo_sapiens \
        --species_attribute name \
        --dest ~/workdir/mc_human_files \
        --build_bowtie

Note that this is all provided as a convenience -- you can download the sequence and annotation files and/or build
the Bowtie index manually.


prep-sequences
~~~~~~~~~~~~~~~
::

    Usage: merrycrispr prep-sequences [OPTIONS]

      Generate target sequences to search for spacers

      Generally, a utility for extracting sequences from a FASTA file for a
      given GTF annotation.

    Options:
      -l, --library_type TEXT  target library type. Accepted values are
                               'knockout', 'repressor/activator', 'excision', or
                               'Cas13'
      -g, --gtf TEXT           input GTF/GFF file
      -f, --fasta TEXT         FASTA sequence file
      -o, --output TEXT        output file
      -t, --feature_type TEXT  feature type
      -n, --gene_name TEXT     gene(s) to extract
      --bound INTEGER          Retrieve a given number of bases on either side of
                               the feature instead of the sequence corresponding
                               to a feature
      --show_features          Scan an annotation file to identify the features present
      --show_genes             Scan an annotation file to identify the genes present
      --show_geneids           Scan an annotation file to identify the geneIDs present
      --help                   Show this message and exit.

The `prep-sequences` subcommand is used to match sequences in a FASTA file to those described in a GTF or GFF annotation
file for downstream use by the `create-library` command.  `prep-sequences` can, for example, be used to select one or
more particular genes, a certain biotype, or regions surrounding a gene.

Usage is fairly straightforward: give `prep-sequences` a FASTA and GTF, tell it the library_type, and give it the name
of an output file to write to.  The `--library_type` options are setup to produce spacers that are appropriate for the
desired function - the "knockout" option extracts exonic sequences; the "repressor/activator" option extracts sequences
-/+ 100 nucleotides surrounding the TSS (or, more precisely, the start of exon 1); the excision option excises a 100
nucleotide block upstream and 100 nucleotides downstream of a feature ("gene", by default), though that can be changed
by passing a number to the `--bound` argument; and "Cas13" returns sequences matching the processed mRNA. Note that the
`--gene_name` argument can take multiple values with just a space between the gene names.

For example, to extract all of the exons from the most current human assembly for use in a knockout library::

    merrycrispr prep-sequences --library_type knockout \
        --gtf ~/workdir/mc_human_files/Homo_sapiens.GRCh38.97.gtf.gz \
        --fasta ~/workdir/mc_human_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        --output ~/workdir/mc_human_files/human_exons.fa



create-library
~~~~~~~~~~~~~~~
::

    Usage: merrycrispr create-library [OPTIONS]

      Build a CRISPR library

    Options:
      -i, --input TEXT                Input FASTA file containing sequences to
                                      target.
      -p, --output TEXT               Name of file to write library to (in CSV
                                      format).
      -n, --nuclease TEXT             Cas nuclease to design for.  Current options
                                      include SpCas9, Cpf1, Cas13a, Csc2, SaCas9
      -r, --reference TEXT            Path to the directory containing the
                                      appropriate Bowtie index.
      --largeindex BOOLEAN
      --rule_set TEXT                 On-target score rule set to use. Current
                                      options include '1', '2', and 'Azimuth'
      -on, --ontarget_score_threshold INTEGER
                                      Spacers with an on-target score below this
                                      will be ignored.
      -off, --offtarget_score_threshold INTEGER
                                      Spacers with an off-target score below this
                                      will be ignored.
      --offtarget_count_threshold INTEGER
                                      Spacers with more than this many off-targets
                                      will be ignored.
      --spacers_per_feature INTEGER   Number of spacers to find for each feature.
      --paired BOOLEAN                Should spacers be designed to work as pairs
                                      (e.g. for excision)?
      --number_upstream_spacers INTEGER
                                      If designing paired spacers, number of
                                      spacers to design that target upstream of
                                      the feature.
      --number_downstream_spacers INTEGER
                                      If designing paired spacers, number of
                                      spacers to design that target downstream of
                                      the feature.
      --min_paired_distance INTEGER   If designing paired spacers, minimum space
                                      required between the up- and downstream
                                      spacers.
      -c, --cores INTEGER             Number of processors to use. By default,
                                      will use all available.
      --help                          Show this message and exit.

The `create-library` subcommand takes the sequences from `prep-sequences` and returns a series of protospacers according
to the given specifications. At current, `create-library` is able to find spacers for SpCas9, Cpf1, Cas13a, Csc2,
and SaCas9, though the ability to add any nuclease is forthcoming.

At current, there are only on-target and off-target scoring functions for SpCas9-based spacers, mainly because I am
unaware of algorithms for the other nucleases.  For Cas9, the on-target scoring is based on the rules from the papers
`Donch et al. Nat Biotechnol. 2014 <https://doi.org/10.1038/nbt.3026>`_ (Rule Set 1) or
`Donch et al. Nat Biotech 2006 <https://doi.org/10.1038/nbt.3437>`_ (Rule Set 2), with the latter as implemented by my
heavily updated version of `Azimuth <https://github.com/milescsmith/Azimuth/>`_. Unless you have a good reason to not do
so, Rules Set 2 should be used. Off-target scoring is a simple algorithm from
`Hsu el al. Nat Biotechnol. 2014 <https://doi.org/10.1038/nbt.2647>`_.  Scores range from 0 to 100, with higher numbers
being more desirable.

Following from the example in `prep-sequences`_ above, the following would find 6 spacers per gene for Cas9 with on- and
off-target scores above 50 using Rule Set 2::

    merrycrispr create-library --input ~/workdir/mc_human_files/human_exons.fa \
        --output ~/workdir/mc_human_files/human_knockout.csv \
        --nuclease SpCas9 \
        --reference ~/workdir/mc_human_files/Homo_sapiens \
        --rule_set 2 \
        --ontarget_score_threshold 50
        --offtarget_score_threshold 50
        --spacers_per_feature 6

Of the options, probably the `--reference` is the easiest to mistake how to use.  Bowtie creates a series of 4-7 files
with the suffix `.ewbt` (or, it the genome is exceptionally large, `.ewbtl`, thus the `--largeindex` flag). So, for
instance, when `new-species` above created the index, the files produced were::

    .
    ├── Homo_sapiens.1.ebwt
    ├── Homo_sapiens.2.ebwt
    ├── Homo_sapiens.3.ebwt
    ├── Homo_sapiens.4.ebwt
    ├── Homo_sapiens.rev.1.ebwt
    └── Homo_sapiens.rev.2.ebwt

When specifying the index to use, Bowtie needs to know the common prefix - that is, everything up to the `.1.ebwt`.


Docker
-------
While the method of installation via :ref:`installation:Docker` is simpler, its self-contained nature means that
by default a Docker container cannot see files outside the container.  To make it work, we need to run the container
with a "bind mount" to an existing directory outside the container - this serves as a link between a real folder and
one inside the container.

To run a container with such a mount::

    docker run -v OUTSIDE_PATH:INSIDE_PATH IMAGE COMMAND [COMMAND OPTIONS]

So to run the `new-species`_ command from above::

    docker run -v "$(pwd)"/:$HOME/workspace merrycrispr merrycrispr new-species \
        --species_value homo_sapiens \
        --species_attribute name \
        --dest $HOME/workspace/

The bind mount flag `-v "$(pwd)"/:$HOME/workspace` flag above instructs Docker to set up a link from the *p*resent
*w*orking *d*irectory to a directory named "workspace"in the user's home folder inside the docker. Note that later in
the command, the location passed to the `--dest` argument is *inside* the container - this is fine as the link is
two-way, meaning that any files written to "$HOME/workspace" are also written to your local $PWD/workspace directory.

Also note that the double `merrycrispr merrycrispr` is not a typo - the first is to create a container from the
merrycrispr image, the second is to run the merrycrispr command inside the new container.