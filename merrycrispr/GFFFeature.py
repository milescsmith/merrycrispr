class GFFFeature(object):
    """ Class used to extract and format information described in a GFF3 file and found in and FASTA file
    Requires a string read in from the GFF3 file for initialization, adding the corresponding sequence later.
    The overridden __str__ method is specific for creating output for the MerryCRISPR script.
    """

    __slots__ = ['record','attributes','seqname','source','feature','start','end','strand','name','parent','geneID','transcriptID','seq']

    def __init__(self, unprocessed, from_ensembl=True):

        # Convert the string into a dictionary corresponding to the standard GFF3 features
        keys = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        values = unprocessed.strip('\n').split('\t')
        self.record = dict(zip(keys, values))
        self.seqname = self.record['seqname']
        self.source = self.record['source']
        self.feature = self.record['feature']
        self.start = int(self.record['start'])
        self.end = int(self.record['end'])
        # self.score = record['score']
        self.strand = self.record['strand']
        # self.frame = record['frame']
        self.attributes = {item.split('=')[0]: item.split('=')[1] for item in self.record['attribute'].strip(';').split(';')}

        # The below are hacks to make this class work with the GCF_000223135.1_CriGri_1.0.gff
        # and gencode.v20.annotation.gff3 files and anything fron Ensembl
        # (Ensembl's crap being why I have to do the whole name from parent lookup stuff)
        # If you find the below to be a real shitshow, please feel free to either clean it up or yell at Ensembl!

        # class required members:
        # name
        # parent
        # geneID
        # transcriptID

        if from_ensembl:
            gene_types = ['gene', 'miRNA_gene', 'rRNA_gene', 'snRNA_gene', 'snoRNA_gene', 'RNA']
            transcript_types = ['transcript', 'miRNA', 'rRNA', 'snRNA', 'snoRNA', 'lincRNA', 'processed_pseudogene']
            if self.feature == 'chromosome' or self.feature == 'supercontig':
                self.name = self.attributes['ID']
                self.parent = None
                self.geneID = None
                self.transcriptID = None
            elif self.feature in gene_types:
                if 'Name' in self.attributes.keys():
                    self.name = self.attributes['Name']
                else:
                    self.name = 'uncharacterized'
                self.geneID = self.attributes['gene_id']
                self.transcriptID = None
                self.parent = None
            elif self.feature in transcript_types:
                self.name = None
                self.transcriptID = self.attributes['transcript_id']
                self.parent = self.attributes['Parent'].split(':')[1]
                self.geneID = self.parent
            elif self.feature == 'lincRNA_gene':
                self.name = 'unnamed lincRNA'
                self.parent = None
                self.geneID = self.attributes['gene_id']
            elif self.feature == 'exon':
                self.name = None
                if self.source == 'lncipedia.org':
                # because, of course, the lncipedia has be to different...
                    self.parent = self.attributes['Parent'].split(':')[0][:-2]
                    self.name = self.parent
                else:
                    self.parent = self.attributes['Parent'].split(':')[1]
                self.transcriptID = self.parent
                if 'rank' in self.attributes.keys():
                    self.exonrank = self.attributes['rank']
                self.geneID = None
            elif self.feature == 'CDS' or self.feature == 'five_prime_UTR' or self.feature == 'three_prime_UTR':
                self.parent = self.attributes['Parent'].split(':')[1]
                self.name = None
                self.transcriptID = self.parent
                self.geneID = None
            elif self.feature == 'mt_gene' or self.feature == 'pseudogene':
                if 'gene_id' in self.attributes.keys():
                    self.name = 'none'
                    self.geneID = self.attributes['gene_id']
                    self.transcriptID = None
                    self.parent = None
                elif 'transcript_id' in self.attributes.keys():
                    self.name = None
                    self.transcriptID = self.attributes['transcript_id']
                    self.parent = self.attributes['Parent'].split(':')[1]
                    self.feature = 'pseudogene_xcript'
                    self.geneID = self.parent
            else:  # look for Parent, gene_id, transcript_id, name.  This should (read: needs) to handle any of the cases I haven't hardwired in above
                if 'Name' in self.attributes.keys() and 'gene_id' in self.attributes.keys():
                    self.name = self.attributes['Name']
                    self.geneID = self.attributes['gene_id']
                    self.transcriptID = None
                    self.parent = None
                elif 'gene_id' in self.attributes.keys() and not 'Name' in self.attributes.keys():
                    self.geneID = self.attributes['gene_id']
                    self.name = 'none'
                elif 'transcript_id' in self.attributes.keys() and not 'gene_id' in self.attributes.keys():
                    self.transcriptID = self.attributes['transcript_id']
                    self.name = None
                else:
                    self.transcriptID = None
                    self.geneID = None
                    self.name = None

                if 'Parent' in self.attributes.keys():
                    self.parent = self.attributes['Parent'].split(':')[1]
                    if self.attributes['Parent'].split(':')[0] == 'gene':
                        self.geneID = self.parent
                    elif self.attributes['Parent'].split(':')[0] == 'transcript':
                        self.transcriptID = self.parent
                else:
                    self.parent = None
                    # if there is a name, set name.  we probably also have a gene_id, so set that
                    # if there isn't a name, see if there is a gene_id.  set that and then set the name to 'none'
                    # else, is there a transcript_id? Set that and see if there is a parent
        else:
            try:
                if self.attributes['gene_id']:
                    self.geneID = self.attributes['gene_id']
                elif self.attributes['Dbxref']:
                    temp = {x.split(':')[0]: x.split(':')[1] for x in self.attributes['Dbxref'].split(',')}
                    self.geneID = temp['GeneID']
            except:
                self.geneID = 'Unknown'

            try:
                if self.attributes['transcript_id']:
                    self.transcriptID = self.attributes['transcript_id']
                elif self.attributes['Dbxref']:
                    temp = {x.split(':')[0]: x.split(':')[1] for x in self.attributes['Dbxref'].split(',')}
                    self.transcriptID = temp['Genbank']
            except:
                self.transcriptID = 'Unknown'

            if 'gene_name' in self.attributes.keys():
                self.name = self.attributes['gene_name']
            elif 'gene' in self.attributes.keys():
                self.name = self.attributes['gene']
            elif 'name' in self.attributes.keys():
                self.name = self.attributes['name']
            elif 'Name' in self.attributes.keys():
                self.name = self.attributes['Name']
            else:
                self.name = None

        self.seq = ''

        # We need to save on memory at every point we can.
        # This actually frees up a good chunk of change right here
        del(self.record)
        del(self.attributes)

    def __repr__(self):
        return 'GeneID: {self.geneID}, GeneName: {self.name}, Transcript: {self.transcriptID}, Chromosome: {self.seqname}, ' \
               'Sequence: {self.seq}, Start: {self.start}, End: {self.end}, Strand: {self.strand}, ' \
               'Type: {self.feature}, Parent: {self.parent}'.format(self=self)

    def __str__(self):
        ''' MerryCRISPR requires a FASTA file with each sequence containing a header in the format:
        GENEID | TRANSCRIPTID | GENENAME | EXON RANK | CONSTITUTIVE EXON | 5' UTR END | 3' UTR STOP | EXON START | EXON END
        '''
        return '>{self.geneID}|{self.transcriptID}|{self.name}|||||{self.start}|{self.end}'.format(self=self)

    def add_seq(self, sequence):
        self.seq = sequence.upper()

    def find_name(self, name_list):
        ''' an obj should have either a gene_id and a name, a transcript_ID and a Parent, or just a Parent
         if obj has a name, skip it
        '''
        if self.name:
            pass
        # if an obj has a gene_id but no name, use Parent to get the name
        elif self.geneID:
            self.name = name_list[self.parent]
        # if an obj has a transcript_id but no gene_id, use Parent to get gene_id and then name
        elif self.transcriptID:
            self.geneID = name_list[self.transcriptID]
            self.name = name_list[self.geneID]