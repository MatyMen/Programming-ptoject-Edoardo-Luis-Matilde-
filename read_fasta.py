import pandas as pd
class Dataset:
    def __init__(self,Data=None,DNA_RNA=None):
        self.data = Data
        self.DNA_RNA = DNA_RNA

    def readfasta_set(self, fasta: dir, DNA_RNA: bool = True):
        f = open(fasta, 'r')
        name = f.readline().strip()
        bases = f.read()
        bases = bases.replace('\n', '')
        mult_genomes = bases.find('>')
        if mult_genomes != -1:
            bases = bases[:mult_genomes]
        bases = bases.upper()
        f.close()
        Data = pd.Series(list(bases), name=name)

        self.data = Data
        self.DNA_RNA = DNA_RNA

    def get_data(self):
        return self.data
    def get_DNA_RNA(self):
        return self.DNA_RNA

