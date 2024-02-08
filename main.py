import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



class GeneticEntities:
    pass


class Sequences(GeneticEntities):

    def __init__(self, Dataset, DNA_RNA:bool=True): #Dataset:Series  #if dna_rna true produce a sequence that is binded to nucleotides, if false binds to aa

        self.name = Dataset.name
        self.strand = Dataset     #could also store the name just as a attribute and turn it to a series, easier( in cod dont put all those recalls)
                       #return list of all the objects of organic elements present
        self.blocks_letters = [i for i in self.strand.unique().tolist()]  #transform the list of building blocks into string, so we can easily associate(? can use this later?)
        if DNA_RNA:
            self.building_blocks = [Nucleotides(i) for i in self.blocks_letters]        # in this way each element of sequence instanciated has its own organic elements instance
        else:
            self.building_blocks = [AA(i) for i in self.blocks_letters]
        self.strand = self.strand.replace(self.blocks_letters, self.building_blocks)  #and replace the elements of dataset, replaceworks here only bc is a dataframe
        self.blocks=dict(zip(self.blocks_letters, self.building_blocks)) #this is the one to be used, could eliminate self.buildingblocks and self.blocksletter
        #this binds the creation of a sequence with a specific aa

    def mantein_organics_elements(self):
        for k,v in self.blocks.items():
            v.set_base(k)


    def get_name(self)->str:
        return self.name


    def transcription(self): #-> mRNA or aachain #DNA_instance.transcription()-->mRNA_instance the whole no start or end here, the others doesn't change
        if isinstance(self, DNA):
            transc = self.strand
            self.blocks['T'].set_base('U') #WE ARE SETTING IT AS IF THE NEGATIVE STRAND IS USED AS TEMPLATE FOR THE POSITIVE ONE
            RNA =mRNA(transc.astype(str))
            self.mantein_organics_elements()
            return RNA
        else:
            return self


    def get_strand(self): #->Series
        return self.strand

    @classmethod
    def produce_graph(cls, ser): #:series ->plt
        plt.cla()
        ser.value_counts().plot(kind='bar')  #use this for graph
        plt.xlabel("organic elements")
        plt.ylabel("frequency")
        plt.title(f'frequency of the strand')
        return plt

    @classmethod

    def gen_data(cls,ser)->list: #:series
        counts = ser.value_counts()
        n = ser.count()
        frequence = counts / n *100
        counts_list= counts.to_string().split('\n')
        counts_list[0] = 'Organic element'
        counts_list[1] += ' Max'
        counts_list[-1] += ' Min'
        freq_list= frequence.to_string().split('\n')
        freq_list[0]='Frequence'

        counts_list.append(f'count of all elements = {n}')
        freq_list.append(f'tot freq = 100')
        return [counts_list, freq_list]

    @classmethod
    def turn_in_str(cls, ser):  # uses only series
        strand=''
        for i in ser :
            strand=strand + str(i)
        return strand


class DNA(Sequences):

    def produce_negative_strand(self): #->Series

        self.blocks['T'].set_base('A')
        self.blocks['A'].set_base('T')
        self.blocks['C'].set_base('G')
        self.blocks['G'].set_base('C')
        ser = self.strand.iloc[::-1].reset_index(drop=True).astype(str)
        self.mantein_organics_elements()
        return ser


class mRNA(Sequences):

    def produce_negative_strandRNA(self): #-> Series
        self.blocks['U'].set_base('A')
        self.blocks['A'].set_base('U')
        self.blocks['C'].set_base('G')
        self.blocks['G'].set_base('C')
        ser = self.strand.iloc[::-1].reset_index(drop=True).astype(str)
        self.mantein_organics_elements()
        return ser

    def translation(self)->dict:
        codons_series_str = self.strand.astype(str)
        codons = codons_series_str + codons_series_str.iloc[1::].reset_index(drop=True) + codons_series_str.iloc[2::].reset_index(drop=True)
        codons = codons.dropna()
        negstrand = self.produce_negative_strandRNA()
        neg_codons_series_str = negstrand.astype(str)
        codonsneg = neg_codons_series_str + neg_codons_series_str.iloc[1::].reset_index(drop=True) + neg_codons_series_str.iloc[2::].reset_index(drop=True)
        codonsneg = codonsneg.dropna()
        # if N present replace with ?
        GeneticCode = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'UAA': '.', 'UAG': '.', 'UGA': '.'}
        codons = codons.map(GeneticCode).fillna('?')
        codonsneg = codonsneg.map(GeneticCode).fillna('?')
        rf_1= AAChain(codonsneg.iloc[::3].reset_index(drop=True), DNA_RNA=False)
        rf_2= AAChain(codonsneg.iloc[1::3].reset_index(drop=True), DNA_RNA=False)
        rf_3= AAChain(codonsneg.iloc[2::3].reset_index(drop=True), DNA_RNA=False)
        rf1= AAChain(codons.iloc[::3].reset_index(drop=True), DNA_RNA=False)
        rf2= AAChain(codons.iloc[1::3].reset_index(drop=True), DNA_RNA=False)
        rf3= AAChain(codons.iloc[2::3].reset_index(drop=True), DNA_RNA=False)
        return {-3:rf_3,-2:rf_2,-1:rf_1,1:rf1,2:rf2,3:rf3}        # we will repeat this for neg strand, we will access the single aa by an input(indexing)


class AAChain(Sequences):  #has the whole strand of aa,so will have another attribute array that coins all the possible proteins(oligos)

    def __init__(self, Dataset, DNA_RNA=False):
        super().__init__(Dataset, DNA_RNA)
        self.starts = (np.where(self.strand == self.blocks['M']))[0]
        self.stops = (np.where(self.strand == self.blocks['.']))[0]
        stop_index = 0
        list_p = []
        list_len = []
        for i in self.starts:
            while stop_index < len(self.stops) and i > self.stops[stop_index]:
                stop_index += 1
            if stop_index < len(self.stops):
                list_p.append(self.strand.iloc[i:self.stops[stop_index]])
                list_len.append(self.stops[stop_index]-i)
        self.aachains = pd.DataFrame({'aachains': list_p, 'length': list_len}).sort_values(by='length').reset_index(drop=True)

    def get_oligos(self): #->Dataset
        ind = (np.where(self.aachains['length'] > 19))[0][0]  #FINDS THE FIRST element that has 20 aa
        return self.aachains.iloc[0:ind]

    def get_proteins(self): #->Dataset
        ind = np.where(self.aachains['length'] > 19)[0][0]
        return self.aachains.iloc[ind:]

    def get_single_aachain(self, ind:int): #->Series
        return self.aachains['aachains'][ind]

    def get_single_prot_len(self, ind:int)->int:
        return self.aachains['length'][ind]


class OrganicElements(GeneticEntities):

    def __init__(self, letter: str):
        self.letter = letter

    def get_letter(self)->str:
        return self.letter

    def set_base(self, base:str):
        self.letter = base

    def __str__(self):
        return self.letter

    def __lt__(self, other):  #< lower than used for count_values() function in dataframe
        return self.letter < self.letter

class Nucleotides(OrganicElements):

    pass

class AA(OrganicElements):  #associate to the single letter
    pass


