import json


def Complement(nucleotide):
    dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'W':'M', 'M':'W'}
    return dict[nucleotide]

def list_to_string(list):
    s = ""
    for e in list:
        s += e
    return s

def replace_sub_lists(dna_sequence, origen, to):
        first = dna_sequence.find_alignment(origen)
        if first == -1:
            return
        dna_sequence.get_nucleotide()[first:first + len(origen)] = to
        replace_sub_lists(dna_sequence[first:])
        return


class DNASequence:
    def __init__(self, nucleotides):
        self.m_nucleotides = nucleotides.copy()
    
    def get_sequence(self):
        return self.m_nucleotides
    
    def get_length(self):
        return len(self.m_nucleotides)
    
    def get_complement(self):
        return [Complement(nuc) for nuc in self.m_nucleotides]
    
    def get_nucleotide(self, index):
        return self.m_nucleotides[index]
    
    def find_alignment(self, seq):
        return list_to_string(self.m_nucleotides).find(seq)
    
    def replace_sequence(self, seq):
        self.m_nucleotides = seq.copy()
    
    
class Enzyme:
    def __init__(self):#do we need it?
        return
    def process(self, dna_sequence):
        return


class Mutase(Enzyme):
    def __init__(self, freq):
        self.m_freq = freq
    
    def process(self, dna_sequence):
        dna_sequence.get_sequence[0::self.m_freq] = dna_sequence.get_complement[0::self.m_freq].copy()
        assert dna_sequence.m_nucleotides[0::self.m_freq] == dna_sequence.get_complement[0::self.m_freq]

class Polymerase(Enzyme):
    def process(self, dna_sequence):
        return dna_sequence.get_complement()
    

class CRISPR(Enzyme):
    def __init__(self, seq):
        self.m_seq = seq.copy()

    def process(self, dna_sequence):
        replace_sub_lists(dna_sequence, self.m_seq, ['W'])


class CRISPR_Cas9(CRISPR):
    def __init__(self, seq, new_seq):
        super().__init__(seq)
    
    def process(self, dna_sequence):
        super().process(dna_sequence)
        to_replace = [nuc for nuc in self.m_seq]
        replace_sub_lists(dna_sequence, ['W'], to_replace)
        replace_sub_lists(dna_sequence, ['M'], DNASequence(to_replace).get_complement())

def processData(dir_path):
    dna_name = dir_path + 'DNA.json'
    with open(dna_name, 'r') as file:
        loaded = json.load(file)

    for e in loaded.values():
        e = DNASequence([nuc for nuc in e])#filles wrong
    
    protocol_name = dir_path + 'protocol.txt'
    with open(protocol_name, 'r') as protocol_file:
        protocol = protocol_file.readlines()

    enzymes = []

    for e in protocol:
        arguments = e.split()
        if(arguments[])