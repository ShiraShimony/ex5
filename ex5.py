import json
import sys

def Complement(nucleotide):
    dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'W':'M', 'M':'W'}
    return dict[nucleotide]

def list_to_string(list):
    s = ""
    for e in list:
        s += e
    return s




class DNASequence:
    def __init__(self, nucleotides):
        self.m_nucleotides = nucleotides.copy()
    
    def get_sequence(self):
        return self.m_nucleotides.copy()
    
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

def replace_sub_lists(dna_sequence, origen, to):
        first = dna_sequence.find_alignment(origen)
        if first == -1:
            return
        dna_sequence.m_nucleotides[first:first + len(origen)] = to
        replace_sub_lists(dna_sequence, origen, to)
        return
    
class Enzyme:
    def __init__(self):#do we need it?
        return
    def process(self, dna_sequence):
        return


class Mutase(Enzyme):
    def __init__(self, freq):
        self.m_freq = freq
    
    def process(self, dna_sequence):
        dna_sequence.m_nucleotides[0::int(self.m_freq)] = dna_sequence.get_complement()[0::int(self.m_freq)]

class Polymerase(Enzyme):
    def process(self, dna_sequence):
        return dna_sequence.get_complement()
    

class CRISPR(Enzyme):
    def __init__(self, seq):
        self.m_seq = seq

    def process(self, dna_sequence):
        replace_sub_lists(dna_sequence, self.m_seq, ['W'])


class CRISPR_Cas9(CRISPR):
    def __init__(self, seq, new_seq):
        super().__init__(seq)
        self.m_new_seq = new_seq
    
    def process(self, dna_sequence):
        super().process(dna_sequence)
        to_replace = [nuc for nuc in self.m_new_seq]
        replace_sub_lists(dna_sequence, 'W', to_replace)
        replace_sub_lists(dna_sequence, 'M', DNASequence(to_replace).get_complement())

def processData(dir_path):
    name_index = 0
    enzyme_index = 1
    first_argument = 2
    second_argument = 3

    dna_name = dir_path + '/DNA.json'
    with open(dna_name, 'r') as file:
        loaded = json.load(file)
    
    protocol_name = dir_path + '/protocol.txt'
    with open(protocol_name, 'r') as protocol_file:
        protocol = protocol_file.readlines()


    for e in protocol:
        arguments = e.strip('\n').split(' ')
        loaded[arguments[name_index]] =  DNASequence([nuc for nuc in loaded[arguments[name_index]]])
        enzyme = Polymerase()
        if arguments[enzyme_index] == "Polymerase":
            enzyme = Polymerase()
        if arguments[enzyme_index] == "Mutase":
            enzyme = Mutase(arguments[first_argument])
        if arguments[enzyme_index] == "CRISPR":
            enzyme = CRISPR(arguments[first_argument])
        if arguments[enzyme_index] == "CRISPR/Cas9":
            enzyme = CRISPR_Cas9(arguments[first_argument], arguments[second_argument])
        enzyme.process(loaded[arguments[name_index]])
        loaded[arguments[name_index]] = list_to_string(loaded[arguments[name_index]].get_sequence())

    name = dir_path + '/ModifiedDNA.json'
    with open(name, 'w') as file:
        json.dump(loaded, file, indent=4)

if __name__ == "__main__":
    processData(sys.argv[1])