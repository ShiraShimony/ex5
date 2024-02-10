import json
import sys

def Complement(nucleotide):
    """
    Returns the complement of a given DNA nucleotide.

    Parameters:
    - nucleotide (str): The input DNA nucleotide.

    Returns:
    - str: The complement of the input nucleotide.
    """
    dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'W':'M', 'M':'W'}
    return dict[nucleotide]

def list_to_string(lst):
    """
    Converts a list of characters into a single string.

    Parameters:
    - lst (list): List of characters.

    Returns:
    - str: A string containing all characters from the input list concatenated.
    """
    return ''.join(lst)

class DNASequence:
    """
    Represents a DNA sequence.

    """
    def __init__(self, nucleotides):
        """
        Initializes a DNASequence object with a list of nucleotides.

        Parameters:
        - nucleotides (list/str/DNASequence): Initial nucleotide sequence.
        """
        if type(nucleotides) == DNASequence:
            self.nucleotides = nucleotides.nucleotides.copy()
        else:
            self.nucleotides = [nuc for nuc in nucleotides]
    
    def get_sequence(self):
        """
        Returns a copy of the DNA sequence.

        Returns:
        - list: A copy of the list of DNA nucleotides.
        """
        return self.nucleotides.copy()
    
    def get_length(self):
        """
        Returns the length of the DNA sequence.

        Returns:
        - int: Length of the DNA sequence.
        """
        return len(self.nucleotides)
    
    def get_complement(self):
        """
        Returns the complement of the DNA sequence.

        Returns:
        - list: List containing the complement the sequence.
        """
        return [Complement(nuc) for nuc in self.nucleotides]
    
    def get_nucleotide(self, index):
        """
        Returns the nucleotide at the specified index.

        Parameters:
        - index (int): Index of the nucleotide in the sequence.

        Returns:
        - str: The nucleotide at the specified index.
        """
        return self.nucleotides[index]
    
    def find_alignment(self, seq):
        """
        Finds the index of the first occurrence of a subsequence in the DNA sequence.

        Parameters:
        - seq (str): Subsequence to search for.

        Returns:
        - int: Index of the first occurrence or -1 if not found.
        """
        return list_to_string(self.nucleotides).find(seq)
    
    def replace_sequence(self, seq):
        """
        Replaces the DNA sequence with a new sequence.

        Parameters:
        - seq (list/str/DNASequence): New DNA sequence.
        """
        self.nucleotides = seq.copy()

def replace_sub_lists(dna_sequence, origin, to):
    """
    Recursively replaces sublists in a DNA sequence.

    Parameters:
    - dna_sequence (DNASequence): The DNA sequence to be modified.
    - origin (str/list): Subsequence to be replaced.
    - to (str/list): New sequence for replacement.
    """
    first = dna_sequence.find_alignment(origin)
    tmp = dna_sequence.nucleotides.copy()
    if first == -1:
        return
    dna_sequence.nucleotides[first:first + len(origin)] = to
    assert dna_sequence.nucleotides == tmp[:first] + to + tmp[first + len(origin):]
    replace_sub_lists(dna_sequence, origin, to)
    return

class Enzyme:
    """
    Represents a generic enzyme.

    Methods:
    - process(dna_sequence): Abstract method for processing a DNA sequence using the enzyme.
    """
    def __init__(self):
        """
        Initializes an Enzyme object. (Not used in the current implementation)
        """
        return
    
    def process(self, dna_sequence):
        """
        Abstract method for processing a DNA sequence using the enzyme.

        Parameters:
        - dna_sequence (DNASequence): The DNA sequence to be processed.
        """
        raise NotImplementedError("Subclasses must implement this method")

class Mutase(Enzyme):
    """
    Represents a Mutase enzyme, a subclass of Enzyme.

    Attributes:
    - freq (int): Frequency of mutation.

    Methods:
    - process(dna_sequence): Mutates the DNA sequence based on the specified frequency.
    """
    def __init__(self, freq):
        """
        Initializes a Mutase enzyme with a given frequency.

        Parameters:
        - freq (int): Frequency of mutation.
        """
        self.freq = freq
    
    def process(self, dna_sequence):
        """
        Mutates the DNA sequence based on the specified frequency.

        Parameters:
        - dna_sequence (DNASequence): The DNA sequence to be mutated.

        Returns:
        - DNASequence: The mutated DNA sequence.
        """
        dna_sequence.nucleotides[(self.freq - 1)::int(self.freq)] = dna_sequence.get_complement()[(self.freq - 1)::int(self.freq)]
        return dna_sequence

class Polymerase(Enzyme):
    """
    Represents a Polymerase enzyme, a subclass of Enzyme.

    Methods:
    - process(dna_sequence): Generates the complement of the DNA sequence using the Polymerase enzyme.
    """
    def process(self, dna_sequence):
        """
        Generates the complement of the DNA sequence using the Polymerase enzyme.

        Parameters:
        - dna_sequence (DNASequence): The DNA sequence to be complemented.

        Returns:
        - DNASequence: A new DNASequence object representing the complement.
        """
        return DNASequence(dna_sequence.get_complement())

class CRISPR(Enzyme):
    """
    Represents a CRISPR enzyme, a subclass of Enzyme.

    Attributes:
    - seq (str): The sequence to be targeted for replacement.

    Methods:
    - process(dna_sequence): Replaces occurrences of the specified sequence in the DNA sequence.
    """
    def __init__(self, seq):
        """
        Initializes a CRISPR enzyme with a specified sequence.

        Parameters:
        - seq (str): The sequence to be targeted for replacement.
        """
        self.seq = seq

    def process(self, dna_sequence):
        """
        Replaces occurrences of the specified sequence in the DNA sequence.

        Parameters:
        - dna_sequence (DNASequence): The DNA sequence to be processed.

        Returns:
        - DNASequence: The modified DNA sequence.
        """
        replace_sub_lists(dna_sequence, list_to_string(self.seq), ['W'])
        return dna_sequence

class CRISPR_Cas9(CRISPR):
    """
    Represents a CRISPR/Cas9 enzyme, a subclass of CRISPR.

    Attributes:
    - seq (str): The target sequence to be replaced.
    - new_seq (str/list/DNASequence): The new sequence for replacement.

    Methods:
    - process(dna_sequence): Executes the CRISPR_Cas9 modification on the DNA sequence.
    """
    def __init__(self, seq, new_seq):
        """
        Initializes a CRISPR/Cas9 enzyme with a target sequence and a new sequence for replacement.

        Parameters:
        - seq (str): The target sequence to be replaced.
        - new_seq (str/list/DNASequence): The new sequence for replacement.
        """
        super().__init__(seq)
        self.new_seq = new_seq
    
    def process(self, dna_sequence):
        """
        Executes the CRISPR_Cas9 modification on the DNA sequence.

        Parameters:
        - dna_sequence (DNASequence): The DNA sequence to be modified.

        Returns:
        - DNASequence: The modified DNA sequence.
        """
        super().process(dna_sequence)
        to_replace = [nuc for nuc in self.new_seq]
        replace_sub_lists(dna_sequence, 'W', to_replace)
        replace_sub_lists(dna_sequence, 'M', DNASequence(to_replace).get_complement())
        return dna_sequence

def processData(dir_path):
    """
    Reads DNA data and enzyme protocol from files, processes the DNA sequences, and saves the modified sequences to a new file.

    Parameters:
    - dir_path (str): Directory path containing DNA data and protocol files.
    """
    name_index = 0
    enzyme_index = 1
    first_argument = 2
    second_argument = 3

    # Load DNA data from JSON file
    dna_name = dir_path + '/DNA.json'
    with open(dna_name, 'r') as file:
        loaded = json.load(file)
    
    # Load enzyme protocol from text file
    protocol_name = dir_path + '/protocol.txt'
    with open(protocol_name, 'r') as protocol_file:
        protocol = protocol_file.readlines()

    # Process each line in the protocol
    for e in protocol:
        arguments = e.strip('\n').split(' ')
        loaded[arguments[name_index]] =  DNASequence([nuc for nuc in loaded[arguments[name_index]]])
        enzyme = Polymerase()  # Default enzyme is Polymerase
        if arguments[enzyme_index] == "Polymerase":
            enzyme = Polymerase()
        if arguments[enzyme_index] == "Mutase":
            enzyme = Mutase(arguments[first_argument])
        if arguments[enzyme_index] == "CRISPR":
            enzyme = CRISPR(arguments[first_argument])
        if arguments[enzyme_index] == "CRISPR/Cas9":
            enzyme = CRISPR_Cas9(arguments[first_argument], arguments[second_argument])
        loaded[arguments[name_index]] = list_to_string(enzyme.process(loaded[arguments[name_index]]).get_sequence())

    # Save modified DNA sequences to a new JSON file
    name = dir_path + '/ModifiedDNA.json'
    with open(name, 'w') as file:
        json.dump(loaded, file, indent=4)

# Example usage:
# processData('/path/to/data_directory')
