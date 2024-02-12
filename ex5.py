import json
import os

def Complement(nucleotide):
    # Arguments:
    #   * nucleotide - the nucleotide
    # returns: the complement nucleotide
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'W':
        return 'M'
    if nucleotide == 'M':
        return 'W'

class DNASequence:
    def __init__(self, nucleotides):
        # c'tor to DNASequence
        # Arguments:
        #   * nucleotides - a string or a list of nucleotides
        self.nucleotides = [nucleotide for nucleotide in nucleotides]

    def get_sequence(self):
        # return the nucleotides of the DNA
        return self.nucleotides
    
    def get_length(self):
        # return the length of the sequence
        return len(self.nucleotides)
    
    def get_complement(self):
        # return the complement sequence
        return [Complement(nucleotide) for nucleotide in self.nucleotides] #1st one-liner
    
    def get_nucleotide(self, index):
        # return the nucleotide in the given index
        return self.nucleotides[index]
    
    def find_alignment(self, seq):
        # return the index of the beginning of the given sequence
        seq_len = len(seq)
        i = 0
        while i != self.get_length():
            print(i, self.get_length())
            if self.nucleotides[i:i+seq_len] == list(seq):
                return i
            i += 1
        return -1 # if seq doesn't appear in the DNA
    
    def replace_sequence(self, seq):
        # replaces the sequence of the DNA with a new sequence "seq" 
        self.nucleotides = seq.copy() #2nd one-liner

class Enzyme:
    def __init__(self):
        # empty c'tor
        pass

    def process(self, dna_sequence):
        return # Basic enzyme does nothing
    
class Polymerase(Enzyme):
    def __init__(self):
        super().__init__()
    
    def process(self, dna_sequence):
        # produces a new DNASequence with the complement sequence
        return DNASequence(dna_sequence.get_complement())

class Mutase(Enzyme):
    def __init__(self, freq):
        super().__init__()
        self.freq = freq
    
    def process(self, dna_sequence):
        # replaces every "freq"-th nucleotide in the sequence with its complement
        res_sequence = dna_sequence.get_sequence()
        for i in range(self.freq - 1, len(res_sequence), self.freq):
            res_sequence[i] = Complement(res_sequence[i])
        dna_sequence.replace_sequence(res_sequence)
        return dna_sequence

class CRISPR(Enzyme):
    def __init__(self, seq):
        super().__init__()
        self.seq = seq

    def process(self, dna_sequence):
        # replaces every apearance of "seq" with 'W' base
        i = dna_sequence.find_alignment(self.seq)
        while i != -1:
            dna_sequence.replace_sequence(
                dna_sequence.get_sequence()[:i] +
                ['W'] +
                dna_sequence.get_sequence()[i+len(self.seq):]
            )
            i = dna_sequence.find_alignment(self.seq)
        return dna_sequence

class CRISPR_Cas9(CRISPR):
    def __init__(self, seq, new_seq):
        super().__init__(seq)
        self.new_seq = new_seq

    def process(self, dna_sequence):
        # replaces every apearance of "seq" with "new_seq"
        dna_sequence = super().process(dna_sequence)
        # Replace every W with new_seq
        i = dna_sequence.find_alignment('W')
        while i != -1:
            dna_sequence.replace_sequence(
                dna_sequence.get_sequence()[:i] +
                DNASequence(self.new_seq).get_sequence() +
                dna_sequence.get_sequence()[i+len(self.seq):]
            )
            i = dna_sequence.find_alignment('W')
        # Replace every M with the complement of new_seq
        j = dna_sequence.find_alignment('M')
        while j != -1:
            dna_sequence.replace_sequence(
                dna_sequence.get_sequence()[:j] +
                DNASequence(self.new_seq).get_complement() +
                dna_sequence.get_sequence()[i+len(self.seq):]
            )
            j = dna_sequence.find_alignment('M')
        return dna_sequence

def processData(dir_path):
    # store the sequences from DNA.json
    with open(os.path.join(dir_path,'DNA.json'), 'r') as f:
        dna_dict = json.load(f)
    for dna_name, dna_seq in dna_dict.items():
        dna_dict[dna_name] = DNASequence(dna_seq)
    # performing a certein process on a certein DNA with given parameters read by file, then writing the new DNA dictionary to a new file
    protocol_list = open(os.path.join(dir_path,'protocol.txt'), 'r')
    for line in protocol_list:
        enzyme_command = line.split()
        if len(enzyme_command) == 0:
            continue
        if enzyme_command[1] == 'Polymerase':
            selected_enzyme = Polymerase()
        if enzyme_command[1] == 'Mutase':
            selected_enzyme = Mutase(int(enzyme_command[2]))
        if enzyme_command[1] == 'CRISPR':
            selected_enzyme = CRISPR(enzyme_command[2])
        if enzyme_command[1] == 'CRISPR/Cas9':
            selected_enzyme = CRISPR_Cas9(enzyme_command[2], enzyme_command[3])

        if enzyme_command[0] in dna_dict:
            dna_dict[enzyme_command[0]] = selected_enzyme.process(dna_dict[enzyme_command[0]])
    with open(os.path.join(dir_path,'ModifiedDNA.json'), 'w') as f:
        for dna_name, dna_seq in dna_dict.items():
            dna_dict[dna_name] = (dna_dict[dna_name]).get_sequence()
        json.dump(dna_dict, f)
    protocol_list.close()

def main(dir_path):
    processData(dir_path)

# ensuring that main() is only called when the script is run directly
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        sys.exit(1)
    main(sys.argv[1])
