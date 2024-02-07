def Complement(nucleotide):
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'

class DNASequence:
    def __init__(self, nucleotides):
        self.nucleotides = nucleotides.copy()

    def get_sequence(self):
        return self.nucleotides
    
    def get_length(self):
        return len(self.nucleotides)
    
    def get_complement(self):
        return [Complement(nucleotide) for nucleotide in self.nucleotides] #1st one-liner
    
    def get_nucleotide(self, index):
        return self.nucleotides[index]
    
    def find_alignment(self, seq):
        nucleotide_str = ""
        for nucleotide in self.nucleotides:
            nucleotide_str += nucleotide
        return nucleotide_str.find(seq)
    
    def replace_sequence(self, seq):
        self.nucleotides = seq.copy() #2nd one-liner

class Enzyme:
    def process(self, dna_sequence):
        return dna_sequence # Basic enzyme does nothing
    
class Polymerase(Enzyme):
    def process(self, dna_sequence):
        return dna_sequence.get_complement()

class Mutase(Enzyme):
    def __init__(self, freq):
        self.freq = freq
    
    def process(self, dna_sequence):
        res_sequence = dna_sequence.get_sequence()
        for i in range(self.freq-1, res_sequence.get_length(), self.freq):
            res_sequence[i] = Complement(res_sequence[i])
        dna_sequence.replace_sequence(res_sequence)

class CRISPR:
    def __init__(self, seq):
        self.seq = seq

    def process(self, dna_sequence):
        i = dna_sequence.find_alignment(self.seq)
        while i != -1:
            dna_sequence.replace_sequence(
                dna_sequence.get_sequence()[:i] + ['W'] + dna_sequence.get_sequence()[i+len(self.seq)]:
            )
            i = dna_sequence.find_alignment(self.seq)

class CRISPR_Cas9(CRISPR):
    def __init__():