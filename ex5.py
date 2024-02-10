def Complement(nucleotide):
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
        self.nucleotides = [nucleotide for nucleotide in nucleotides]

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
        return # Basic enzyme does nothing
    
class Polymerase(Enzyme):
    def process(self, dna_sequence):
        return DNASequence(dna_sequence.get_complement())

class Mutase(Enzyme):
    def __init__(self, freq):
        self.freq = freq
    
    def process(self, dna_sequence):
        res_sequence = dna_sequence.get_sequence()
        for i in range(self.freq - 1, len(res_sequence), self.freq):
            res_sequence[i] = Complement(res_sequence[i])
        dna_sequence.replace_sequence(res_sequence)
        return dna_sequence

class CRISPR:
    def __init__(self, seq):
        self.seq = seq

    def process(self, dna_sequence):
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
        dna_sequence = super().process(dna_sequence)
        # Replace every W with new_seq
        i = dna_sequence.find_alignment("W")
        while i != -1:
            dna_sequence.replace_sequence(
                dna_sequence.get_sequence()[:i] +
                DNASequence(self.new_seq).get_sequence() +
                dna_sequence.get_sequence()[i+len(self.seq):]
            )
            i = dna_sequence.find_alignment("W")
        # Replace every M with the complement of new_seq
        j = dna_sequence.find_alignment("M")
        while j != -1:
            dna_sequence.replace_sequence(
                dna_sequence.get_sequence()[:j] +
                DNASequence(self.new_seq).get_complement() +
                dna_sequence.get_sequence()[i+len(self.seq):]
            )
            j = dna_sequence.find_alignment("W")
        return dna_sequence