from functools import reduce

# from collections import namedtuple
# EpitopeEntry = namedtuple('EpitopeEntry', ['aa', 'idx'])

_aa_polarity = {
    # 0 -- hydrophobic
    # 1 -- polar
    # 2 -- positive
    # 3 -- negative
    'A' : 0,
    'C' : 0,
    'D' : 3,
    'E' : 3,
    'F' : 0,
    'G' : 0,
    'H' : 2,
    'I' : 0,
    'K' : 2,
    'L' : 0,
    'M' : 0,
    'N' : 1,
    'P' : 0,
    'Q' : 1,
    'R' : 2,
    'S' : 1,
    'T' : 1,
    'V' : 0,
    'W' : 0,
    'Y' : 0
}

_aa_weights={
        'A':89.1, 
        'R':174.2,
        'N':132.1,
        'D':133.1,
        'C':121.2,
        'E':147.1,
        'Q':146.2,
        'G':75.1,
        'H':155.2,
        'I':131.2,
        'L':131.2,
        'K':146.2,
        'M':149.2,
        'F':165.2,
        'P':115.1,
        'S':105.1,
        'T':119.1,
        'W':204.0,
        'Y':181.2,
        'V':117.1
}

_aa_aromatic = {
    'A' : 0,
    'C' : 0,
    'D' : 0,
    'E' : 0,
    'F' : 1,
    'G' : 0,
    'H' : 0,
    'I' : 0,
    'K' : 0,
    'L' : 0,
    'M' : 0,
    'N' : 0,
    'P' : 0,
    'Q' : 0,
    'R' : 0,
    'S' : 0,
    'T' : 0,
    'V' : 0,
    'W' : 1,
    'Y' : 1
}

_aa_hydroxilic = {
    'A' : 0,
    'C' : 0,
    'D' : 0,
    'E' : 0,
    'F' : 0,
    'G' : 0,
    'H' : 0,
    'I' : 0,
    'K' : 0,
    'L' : 0,
    'M' : 0,
    'N' : 0,
    'P' : 0,
    'Q' : 0,
    'R' : 0,
    'S' : 1,
    'T' : 1,
    'V' : 0,
    'W' : 0,
    'Y' : 0
}

_aa_amidic = {
    'A' : 0,
    'C' : 0,
    'D' : 0,
    'E' : 0,
    'F' : 0,
    'G' : 0,
    'H' : 0,
    'I' : 0,
    'K' : 0,
    'L' : 0,
    'M' : 0,
    'N' : 1,
    'P' : 0,
    'Q' : 1,
    'R' : 0,
    'S' : 0,
    'T' : 0,
    'V' : 0,
    'W' : 0,
    'Y' : 0
}

_aa_aliphatic = {
    'A' : 1,
    'C' : 0,
    'D' : 0,
    'E' : 0,
    'F' : 0,
    'G' : 1,
    'H' : 0,
    'I' : 1,
    'K' : 0,
    'L' : 1,
    'M' : 0,
    'N' : 0,
    'P' : 1,
    'Q' : 0,
    'R' : 0,
    'S' : 0,
    'T' : 0,
    'V' : 1,
    'W' : 0,
    'Y' : 0
}

_aa_essential = {
    'A' : 0,
    'C' : 0,
    'D' : 0,
    'E' : 0,
    'F' : 1,
    'G' : 0,
    'H' : 1,
    'I' : 1,
    'K' : 1,
    'L' : 1,
    'M' : 1,
    'N' : 0,
    'P' : 0,
    'Q' : 0,
    'R' : 0,
    'S' : 0,
    'T' : 1,
    'V' : 1,
    'W' : 1,
    'Y' : 0
}

_aa_hydrophobicity = {
    'A' : 1.8,
    'C' : 2.5,
    'D' : -3.5,
    'E' : -3.5,
    'F' : 2.8,
    'G' : -0.4,
    'H' : -3.2,
    'I' : 4.5,
    'K' : -3.9,
    'L' : 3.8,
    'M' : 1.9,
    'N' : -3.5,
    'P' : -1.6,
    'Q' : -3.5,
    'R' : -4.5,
    'S' : -0.8,
    'T' : -0.7,
    'V' : 4.2,
    'W' : -0.9,
    'Y' : -1.3
}

_aa_hydration_energy = {
    'A' : 2.06,
    'C' : -0.61,
    'D' : 0.82,
    'E' : 1.04,
    'F' : -0.33,
    'G' : 1.48,
    'H' : -7.6,
    'I' : 3.07,
    'K' : -2.8,
    'L' : 3.01,
    'M' : 1.43,
    'N' : -4.43,
    'P' : -1.33,
    'Q' : -3.94,
    'R' : -9.91,
    'S' : -5.79,
    'T' : -4.18,
    'V' : 2.79,
    'W' : -3.49,
    'Y' : -7.13
}

_aa_volume_sidechains = {
    'A' : 27.5,
    'C' : 44.6,
    'D' : 40,
    'E' : 62,
    'F' : 115.5,
    'G' : 0,
    'H' : 79,
    'I' : 93.5,
    'K' : 100,
    'L' : 93.5,
    'M' : 94.1,
    'N' : 58.7,
    'P' : 41.9,
    'Q' : 80.7,
    'R' : 105,
    'S' : 29.3,
    'T' : 51.3,
    'V' : 71.5,
    'W' : 145.5,
    'Y' : 117.3    
}

_aa_polarity_quant = {
    'A' : 8.1,
    'C' : 5.5,
    'D' : 13,
    'E' : 12.3,
    'F' : 5.2,
    'G' : 9,
    'H' : 10.4,
    'I' : 5.2,
    'K' : 11.3,
    'L' : 4.9,
    'M' : 5.7,
    'N' : 11.6,
    'P' : 8,
    'Q' : 10.5,
    'R' : 10.5,
    'S' : 9.2,
    'T' : 8.6,
    'V' : 5.9,
    'W' : 5.4,
    'Y' : 6.2 
}

_aa_sasa = {
    'A' : 1.181,
    'C' : 1.461,
    'D' : 1.587,
    'E' : 1.862,
    'F' : 2.228,
    'G' : 0.881,
    'H' : 2.025,
    'I' : 1.81,
    'K' : 2.258,
    'L' : 1.931,
    'M' : 2.034,
    'N' : 1.655,
    'P' : 1.468,
    'Q' : 1.932,
    'R' : 1.932,
    'S' : 1.298,
    'T' : 1.525,
    'V' : 1.645,
    'W' : 2.663,
    'Y' : 2.368 
}

class Analyze:

    epitope_sequence = None
    epitope_split = None
    total_weight = None
    num_sulphur = None
    _hydrophobicity_thresh = 2
    _volume_thresh = 70

    def __init__(self, epitope_sequence):
        self.epitope_sequence = epitope_sequence

        # calculating the epitope split for convenience
        self.epitope_split = [{"aa":entry.strip()[0], "idx":int(entry.strip()[1:])}
                                            for entry in epitope_sequence.split(',')]


    def get_num_sulphur(self):

        """
            getting the number of sulphurs which can create
            a disulphurous bond.
        """
        
        if self.num_sulphur is not None:
            return self.num_sulphur

        self.num_sulpurs = 0
        for entry in self.epitope_split:
            self.num_sulpurs += 1 if entry['aa'] == 'C' else 0
        
        return self.num_sulpurs

    def get_num_aa(self, aa):

        """
            getting the number of a particular amino acid from
            the code given.
        """

        num_aa = 0
        for entry in self.epitope_split:
            num_aa += 1 if entry['aa'] == aa else 0
        
        return num_aa
    
    def get_total_weight(self):
        
        """
            calculate the total moklecular weight of the epitope
        """
        if self.total_weight is not None:
            # do nothing
            
            return self.total_weight
        else:
            self.total_weight = 0
            for entry in self.epitope_split:
                self.total_weight += _aa_weights[ entry['aa'] ]
        
        return self.total_weight

    def _filter_sum(self, condition):

        sum = 0
        for entry in self.epitope_split:
            # print("AA: ", entry['aa'], " COND: ", condition(entry['aa']))
            sum += 1 if condition(entry['aa']) else 0
        
        return sum

    def get_polar_non_positive_non_negative(self):
        
        """
            gets the number of polar residues which are not positive 
            or negative (in the sequence)
        """ 

        condition = lambda aa: True if _aa_polarity[aa] == 1 else False
        return self._filter_sum(condition=condition)

    def get_polar(self):
        
        """
            gets the number of polar residues in the sequence
        """ 

        condition = lambda aa: True if _aa_polarity[aa] == 1 or _aa_polarity[aa] == 2 or _aa_polarity[aa] == 3 else False
        return self._filter_sum(condition=condition)

    def get_positive(self):

        """
            gets the number of positively charged residues in the sequence
        """ 

        condition = lambda aa: True if _aa_polarity[aa] == 2 else False
        return self._filter_sum(condition=condition)

    def get_negative(self):

        """
            gets the number of negatively charged residues in the sequence
        """

        condition = lambda aa: True if _aa_polarity[aa] == 3 else False
        return self._filter_sum(condition=condition)

    def get_hydrophobic(self):

        """
            gets the number of hydrophobic residues in the sequence
        """

        condition = lambda aa: True if _aa_polarity[aa] == 0 else False
        return self._filter_sum(condition=condition)
    
    def get_aromatic(self):

        """
            gets the number of aromatic residues in the sequence
        """

        condition = lambda aa: True if _aa_aromatic[aa] == 1 else False
        return self._filter_sum(condition=condition)
    
    def get_non_aromatic(self):

        """
            gets the number of non-aromatic residues in the sequence
        """

        condition = lambda aa: True if _aa_aromatic[aa] == 0 else False
        return self._filter_sum(condition=condition)
    
    def get_hydroxylic(self):

        """
            gets the number of hydroxylic residues in the sequence
        """

        condition = lambda aa: True if _aa_hydroxilic[aa] == 1 else False
        return self._filter_sum(condition=condition)
    
    def get_amidic(self):

        """
            gets the number of amidic groups in the sequence
        """

        condition = lambda aa: True if _aa_amidic[aa] == 1 else False
        return self._filter_sum(condition=condition)
    
    def get_aliphatic(self):

        """
            gets the number of aliphatic in the sequence
        """

        condition = lambda aa: True if _aa_aliphatic[aa] == 1 else False
        return self._filter_sum(condition=condition)

    def get_essential(self):

        """
            gets the number of essential residues in the sequence
        """

        condition = lambda aa: True if _aa_essential[aa] == 1 else False
        return self._filter_sum(condition=condition)

    def get_hydrophobic(self):

        """
            gets the number of hydrophobic residues in the sequence
        """

        condition = lambda aa: True if _aa_hydrophobicity[aa] > self._hydrophobicity_thresh else False
        return self._filter_sum(condition=condition)
    
    # do I calculate total hydrophobicity?

    def __aggregate_calc(self, __dict_to_map):

        return reduce(lambda x, y: x + y, map(lambda x: __dict_to_map[x['aa']], self.epitope_split))
    
    def get_total_hydrophobicity(self):

        """
            gets the hydrophobicity of the sequence
        """

        return self.__aggregate_calc(_aa_hydrophobicity)
    
    def get_total_hydration_energy(self):

        """
            gets the hydration-energy of the sequence
        """

        return self.__aggregate_calc(_aa_hydration_energy)
    
    def get_large_chains(self):

        """
            gets the number of large residues in the sequence
        """        

        condition = lambda aa: True if _aa_volume_sidechains[aa] > self._volume_thresh else False
        return self._filter_sum(condition=condition)
    
    def get_total_sidechain_volume(self):

        """
            gets the seidechain volume of the sequence
        """

        return self.__aggregate_calc(_aa_volume_sidechains)
    
    def get_total_polarity(self):

        """
            gets the total polarity of the sequence
        """

        return self.__aggregate_calc(_aa_polarity_quant)
    
    def get_total_sasa(self):

        """
            gets the total SASA of the sequence
        """

        return self.__aggregate_calc(_aa_sasa)

    # GRAVY

    # Flexibility

    # Secondary Structurr Fraction

    # Exposed Surface

if __name__ == '__main__':

    trial_seq = 'W155, S156, I157, Q158, N159, C203, Y204, D205, M206, K207, T208, T209, C210'
    analy = Analyze(trial_seq)
    
    print("NUM_SULPH: ", analy.get_num_sulphur())
    print("GET_AA: ", analy.get_num_aa('T'))
    print("TOTAL_WT: ", analy.get_total_weight())
    
    print("HP: ", analy.get_hydrophobic())
    print("POLAR: ", analy.get_polar())
    print("ONLY_POLAR: ", analy.get_polar_non_positive_non_negative())
    print("POSITIVE: ", analy.get_positive())

    print("HYDROXYLIC: ", analy.get_hydroxylic())
    print("AMIDIC: ", analy.get_amidic())
    print("AROMATIC: ", analy.get_aromatic())
    print("NON_AROMATIC: ", analy.get_non_aromatic())
    print("ALIPHATIC: ", analy.get_aliphatic())
    print("ESSENTIAL: ", analy.get_essential())
    print("SASA: ", analy.get_total_sasa())