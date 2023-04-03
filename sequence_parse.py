from functools import reduce
import typing
import json

import pandas as pd

from aadict import *

# from collections import namedtuple
# EpitopeEntry = namedtuple('EpitopeEntry', ['aa', 'idx'])

# exceptions
EXCEP = '[EXCEP] An exception occured:'
FILE_NOT_OPENED = '[EXCEP] Unable to open file, aborting'
PANDAS_TO_CSV = '[EXCEP] Pandas to_csv failed'

# abort on exception
def _abort_on_exception(excep, msg):

    """
        util to print exception, and abort
    """
    print(EXCEP, excep)
    print(msg)
            
    # aborting
    exit()

class AAIndex:

    epitope_sequence: list = None
    epitope_split: list = None
    aaidx: dict = None
    frame: pd.DataFrame = None

    def __init__(self, epitope_sequence: list):
        self.epitope_sequence = epitope_sequence
        self.epitope_split = [ [{"aa":entry.strip()[0], "idx":int(entry.strip()[1:])}
                                                                        for entry in _sequence.split(',')]
                                                                            for _sequence in epitope_sequence ]

        # get the data read
        with open('data/aaindex1.json', 'r') as aaindex_handle:
            self.aaidx = json.load(aaindex_handle)
            self.aaidx = {key: value['values'] for key, value in self.aaidx.items()}

    def __aggregate_calc(self, __dict_to_map: dict, _epitope_split: list):

        return reduce(lambda x, y: x + y, map(lambda x: __dict_to_map[x['aa']], _epitope_split))

    def __call__(self, fname='data/epitope.aaindex.csv'):

        _frame = [ [ self.__aggregate_calc(_dict, e)
                                for _dict in self.aaidx.values() ]
                                            for e in self.epitope_split ]
        self.frame = pd.DataFrame(_frame, columns = self.aaidx.keys())

        try:
            self.frame.to_csv(fname)  
        except:
            _abort_on_exception(PANDAS_TO_CSV)
            return -1

        return 0

class CustomIndex:

    epitope_sequence: list = None
    epitope_split: list = None
    frame: pd.DataFrame = None
    _hydrophobicity_thresh = 2
    _volume_thresh = 70

    def __init__(self, epitope_sequence: list):
        self.epitope_sequence = epitope_sequence
        self.epitope_split = [ [{"aa":entry.strip()[0], "idx":int(entry.strip()[1:])}
                                                                        for entry in _sequence.split(',')]
                                                                            for _sequence in epitope_sequence ]

    def _filter_sum(self, condition):

        # sum = 0
        # for entry in self.epitope_split:
        #     # print("AA: ", entry['aa'], " COND: ", condition(entry['aa']))
        #     sum += 1 if condition(entry['aa']) else 0
    
        return [reduce(lambda x, y: x + y, map(lambda x: 1 if condition(x['aa']) else 0, epitope), 0) for epitope in self.epitope_split]

    def get_num_sulphur(self):

        """
            getting the number of sulphurs which can create
            a disulphurous bond.
        """

        condition = lambda aa: True if aa == 'C' else False
        return 'num_sulphur', self._filter_sum(condition)

    def get_num_aa(self, aa):

        """
            getting the number of a particular amino acid from
            the code given.
        """

        condition = lambda x: True if x == aa else False
        return 'num_aa', self._filter_sum(condition=condition)

    def get_polar_non_positive_non_negative(self):
        
        """
            gets the number of polar residues which are not positive 
            or negative (in the sequence)
        """ 

        condition = lambda aa: True if aa_polarity[aa] == 1 else False
        return 'num_polar_non_positive_non_negative', self._filter_sum(condition=condition)

    def get_polar(self):
        
        """
            gets the number of polar residues in the sequence
        """ 

        condition = lambda aa: True if aa_polarity[aa] == 1 or aa_polarity[aa] == 2 or aa_polarity[aa] == 3 else False
        return 'num_polar', self._filter_sum(condition=condition)

    def get_positive(self):

        """
            gets the number of positively charged residues in the sequence
        """ 

        condition = lambda aa: True if aa_polarity[aa] == 2 else False
        return 'num_positive', self._filter_sum(condition=condition)

    def get_negative(self):

        """
            gets the number of negatively charged residues in the sequence
        """

        condition = lambda aa: True if aa_polarity[aa] == 3 else False
        return 'num_negative', self._filter_sum(condition=condition)

    def get_hydrophobic(self):

        """
            gets the number of hydrophobic residues in the sequence
        """

        condition = lambda aa: True if aa_polarity[aa] == 0 else False
        return 'num_hydrophobic', self._filter_sum(condition=condition)
    
    def get_aromatic(self):

        """
            gets the number of aromatic residues in the sequence
        """

        condition = lambda aa: True if aa_aromatic[aa] == 1 else False
        return 'num_aromatic', self._filter_sum(condition=condition)
    
    def get_non_aromatic(self):

        """
            gets the number of non-aromatic residues in the sequence
        """

        condition = lambda aa: True if aa_aromatic[aa] == 0 else False
        return 'num_non_aromatic', self._filter_sum(condition=condition)
    
    def get_hydroxylic(self):

        """
            gets the number of hydroxylic residues in the sequence
        """

        condition = lambda aa: True if aa_hydroxilic[aa] == 1 else False
        return 'num_hydroxylic', self._filter_sum(condition=condition)
    
    def get_amidic(self):

        """
            gets the number of amidic groups in the sequence
        """

        condition = lambda aa: True if aa_amidic[aa] == 1 else False
        return 'num_amidic', self._filter_sum(condition=condition)
    
    def get_aliphatic(self):

        """
            gets the number of aliphatic in the sequence
        """

        condition = lambda aa: True if aa_aliphatic[aa] == 1 else False
        return 'num_aliphatic', self._filter_sum(condition=condition)

    def get_essential(self):

        """
            gets the number of essential residues in the sequence
        """

        condition = lambda aa: True if aa_essential[aa] == 1 else False
        return 'num_essential', self._filter_sum(condition=condition)

    def get_hydrophobic(self):

        """
            gets the number of hydrophobic residues in the sequence
        """

        condition = lambda aa: True if aa_hydrophobicity[aa] > self._hydrophobicity_thresh else False
        return 'num_hydrophobic', self._filter_sum(condition=condition)
    
    # do I calculate total hydrophobicity?

    def __aggregate_calc(self, __dict_to_map):

        return [reduce(lambda x, y: x + y, map(lambda x: __dict_to_map[x['aa']], epitope), 0) for epitope in self.epitope_split]

    def get_total_weight(self):
        
        """
            calculate the total moklecular weight of the epitope
        """
        
        return 'weight', self.__aggregate_calc(aa_weights)

    def get_total_hydrophobicity(self):

        """
            gets the hydrophobicity of the sequence
        """

        return 'hydrophobicity', self.__aggregate_calc(aa_hydrophobicity)
    
    def get_total_hydration_energy(self):

        """
            gets the hydration-energy of the sequence
        """

        return 'hydration_energy', self.__aggregate_calc(aa_hydration_energy)
    
    def get_large_chains(self):

        """
            gets the number of large residues in the sequence
        """        

        condition = lambda aa: True if aa_volume_sidechains[aa] > self._volume_thresh else False
        return 'num_large_chains', self._filter_sum(condition=condition)
    
    def get_total_sidechain_volume(self):

        """
            gets the seidechain volume of the sequence
        """

        return 'sidechain_volume', self.__aggregate_calc(aa_volume_sidechains)
    
    def get_total_polarity(self):

        """
            gets the total polarity of the sequence
        """

        return 'polarity', self.__aggregate_calc(aa_polarity_quant)
    
    def get_total_sasa(self):

        """
            gets the total SASA of the sequence
        """

        return 'sasa', self.__aggregate_calc(aa_sasa)

    # GRAVY

    # Flexibility

    # Secondary Structurr Fraction

    # Exposed Surface


    # def __aaidx_extract(self, feature):

    #     """
    #         util to aggregate a given feature from `aaindex`
    #     """

    #     _record = prop_dict[feature].values
    #     # print("[LOG] feature:", feature, " _record:", _record)
    #     return feature, self.__aggregate_calc(_record)

    # def get_total_aaindex(self):

    #     """
    #         extracts the value of all the
    #     """

    #     return [ self.__aaidx_extract( prop ) for prop in property_codes ]
    
    def __call__(self, fname='data/epitope.cstmidx.csv') -> int:

        """
            creates a dataframe containing all the features
            as different columns is generated for the given data
        """

        _data = [getattr(self, mthd)() 
                            for mthd in self.__dir__()
                                            if callable( getattr(self, mthd) ) and
                                                mthd.startswith('get_') and 
                                                    not mthd.endswith('get_num_aa')]
        _data = {pt[0]: pt[1] for pt in _data}

        self.frame = pd.DataFrame(_data)
        try:
            self.frame.to_csv(fname)
        except:
            _abort_on_exception(PANDAS_TO_CSV)
            return -1

        return 0

if __name__ == '__main__':

    trial_seq = ['W155, S156, I157, Q158, N159, C203, Y204, D205, M206, K207, T208, T209, C210','E53, R98, E107', 'E170, E172', 'L87, L88, V90, R91, S92, E132', 'W155, S156, I157, Q158, N159, C203, Y204, D205, M206, K207, T208, T209, C210', 'D429, R441, D454', 'D429, R441, D454, D463', 'D429, R441, E452, D454', 'D454', 'P462', 'G104, E126, W231', 'G104, G106, E126, W231', 'G104, G106, L107, W231', 'G106, L107', 'T76, G104, G106, L107, E126, W231', 'T76, G106, L107, W231', 'G302', 'I126', 'K179', 'Q52, I126, K136, S275', 'I258', 'N245', 'P199', 'E472, K514', 'G138, R146', 'G138, T144, S145, R146', 'S145, D153, T197, N201, T206, L234', 'S145, T147, K165, P170, T197, T206', 'K307, T330', 'Y302, S306, K307, A308, F309, T330, G331, T332, D333, A365, T366, A367, N368, G389, E390, Q391', 'L271', 'A10', 'A10, A46', 'A10, A46, A95', 'A10, A46, E51, K91', 'A10, E11, Y12, G33, K34, R35, E36, E51, W88, A95', 'A10, Y12, G33, K34, R35, E36, E51, W88, K91, A95', 'A46', 'A198', 'D63', 'G144', 'G158', 'N54', 'S145', 'S157', 'S193', 'S199', 'S205', 'A198', 'D63', 'G129']
    analy = CustomIndex(trial_seq)
    
    print('CALL: ', analy())
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

    lst_seq = ['E53, R98, E107', 'E170, E172', 'L87, L88, V90, R91, S92, E132', 'W155, S156, I157, Q158, N159, C203, Y204, D205, M206, K207, T208, T209, C210', 'D429, R441, D454', 'D429, R441, D454, D463', 'D429, R441, E452, D454', 'D454', 'P462', 'G104, E126, W231', 'G104, G106, E126, W231', 'G104, G106, L107, W231', 'G106, L107', 'T76, G104, G106, L107, E126, W231', 'T76, G106, L107, W231', 'G302', 'I126', 'K179', 'Q52, I126, K136, S275', 'I258', 'N245', 'P199', 'E472, K514', 'G138, R146', 'G138, T144, S145, R146', 'S145, D153, T197, N201, T206, L234', 'S145, T147, K165, P170, T197, T206', 'K307, T330', 'Y302, S306, K307, A308, F309, T330, G331, T332, D333, A365, T366, A367, N368, G389, E390, Q391', 'L271', 'A10', 'A10, A46', 'A10, A46, A95', 'A10, A46, E51, K91', 'A10, E11, Y12, G33, K34, R35, E36, E51, W88, A95', 'A10, Y12, G33, K34, R35, E36, E51, W88, K91, A95', 'A46', 'A198', 'D63', 'G144', 'G158', 'N54', 'S145', 'S157', 'S193', 'S199', 'S205', 'A198', 'D63', 'G129']
    aai = AAIndex(lst_seq)
    aai()