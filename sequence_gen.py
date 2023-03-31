from random import randint

class RandomGenerator:
  
  base_sequence = None
  base_seq_size = None
  
  def __init__(self, base_sequence):
    """
      adds the sequence as a list of AAs
    """
    self.base_sequence = [{'aa': aa[1], 'idx': aa[0]} for aa in list(enumerate(base_sequence))]
    self.base_seq_size = len(base_sequence)

  def _check(self, seq):
    """
      checks if the thresholds match
    """

    # TODO

    return True

  def __call__(self, n):
    """
      generate with n% as the number of residues
    """
    num_pos_to_select = int( self.base_seq_size * (n / 100) )
    random_sequence = list( map(lambda x: self.base_sequence[x], [randint(0, self.base_seq_size - 1) for i in range(num_pos_to_select)]) )
    
    return random_sequence if self._check(random_sequence) else self.__call__(n)
    
if __name__ == '__main__':
  seq = 'AAGGBCQER'
  randgen = RandomGenerator(seq)
  print([randgen(40) for i in range(10)])