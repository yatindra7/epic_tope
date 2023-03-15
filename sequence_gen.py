from random import randint

class RandomGenerator:
  
  base_sequence = None
  base_seq_size = None
  
  def __init__(self, base_sequence):
    """
      adds the sequence as a list of AAs
    """
    self.base_sequence = list(enumerate(base_sequence))
    self.base_seq_size = len(base_sequence)

  def sequence(self, n):
    """
      generate with n% as the number of residues
    """
    num_pos_to_select = int( base_seq_size * (n / 100) )
    return map(lambda x: {'aa': base_sequence[x], 'idx': x}, [randint(0, base_seq_size - 1) for i in range(num_pos_to_select)])
    
