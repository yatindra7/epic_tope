from sequence_gen import RandomGenerator
from sequence_get import CSV, XLSX
from sequence_parse import AAIndex, CustomIndex

# getting the epitopes data
epitopes = XLSX()

# getting the chikunguniya data
chkngn = CSV(filename='data/disease/Chikungunya.csv')
chkngn = list(chkngn(col=2))

chkngn_aaidx = AAIndex(chkngn)
chkngn_cstmidx = CustomIndex(chkngn)

# getting the dengue data
dengue = CSV(filename='data/disease/Dengue.csv')
dengue = list(dengue(col=2))

dengue_aaidx = AAIndex(dengue)
dengue_cstmidx = CustomIndex(dengue)

# getting the malaria data
malaria = CSV(filename='data/disease/Chikungunya.csv')
malaria = list(malaria(col=2))

malaria_aaidx = AAIndex(malaria)
malaria_cstmidx = CustomIndex(malaria)