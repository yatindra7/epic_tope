from sequence_get import XLSX
from sequence_parse import AAIndex, CustomIndex

xls = XLSX()
list_of_epitope = list(xls())

print([el for el in list_of_epitope if len(el) < 4])

aaidx = AAIndex(list_of_epitope)
custidx = CustomIndex(list_of_epitope)

aaidx()
custidx()