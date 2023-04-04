from sequence_get import CSV, XLSX
from sequence_parse import AAIndex, CustomIndex

xls = XLSX()
#xls = CSV(filename='data/disease/Malaria.csv')
#list_of_epitope = list(xls(col=2))
list_of_epitope = list(xls())

print([el for el in list_of_epitope if len(el) < 4])

aaidx = AAIndex(list_of_epitope)
custidx = CustomIndex(list_of_epitope)

aaidx()
custidx()