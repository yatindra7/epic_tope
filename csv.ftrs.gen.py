from sequence_get import CSV, XLSX
from sequence_parse import AAIndex, CustomIndex

xls = XLSX()
xls = CSV(filename='data/fake.csv')
#list_of_epitope = list(xls(col=2))
list_of_epitope = list(xls())

print([el for el in list_of_epitope if len(el) < 4])

aaidx = AAIndex(list_of_epitope)
custidx = CustomIndex(list_of_epitope)

aaidx(fname='data/fake.aaidx.csv', save=True)
custidx(fname='data/fake.custidx.csv', save=True)