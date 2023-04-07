# epic_tope
An epitope prediction util under development!

## How to use it?

We will be making a Python package out of this soon!

```
from sequence_get import CSV
from sequence_parse import AAIndex, CustomIndex

xls = CSV(filename='data/fake.csv')
list_of_epitope = list(xls())

print([el for el in list_of_epitope if len(el) < 4])

aaidx = AAIndex(list_of_epitope)
custidx = CustomIndex(list_of_epitope)

aaidx(fname='data/fake.aaidx.csv', save=True)
custidx(fname='data/fake.custidx.csv', save=True)
```
the above block shows how we can import from the module:
- `sequence_get`: used for extracting sequences and data
  - `XLS` for xls files
  - `CSV` for csv files
  - `EXPASY` for scraping
- `sequence_parse`: used for calculating features for epitope
  - `AAIndex` using `AAIndex(epitope_list: list)` constructor,
    for a created variable, pass `save` as True, if you wish to save
    the results, and pass the save path
  - `CustomIndex` same as AAIndex
- `sequence_gen`: use for generating random epitopes
  - `RandomGenerator` using `RandomGenerator(protein)`, you can
    create epitopes over an existing protein by calling the created
    object
    ```
    randgen = RandomGenerator()
    val = randgen()
    ```

## notebooks
- `cbp_plot*` notebooks contain plots and analysis of our work
- `cbp_dengue` notebook contains dengue analysis
