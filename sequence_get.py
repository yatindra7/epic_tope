import pandas as pd

"""
    This module will be able to have get facility friom an
    Excel diocument as well as directly from UniProt for 
    important artifacts
"""

# exceptions
EXCEP = '[EXCEP] An exception occured:'
FILE_NOT_OPENED = '[EXCEP] Unable to open file, aborting'

class XLSX:

    # class vars
    filename = 'data/Discontinuous.xlsx'
    file_handle = None
    rows = None
    cols = None

    def __init__(self, filename = 'Discontinuous.xlsx'):

        self.filename = filename

        # handling the incorrect filename
        try:
            self.file_handle = pd.read_excel(filename)
        except Exception as ex:
            print(EXCEP, ex)
            print(FILE_NOT_OPENED)
            
            # aborting
            exit()
        
        self.rows = len(self.file_handle.index)
        self.cols = len(self.file_handle.columns)

    def __call__(self):

        """
            reads an xlsx file of the format provided in week-1
            `TBC`
        """

        row_get = self.file_handle.iloc
        return (row_get[row][0]
                    for row in range(1, self.rows) )


if __name__ == '__main__':
    xls = XLSX()

    print("Type: ")
    print(type(xls()))
    print("An overview: ")
    print(list(xls())[:50])
    