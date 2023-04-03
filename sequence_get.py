import pandas as pd
import re
import time

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager

"""
    This module will be able to have get facility friom an
    Excel diocument as well as directly from UniProt for 
    important artifacts
"""

# exceptions
EXCEP = '[EXCEP] An exception occured:'
FILE_NOT_OPENED = '[EXCEP] Unable to open file, aborting'

# abort on exception
def _abort_on_exception(excep, msg):

    """
        util to print exception, and abort
    """
    print(EXCEP, excep)
    print(msg)
            
    # aborting
    exit()

class XLSX:

    # class vars
    filename = 'data/Discontinuous.xlsx'
    file_handle = None
    rows = None
    cols = None

    def __init__(self, filename = 'data/Discontinuous.xlsx'):

        self.filename = filename

        # handling the incorrect filename
        try:
            self.file_handle = pd.read_excel(filename)
        except Exception as ex:
            _abort_on_exception(ex, FILE_NOT_OPENED)
        
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

class EXPASY:

    # class vars
    filename = 'data/Discontinuous.xlsx'
    file_handle = None
    rows = None
    cols = None

    def __init__(self, filename = 'data/Discontinuous.xlsx'):

        self.filename = filename

        # handling the incorrect filename
        try:
            self.file_handle = pd.read_excel(filename)
        except Exception as ex:
            _abort_on_exception(ex, FILE_NOT_OPENED)
        
        self.rows = len(self.file_handle.index)
        self.cols = len(self.file_handle.columns)

    def __extract_data(self, paragraph, i):
        
        # Extract the desired information using regular expressions
        aa_count = int(re.search(r"Number of amino acids: (\d+)", paragraph).group(1))
        
        mw = float(re.search(r"Molecular weight: ([\d.]+)", paragraph).group(1))
        pi = float(re.search(r"Theoretical pI: ([\d.]+)", paragraph).group(1))

        # Print the extracted information
        print(f"Number of amino acids: {aa_count}")
        print(f"Molecular weight: {mw}")
        print(f"Theoretical pI: {pi}")

        colheads = ['mw', 'aa_count', 'pi']
        coldata = [mw, aa_count, pi]

        # Extract the amino acid composition using regular expressions
        aa_composition = {}
        aa_pattern = r"([A-Z][a-z]{2}) \(([A-Z])\) +(\d+) +([\d.]+)%"
        aa_matches = re.findall(aa_pattern, paragraph)
        for aa in aa_matches:
            aa_code = aa[1]
            aa_count = int(aa[2])
            aa_percentage = float(aa[3])
            aa_composition[aa_code] = {"percentage": aa_percentage, "count": aa_count}
            colheads.append(aa_code)
            coldata.append(aa_percentage)


        # Print the amino acid composition
        print("Amino acid composition:")
        for aa in aa_composition:
            print(f"{aa}: percentage = {aa_composition[aa]['percentage']}%")

        neg_charge = int(aa_composition['D']['count']) + int(aa_composition['E']['count'])
        pos_charge = int(aa_composition['R']['count']) + int(aa_composition['K']['count'])

        coldata.append(neg_charge)
        coldata.append(pos_charge)

        colheads.append('total_neg_res')
        # coldata.append(total_neg_charged_res)
        colheads.append('total_pos_res')
        # coldata.append(total_pos_charged_res)

        # Print the total number of negatively charged and positively charged residues
        print(f"Total number of negatively charged residues (Asp + Glu): {neg_charge}")
        print(f"Total number of positively charged residues (Arg + Lys): {pos_charge}")

        # Extract the important data using regular expressions
        atom_composition = {}
        atom_pattern = r"([A-Z][a-z]*) +([A-Z]+) +(\d+)"
        atom_matches = re.findall(atom_pattern, paragraph)
        for atom in atom_matches:
            atom_name = atom[0]
            atom_symbol = atom[1]
            atom_count = int(atom[2])
            atom_composition[atom_name] = {"symbol": atom_symbol, "count": atom_count}
            colheads.append(atom_name+' count')
            coldata.append(atom_count)

        # Print the atomic composition
        print(f"Atomic composition:")
        for atom in atom_composition:
            print(f"{atom_composition[atom]['symbol']} {atom_composition[atom]['count']}")

        # Extract the important data using regular expressions
        total_atoms_pattern = r"Total number of atoms: (\d+)"
        total_atoms_matches = re.findall(total_atoms_pattern, paragraph)
        total_atoms = int(total_atoms_matches[0])

        colheads.append('total_atoms')
        coldata.append(total_atoms)

        # Print the total number of atoms
        print(f"Total number of atoms: {total_atoms}")

        pattern = r"Ext\.? coeff(?:icient)?\s+(\d+\.\d+)\s+"

        # Search for matches in the text
        ext_coeff_matches = re.findall(pattern, paragraph)

        colheads.append('ext_coeff')
        
        if ext_coeff_matches:
            print("Extinction coefficients found:")
            for match in ext_coeff_matches:
                print(match)
                coldata.append(match)
        else:
            coldata.append('NIL')
            print("No extinction coefficients found.")

        mammalian_reticulocytes_pattern = r"(\d+(?:\.\d+)?) hours \(mammalian reticulocytes, in vitro\)"
        mammalian_reticulocytes_match = re.search(mammalian_reticulocytes_pattern, paragraph)

        colheads.append('half_life')
        if mammalian_reticulocytes_match:
            mammalian_reticulocytes_half_life = float(mammalian_reticulocytes_match.group(1))
            print("The half-life in mammalian reticulocytes is:", mammalian_reticulocytes_half_life, "hours.")
            coldata.append(mammalian_reticulocytes_half_life)
        else:
            coldata.append('NIL')
            print("No half-life information found for mammalian reticulocytes.")

        instability_regex = r"instability index \(II\) is computed to be ([\d\.]+)"
        aliphatic_regex = r"Aliphatic index: ([\d\.]+)"
        gravy_regex = r"Grand average of hydropathicity \(GRAVY\): ([\d\.-]+)"

        # Extract values using regular expressions
        instability_match = re.search(instability_regex, paragraph)
        aliphatic_match = re.search(aliphatic_regex, paragraph)
        gravy_match = re.search(gravy_regex, paragraph)
        
        colheads.append('instability_index')
        # Print results
        if instability_match:
            instability_index = float(instability_match.group(1))
            coldata.append(instability_index)
            print(f"Instability index: {instability_index}")
        else:
            coldata.append('NIL')
            print("No instability index found.")

        colheads.append('aliphatic_index')
        if aliphatic_match:
            aliphatic_index = float(aliphatic_match.group(1))
            print(f"Aliphatic index: {aliphatic_index}")
            coldata.append(aliphatic_index)
        else:
            coldata.append('NIL')
            print("No Aliphatic index found.")
        
        colheads.append('GRAVY')
        if gravy_match:
            gravy_match = float(gravy_match.group(1))
            coldata.append(gravy_match)
            print(f"Gravy Match: {gravy_match}")
        else:
            coldata.append('NIL')
            print("No Gravy Match found.")
        
        self.file_handle.loc[i, colheads] = coldata
    
    def __call__(self):

        """
            retrieves all data, from the web
        """

        with webdriver.Chrome(ChromeDriverManager().install()) as driver:
            for idx, value in enumerate(self.file_handle.values):
                
                print(f"starting for {value}")
                # extracting the sequence
                sequence = value[0]
                
                # setting up webdriver
                wait = WebDriverWait(driver, 10)
                driver.get('https://web.expasy.org/protparam/')

                seq_field = wait.until(EC.presence_of_element_located((By.XPATH, '//*[@id="sib_body"]/form/textarea')))
                seq_field.send_keys(sequence)

                submit_btn = wait.until(EC.presence_of_element_located((By.XPATH, '//*[@id="sib_body"]/form/p[1]/input[2]')))
                submit_btn.click()

                data_div = wait.until(EC.presence_of_element_located((By.XPATH, '//*[@id="sib_body"]/pre[2]')))

                print(f"retrieved, sleeping")
                time.sleep(10)

                print(f"extracting, done sleeping")
                self.__extract_data(data_div.text , idx)

        row_get = self.file_handle.iloc
        return (row_get[row][0]
                    for row in range(1, self.rows) )

if __name__ == '__main__':
    xls = XLSX()

    print("Type: ")
    print(type(xls()))
    print("An overview: ")
    print(list(xls())[:50])

    # exp = EXPASY()
    # print("Type: ")
    # ret = exp()
    # print(type(ret))
    # print("An overview: ")
    # print(list(ret)[:50])
