"""This module contains a number of simple, but useful functons for handling RODEO input/output.

"""
import csv

def split_file(file_path, n = 1000):
    """Split a file into chunks with a given number of lines in each.
    
    Resulting files are stored in the same directory as input.
    Useful to split files for RODEO since it doesn't accept more than 1000 accession numbers. 

    Parameters
    ----------
    file_path : str
        Path to the file to split.
    n : int, optional
        The number of lines in each chunk.
    
    """
    with open(file_path, 'r') as infile:
        inlist = [x[0] for x in csv.reader(infile)]
    
    base_name = file_path.rstrip('.txt')
    
    for i in range(len(inlist)/n + 1):
        name = base_name + '_' + str(i * n) + '-' + str((i + 1) * n) + '.txt'
        print name
        with open(name, 'w') as outf:
            for j in inlist[i * n : (i + 1) * n]:
                outf.write(j + '\n')