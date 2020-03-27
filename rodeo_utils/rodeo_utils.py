"""Main module of the package.
Provides RODEO output iterator and a class to store information for a single query.
Also contains a number of simple, but useful functons. 

"""
import os
import csv
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

class RodeoOutput(object):
    """Object that stores rodeo output for a single query as attributes.

    """    
    def __init__(self, table):
        self.query = table[0][0]
        """str: Accession number of the query protein."""
        self.table = table
        """list of list: Corresponding part of rodeo output main_co_occur.csv table. Contains rows of the table as lists."""
        # TODO: Add attributes: lists for left and rigth query neighbours.
    
    def table_proccessing(self, bg_domains, n):
        # TODO: Split into two methods. 
        self.biosynthetic_genes = []
        """list of str: List of proteins (accession numbers) that contain domains from a given list"""
                
        operon_buffer = []
        operon_domain_buffer = []
        
        prev_end = 0
        prev_strand = ''
        
        for row in self.table:

            if len(row) > 7:
                domain = row[7]
                if row[7] in bg_domains:
                    self.biosynthetic_genes.append(row[3])
            else:
                domain = 'no_match'

            start = min(int(row[4]), int(row[5]))
            if (row[6] == prev_strand) and (start - prev_end < n):
                operon_buffer.append(row[3])
                operon_domain_buffer.append(domain)
            else:
                if self.query in operon_buffer:
                    self.operon_accs = operon_buffer
                    """list of str: Group of protein accession numbers that includes query. Genes of these proteins form so called quasioperon: they are adjacent and located on the same strand."""
                    self.operon_domains = operon_domain_buffer
                    """list of str: List of Pfam or TIGRFAMs IDs corresponding to proteins from quasioperon."""
                operon_buffer = [row[3]]
                operon_domain_buffer = [domain]        
            
            prev_end = max(int(row[4]), int(row[5]))
            prev_strand = row[6]
                     
        if self.query in operon_buffer:
            self.operon_accs = operon_buffer
            self.operon_domains = operon_domain_buffer

class RodeoDirError(Exception):
    """Exception to raise if rodeo_output_iterator failed to infer the structure of the given directory.

    """
    def __init__(self, rod_dir):
        self.rod_dir = rod_dir
    def __str__(self):
        return 'Failed to infer the structure of %s directory. It must be either RODEO output or rodout/ folder of the RIPPER output.' % self.rod_dir

def rodeo_output_iterator(rod_dir):
    """Iterate over RODEO output.

    Parameters
    ----------
    rod_dir : str
        Path to RODEO output.

    Yields
    ------
    RodeoOutput
        An object with the results for a single query.

    Raises
    --------
    RodeoDirError
        If the given directory is neither RODEO or RIPPER output

    """
    rod_dir_type = None
    if 'main_co_occur.csv' in os.listdir(rod_dir):
        rod_dir_type = 'RODEO'
    else:
        for folder in os.listdir(rod_dir):
            if os.path.isdir('%s/%s' % (rod_dir, folder)):
                if 'main_co_occur.csv' in os.listdir('%s/%s' % (rod_dir, folder)):
                    rod_dir_type = 'RIPPER'
                    break
        if rod_dir_type is None:
            raise RodeoDirError(rod_dir)
    
    if rod_dir_type == 'RODEO':
        # It woulod be easier to use Pandas, but this way is more memory-efficient.
        # It may be crucial when working with large datasets
        with open('%s/main_co_occur.csv' % rod_dir) as infile:
            infile.next()
            prev_seed = None
            table = []
            for row in csv.reader(infile):
                
                if prev_seed is None:
                    prev_seed = row[0]
                
                if row[0] == prev_seed:
                    table.append(row)
                
                else:
                    yield RodeoOutput(table)
                    prev_seed = row[0]
                    table = [row]
                
            yield RodeoOutput(table)
            
    elif rod_dir_type == 'RIPPER':
        for folder in os.listdir(rod_dir):
            if 'main_co_occur.csv' in os.listdir('%s/%s' % (rod_dir, folder)):
                with open('%s/%s/main_co_occur.csv' % (rod_dir, folder)) as infile:
                    infile.next()
                    infile_as_list = list(csv.reader(infile))
                    if len(infile_as_list) != 0:
                        yield RodeoOutput(infile_as_list)



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