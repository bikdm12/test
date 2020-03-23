import os
import csv
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

class RodeoOutput(object):
    
    def __init__(self, table):
        
        self.query = table[0][0]
        self.table = table
    
    def table_proccessing(self, bg_domains, n):
    
        self.biosynthetic_genes = []
                
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
                    self.operon_domains = operon_domain_buffer
                
                operon_buffer = [row[3]]
                operon_domain_buffer = [domain]        
            
            prev_end = max(int(row[4]), int(row[5]))
            prev_strand = row[6]
                     
        if self.query in operon_buffer:
            self.operon_accs = operon_buffer
            self.operon_domains = operon_domain_buffer



def rodeo_output_iterator(rod_dir, rod_dir_type):
    
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
    with open(file_path, 'r') as infile:
        inlist = [x[0] for x in csv.reader(infile)]
    
    base_name = file_path.rstrip('.txt')
    
    for i in range(len(inlist)/n + 1):
        name = base_name + '_' + str(i * n) + '-' + str((i + 1) * n) + '.txt'
        print name
        with open(name, 'w') as outf:
            for j in inlist[i * n : (i + 1) * n]:
                outf.write(j + '\n')