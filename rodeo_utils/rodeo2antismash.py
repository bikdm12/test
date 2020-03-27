"""Convert sequence in a common genbank format to genbank output of antiSMASH.

The main function is convert_gbk(...).

"""

import os
import csv
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation



def check_if_border(feature, operon_borders):
    """Return start/end coordinate of the SeqFeature if its accession is the first/last in the given tuple.

    Parameters
    ----------
    feature : SeqFeature
        Directory with input genbank files.
    operon_borders : str
        Directory to store the output.
        
    Returns
    -------
    tuple(str, int)
        The first value is either 'end' or 'start'.
        The second value is the corresponding coordinate.
        
        If SeqFeature's ID is not in the given tuple, returns None. 
    
    """
    prot_id_regexp = re.compile('[A-Z]{2}_[0-9]+\.[0-9]')
        
    start_id = operon_borders[0]
    end_id = operon_borders[-1]
    
    if 'protein_id' in feature.qualifiers:
                    
        if feature.qualifiers['protein_id'][0]  == start_id:
            return ('start', feature.location.start + 1) # genbank is 1-based, python is 0-based
                           
        if feature.qualifiers['protein_id'][0]  == end_id:
            return ('end', int(feature.location.end))
                
    elif 'pseudo' in feature.qualifiers:        
        if 'inference' in feature.qualifiers:
                        
            inference_prot_id_search = prot_id_regexp.search(feature.qualifiers['inference'][0])
                        
            if inference_prot_id_search is not None:
                            
                inference_prot_id = inference_prot_id_search.group(0)
                            
                if inference_prot_id  == start_id:
                    return('start', feature.location.start + 1)
                        
                if inference_prot_id  == end_id:
                    return ('end', int(feature.location.end))
                    
                else:
                    if feature.qualifiers['locus_tag'][0]  == start_id:
                        return ('start', feature.location.start + 1)

                    if feature.qualifiers['locus_tag'][0] == end_id:
                        return ('end', int(feature.location.end))
                    
    return None



def convert_gbk(gb_dir, gb_out_dir, rodeo_output, bg_domains, max_intergenic_distance = 100, product_class = 'thiopeptide'):
    """Convert a common genbank file to the genbank antiSMASH output.
    
    Adds a feature “cluster” with information about the class of the product.
    The coordinates of this feature are boundaries of the group of adjacent genes on the same strand that includes RODEO query.
    Marks genes with given domains as biosynthetic.

    Parameters
    ----------
    gb_dir : str
        Directory with input genbank files.
    gb_out_dir : str
        Directory to store the output.
    rodeo_output: RodeoOutput
        RODEO output to use as a reference.
    bg_domains : list
        List of Pfam or TIGRFAMs IDs for domains that are important for your product biosynthesis.
    max_intergenic_distance : int, optional
        Maximum distance (nt) between genes within the biosynthetic gene cluster (default: 100).
    product_class : string, optional
        A putative class of the final product (default: thiopeptide).
    
    Returns
    -------
    bool
        True if successful, False otherwise.
    
    """
    rodeo_output.table_proccessing(bg_domains, max_intergenic_distance)
    operon_border_accs = (rodeo_output.operon_accs[0], rodeo_output.operon_accs[-1])
    biosynthetic_genes = rodeo_output.biosynthetic_genes
    
    contig_edge = False
    prot_id = rodeo_output.query
    try:    
        genbank = SeqIO.parse('%s%s.gbk' % (gb_dir, prot_id), 'genbank')
        for record in genbank: # Every file is expected to contain only one record
            
            cluster_coords = OrderedDict([('start', 1), ('end', len(record))])
            
            for feature in record.features:
                if feature.type == 'CDS':
                    
                    border_check = check_if_border(feature, operon_border_accs)
                    if border_check is not None:
                        cluster_coords[border_check[0]] = border_check[1]

                    if 'protein_id' in feature.qualifiers:
                        if feature.qualifiers['protein_id'][0] in biosynthetic_genes:
                            feature.qualifiers['sec_met'] = ['Kind: biosynthetic']
            
            start, end = cluster_coords.values()
            cluster_location = FeatureLocation(start, end)
            cluster_qualifiers = OrderedDict([('contig_edge', str(contig_edge)), ('product', product_class)])
            cluster = SeqFeature(location = cluster_location, type = 'cluster', qualifiers = cluster_qualifiers)
            record.features = [cluster] + record.features
            
            SeqIO.write(record, '%s%s.gbk' % (gb_out_dir, prot_id), 'genbank')
            return True
    
    except Exception as e:
        print e
        return False