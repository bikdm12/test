import os
import csv
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation



def check_if_border(feature, operon_borders):
    
    prot_id_regexp = re.compile('[A-Z]{2}_[0-9]+\.[0-9]')
        
    start_id, end_id = operon_borders
    
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



def convert_gbk(gb_dir, gb_out_dir, rodeo_output, bg_domains, n, product_class):
    
    rodeo_output.table_proccessing(bg_domains, n)
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
    
    except Exception as e:
        print e