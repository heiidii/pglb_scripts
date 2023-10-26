
from Bio import AlignIO
from Bio.SeqUtils import seq1
import os
import numpy as np


def write_sequences(path, chain_name='A', outfile='fasta.fasta'):
  '''
  '''
  import glob
  from Bio.PDB import PDBParser
  pdb_files = [t for t in glob.glob(path+'/*.pdb')]
  
  dict_files = {}
  dict_seqs = {}
  for pdb_file in pdb_files:
    id = os.path.basename(pdb_file)[:4]
    dict_files[id] = pdb_file
    parser = PDBParser()
    structure = parser.get_structure(id, pdb_file)
    
    for chain in structure.get_chains():
        if chain_name == chain.id:
          seq = seq1(''.join([residue.resname for residue in chain\
                      if 'CA' in residue]))
          print(seq)
          dict_seqs[id] = seq

  outstr = ['>{}_{}\n{}\n'.format(key, chain_name, dict_seqs[key])
            for key in dict_seqs]
  with open(outfile, 'w') as f:
    f.write('\n'.join(outstr))

def write_readable_lists(positions_dict,filename):
  outf=open(filename,'w')
  for key in positions_dict:
    outf.write('%s\t%s\n' %(key,','.join(positions_dict[key])))
  outf.close()

def write_readable_lists_pdbnumbering(positions_dict, filename):
 outf=open(filename,'w')
 for key in positions_dict:
  pdb_numbering = ['{}{}'.format(t[0], int(t[1:])) for t in positions_dict[key]]
  outf.write('%s,%s\n' %(key,','.join(pdb_numbering)))
 outf.close()

def column_to_sequence(n_col,align,name, include_aa=False):
  align_array = np.array([list(rec) for rec in align], np.character)
  rows= [i for i in range(len(align)) if str(align[i].id)==name]
  assert(len(rows)==1)
  row=rows[0]
  seq_upto_column = [t for t in align[row,:(n_col+1)] if str(t).isalpha()]
  n_seq = len(seq_upto_column)
  if include_aa:
    return n_seq
  return n_seq

def sequence_to_column(align,ref="Dmar_A",lookup=375, include_aa=False):
  align_array = np.array([list(rec) for rec in align], np.character)
  print(align_array.shape)
  rows= [i for i in range(len(align)) if str(align[i].id)==ref]
  assert(len(rows)==1)
  row=rows[0]
  #print('lookup', lookup,row)
  count = 0
  column = -9999
  aa_column = 'X'
  for record in align:
    #print('Rec ',record.id)
    if str(record.id)==ref:
      for i in range(0,align_array.shape[1]):
        if str(align[row,i]).isalpha():
          if count==(lookup-1):
            column=i
            aa_column=align[row,i]
            break
          count+=1
      #print('dat', column,align[row,column-3:column+3],align[row,column])
  if include_aa:
    #print('AAref', aa_column)
    return column, aa_column
  return column

def get_positions(alignment_file='sequences.aln',
                  lookup_pos=[133,136,418,422,367,371,374,375],
                  ref='Dmar_A', include_aa=False):
 columns=[]
 align = AlignIO.read(alignment_file, "clustal")
 for lookup in lookup_pos:
  if include_aa:
    column, aligned_aa = sequence_to_column(align,ref=ref,lookup=lookup, include_aa=include_aa)
  else:
    column = sequence_to_column(align,ref=ref,lookup=lookup)
  columns.append(column) 

 positions_dict={}
 for row_ in range(len(align)):
  positions=[]
  for lookup,column in zip(lookup_pos,columns):
    n_seq = column_to_sequence(column,align,str(align[row_].id),include_aa=include_aa)
    aa_seq = align[row_, column]
    if include_aa:
      positions.append('{}{}'.format(aa_seq, n_seq))
    else:
      positions.append(str(n_seq))
  positions_dict[str(align[row_].id)]=positions
 return positions_dict

def get_positions_m2(alignment_file, out_dir='.',
                      lookup_pos=[133,136,137,418,422,367,371,374,375],
                      suffix=''): #Q/D
 dict_ =  get_positions(alignment_file,
                        lookup_pos=lookup_pos,
                        ref='Dmar_A')
 import pickle
 pickle.dump(dict_,open('{}/residue_positions_m2{}.p'.format(out_dir, suffix),'wb')) 
 import json
 json.dump(dict_,open('{}/residue_positions_m2{}.json'.format(out_dir, suffix),'w'))
 write_readable_lists(dict_,'{}/aligned_positions_m2{}.txt'.format(out_dir, suffix))

def get_positions_m1(): #Q/D
 dict_ =  get_positions(lookup_pos=[52,53,54,136,137,139,142,476,478,479,481],ref='Dmar_A')
 import pickle
 pickle.dump(dict_,open('residue_positions_m1.p','wb'))
 import json
 json.dump(dict_,open('residue_positions_m1.json','w'))
 write_readable_lists(dict_,'aligned_positions_m1.txt')

def get_positions_el5():
 dict_c_term = get_positions(lookup_pos=range(305,326),ref='lari')
 dict_n_term = get_positions(lookup_pos=range(280,305),ref='lari')
 import pickle
 pickle.dump(dict_c_term,open('residue_positions_el5loop_c.p','wb'))
 pickle.dump(dict_n_term,open('residue_positions_el5loop_n.p','wb'))
 import json
 json.dump(dict_c_term,open('residue_positions_el5loop_c.json','w'))
 json.dump(dict_n_term,open('residue_positions_el5loop_n.json','w'))
 write_readable_lists(dict_n_term,'aligned_positions_cl5loop_n.txt')
 write_readable_lists(dict_c_term,'aligned_positions_el5loop_c.txt')


def get_positions_selection(alignment_file,
                            outdir='./',
                            outfile_basename='aligned_positions_selection',
                            ref='Clar'):
  positions_pdbnum = [54, 56, 154, 156, 319, 196, 375, 468, 572,
			463, 464, 465, 482, 483, 484, 331]
  lookup_pos = [t-1 for t in positions_pdbnum]
  dict_pos =  get_positions(alignment_file,
                        lookup_pos=lookup_pos,
                        ref=ref, include_aa=True)
  write_readable_lists_pdbnumbering(dict_pos,f'{outdir}/{outfile_basename}_pdbnum.csv')
  dict_pos =  get_positions(alignment_file,
                        lookup_pos=lookup_pos,
                        ref=ref)
  write_readable_lists(dict_pos, f'{outdir}/{outfile_basename}.txt')


def get_positions_selection_for_TIXE(alignment_file,
                                     outdir='./',
                                     outfile_basename='aligned_positions_selection_tixe',
                                     ref='Clar'
                                     ):
  positions_pdbnum = [316, 317, 318, 319]
  lookup_pos = [t-1 for t in positions_pdbnum]
  dict_pos =  get_positions(alignment_file,
                        lookup_pos=lookup_pos,
                        ref=ref, include_aa=True)
  write_readable_lists_pdbnumbering(dict_pos,f'{outdir}/{outfile_basename}_pdbnum.csv')
  dict_pos =  get_positions(alignment_file,
                        lookup_pos=lookup_pos,
                        ref=ref)
  write_readable_lists(dict_pos, f'{outdir}/{outfile_basename}.txt')  


if __name__=='__main__':
 
 outdir='aligned_positions'
 os.makedirs(outdir, exist_ok=True)
 get_positions_selection(alignment_file='Sequences_all_for_alignment.aln',
                         outdir=outdir,
                         ref='Campylobacter_lari')