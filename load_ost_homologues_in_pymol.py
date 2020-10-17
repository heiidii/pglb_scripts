import os
import sys
import glob
#print(sys.executable)
#print(sys.path)
import pickle
import pandas as pd
from pymol import cmd

colors={'DmPglB':'lightblue','DgPglB':'argon',
        'DvPglB':'lightorange','DdPglB':'paleyellow',
        'ClPglB':'lightpink','CjPglB':'hotpink'}
dict_names={'DmPglB':'marinus','DvPglB':'vulgaris','DgPglB':'gigas','DdPglB':'alaskensis',
              'CjPglB':'jejuni','ClPglB':'lari'}
color_ref='paleyellow'
color_surface = 'white'
color_cartoon = 'scandium'
sele_Thr_glyc = "/%s//P/THR`3"
surface_transparency=0.3

def show_selected_residues(curprot,select_residues,rep='sticks'):
  res_selname = '%s_%s' %(curprot,key)
  select='%s and ( resi '%curprot + ','.join(select_residues[key])+')'          print(res_selname,select)
            cmd.select(res_selname,select)
            #doesnt work - turns on all sidechains
            #select_sidechain = 'sc. ( %s )' %res_selname
            #print(select_sidechain)
            cmd.show(rep,res_selname)
            cmd.hide('sticks','hydrogen')
            cmd.color(color_current,res_selname)
            cmd.util.cnc(res_selname) 

def open_from_df(df,top_n=1,pdb_column='description',prefix='./',ref_argname=None,curprot=None,select_residues={},color_current='white'):
  print(df.columns)
  for i in range(0,top_n): #topN pdbs
    pdb=df[pdb_column][i]
    print(pdb)
    pdbfile = glob.glob(os.path.join(prefix,pdb)+'*')[0]
    if os.path.exists(pdbfile):
        cmd.load(pdbfile)
        argname_base = pdbfile.split('/')[-1].split(".pdb")[0]
        if curprot is None:
            curprot = 'prot_%d' %i
        cmd.select(argname_base + ' and polymer')
        cmd.create(curprot,"sele")
        cmd.show_as('cartoon',curprot)
        cmd.set('transparency',surface_transparency,curprot)
        cmd.disable(argname_base)
        if not ref_argname is None:
            cmd.align(curprot,ref_argname)
        cmd.color(color_surface,curprot)
        cmd.set('sphere_scale',0.36)
        cmd.color('scandium','resname MN')
        if select_residues:
          for key in select_residues:
            res_selname = '%s_%s' %(curprot,key)
            select='%s and ( resi '%curprot + ','.join(select_residues[key])+')'
            print(res_selname,select)
            cmd.select(res_selname,select)
            #doesnt work - turns on all sidechains
            #select_sidechain = 'sc. ( %s )' %res_selname
            #print(select_sidechain)
            cmd.show('sticks',res_selname)
            cmd.hide('sticks','hydrogen')
            cmd.color(color_current,res_selname)
            cmd.util.cnc(res_selname)
    else:
        print('Cannot find file: top %d: %s' %(top_n,pdbfile))


def generate_session_homologues(preprefix,dirname,enzymes,dict_positions,session_name,by='total_score',suffix='_enzyme'):
  ref_argname=None #enzyme at position 0 in list enzymes will be used as reference for alignment
  for ie,name in enumerate(enzymes):
    print('Processing %s' %name)
    prefix='%s/%s/%s' %(preprefix,name,dirpath)
    scorefile='%s/output_enzyme/score.sc' %(prefix)
    pickled_path = os.path.join(prefix,'pickled_files')
    csvfile = os.path.join(pickled_path,'df_sorted-by_%s%s.csv' %(by,suffix))
    df = pd.read_csv(csvfile) #pd.read_pickle(serialized_db) 
    select_residues={}
    if dict_names[name] in dict_positions:
     select_residues['m2']=dict_positions[dict_names[name] ]
    color_current=colors[name]
    open_from_df(df,prefix=prefix,curprot=name,ref_argname=ref_argname,select_residues=select_residues,color_current=color_current)
    if ie==0:
     ref_argname=name
  cmd.save(session_name)

    
if __name__=='__main__':
  dirpath='01_prepare_structure'
  preprefix='/home/saipooja/Glycosyltransferases'
  enzymes=['DmPglB','DgPglB','DdPglB','CjPglB','ClPglB']#,'DvPglB']
  dict_positions_m1 = pickle.load(open('OST_sequences/residue_positions_m1.p','rb'))
  dict_positions_m2 = pickle.load(open('OST_sequences/residue_positions_m2.p','rb'))
  #generate_session_homologues(preprefix,dirpath,enzymes,dict_positions_m1,'session_m1.pse')
  generate_session_homologues(preprefix,dirpath,enzymes,dict_positions_m2,'session_m2.pse')
