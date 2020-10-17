import os,sys
import pandas as pd
import pickle
import glob
import re
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
#from HelperPlotting import *

listres = ['L','V','A','I','M','P','G','F','W','Y','S','T','N','Q','H','K','R','D','E']
cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
#os.system('mkdir -p PickledFiles')
top_pdbs_dir = 'top_pdbs'
os.system('mkdir -p %s' %top_pdbs_dir)
figures_dir = 'figures'
pdb_path = 'pdbs_initial_relax_constrained'
columns = ['total_score',
          'description']


def fixname(entry,pdb_path_local=None):
  basepdbname = '/'.join(entry.split('/')[-1:])
  if pdb_path_local is None:
    return '/'.join([pdb_path,basepdbname])    
  else:
    return '/'.join([pdb_path_local,basepdbname])



def scorefile_to_dataframe(scorefile):
  df = pd.read_fwf(scorefile,skiprows=1)
  print(df)
  df_clean =  pd.DataFrame()
  
  for column in columns:
    if column in df:
      df_clean[column]= df[column]
  return df_clean

def generateplots(df,outfile,plottype='box',hue=None,show=False):
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize=(18,12))
  for i,column in enumerate(columns[:12]):
    if column=='description': continue
    ax1 = fig.add_subplot(3,4,i+1)
    if plottype=='box':
      sns.boxplot( column , data=df , ax=ax1 )
    if plottype=='swarm':
      im = sns.swarmplot( x=df[column]  , ax=ax1 , y=[""]*len(df) , palette= 'Reds_r')
      #ax1.get_legend().remove()
    if plottype=='violin':
      sns.violinplot( column , data=df , ax=ax1 )
    ax1.tick_params(labelsize=20, grid_linewidth=1.0, grid_linestyle='-', grid_color='black',grid_alpha=0.6)

  plt.tight_layout()
  plt.savefig(outfile,transparent='True',dpi=fig.dpi)
  if show:
    plt.show()
  plt.close()


def update_dataframe_pdb_and_fastas(df,top_pdbfiles,pdb_path_local=None):
  cleanpdbcolumn=[]

  for i,entry in enumerate(list(df['description'])):
    #def fixname(entry):
    #  return '/'.join(entry.split('/')[-2:])
    entry = fixname(entry,pdb_path_local)
    foundentry = glob.glob(entry+'*')
    
    if len(foundentry)==1:
      cleanpdbcolumn.append(foundentry[0])
    else:
      print(entry + ' not found\n')
    
  df['description_clean']=cleanpdbcolumn
  return df  

def write_files_cumulative(df,outfile_handle,top_N=10,prefix='',field='description_clean'):
  for i,pdb in enumerate(list(df[field])[:top_N]):
    outfile_handle.write(prefix+pdb+'\n')

def write_pdbfiles(df,outfile,top_N=10):
  outf=open(outfile,'w')
  write_files_cumulative(df,outf,top_N=top_N,field='description_clean')
  outf.close()

def write_pymol_session(df,outfile,top_N=10,prefix=''): 
  outf=open(outfile,'w')
  outf.write('pymol ')
  for i,pdb in enumerate(list(df['description_clean'])[:top_N]):
    outf.write(os.path.join(prefix,pdb)+'\t')
  outf.write('& \n')
  outf.close()
  os.chmod(outfile, 0o775)

def serialize_data(df,infile,outname,dict_metadata={}):
  data_={}
  data_['dataframe'] = df
  data_['source']=infile
  data_['metadata']=''
  data_['script']=sys.argv[0]
  data_['path']=os.getcwd()
  if dict_metadata:
    for key in dict_metadata:
      data_[key]=dict_metadata[key]
  #pickle.dump(data_,open(outname,'wb'))
  #df.to_pickle(outname)
  df.to_csv(outname)
  return outname


def scorefile_to_df_file(infile,outfile='dataframe_scorefile.p',N=50):
  df = scorefile_to_dataframe(infile)
  df_totalscore_sorted = df.sort_values(by='total_score')
  df_totalscore_sorted.reset_index(drop=True,inplace=True)
  serialize_data(df_totalscore_sorted,infile,outname=outfile)
  

def process_scorefile(infile,outfile_p='dataframe_scorefile.p',pdb_path_local=None,suffix='',outfile_handle=None,prefix='',by='total_score',top_N=500,secondary_filter=''):
  df = scorefile_to_dataframe(infile)

  graphtypes=['violin','swarm']
  path = os.path.join(prefix,figures_dir)
  if not os.path.exists(path):
    os.mkdir(path)
  for graphtype in graphtypes:
    outfile = os.path.join(path, 'alldata_%s%s.png' %(graphtype,suffix))
    generateplots(df,outfile,graphtype)

  df_totalscore_sorted = df.sort_values(by=by)
  df_totalscore_sorted.reset_index(drop=True,inplace=True)

  df_totalscore_sorted_topN = df_totalscore_sorted[:top_N]
  df_totalscore_sorted_topN.reset_index(drop=True,inplace=True)

  top_pdbfiles = df_totalscore_sorted_topN['description']
  
  df_temp = update_dataframe_pdb_and_fastas(df_totalscore_sorted_topN,top_pdbfiles,pdb_path_local)

  del df_totalscore_sorted_topN
  df_totalscore_sorted_topN = df_temp
  del df_temp

  if secondary_filter != '':
     if secondary_filter == 'distance_catalysis_HWUO1B-THR7N':
       df_filter = df_totalscore_sorted_topN[ df_totalscore_sorted_topN[secondary_filter] < 5.5 ]
       del df_totalscore_sorted_topN
       df_totalscore_sorted_topN = df_filter

  dfs = [df_totalscore_sorted_topN]#, df_complex_normalized_sorted_topN]
  names = [ 'sorted_totalscore']#,'sorted_complexnormalized']
  for curdf, name in zip(dfs,names):

    #outfile = os.path.join(path, '/'+name+'.png')
    #generateplots(curdf,outfile,plottype='swarm',hue='total_score')
    path_pdb = os.path.join(prefix,top_pdbs_dir)
    if not os.path.exists(path_pdb):
      os.mkdir(path_pdb)
    for j in [5,10,100]:
      outfile = os.path.join(path_pdb,'toppdbs_top%03d_%s%s.txt' %(j,name,suffix))
      write_pdbfiles(df_totalscore_sorted_topN,outfile,top_N=j)
      outfile = os.path.join(path_pdb,'run_pymol_top%03d_%s%s' %(j,name,suffix))
      write_pymol_session(df_totalscore_sorted_topN,outfile,top_N=j,prefix='')
  
  if not outfile_handle is None:    
    write_files_cumulative(df_totalscore_sorted_topN,outfile_handle,top_N=1,prefix=prefix)
  #pickled_path = os.path.join(prefix,'pickled_files')
  #if not os.path.exists(pickled_path):
  #  os.mkdir(pickled_path)
  #outname = os.path.join(pickled_path,'df_sorted-by_%s_%s.p' %(by,suffix))
  #serialized_db = serialize_data(df_totalscore_sorted,infile,outname=outfile_p)
  return df_totalscore_sorted

def process_scorefiles_with_pattern(pattern,prefix='./',by='total_score',secondary_filter='',suffix='',top_pdb_paths=1):
  files = glob.glob(pattern)
  files.sort()
  path=os.path.join(prefix,top_pdbs_dir)
  if not os.path.exists(path):
    os.mkdir(path)
  outfile_handle = open(os.path.join(path,'top1_pdbs%s.txt' %suffix),'w')
  for infile in files:
    pdb_path_local = infile.split('score.sc')[0]
    if secondary_filter != '':
       suffix += '_%s_%s' %(by,secondary_filter)
    else:
       suffix += '_%s' %by #infile.split('/')[0].split('FastRelax')[1]
    outfile='dataframe%s.p' %(suffix)
    df_totalscore_sorted = process_scorefile(infile,outfile_p=outfile,pdb_path_local= pdb_path_local,outfile_handle=outfile_handle,suffix=suffix,prefix=prefix,by=by,secondary_filter=secondary_filter)
  outfile_handle.close()
  return df_totalscore_sorted

def generate_scatterplots_2d(pattern,top_N=200,filter_=True,filter_by='',xfield='x',yfield='y',zfield='total_score',xlim=[],ylim=[]):
  files = glob.glob(pattern)
  files.sort()
  for infile in files:
    pdb_path_local = infile.split('score.sc')[0]
    df = scorefile_to_dataframe(infile)
    if filter_:
      df_totalscore_sorted = df.sort_values(by=filter_by)
      df_totalscore_sorted.reset_index(drop=True,inplace=True)
      df_totalscore_sorted_topN = df_totalscore_sorted[:top_N]
      df_totalscore_sorted_topN.reset_index(drop=True,inplace=True)
      del df
      df = df_totalscore_sorted_topN

    print(xfield,yfield,zfield)
    print(len(list(df[xfield])),len(list(df[yfield])),len(list(df[zfield])))
    suffix+=xfield+'_vs_'+yfield+'_z'+zfield
    fig,ax = plt.subplots(figsize=(4,3))
    s=2
    ax = sns.scatterplot(x=df[xfield],y=df[yfield],data=df,hue=df[zfield],s=s,edgecolor=None,linewidths=0.1,palette="BuPu_r")
    if len(xlim) >0:
        ax.set(xlim=xlim)
    if len(ylim) >0:
        ax.set(ylim=ylim)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.tight_layout()
    os.system('mkdir -p results/scatter')
    if filter_:
      outfile = 'results/scatter/scatterplot_%s_topN%d.png' %(suffix,top_N)
    else:
      outfile = 'results/scatter/scatterplot_%s.png' %suffix
    plt.savefig(outfile,transparent=True,dpi=600)
    #plt.show() 
    plt.close()

def generate_scatterplots(pattern,top_N=200,filter_=True,xfield='bb_heavy_rmsd',yfield='total_score',xlim=[],ylim=[]):
  files = glob.glob(pattern)
  files.sort()
  for infile in files:
    pdb_path_local = infile.split('score.sc')[0]
    df = scorefile_to_dataframe(infile)
    if filter_:
      df_totalscore_sorted = df.sort_values(by=yfield)
      df_totalscore_sorted.reset_index(drop=True,inplace=True)
      df_totalscore_sorted_topN = df_totalscore_sorted[:top_N]
      df_totalscore_sorted_topN.reset_index(drop=True,inplace=True)
      del df
      df = df_totalscore_sorted_topN
    
    suffix=xfield+'_vs_'+yfield
    fig,ax = plt.subplots(figsize=(4,3))
    s=2
    x = sns.scatterplot(x=df[xfield],y=df[yfield],data=df,s=s,color='orchid',edgecolor=None,linewidths=0.1)
    if len(xlim) >0:
        ax.set(xlim=xlim)
    if len(ylim) >0:
        ax.set(ylim=ylim)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.tight_layout()
    os.system('mkdir -p results/scatter')  
    if filter_:
      outfile = 'results/scatter/scatterplot_%s_topN%d.png' %(suffix,top_N)
    else:
      outfile = 'results/scatter/scatterplot_%s.png' %suffix
    plt.savefig(outfile,transparent=True,dpi=600)
    #plt.show() 
    plt.close()

if __name__ == '__main__':
  dirpath='01_prepare_structure'
  enzymes=['DvPglB','DgPglB','DdPglB','CjPglB','DmPglB','ClPglB']
  for name in enzymes:
    print('Processing %s' %name)
    prefix='%s/%s' %(name,dirpath)
    scorefile='%s/output_enzyme/score.sc' %(prefix)
    suffix='_enzyme'
    by='total_score'
    pickled_path = os.path.join(prefix,'pickled_files')
    if not os.path.exists(pickled_path):
        os.mkdir(pickled_path)
    outname = os.path.join(pickled_path,'df_sorted-by_%s%s.csv' %(by,suffix))
    if not os.path.exists(outname):
      df = process_scorefiles_with_pattern(scorefile,by=by,suffix=suffix,prefix=prefix)
      serialized_db = serialize_data(df,scorefile,outname=outname, dict_metadata=dict(name=name))
    print(outname)
      
