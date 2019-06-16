import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import json
import requests
from collections import defaultdict, Counter
import numpy as np
import scipy.stats as ss
from itertools import product

   
     
def scrubInput(run_name,local_alignment=True,similarity=True,domain='SCOP'):

     rows_list = []

     domain_dict=readDomains(domain)
     #homomer_domains=getUniqueHomodimerDomains(100)

          
     for result in glob.glob('/scratch/asl47/PDB/results/{}/*.results'.format(run_name)):
          d={'id':result.split('.')[0].split('/')[-1].lower()}
          raw=[line.rstrip() for line in open(result)]
          
          try:
               d['homomer']=bool(raw[0]=='1')
                    
               if len(raw)==1:
                    print("Empty ID on {}".format(d['id']))
                    continue
               d['chains'] = tuple(sorted(raw[1].split()))
               
               if any('#' in line for line in raw):
                    homologies=tuple(filter(lambda element: '#' in element,raw))
                    d['homology']=float(homologies[2*local_alignment+similarity].split()[-1][:-1]) 
                    d['H_pred']='yes' if d['homology']>=30 else 'no'

               if any('Total' in line for line in raw):
                    bsa=tuple(filter(lambda element: 'Total' in element,raw))
                    d['BSA']=((float(bsa[2].split()[-1])+float(bsa[1].split()[-1]))-float(bsa[0].split()[-1]))/2
                    
               if d['id'][:4] in domain_dict:
                    d['domain']=domainMatch(domain_dict[d['id'][:4]],*d['chains'])
                    if d['domain']=='full':
                         d['arch']=domain_dict[d['id'][:4]][d['chains'][0]]
                    else:
                         d['arch']=None
                         #or d['domain']=='partial':

                         #d['d_group']=any([tuple(domain_dict[d['id'][:4]][c]) in homomer_domains for c in d['chains']])
                    #elif 'Homodimer' in path:
                    #     d['d_group']=None
                    #else:
                    #     d['d_group']='none'
                    
                         
                         
               else:
                    d['domain']='unknown'

               if any('TM-score' in line for line in raw):
                    aligns=tuple(filter(lambda element: 'TM-score' in element,raw))
                    d['TM']=float(aligns[0].split()[1])
                    d['S_pred']= 'no' if d['TM']<=.2 else 'yes' if d['TM']>=.5 else 'maybe'
                    

          except Exception as e:
               print("Exception",e)
               print("except on ", d['id'], "len ",len(raw))
               continue

          rows_list.append(d)
     return pd.DataFrame(rows_list)
          

def plotData(data,stat_func=''):
     #g=sns.regplot(x="TM", y="BSA", data=data)
     
                   
     #g = sns.violinplot(y="BSA",data=low)#xlim=(0,1))
     #plt.figure()
     #g3 = sns.violinplot(y="BSA",data=high)#xlim=(0,1))
     #g2 = sns.jointplot(x="TM", y="BSA",data=data,xlim=(0,1))
     sns.relplot(x='TM', y='BSA', size="BSA",sizes=(40, 400),hue='domain', alpha=.75, height=6, data=data)
     plt.figure()
     ax = sns.violinplot(x="domain", y="BSA",hue='d_group', data=data, palette="muted",split=False,  scale="width",scale_hue=True,inner="quartile",bw=.1)
     
     plt.show(block=False)

     low=data.loc[data['S_pred']=='no']
     high=data.loc[data['S_pred']=='yes']

     for val, dic in zip(('High','Low'),(high,low)):
          print("{}: ({})".format(val,len(dic))," = ", np.nanmedian(dic['BSA']))
     
     
     print("\np-value: {}\n".format(ss.mannwhitneyu(high['BSA'],low['BSA'],alternative='greater')[1]))
     for overlap in ('full','partial','none','NA'):
          print(overlap,"({})".format(len(data.loc[data['domain']==overlap])),np.nanmedian([float(i) for i in data.loc[data['domain']==overlap]['BSA']]))

     domain_yes=data.loc[(data['domain']=='full') | (data['domain']=='partial')]
     domain_no=data.loc[data['domain']=='none']
     #data.loc[~(data['domain']=='full') & ~(data['domain']=='partial')]
     print("\np-value: {}".format(ss.mannwhitneyu(domain_yes['BSA'],domain_no['BSA'],alternative='greater')[1]))
     print('\n CLES: ',commonLanguageES(domain_no['BSA'],domain_yes['BSA']))
     print(np.nanmedian(domain_yes['BSA']),np.nanmedian(domain_no['BSA']))

##find fraction of pairs where the higher element is actually higher than the lower element
def commonLanguageES(lows, highs):
     assert len(lows) and len(highs), 'invalid arguments'

     return sum(h>l for l,h in product(lows,highs))/(len(lows)*len(highs))


def loadLevy():
     data=[]
     for line in open('../../../../Downloads/levydata'):
          vs=line.split()
          data.append((vs[0],float(vs[3]),*[int(i) for i in vs[6:9]]))
     return data

def plotLevy(df,plo='BSA'):
     plt.figure()
     if plo == 'BSA':
          plt.hist(df['BSA'],bins=100)   
          plt.show(block=False)

     elif plo == 'Hom':
          ax = sns.violinplot(x="homol", y="BSA", data=df.loc[df['iden']==0], palette="Set2",  scale="count",inner="quartile",bw=.1)
          plt.show(block=False)
     data=df.loc[df['iden']==0]
     #return data
     print(np.median(data.loc[data['homol']==1]['BSA']),np.median(data.loc[data['homol']==0]['BSA']))
     print(stats.mannwhitneyu(data.loc[data['homol']==1]['BSA'],data.loc[data['homol']==0]['BSA'],alternative='greater'))

       
                    

     
def correspondingHomodimers(heteromerics, homomerics):

     domain_groupings={}
     for idx,data in enumerate((heteromerics,homomerics)):
          for _,row in data.iterrows():
               if row['arch'] is None:
                    continue
               domain=' '.join(row['arch'])
               if domain not in domain_groupings:
                    domain_groupings[domain]=([],[])
               else:
                    domain_groupings[domain][idx].append('{}_{}_{}'.format(row['id'][:4],*row['chains']))

     non_trivial={domains: pdb_pairs for domains,pdb_pairs in domain_groupings.items() if all(pdb_pairs)}
     with open('domain_groups.json', 'w') as file_:
          file_.write(json.dumps(non_trivial))
                
     return non_trivial


def loadDict(run):
     return json.load(open('full_run_{}.dict'.format(run),'r'))

def plotData(data):
     clean_data = np.log10(list(filter(lambda x: isinstance(x,float),data.values())))
     f,ax = plt.subplots()
     
     ax.hist(clean_data,range=(-20,0),bins=300)
     ax.set_yscale('log')
     ax.axvline(np.log10(.05),c='r',ls='--',lw=2)
     
     plt.show(block=False)
     
