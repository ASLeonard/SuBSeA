import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import json
import requests
from collections import defaultdict, Counter
import numpy as np
from itertools import product
from scipy.stats import linregress, ks_2samp, anderson_ksamp, mannwhitneyu, epps_singleton_2samp

from periodic import loadCSV,overlaps

from domains import readDomains, domainMatch

def makeCSV(df,domain):
    domain_dict=readDomains(domain)

    new_rows = []

    for _, row in df.iterrows():
        domain_info = ';'.join([chain+':{}'.format(tuple(domain_dict[row['id'][:4]][chain])) if chain in domain_dict[row['id'][:4]] else '' for chain in row['chains']] )
        
        new_rows.append({'PDB_id':row['id'][:4], 'interfaces':'-'.join(sorted(row['chains'])), 'domains':domain_info, 'BSAs':round(row['BSA']) if not pd.isna(row['BSA']) else ''})

    return pd.DataFrame(new_rows)

def writeCSV(df,fname):
    df.to_csv(fname,index=False,columns=['PDB_id','interfaces','domains','BSAs'])

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
                #print("Empty ID on {}".format(d['id']))
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

def convertH(df):
    rows = []
    for _,row in df.iterrows():
        for interface in row['interfaces']:
            rows.append({'id':'a','shared':domainMatch(row['domains'],*interface.split('-')),'BSA':row['BSAs'][interface]})

    return pd.DataFrame(rows)

def plotDataX(data,stat_func=mannwhitneyu):
    #g=sns.regplot(x="TM", y="BSA", data=data)                  
    #g = sns.violinplot(y="BSA",data=low)#xlim=(0,1))
    #plt.figure()
    #g3 = sns.violinplot(y="BSA",data=high)#xlim=(0,1))
    #g2 = sns.jointplot(x="TM", y="BSA",data=data,xlim=(0,1))
    #sns.relplot(x='TM', y='BSA', size="BSA",sizes=(40, 400),hue='domain', alpha=.75, height=6, data=data)
    plt.figure()
    ax = sns.violinplot(x="shared", y="BSA", data=data, palette="muted",split=False,  scale="count",scale_hue=True,inner="quartile",cut=0,bw=.3)

    
    plt.show(block=False)
    low = data.loc[(data['shared']=='none')]# | (data['shared']=='NA')]
    high = data.loc[(data['shared']=='full') | (data['shared']=='partial')]
    homodimer_ref = data.loc[(data['shared']=='hom')]
    for val, dic in zip(('High','Low','Ref'),(high,low,homodimer_ref)):
        print("{}: ({})".format(val,len(dic))," = ", np.nanmedian(dic['BSA']))
     
    print("\np-value: {}\n".format(stat_func(high['BSA'],low['BSA'],alternative='greater',distribution='normal')[1]))
    print("\np-value: {}\n".format(stat_func(homodimer_ref['BSA'],low['BSA'],distribution='normal')[1]))
    print("\np-value: {}\n".format(stat_func(homodimer_ref['BSA'],high['BSA'],distribution='normal')[1]))
    
    #for overlap in ('full','partial','none'): #'NA'
    #     print(overlap,"({})".format(len(data.loc[data['shared']==overlap])),np.nanmedian([float(i) for i in data.loc[data['shared']==overlap]['BSA']]))
    
    #domain_yes=data.loc[(data['shared']=='full') | (data['shared']=='partial')]
    #domain_no=data.loc[data['shared']=='none']
    #data.loc[~(data['domain']=='full') & ~(data['domain']=='partial')]
    #print("\np-value: {}".format(mannwhitneyu(domain_yes['BSA'],domain_no['BSA'],alternative='greater')[1]))
    
    #print('\n CLES: ',commonLanguageES(low['BSA'],high['BSA']))
    #print('\n CLES: ',commonLanguageES(high['BSA'],homodimer_ref['BSA']))
    #print('\n CLES: ',commonLanguageES(low['BSA'],homodimer_ref['BSA']))

    #print(np.nanmedian(domain_yes['BSA']),np.nanmedian(domain_no['BSA']))


def plotHeteromerBSA(df):

    domain_yes=df.loc[(df['shared']=='full') | (df['shared']=='partial')]
    domain_no=df.loc[(df['shared']=='none') | (df['shared']=='NA')]

    f, ax = plt.subplots()
    
    for overlap in ('full','partial','none','NA','hom'):
        sns.kdeplot(ax=ax,data=df.loc[df['shared']==overlap]['BSA'],label=overlap)

    f, ax = plt.subplots()
    
    for dd in (domain_yes,domain_no):
        sns.kdeplot(ax=ax,data=dd['BSA'])

    plt.show(0)
        
    

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

def loadDict(file_ID):
    return json.load(open(f'{file_ID}_comparison.dict','r'))

def loadND(file_ID):
    return np.fromfile(open(f'{file_ID}_comparison.ND','r')).reshape(-1,2)

def loadDF(file_ID,json_save = False):
    raw_data = (loadDict if json_save else loadND)(file_ID) 
    rows = []
    for results in (raw_data.values() if json_save else raw_data):
        pval, similarity = results
        rows.append({'pval':pval or 1,'similarity':similarity,'sg':int(similarity//10)})
        
    df = pd.DataFrame(rows)
    df['pval'] = np.log10(df['pval'])
    return df

def getFrac(data,key):
    vals = list(data.values())
    print('{:.3f}'.format(vals.count(key)/len(vals))) 

def plotData(datas,ax=None,stat_func=ks_2samp,merge_nones=True):
    if isinstance(datas,list) and isinstance(datas[0],dict):
        labels = ['1']*len(datas)
    else:
        labels = datas
        datas = [loadDict(d) for d in labels]
    if merge_nones:
        cleaned_datas = [np.log10([val[0] or 1 for val in data.values() if (val!='error' and val[1]<300)]) for data in datas]
    else:
        cleaned_datas = [np.log10(list(filter(lambda x: isinstance(x[0],float),data.values()))) for data in datas]
    #

    main_range=(-25,0)   
    
    if not ax:
        f,ax = plt.subplots()

    def logSlope(counts,bins):
        slope_ROI=slice(round(.5*len(counts)),round(.85*(len(counts))))
        counts_ROI=np.log10(counts[slope_ROI])
        bins_ROI=(bins[slope_ROI]+bins[slope_ROI.start+1:slope_ROI.stop+1])/2
        y,b=linregress(bins_ROI,counts_ROI)[:2]
        print('Slope: ',y)
        ax.plot(bins[slope_ROI],10**(bins[slope_ROI]*y+b),ls='--',lw=3,c=col2)

    for clean_data,label,(col1,col2) in zip(cleaned_datas,labels,(('royalblue','skyblue'),('orangered','coral'),('g','g'),('m','m'),('k','k'))):
        print(len(clean_data), sum(clean_data<np.log10(.05))/len(clean_data))
        counts, bins, _ = ax.hist(clean_data,range=main_range,bins=100,density=True,histtype='step',color=col1,label=label)
        logSlope(counts,bins)
        
        outlier_low = sum(clean_data<main_range[0])/len(clean_data)
        ax.bar(bins[0],outlier_low,width=-5*(bins[1]-bins[0]),align='edge',edgecolor='w',hatch='////',facecolor=col1,alpha=0.3)

       
          
    ax.set_yscale('log')
    ax.axvline(np.log10(.05),c='darkgrey',ls='--',lw=5)

    plt.legend()
    plt.show(block=False)
    #print(stat_func(*cleaned_datas))
    return ax
     
def plotSNS(df):
    sns.jointplot(data=df,x='pval',y='similarity',kind='hex',joint_kws={'gridsize':(1000,100),'bins':'log'})
    plt.show(block=False)

def plotSNS2(df):
    grid = sns.JointGrid(x='pval', y='similarity', data=df)

    g = grid.plot_joint(plt.hexbin,gridsize=(300,50),bins='log',cmap='RdGy',mincnt=1,extent=(-25,0,0,100))
    sns.kdeplot(df['pval'], ax=g.ax_marg_x, legend=False,clip=(-10,0))
    g.ax_marg_x.set_yscale('log')
    g.ax_marg_x.set_ylim(1e-6,1)
    sns.kdeplot(df['similarity'], ax=g.ax_marg_y, vertical=True, legend=False)
    plt.show(block=False)

def plotCats(df,ax=None,ls='-'):
    if ax is None:
        f, ax =plt.subplots()
    cmap = plt.get_cmap('tab20')
    for i in range(0,3):
        if len(df.loc[df['sg']==i]['pval']) < 10:
            continue
        print(len(df.loc[df['sg']==i]['pval']))
        sns.distplot(a=df.loc[df['sg']==i]['pval'],bins=np.linspace(-50,0,51),ax=ax,norm_hist=True,color=cmap(i/11),label=i,kde=False,kde_kws={'cut':0,'kernel':'epa','bw':.1,'ls':ls},hist_kws={'histtype':'step','alpha':1,'lw':2,'ls':ls})#kde_kws={'ls':ls,'alpha':1})
    plt.yscale('log',nonposy='mask')
    plt.legend()
    plt.show(block=False)
    return ax


