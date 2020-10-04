import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import json
from collections import defaultdict
import numpy as np
from itertools import product

from scipy.stats import linregress,ks_2samp,brunnermunzel,spearmanr
from SuBSeA.utility import loadCSV
from SuBSeA.domains import readDomains, duplicateIntersection

from matplotlib.colors import Normalize
 


def groupHeteromericBSAs(threshold=100,filter_immunoglobulins=False,domain_type='CATH'):

    def domainMatchFractional(domains,c1,c2):
        return len(duplicateIntersection(domains[c1],domains[c2])) / max(len(domains[c1]),len(domains[c2]))

    df = loadCSV(f'Heteromeric_complexes_{domain_type}_{threshold}.csv')

    rows = []
    for _,row in df.iterrows():
        for interface in row['interfaces']:
            rows.append({'id':f'{row["PDB_id"]}_{interface.replace("-","_")}', 'shared':domainMatchFractional(row['domains'],*interface.split('-')),
    'BSA':row['BSAs'][interface]})
            if filter_immunoglobulins and rows[-1]['shared']>=0 and ('2.60.40.10') in row['domains'][interface[0]]:
                del rows[-1]

    bsa_df = pd.DataFrame(rows)
    
    bsa_df ['ancestry'] = ['homologous' if x>0 else 'analogous' for x in bsa_df['shared']]
    
    bsa_df['filter'] = threshold
    return bsa_df

## plots violinplot for heteromeric BSA data give by "groupHeteromericBSAs" function.
## stat_func can take on most scipy.stats 2 dataset tests
def plotGroupedHeteromericBSAs(data,stat_func=brunnermunzel):
    ## make the figure
    plt.figure()
    data = data[(data.ancestry=='homologous') | (data.ancestry=='analogous')]
    sns.violinplot(x="ancestry", y="BSA", hue='ancestry',data=data, palette="muted",scale_hue=False,inner="quartile",cut=0,bw=.3,split='True')
    plt.yscale('log')
    plt.show(block=False)

    ## crunch the stats
    analogous = data.loc[data['ancestry']=='analogous']
    homologous = data.loc[data['ancestry']=='homologous']
    
    for val, data in zip(('Homologous','Analogous'),(homologous,analogous)):
        print(f'{val} median = {np.nanmedian(data["BSA"]):.0f} ({len(data)} points)')

    stat_args = defaultdict(dict,{brunnermunzel:{"alternative":"greater","distribution":"normal"}})
    
    print(f'\np-value: {stat_func(homologous["BSA"],analogous["BSA"],**stat_args[stat_func])[1]:.3e}')
    print(f'CLES: {commonLanguageES(analogous["BSA"],homologous["BSA"]):.3f}')
      

##find fraction of pairs where the higher element is actually higher than the lower element
def commonLanguageES(lows, highs):
    assert len(lows) and len(highs), 'invalid arguments'
    return sum(h>l for l,h in product(lows,highs))/(len(lows)*len(highs))


def makeCSV(df,domain):
    domain_dict=readDomains(domain)
    new_rows = []

    for _, row in df.iterrows():
        domain_info = ';'.join([chain+':{}'.format(tuple(domain_dict[row['id'][:4]][chain])) if chain in domain_dict[row['id'][:4]] else ''
            for chain in row['chains']] )
        
        new_rows.append({'PDB_id':row['id'][:4], 'interfaces':'-'.join(sorted(row['chains'])),
            'domains':domain_info, 'BSAs':round(row['BSA']) if not pd.isna(row['BSA']) else ''})

    return pd.DataFrame(new_rows)
     
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
    return np.fromfile(open(f'{file_ID}_comparison.ND','r')).reshape(-1,3)

def loadDF(file_ID,json_save = False,csv_save=False):
    if csv_save:
        df = pd.read_csv(f'{file_ID}_comparison.csv')
    else:
        raw_data = (loadDict if json_save else loadND)(file_ID)
        rows = []
        for results in (raw_data.values() if json_save else raw_data):
            if results == 'error':
                continue
            pval, N_hits, similarity = results
            rows.append({'pval':pval or 1,'similarity':similarity,'sg':int(similarity//5),'hits':N_hits})
        
        df = pd.DataFrame(rows)
        
    for pval in ('pval_F','pval_S','pval_T','pval_F2','pval_S2','pval_T2'):
            df[pval] = -1*np.log10(df[pval])
            print(pval)

    return df

def splitData(df,thresh=80):
    f80 = df.loc[(df['code']=='MUT') & (df['similarity']<thresh)]
    M80 = df.loc[(df['code']=='MPA') & (df['similarity']<thresh)]
    R80 = df.loc[(df['code']=='DNO') & (df['similarity']<thresh)]
    return f80,M80,R80

def loadALL(sample=None,rscratch=True):
    rscratch = '/rscratch/asl47/PDB_results/' if rscratch else ''

    DFs = [loadDF(f'{rscratch}{NAME}',csv_save=True) for NAME in ('Trial_70','Trial_90')]

    if sample:
        for i,df in enumerate(DFs):
            DFs[i] = df.sample(sample)
        

    for i, df  in enumerate(DFs):
        df['gaps'] = df['align_length']-df['overlap']
        df['sg'] =  df['similarity']//5
        df['norm_OVR'] = df['overlap']/df['align_length']
        df['norm_SCR'] = df['score']/df['align_length']
        
        for pval in ('pval_F','pval_S','pval_T','pval_F2','pval_S2','pval_T2'):
            df[pval] = -1*np.log10(df[pval])
        df['split'] = df['pval_S2']-df['pval_S']
        DFs[i]=df
    #print(df)

    
    return DFs#pd.concat(DFs,ignore_index=True)
    


def getFrac(data,key):
    vals = list(data.values())
    print('{:.3f}'.format(vals.count(key)/len(vals)))

def plotData(datas,ax=None,stat_func=ks_2samp,merge_nones=True):
    #if isinstance(datas,list) and not isinstance(datas[0],str):
    labels = ['1']*len(datas)
    cleaned_datas= datas
    main_range=(-25,0)
    
    if not ax:
        _,ax = plt.subplots()

    def logSlope(counts,bins):
        slope_ROI=slice(round(.5*len(counts)),round(.85*(len(counts))))
        counts_ROI=np.log10(counts[slope_ROI])
        bins_ROI=(bins[slope_ROI]+bins[slope_ROI.start+1:slope_ROI.stop+1])/2
        y,b=linregress(bins_ROI,counts_ROI)[:2]
        print('Slope: ',y)
        ax.plot(bins[slope_ROI],10**(bins[slope_ROI]*y+b),ls='--',lw=3,c=col2)

    for clean_data,label,(col1,col2) in zip(cleaned_datas,labels,(('royalblue','skyblue'),('orangered','coral'),('g','g'),('m','m'),('k','k'))):
        clean_data=np.log10(clean_data[:,0])
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

    
def plotSNS2(df,X_C='pval_F',Y_C='similarity',code=None):

    if code:
        df = df[df.code==code]
    df = df[df.similarity<95] #(df.norm_OVR>.25) & (df.pval_S2<-1)]
     
    extent_codes = {'pval_F':(30,0),'pval_S':(0,50),'pval_T':(-80,0),'similarity':(0,100),'norm_OVR':(0,1),'norm_SCR':(0,3)}
    grid = sns.JointGrid(x=X_C, y=Y_C, data=df)
    print(min(df[X_C]))

   

    g = grid.plot_joint(plt.hexbin,gridsize=(100,100),bins='log',cmap='cividis',mincnt=1,extent=extent_codes[X_C]+extent_codes[Y_C])

    g = g.plot_marginals(sns.distplot, kde=True, color="m",kde_kws={'cut':0,'kernel':'epa','bw':.05})
    
    g.ax_marg_x.set_yscale('log')
    g = g.annotate(spearmanr,loc="center left")
    plt.colorbar()

    plt.show(block=False)

def plotCats(df,ax=None,ls='-',cmap='green',pval_type=''):
    N=6
    if ax is None:
        _, ax =plt.subplots(3,2)
    ax= ax.flatten()
    #cmap = plt.get_cmap(cmap)
    for i in range(0,N):
        if len(df.loc[df['sg']==i][f'pval{pval_type}']) < 10:
            continue
        print(len(df.loc[df['sg']==i][f'pval{pval_type}']))

        color = cmap#cmap(i/(2*N))
        
        sns.distplot(a=np.clip(df.loc[df['sg']==i][f'pval{pval_type}'],-10,0),bins=np.linspace(-10,0,101),
            ax=ax[i],norm_hist=True,color=color,label=i,kde=False,
            kde_kws={'cut':0,'kernel':'epa','ls':ls},hist_kws={'histtype':'step','alpha':1,'lw':2})#kde_kws={'ls':ls,'alpha':1})
        CfD(df.loc[df['sg']==i][f'pval{pval_type}'],ax[i],color,'--')
        
        #fit_p = pareto.fit(-1*np.array(df.loc[df['sg']==i]['pval']))
        #plt.plot(np.arange(-25,0),pareto(*fit_p).pdf(np.arange(1,26)),'-.')
        ax[i].set_yscale('log',nonposy='mask')
        ax[i].text(.5,.8,f'{5*i} - {5*(i+1)} % similarity',va='center',ha='center',transform=ax[i].transAxes)
    #plt.legend()
    plt.show(block=False)
    return ax

def CfD(data,ax,c,ls):
    xs = np.linspace(0, 1, len(data), endpoint=False)
    data_sorted = np.sort(data)
    ax.plot(data_sorted,xs,c=c,lw=1,ls=ls,alpha=0.75)


def plotGrid(df):
    df = df.loc[df['sg']<6]
    df = df.loc[df['hits']>0]
    df = df.loc[df['hits']<=8]

    g = sns.FacetGrid(df, row="sg", col="hits",hue='match', margin_titles=True,sharex=True,sharey=True) 
    g.map(plt.hist, "pval_T", bins=np.linspace(-15,0,201),density=0,histtype='step',alpha=0.75).add_legend()
    g.set(yscale = 'log',ylim=[1,1e5])#1e-5,1])
    #return g
    plt.show(0)

def hexbin(x, y, color, **kwargs):
    cmap = plt.get_cmap('cividis')
    #cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, cmap=cmap, **kwargs)
    plt.text(-10,.1,len(x),ha='left',va='bottom')
    plt.colorbar()
    print(spearmanr(x,y))
    print(np.median(x),np.mean(x))

def plotNormOVR(df,sigma=-np.log10(.05),x_var='norm_OVR',stat_func=brunnermunzel):
    df['sig']=df.pval_S>=sigma
    df.loc[df.code=='MPA','code'] = 'MUT'

    print(f'Statistics for x-var: {x_var}')
    for domain_match in ('MUT','DNO'):
        positive_group = df.loc[(df.code==domain_match) & (df.sig)]
        negative_group = df.loc[(df.code==domain_match) & (not df.sig)]
        print(f'\tDomain overlap: {domain_match} = {stat_func(negative_group[x_var],positive_group[x_var])[1]}')
        print(f'\t  Positive mean: {np.mean(positive_group[x_var]):.2f} ({np.std(positive_group[x_var]):.2f})')
        print(f'\t  Negative mean: {np.mean(negative_group[x_var]):.2f} ({np.std(negative_group[x_var]):.2f})')
        #print(commonLanguageES(negative_group[x_var],positive_group[x_var]))

    sns.violinplot('code',x_var,data=df,hue='sig',cut=0,orient='v',inner='quartile',scale_hue=True,scale='area',split=True)
    plt.show(0)

def plotHeteromericConfidenceDecay(df,x_var = 'pval_S',sigma=-1,sim_thresh=95,MIN_SF=0.01):
    ## filter out "overly" similar protein subunits
    df = df.loc[df['similarity'] <= sim_thresh]

    ## filter out alignments below significance threshold based on variance
    var = np.var(df[x_var])
    print(f'Variance in x_var ({x_var}) is {var:.2f}')
    val_cut = sigma if isinstance(sigma,float) else sigma*var
    df = df[df[x_var]>val_cut]

    plt.figure()
    cdc= {'MUT':('H','darkorange'),'MPA':('P','forestgreen'),'DNO':('p','royalblue')}
    df.loc[df.code=='MPA','code'] = 'MUT'

    P_STARS = np.linspace(0,40,501)

    for code in ('MUT','MPA','DNO'):
        df_c = df[df.code==code][x_var]
        if len(df_c)==0:
            continue
        print(f'Heteromeric grouping: {code}, {len(df_c)} entries.')

        cdf = np.array([np.sum(df_c>=p_star) for p_star in P_STARS])

        np.seterr(divide='ignore')
        cdf = np.log10(cdf/np.nanmax(cdf))
        np.seterr(divide='warn')

        #cut off data to prevent outliers dominating
        cdf = cdf[cdf>=np.log10(MIN_SF)]
        #similarly trim data so fit is only in statistically strong region
        cdf_sl = cdf[cdf>=np.log10(MIN_SF)]

        plt.plot(P_STARS[:len(cdf)],cdf,c=cdc[code][1],marker=cdc[code][0],ms=20,mfc='None',mew=3,lw=3,ls='--',markevery=[-1])

        slope, iner, *var = linregress(P_STARS[0:len(cdf_sl)],cdf_sl[0:])
        plt.plot(P_STARS[0:len(cdf)],slope*P_STARS[0:len(cdf)]+iner,':',c=cdc[code][1],lw=2)

        print(f'Power law fit: {slope:.3f}, initial drop is {cdf[1]*100:.1f}%')

    plt.show(block=False)

def plotGap(df,sigma=3,X_C='similarity',Y_C='gapX',H_C='norm_OVR'):
    df = df.loc[df['similarity'] < 95]

    var = np.var(df[df.pval_S2>0]['pval_S'])

    df = df.loc[(df['norm_OVR'] > -.01) & (df['pval_S']>1*sigma*var)]
    #df_scaled.pval_S = df_scaled.pval_S*-1
    #df_scaled.pval_S2 = df_scaled.pval_S2*-1
    #df = df_scaled
    df['gapX'] = df.norm_OVR - df.similarity/100

    g = sns.FacetGrid(df,col="code", height=4,col_wrap=2,col_order=['MUT','MPA','DNO'])
    def fcc(x,y,c,**kwargs):
        kwargs.pop('color')
        plt.scatter(x,y,c=c,norm=Normalize(0,15),cmap='plasma',alpha=0.6,**kwargs)
        #plt.plot([0,100],[0,1],'k--')

    g.map(fcc, X_C,Y_C,H_C)
    #g.map(sns.kdeplot,'pval_S2','norm_OVR')
    plt.show(0)

def biaxial(cath,scop):
    cath = loadDF('Data_CATH_90',csv_save=True)
    scop = loadDF('Data_SCOP_90',csv_save=True)
    
    mut = set(cath['id']).intersection(set(scop['id']))
    cath_only = set(cath['id']).difference(set(scop['id']))
    scop_only = set(scop['id']).difference(set[cath['id']))

