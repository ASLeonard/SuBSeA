import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import json
from collections import defaultdict
import numpy as np
from itertools import product

from scipy.stats import linregress, ks_2samp, anderson_ksamp, mannwhitneyu, epps_singleton_2samp, brunnermunzel

from utility import loadCSV
from domains import readDomains, domainMatch, duplicateIntersection


def groupHeteromericBSAs(threshold=100,filter_immunoglobulins=False):

    def domainMatchFractional(domains,c1,c2):
        return len(duplicateIntersection(domains[c1],domains[c2])) / max(len(domains[c1]),len(domains[c2]))

    df = loadCSV(f'Heteromeric_complexes_{threshold}.csv')
    #bsa_het = convertH(df_het,filter_immunoglobulins)

    rows = []
    for _,row in df.iterrows():
        for interface in row['interfaces']:
            rows.append({'id':f'{row["PDB_id"]}_{interface.replace("-","_")}','shared':domainMatchFractional(row['domains'],*interface.split('-')),'BSA':row['BSAs'][interface]})
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
    ax = sns.violinplot(x="ancestry", y="BSA", hue='ancestry',data=data, palette="muted",scale_hue=False,inner="quartile",cut=0,bw=.3,split='True')
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
        domain_info = ';'.join([chain+':{}'.format(tuple(domain_dict[row['id'][:4]][chain])) if chain in domain_dict[row['id'][:4]] else '' for chain in row['chains']] )
        
        new_rows.append({'PDB_id':row['id'][:4], 'interfaces':'-'.join(sorted(row['chains'])), 'domains':domain_info, 'BSAs':round(row['BSA']) if not pd.isna(row['BSA']) else ''})

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
        return pd.read_csv(f'{file_ID}_comparison.csv')
    raw_data = (loadDict if json_save else loadND)(file_ID) 
    rows = []
    for results in (raw_data.values() if json_save else raw_data):
        if results == 'error':
            continue
        pval, N_hits, similarity = results
        rows.append({'pval':pval or 1,'similarity':similarity,'sg':int(similarity//5),'hits':N_hits})
        
    df = pd.DataFrame(rows)
    df['pval'] = np.log10(df['pval'])
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
        for i in range(len(DFs)):
            DFs[i] = DFs[i].sample(sample)
        

    for i in range(len(DFs)):
        df=DFs[i]
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
    #else:
    #    labels = datas
    #    datas = [loadDict(d) for d in labels]
    #if merge_nones:
    #    cleaned_datas = [np.log10([val[0] or 1 for val in data.values() if (val!='error' and val[1]<300)]) for data in datas]
    #else:
    #    cleaned_datas = [np.log10(list(filter(lambda x: isinstance(x[0],float),data.values()))) for data in datas]
    #
    cleaned_datas= datas
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


from scipy.stats import brunnermunzel, pareto, pearsonr, spearmanr,kendalltau,theilslopes
def plotCats(df,ax=None,ls='-',cmap='green',pval_type=''):
    N=6
    if ax is None:
        f, ax =plt.subplots(3,2)
    ax= ax.flatten()
    #cmap = plt.get_cmap(cmap)
    for i in range(0,N):
        if len(df.loc[df['sg']==i][f'pval{pval_type}']) < 10:
            continue
        print(len(df.loc[df['sg']==i][f'pval{pval_type}']))

        color = cmap#cmap(i/(2*N))
        
        sns.distplot(a=np.clip(df.loc[df['sg']==i][f'pval{pval_type}'],-10,0),bins=np.linspace(-10,0,101),ax=ax[i],norm_hist=True,color=color,label=i,kde=False,kde_kws={'cut':0,'kernel':'epa','ls':ls},hist_kws={'histtype':'step','alpha':1,'lw':2})#kde_kws={'ls':ls,'alpha':1})
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

#ff= pd.concat([mp,rp],ignore_index=1)
def plotGrid(df):
    df = df.loc[df['sg']<6]
    df = df.loc[df['hits']>0]
    df = df.loc[df['hits']<=8]
    
    #print(len(df))

    #match_C = sum(df['match']=='match')
    #random_C = sum(df['match']=='random')
    
    #drop_samp = np.random.choice(np.arange(match_C,len(df)),size=random_C-match_C,replace=False)
    #print(np.mean(drop_samp))
    
    #df.drop(df.index[drop_samp],inplace=True)



    g = sns.FacetGrid(df, row="sg", col="hits",hue='match', margin_titles=True,sharex=True,sharey=True)
        
    g.map(plt.hist, "pval_T", bins=np.linspace(-15,0,201),density=0,histtype='step',alpha=0.75).add_legend()
    g.set(yscale = 'log',ylim=[1,1e5])#1e-5,1])
    #return g
    plt.show(0)


def explore(df,filter_fails=True):
    
    list_of_cmaps=['Blues_r','Purples_r','Greens_r','Reds_r']
    if filter_fails:
        df = df.loc[df['hits']>0]
        df = df.loc[df['similarity']<1000]
        
    df['glow'] = df['overlap']/df['align_length']
    df['sco'] = df['score']/df['align_length']
    df['similarity'] = df['similarity']/100
    g = sns.PairGrid(df,hue_order=['M','MP','R','RP'],hue_kws={"cmap":list_of_cmaps},hue='match',vars=['sco','glow','similarity','pval_F'])#'pval_F','pval_T','pval_S'])#['score','glow','pval_F','pval_S','pval_T','gaps'])
    g.map_upper(plt.scatter,alpha=0.5,s=6)
    #g.map_diag(sns.kdeplot, lw=2,alpha=0.8)
    #g.map_lower(plt.hexbin)
    #g.map_lower(sns.kdeplot,n_levels=10,gridsize=100,bw=.05)
    g.add_legend()
    
    plt.show(0)

from matplotlib.colors import LogNorm, Normalize
def hexbin(x, y, color, **kwargs):
    cmap = plt.get_cmap('cividis')
    #cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, cmap=cmap, **kwargs)
    plt.text(-10,.1,len(x),ha='left',va='bottom')
    plt.colorbar()
    print(spearmanr(x,y))
    print(np.median(x),np.mean(x))

from scipy.stats import expon

def hexIT(df,X_C='pval_S',Y_C='norm_OVR',sim_thresh=95,sigma=-1*np.log10(.05)):
    extent_codes = {'pval_F':(0,20),'pval_F2':(0,20),'pval_S':(0,20),'pval_S2':(0,20),'pval_T':(0,20),'pval_T2':(0,20),'similarity':(0,100),'norm_OVR':(0,1),'norm_SCR':(0,3),'split':(0,3)}

    df = df.loc[df['similarity'] <= sim_thresh]

    var = np.var(df[X_C])#[df[X_C]>0][X_C])
    print(var)

    val_cut = sigma if isinstance(sigma,float) else sigma*var
    df = df[df[X_C]>val_cut]
    
    df_scaled = df.loc[(df['norm_OVR'] > -.01) & (df[X_C]>sigma*var)]
    #df.pval_S2 = df.pval_S2*-1
    #print(len(df_scaled))
    
    #fig=plt.figure()

  
    
    ax=sns.lmplot(x=X_C, y=Y_C, hue="code",hue_order=['DNO','MPA','MUT'], data=df, markers=["o", "P",'d'], palette={'MUT':'darkorange','DNO':'royalblue','MPA':'forestgreen'},scatter_kws={'alpha':.75,'s':150,'facecolor':'None','lw':3},robust=False,truncate=True)
    plt.plot([-np.log10(.05)]*2,[-1,2])
    #plt.plot([0,20],[0,-20],'k--',lw=3)
    #g.map(sns.residplot,data=df,x='pval_S2',y='norm_OVR',robust=True)

def plotHeteromericConfidenceDecay(df,x_var = 'pval_S',sigma=-1,sim_thresh=90):

    df = df.loc[df['similarity'] <= sim_thresh]
    var = np.var(df[x_var])
    print(f'Variance in x_var is {var}')
    val_cut = sigma if isinstance(sigma,float) else sigma*var    
    df = df[df[x_var]>val_cut]
    
    
    plt.figure()
    cdc= {'MUT':('H','darkorange'),'MPA':('P','forestgreen'),'DNO':('p','royalblue')}
    df.loc[df.code=='MPA',1] = 'MUT'
    P_STARS = np.linspace(0,15,101)
    
    for ind,code in enumerate(('MUT','MPA','DNO')):
        df_c = df[df.code==code][x_var]
        if len(df_c)==0:
            continue
        print(f'Heteromeric grouping: {code}, {len(df_c)} entries.')
        
        cdf = np.array([np.sum(df_c>=p_star) for p_star in P_STARS])
        cdf = np.log10(cdf/cdf.max())
                
        #cut off data to prevent outliers dominating
        cdf = cdf[cdf>=np.log10(.01)]
        #similarly trim data so fit is only in statistically strong region
        cdf_sl = cdf[cdf>=np.log10(.01)]
        
        plt.plot(P_STARS[:len(cdf)],cdf,c=cdc[code][1],marker=cdc[code][0],ms=20,mfc='None',mew=3,lw=3,ls='--',markevery=[-1])

        slope, iner, *var = linregress(P_STARS[0:len(cdf_sl)],cdf_sl[0:])
        plt.plot(P_STARS[0:len(cdf)],slope*P_STARS[0:len(cdf)]+iner,':')
        
        print(f'Power law fit: {slope:.3f}, initial drop is {cdf[1]*100:.1f}%')

    plt.show(0)

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
