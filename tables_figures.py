#This code has been tested with Python version 3.10.4
import pandas as pd #Tested with version 1.4.2
import numpy as np #Tested with version 1.22.3
from scipy.stats import norm, chi2 #Tested with version 1.8.0
import seaborn as sns #Tested with version 0.11.2
import matplotlib.pyplot as plt #Tested with version 3.5.2


def get_dist_vals(x, s, par=False):
    '''format data for table'''
    sf = '{:.' + str(s) + 'f}'
    if par:
        s1 = '('
        s2 = ')'
    else:
        s1 = ''
        s2 = ''
    return (' & ' + s1 + str(sf.format(x.mean())) + s2 +
            ' & ' + s1 + str(sf.format(x.min())) + s2 +
            ' & ' + s1 + str(sf.format(x.max())) + s2 +
            ' & ' + s1 + str(sf.format(np.quantile(x, .5))) + s2 +
            ' & ' + s1 + str(sf.format(np.quantile(x, .75))) + s2 +
            ' & ' + s1 + str(sf.format(np.quantile(x, .95))) + s2 +
            ' & ' + s1 + str(sf.format(np.quantile(x, .999))) + s2 +
            ' \\\\')


def gen_tab_prep_fig(file, N, s1=1, s2=4):
    '''Prepare data for the figure and output data for the table'''
    df = pd.read_csv(file)
    df = df.reset_index()
    df.drop(columns=('index'), inplace=True)
    df['id'] = df.index
    
    df.rename(columns={'as_cor_alg': 'as_cor_rerand5'}, inplace=True)


    dfl = pd.wide_to_long(df, stubnames=['as_cor_rerand'],
                         i='id', j='type').reset_index()

    return dfl, ('$p_A=1$' + get_dist_vals(df['as_cor_rerand1'], s2) + '\n' +
            '\\addlinespace\n' + 
            '$p_A=0.1$' + get_dist_vals(df['as_cor_rerand2'], s2) + '\n' +
            '\\addlinespace\n' + 
            '$p_A=0.01$' + get_dist_vals(df['as_cor_rerand3'], s2) + '\n' +
            '\\addlinespace\n' + 
            '$p_A=0.001$' + get_dist_vals(df['as_cor_rerand4'], s2) + '\n' +
            '\\addlinespace\n' + 
            'PS-alg' + get_dist_vals(df['as_cor_rerand5'], s2) + '\n' +
            '\\addlinespace\n')


def gen_seaborn_figs(dfn, dfl, outfile):
    '''Generate figure'''
    di1 = {1: '$p_A=1$', 2: '$p_A=.1$', 3: '$p_A=.01$', 4: '$p_A=.001$',
           5: 'PS-alg'}
    di2 = {2: '$p_A=.1$', 3: '$p_A=.01$', 4: '$p_A=.001$', 5: 'PS-alg'}
    
    dfn['lab'] = dfn['type'].replace(di1)
    dfl['lab'] = dfl['type'].replace(di1)
    
    fc = 18
    
    plt.rc('xtick', labelsize=fc-2) 
    plt.rc('ytick', labelsize=fc-2)
    
    fig, axes = plt.subplots(1, 2, figsize=(16,8))
    
    sns.stripplot('lab', 'as_cor_rerand', jitter=.4, data=dfn,
                       color='#F5793A', s=1.5, alpha=1, ax=axes[0])
    axes[0].set_ylim([0.0198, 0.0242])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('')
    axes[0].set_yticks([.02, .021, .022, .023, .024])
    axes[0].axhline(1/(N-1), ls='--')
    
    sns.stripplot('lab', 'as_cor_rerand', jitter=.4, data=dfl,
                       color='#F5793A', s=1.5, alpha=1, ax=axes[1])
    axes[1].set_ylim([0, .18])
    axes[1].set_xlabel('')
    axes[1].set_ylabel('')
    axes[1].set_yticks([0, .05, .1, .15])
    axes[1].axhline(1/(N-1), ls='--')
    
    axes[0].set_title('Normal covariates,' + '\n' +
                        'assignment correlation', fontsize=fc)
    axes[1].set_title('Log-normal covariates,' + '\n' +
                        'assignment correlation', fontsize=fc)

    plt.savefig(outfile, bbox_inches='tight', pad_inches=0)


#Seed is only relevant to exactly reproduce the jitter in the stripplots
np.random.seed(12345)

#The code can run for N=100 as well. The y-axis have to be adjusted in the
#gen_seaborn_figs() function to display the y-axis correctly
N = 50

#Files where results from the simulation in sim_rerand_alg.jl are stored
filename_normal = 'rerand_alg_normal_' + str(N) + '.csv'
filename_lognormal =  'rerand_alg_lognormal_' + str(N) + '.csv'

dfnorm, a1 = gen_tab_prep_fig(filename_normal, N)
dflognorm, a2 = gen_tab_prep_fig(filename_lognormal, N)

#Generate Figure 1 in the paper (and Figure A1 in the supplementary material
#for N=100)
gen_seaborn_figs(dfnorm, dflognorm, 'fig_rerand_alg_' + str(N) + '.pdf')


q = ('\\begin{tabular*}{\\textwidth}{@{\\hskip\\tabcolsep\\extracolsep\\fill}'
     'l*{1}{cccccccc}}\\toprule\n & & & & \\multicolumn{4}{c}{Quantiles}\\\\\n' 
     '\\cmidrule(lr){5-8}& Mean & Min & Max & 0.5 & 0.75 & 0.975 & 0.999\\\\\n' 
     '\\midrule\\addlinespace\\multicolumn{8}{l}{\\textit{Five standard normal'
     ' covariates}}\\\\\n\\addlinespace\n')
w = ('\\addlinespace\\multicolumn{8}{l}{\\textit{Five log-normal covariates}}'
     '\\\\\n\\addlinespace\n')
v = '\\bottomrule\\end{tabular*}'

#Generate table 1 (and Table A2 in the supplementary material for N=100)
out = q + a1 + w + a2 + v
with open('tab_dist_as_cor_' + str(N) + '.tex', 'w') as f:
    f.write(out)
