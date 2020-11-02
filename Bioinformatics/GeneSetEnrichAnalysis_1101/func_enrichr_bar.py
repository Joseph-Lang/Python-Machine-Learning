#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math


# In[2]:


def enrichment_bar(table, xlim1=[0, 50], xlim2=[100, 150], ann=False):
    def Terminology(x, split):
        ls = x.split(split)
        ls.reverse()
        ls[-1] = '(' + ls[-1] + ')'
        term = ' '.join(ls)

        return term
    
    if ann=='GO':
        table['Term'] = table['Term'].apply(Terminology, args='~')
    elif ann=='InterPro':
        split=':'
        table['Term'] = table['Term'].apply(Terminology, args=':') 
    
    cols = ['Term', 'Benjamini']
    table = table[cols]
    table['-log(p.adj)'] = table['Benjamini'].apply(lambda v: -math.log10(v))
    table.sort_values(by='-log(p.adj)', ascending=False, inplace=True)
    table = table[0:15]
    
    fig = plt.figure(figsize=(8, 5))
    gs = fig.add_gridspec(1, 3, width_ratios=(1, 1, 0.1),
                          left=0.001, right=0.9, top=0.99, bottom=0.085,
                          wspace=0.1)

    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax)
    ax3 = fig.add_subplot(gs[0, 2])

    sns.barplot(x='-log(p.adj)', y='Term', data=table, palette='YlOrRd_r', ax=ax)
    sns.barplot(x='-log(p.adj)', y='Term', data=table, palette='YlOrRd_r', ax=ax2)
    yticklabels = table['Term'].to_list()

    ax.set_xlabel('-log(p.adj)', x=1.0, fontsize=18)
    ax.set_ylabel(None)
    ax.set_xlim(xlim1)
    ax.set_yticklabels([])
    ax.tick_params(axis='both', length=0)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    sns.despine(ax=ax)

    ax2.set_xlabel('')
    ax2.set_ylabel(None)
    ax2.set_xlim(xlim2)
    ax2.tick_params(axis='both', length=0)
    ax2.spines['bottom'].set_linewidth(2)
    sns.despine(ax=ax2, left=True)

    # how big to make the diagonal lines in axes coordinates
    d = .015
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d), (-d,+d), **kwargs)
    # switch to the bottom axes
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d,+d), (-d,+d), **kwargs)

    ax3.set_ylim(ax.get_ylim())
    ax3.set_yticks(ax.get_yticks())
    ax3.set_yticklabels(yticklabels, fontsize=18)
    ax3.set_xticklabels([])
    ax3.tick_params(axis='x', length=0)
    ax3.yaxis.tick_right()
    sns.despine(ax=ax3, left=True, bottom=True)

    plt.tight_layout()


# In[3]:


dat = pd.read_table('../mission_10_19/Result2_Dataset/GeneSetA_INTERPRO.txt')
enrichment_bar(dat, xlim1=[0, 45], xlim2=[100, 150], ann='InterPro')

