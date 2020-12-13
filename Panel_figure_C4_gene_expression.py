#!/usr/bin/env python
# coding: utf-8

# In[1]:

import sys
import glob
import pandas as pd
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.pyplot import subplots_adjust
import matplotlib.ticker as mticks
import matplotlib.gridspec as gridspec


# ### Maize mesophyll vs bundle sheath cell expression

# #### Chang et al 2012 data

# In[2]:


def maize_M_BS_chang(ax, orthogroup, orthogroup_fasta_file, zm_colour_dict, targetp_dict):
    maize_genes = get_species_genes(orthogroup_fasta_file, 'Zm0', '_', 1)
    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Zm0' in key:
            orthogroup_location_dict[(key.split('_'))[0]] = targetp_dict[key]

    chang_df = pd.read_csv("data/expression_data/SRA047278_maize_chang_TPMs_concant_gene_models.csv", delim_whitespace=True)
    chang_df = chang_df.loc[chang_df['Name'].isin(maize_genes)]
    chang_df['M_mean_TPM'] = chang_df[['SRR354212_TPM', 'SRR354213_TPM']].mean(axis = 1)
    chang_df['BS_mean_TPM'] = chang_df[['SRR354214_TPM', 'SRR354215_TPM']].mean(axis = 1)

    if plot_highly_expressed:
        chang_df = chang_df[(chang_df['M_mean_TPM'] + chang_df['BS_mean_TPM']) > min_mean_TPM]
        if len(chang_df) == 0:
            chang_df = chang_df.append(pd.Series(0, index=chang_df.columns), ignore_index=True)
            chang_df['Name'] = 'All genes'
            zm_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    gray = 0.1
    for i in maize_genes:
        if i not in zm_colour_dict.keys():
            zm_colour_dict[i] = (gray, gray, gray, 1.0)
            if gray > 0.95:
                gray = 0.1
            gray = gray + 0.05
    rgb_colours = []
    for i, row in chang_df.iterrows():
        rgb_colours.append(zm_colour_dict[row['Name']])
    chang_df['colour'] = rgb_colours

    chang_df['M_summed_mean'] = chang_df['M_mean_TPM'].groupby(chang_df['colour']).transform('sum')
    chang_df['BS_summed_mean'] = chang_df['BS_mean_TPM'].groupby(chang_df['colour']).transform('sum')

    cumval_M = 0
    cumval_BS = 0

    sorted_M = chang_df.set_index('Name').sort_values('M_summed_mean')['M_mean_TPM']
    sorted_BS = chang_df.set_index('Name').sort_values('BS_summed_mean')['BS_mean_TPM']
    colours = chang_df.set_index('Name')['colour']

    for name, value in sorted_M.iteritems():
        ax.bar('M', value, bottom=cumval_M, label=(orthogroup_location_dict[name] + ' ' + name), color=colours[name])
        cumval_M += value

    for name, value in sorted_BS.iteritems():
        ax.bar('BS', value, bottom=cumval_BS, color=colours[name])
        cumval_BS += value

    ax.set_ylabel('TPM')
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# #### Tausta et al 2012 data

# In[3]:


def maize_M_BS_tausta(ax1, ax2, orthogroup, orthogroup_fasta_file, zm_colour_dict, targetp_dict):
    maize_genes = get_species_genes(orthogroup_fasta_file, 'Zm0', '_', 1)
    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Zm0' in key:
            orthogroup_location_dict[(key.split('_'))[0]] = targetp_dict[key]

    tausta_df = pd.read_csv("data/expression_data/SRP035577_merged_TPMs_concated_gene_models.csv", delim_whitespace=True)
    tausta_df = tausta_df.loc[tausta_df['Name'].isin(maize_genes)]

    tausta_df['M_mean_sec_4_TPM'] = tausta_df[['GSM1311350_TPM', 'GSM1311351_TPM']].mean(axis = 1)
    tausta_df['M_mean_sec_9_TPM'] = tausta_df[['GSM1311354_TPM', 'GSM1311355_TPM']].mean(axis = 1)
    tausta_df['M_mean_sec_14_TPM'] = tausta_df[['GSM1311358_TPM', 'GSM1311359_TPM']].mean(axis = 1)
    tausta_df['BS_mean_sec_4_TPM'] = tausta_df[['GSM1311348_TPM', 'GSM1311349_TPM']].mean(axis = 1)
    tausta_df['BS_mean_sec_9_TPM'] = tausta_df[['GSM1311352_TPM', 'GSM1311353_TPM']].mean(axis = 1)
    tausta_df['BS_mean_sec_14_TPM'] = tausta_df[['GSM1311356_TPM', 'GSM1311357_TPM']].mean(axis = 1)

    if plot_highly_expressed:
        tausta_df = tausta_df[((tausta_df['M_mean_sec_4_TPM'] + tausta_df['BS_mean_sec_4_TPM']) > min_mean_TPM) |
        ((tausta_df['M_mean_sec_9_TPM'] + tausta_df['BS_mean_sec_9_TPM']) > min_mean_TPM) |
        ((tausta_df['M_mean_sec_14_TPM']+ tausta_df['BS_mean_sec_14_TPM']) > min_mean_TPM)]

        if len(tausta_df) == 0:
            tausta_df = tausta_df.append(pd.Series(0, index=tausta_df.columns), ignore_index=True)
            tausta_df['Name'] = 'All genes'
            zm_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    for i in maize_genes:
        if i not in zm_colour_dict.keys():
            zm_colour_dict[i] = (0.1, 0.3, 0.5, 1.0)
    rgb_colours = []
    for i, row in tausta_df.iterrows():
        rgb_colours.append(zm_colour_dict[row['Name']])
    tausta_df['colour'] = rgb_colours

    lines = ["-","--","-.",":"]
    tausta_df.sort_values(by=['colour'], inplace=True)
    tausta_df.reset_index(inplace=True)
    maize_tp_name = []
    line_styles = []
    count = 1
    for i, row in tausta_df.iterrows():
        maize_tp_name.append(orthogroup_location_dict[row['Name']] + ' ' + row['Name'])
        if i == 0:
            line_styles.append(lines[0])
            count = 1
        elif (tausta_df.iloc[i]['colour']) == (tausta_df.iloc[i-1]['colour']):
            if count == len(lines):
                count = 0
            line_styles.append(lines[count])
            count+=1
        else:
            line_styles.append(lines[0])
            count = 1

    tausta_df['maize_tp_name'] = maize_tp_name
    tausta_df['line'] = line_styles

    tausta_df = tausta_df.set_index('maize_tp_name')

    tausta_df[['M_mean_sec_4_TPM', 'M_mean_sec_9_TPM', 'M_mean_sec_14_TPM']].transpose().plot(ax=ax1, color=tausta_df['colour'], style=list(tausta_df['line']), linewidth=1.5)
    tausta_df[['BS_mean_sec_4_TPM', 'BS_mean_sec_9_TPM', 'BS_mean_sec_14_TPM']].transpose().plot(ax=ax2, legend=False, color=tausta_df['colour'], style=list(tausta_df['line']), linewidth=1.5)
    N=3
    ax1.set_ylabel('TPM')
    ax1.set_xticklabels(['-1', '+4', '+9'])
    ax2.set_xticklabels(['-1', '+4', '+9'])
    ax1.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))
    ax2.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))

    ax1.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax1.set_title('Maize M (Tausta et al 2014)')
    ax2.set_title('Maize BS (Tausta et al 2014)')

    ax1.set_xlabel('cm from transition zone')
    ax2.set_xlabel('cm from transition zone')


    # Put a legend to the right of the current axis
    ax1.legend(loc='center', bbox_to_anchor=(0.5, 0.75), fontsize=10)


# #### Denton et al 2017 data

# In[4]:


def maize_M_BS_denton(ax1, ax2, orthogroup, orthogroup_fasta_file, zm_colour_dict, targetp_dict):
    maize_genes = get_species_genes(orthogroup_fasta_file, 'Zm0', '_', 1)
    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Zm0' in key:
            orthogroup_location_dict[(key.split('_'))[0]] = targetp_dict[key]

    denton_df = pd.read_csv("data/expression_data/SRP052802_merged_TPMs_concated_gene_models_denton.csv", delim_whitespace=True)
    denton_df = denton_df.loc[denton_df['Name'].isin(maize_genes)]

    denton_df['M_mean_sec_1_TPM'] = denton_df[['SRR2186729_TPM', 'SRR2186713_TPM', 'SRR2186637_TPM']].mean(axis = 1)
    denton_df['M_mean_sec_2_TPM'] = denton_df[['SRR2186732_TPM', 'SRR2186715_TPM', 'SRR2186662_TPM']].mean(axis = 1)
    denton_df['M_mean_sec_3_TPM'] = denton_df[['SRR2186738_TPM', 'SRR2186719_TPM', 'SRR2186705_TPM']].mean(axis = 1)
    denton_df['M_mean_sec_4_TPM'] = denton_df[['SRR2186625_TPM', 'SRR2186624_TPM', 'SRR2186623_TPM']].mean(axis = 1)
    denton_df['M_mean_sec_5_TPM'] = denton_df[['SRR2186622_TPM', 'SRR2186621_TPM', 'SRR2186620_TPM']].mean(axis = 1)
    denton_df['BS_mean_sec_1_TPM'] = denton_df[['SRR2186726_TPM', 'SRR2186711_TPM', 'SRR2186626_TPM']].mean(axis = 1)
    denton_df['BS_mean_sec_2_TPM'] = denton_df[['SRR2186730_TPM', 'SRR2186714_TPM', 'SRR2186654_TPM']].mean(axis = 1)
    denton_df['BS_mean_sec_3_TPM'] = denton_df[['SRR2186734_TPM', 'SRR2186717_TPM', 'SRR2186703_TPM']].mean(axis = 1)
    denton_df['BS_mean_sec_4_TPM'] = denton_df[['SRR2186739_TPM', 'SRR2186720_TPM', 'SRR2186706_TPM']].mean(axis = 1)
    denton_df['BS_mean_sec_5_TPM'] = denton_df[['SRR2186741_TPM', 'SRR2186723_TPM', 'SRR2186708_TPM']].mean(axis = 1)

    if plot_highly_expressed:
        denton_df = denton_df[((denton_df['M_mean_sec_1_TPM'] + denton_df['BS_mean_sec_1_TPM']) > min_mean_TPM) |
        ((denton_df['M_mean_sec_2_TPM'] + denton_df['BS_mean_sec_2_TPM']) > min_mean_TPM) |
        ((denton_df['M_mean_sec_3_TPM'] + denton_df['BS_mean_sec_3_TPM']) > min_mean_TPM) |
        ((denton_df['M_mean_sec_4_TPM'] + denton_df['BS_mean_sec_4_TPM']) > min_mean_TPM) |
        ((denton_df['M_mean_sec_5_TPM'] + denton_df['BS_mean_sec_5_TPM']) > min_mean_TPM)]

        if len(denton_df) == 0:
            denton_df = denton_df.append(pd.Series(0, index=denton_df.columns), ignore_index=True)
            denton_df['Name'] = 'All genes'
            zm_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    for i in maize_genes:
        if i not in zm_colour_dict.keys():
            zm_colour_dict[i] = (0.1, 0.3, 0.5, 1.0)
    rgb_colours = []
    for i, row in denton_df.iterrows():
        rgb_colours.append(zm_colour_dict[row['Name']])
    denton_df['colour'] = rgb_colours

    lines = ["-","--","-.",":"]
    denton_df.sort_values(by=['colour'], inplace=True)
    denton_df.reset_index(inplace=True)
    maize_tp_name = []
    line_styles = []
    count = 1
    for i, row in denton_df.iterrows():
        maize_tp_name.append(orthogroup_location_dict[row['Name']] + ' ' + row['Name'])
        if i == 0:
            line_styles.append(lines[0])
            count = 1
        elif (denton_df.iloc[i]['colour']) == (denton_df.iloc[i-1]['colour']):
            if count == len(lines):
                count = 0
            line_styles.append(lines[count])
            count+=1
        else:
            line_styles.append(lines[0])
            count = 1

    denton_df['maize_tp_name'] = maize_tp_name
    denton_df['line'] = line_styles

    denton_df = denton_df.set_index('maize_tp_name')

    denton_df[['M_mean_sec_5_TPM', 'M_mean_sec_4_TPM', 'M_mean_sec_3_TPM', 'M_mean_sec_2_TPM', 'M_mean_sec_1_TPM']].transpose().plot(ax=ax1, color=denton_df['colour'], style=list(denton_df['line']), linewidth=1.5)
    denton_df[['BS_mean_sec_5_TPM', 'BS_mean_sec_4_TPM', 'BS_mean_sec_3_TPM', 'BS_mean_sec_2_TPM', 'BS_mean_sec_1_TPM']].transpose().plot(ax=ax2, legend=False, color=denton_df['colour'], style=list(denton_df['line']), linewidth=1.5)
    N=5
    ax1.set_ylabel('TPM')
    ax1.set_xticklabels(['0', '4', '8', '12', '16'])
    ax2.set_xticklabels(['0', '4', '8', '12', '16'])
    ax1.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))
    ax2.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))

    ax1.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax1.set_title('Maize M (Denton et al 2017)')
    ax2.set_title('Maize BS (Denton et al 2017)')

    ax1.set_xlabel('cm from ligule')
    ax2.set_xlabel('cm from ligule')


    # Put a legend to the right of the current axis
    #ax2.legend(loc='center left', bbox_to_anchor=(0.5, 0.5), fontsize=10)
    ax1.legend(loc='center', bbox_to_anchor=(0.5, 0.75), fontsize=10)





# ### Other mesophyll vs bundle sheath cell expression datasets

# #### Emms et al 2016 data for Sourghum bicolor

# In[5]:


def sbicolor_M_BS_Oxford(ax, orthogroup, orthogroup_fasta_file, sb_colour_dict, targetp_dict):
    sbicolor_genes = get_species_genes(orthogroup_fasta_file, 'Sobic', '.', 2)

    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Sobic' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 2)] = targetp_dict[key]

    sbicolor_df = pd.read_csv('data/expression_data/ERP013053_Oxford_2016_S.bicolor-M-BS_merged_TPMs_concat_gene_models.csv', delim_whitespace=True)
    sbicolor_df = sbicolor_df.loc[sbicolor_df['Name'].isin(sbicolor_genes)]
    sbicolor_df['M_mean_TPM'] = sbicolor_df[['ERR1109878_TPM', 'ERR1109879_TPM', 'ERR1109880_TPM']].mean(axis=1)
    sbicolor_df['BS_mean_TPM'] = sbicolor_df[['ERR1109875_TPM', 'ERR1109876_TPM', 'ERR1109877_TPM']].mean(axis=1)

    if plot_highly_expressed:
        sbicolor_df = sbicolor_df[((sbicolor_df['M_mean_TPM'] + sbicolor_df['BS_mean_TPM']) > min_mean_TPM)]
        if len(sbicolor_df) == 0:
            sbicolor_df = sbicolor_df.append(pd.Series(0, index=sbicolor_df.columns), ignore_index=True)
            sbicolor_df['Name'] = 'All genes'
            sb_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    red = 0.1
    for i in sbicolor_genes:
        if i not in sb_colour_dict.keys():
            sb_colour_dict[i] = (red, 0, 0, 1.0)
            if red > 0.95:
                red = 0.1
            red += 0.05

    rgb_colours = []
    for i, row in sbicolor_df.iterrows():
        rgb_colours.append(sb_colour_dict[row['Name']])
    sbicolor_df['colour'] = rgb_colours

    sbicolor_df['M_summed_mean'] = sbicolor_df['M_mean_TPM'].groupby(sbicolor_df['colour']).transform('sum')
    sbicolor_df['BS_summed_mean'] = sbicolor_df['BS_mean_TPM'].groupby(sbicolor_df['colour']).transform('sum')

    cumval_M = 0
    cumval_BS = 0

    sorted_M = sbicolor_df.set_index('Name').sort_values('M_summed_mean')['M_mean_TPM']
    sorted_BS = sbicolor_df.set_index('Name').sort_values('BS_summed_mean')['BS_mean_TPM']
    colours = sbicolor_df.set_index('Name')['colour']

    for name, value in sorted_M.iteritems():
        ax.bar('M', value, bottom=cumval_M, label=(orthogroup_location_dict[name] + ' ' + name), color=colours[name])
        cumval_M += value

    for name, value in sorted_BS.iteritems():
        ax.bar('BS', value, bottom=cumval_BS, color=colours[name])
        cumval_BS += value

    ax.set_ylabel('TPM')
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# #### John et al 2014 data for Setaria viridis

# In[6]:


def setaria_M_BS_john(ax, orthogroup, orthogroup_fasta_file, si_colour_dict, targetp_dict):
    sitalica_genes = get_species_genes(orthogroup_fasta_file, 'Seita', '.', 2)

    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Seita' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 2)] = targetp_dict[key]

    john_df = pd.read_csv('data/expression_data/ERP004434_John_2014_merged_TPMs_concated_gene_models.csv', delim_whitespace=True)
    john_df = john_df.loc[john_df['Name'].isin(sitalica_genes)]
    john_df['M_mean_TPM'] = john_df[['ERR385861_TPM', 'ERR385862_TPM', 'ERR385863_TPM']].mean(axis=1)
    john_df['BS_mean_TPM'] = john_df[['ERR385864_TPM', 'ERR385865_TPM', 'ERR385866_TPM']].mean(axis=1)

    if plot_highly_expressed:
        john_df = john_df[((john_df['M_mean_TPM'] + john_df['BS_mean_TPM']) > min_mean_TPM)]
        if len(john_df) == 0:
            john_df = john_df.append(pd.Series(0, index=john_df.columns), ignore_index=True)
            john_df['Name'] = 'All genes'
            si_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    green = 0.1
    for i in sitalica_genes:
        if i not in si_colour_dict.keys():
            si_colour_dict[i] = (0, green, 0, 1.0)
            if green > 0.95:
                green = 0.1
            green += 0.05

    rgb_colours = []
    for i, row in john_df.iterrows():
        rgb_colours.append(si_colour_dict[row['Name']])
    john_df['colour'] = rgb_colours

    john_df['M_summed_mean'] = john_df['M_mean_TPM'].groupby(john_df['colour']).transform('sum')
    john_df['BS_summed_mean'] = john_df['BS_mean_TPM'].groupby(john_df['colour']).transform('sum')

    cumval_M = 0
    cumval_BS = 0

    sorted_M = john_df.set_index('Name').sort_values('M_summed_mean')['M_mean_TPM']
    sorted_BS = john_df.set_index('Name').sort_values('BS_summed_mean')['BS_mean_TPM']
    colours = john_df.set_index('Name')['colour']

    for name, value in sorted_M.iteritems():
        ax.bar('M', value, bottom=cumval_M, label=(orthogroup_location_dict[name] + ' ' + name), color=colours[name])
        cumval_M += value

    for name, value in sorted_BS.iteritems():
        ax.bar('BS', value, bottom=cumval_BS, color=colours[name])
        cumval_BS += value

    ax.set_ylabel('TPM')
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# #### Rao et al 2016 data for Panicum virgatum

# In[7]:


def pvirgatum_M_BS_rao(ax, orthogroup, orthogroup_fasta_file, pv_colour_dict, targetp_dict):
    pvirgatum_genes = get_species_genes(orthogroup_fasta_file, 'Pavir', '.', 2)

    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Pavir' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 2)] = targetp_dict[key]

    rao_df = pd.read_csv('data/expression_data/SAMN040029664_Pvigatum_TPMs_Rao_et_al.csv', delim_whitespace=True)
    rao_df = rao_df.loc[rao_df['Name'].isin(pvirgatum_genes)]
    rao_df['M_mean_TPM'] = rao_df[['SRR3217256_TPM', 'SRR3217257_TPM']].mean(axis=1)
    rao_df['BS_mean_TPM'] = rao_df[['SRR3217892_TPM', 'SRR3217893_TPM']].mean(axis=1)

    if plot_highly_expressed:
        rao_df = rao_df[((rao_df['M_mean_TPM'] + rao_df['BS_mean_TPM']) > min_mean_TPM)]
        if len(rao_df) == 0:
            rao_df = rao_df.append(pd.Series(0, index=rao_df.columns), ignore_index=True)
            rao_df['Name'] = 'All genes'
            pv_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    blue = 0.1
    for i in pvirgatum_genes:
        if i not in pv_colour_dict.keys():
            pv_colour_dict[i] = (0.5, 0.5, blue, 1.0)
            if blue > 0.95:
                blue = 0.1
            blue += 0.01

    rgb_colours = []
    for i, row in rao_df.iterrows():
        rgb_colours.append(pv_colour_dict[row['Name']])
    rao_df['colour'] = rgb_colours

    rao_df['M_summed_mean'] = rao_df['M_mean_TPM'].groupby(rao_df['colour']).transform('sum')
    rao_df['BS_summed_mean'] = rao_df['BS_mean_TPM'].groupby(rao_df['colour']).transform('sum')

    cumval_M = 0
    cumval_BS = 0

    sorted_M = rao_df.set_index('Name').sort_values('M_summed_mean')['M_mean_TPM']
    sorted_BS = rao_df.set_index('Name').sort_values('BS_summed_mean')['BS_mean_TPM']
    colours = rao_df.set_index('Name')['colour']

    for name, value in sorted_M.iteritems():
        ax.bar('M', value, bottom=cumval_M, label=(orthogroup_location_dict[name] + ' ' + name), color=colours[name])
        cumval_M += value

    for name, value in sorted_BS.iteritems():
        ax.bar('BS', value, bottom=cumval_BS, color=colours[name])
        cumval_BS += value

    ax.set_ylabel('TPM')
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# ### Other expression datasets of relevance

# #### van Campen et al 2016 data for rice leaf development

# In[8]:


def rice_dev_vancampen(ax, orthogroup, orthogroup_fasta_file, os_colour_dict, targetp_dict):
    rice_genes = get_species_genes(orthogroup_fasta_file, 'LOC_Os', '.', 1)

    #os_colour_dict = get_species_colour_dict('Osativa', zm_colour_dict, orthogroup, '.', 1)

    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'LOC_Os' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 1)] = targetp_dict[key]

    vancampen_df = pd.read_csv('data/expression_data/SRP062323_Fleming_merged_TPMs_concated_gene_models.csv', delim_whitespace=True)
    vancampen_df = vancampen_df.loc[vancampen_df['Name'].isin(rice_genes)]

    vancampen_df['P3_mean_TPM'] = vancampen_df[['SRR2156305_TPM', 'SRR2156307_TPM']].mean(axis = 1)
    vancampen_df['P4_mean_TPM'] = vancampen_df[['SRR2156309_TPM', 'SRR2156312_TPM']].mean(axis = 1)
    vancampen_df['P5_mean_TPM'] = vancampen_df[['SRR2156314_TPM', 'SRR2156315_TPM']].mean(axis = 1)

    if plot_highly_expressed:
        vancampen_df = vancampen_df[(vancampen_df['P3_mean_TPM'] > min_mean_TPM) |
        (vancampen_df['P4_mean_TPM'] > min_mean_TPM) |
        (vancampen_df['P5_mean_TPM'] > min_mean_TPM) ]

        if len(vancampen_df) == 0:
            vancampen_df = vancampen_df.append(pd.Series(0, index=vancampen_df.columns), ignore_index=True)
            vancampen_df['Name'] = 'All genes'
            os_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    for i in rice_genes:
        if i not in os_colour_dict.keys():
            os_colour_dict[i] = (0.5, 0.5, 0.5, 1.0)
    rgb_colours = []
    for i, row in vancampen_df.iterrows():
        rgb_colours.append(os_colour_dict[row['Name']])
    vancampen_df['colour'] = rgb_colours

    lines = ["-","--","-.",":"]
    vancampen_df.sort_values(by=['colour'], inplace=True)
    vancampen_df.reset_index(inplace=True)
    gene_tp_name = []
    line_styles = []
    count = 1

    for i, row in vancampen_df.iterrows():
        gene_tp_name.append(orthogroup_location_dict[row['Name']] + ' ' + row['Name'])
        if i == 0:
            line_styles.append(lines[0])
            count = 1
        elif (vancampen_df.iloc[i]['colour']) == (vancampen_df.iloc[i-1]['colour']):
            if count == len(lines):
                count = 0
            line_styles.append(lines[count])
            count+=1
        else:
            line_styles.append(lines[0])
            count = 1

    vancampen_df['gene_tp_name'] = gene_tp_name
    vancampen_df['line'] = line_styles

    vancampen_df = vancampen_df.set_index('gene_tp_name')

    vancampen_df[['P3_mean_TPM', 'P4_mean_TPM', 'P5_mean_TPM']].transpose().plot(ax=ax, color=vancampen_df['colour'], style=list(vancampen_df['line']), linewidth=1.5)

    N=3
    ax.set_ylabel('TPM')
    ax.set_xticklabels(['P3', 'P4', 'P5'])
    ax.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))

    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Developmental stage')



    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# #### Woo et al 2016 data for arabidopsis leaf development

# In[9]:


def arabidopsis_dev_woo(ax, orthogroup, orthogroup_fasta_file, targetp_dict, at_colour_dict):
    arabidopsis_genes = get_species_genes(orthogroup_fasta_file, 'AT', '.', 1)
    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'AT' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 1)] = targetp_dict[key]

    woo_df = pd.read_csv('data/expression_data/SRP018034_Woo_merged_TPMs_concated_gene_models.csv', delim_whitespace=True)
    woo_df = woo_df.loc[woo_df['Name'].isin(arabidopsis_genes)]

    woo_df['4_mean_TPM'] = woo_df[['SRR2079771_TPM', 'SRR2079785_TPM']].mean(axis = 1)
    woo_df['6_mean_TPM'] = woo_df[['SRR2079772_TPM', 'SRR2079786_TPM']].mean(axis = 1)
    woo_df['8_mean_TPM'] = woo_df[['SRR2079773_TPM', 'SRR2079787_TPM']].mean(axis = 1)
    woo_df['10_mean_TPM'] = woo_df[['SRR2079774_TPM', 'SRR2079788_TPM']].mean(axis = 1)
    woo_df['12_mean_TPM'] = woo_df[['SRR2079775_TPM', 'SRR2079789_TPM']].mean(axis = 1)
    woo_df['14_mean_TPM'] = woo_df[['SRR2079776_TPM', 'SRR2079790_TPM']].mean(axis = 1)
    woo_df['16_mean_TPM'] = woo_df[['SRR2079777_TPM', 'SRR2079791_TPM']].mean(axis = 1)
    woo_df['18_mean_TPM'] = woo_df[['SRR2079778_TPM', 'SRR2079792_TPM']].mean(axis = 1)
    woo_df['20_mean_TPM'] = woo_df[['SRR2079779_TPM', 'SRR2079793_TPM']].mean(axis = 1)
    woo_df['22_mean_TPM'] = woo_df[['SRR2079780_TPM', 'SRR2079794_TPM']].mean(axis = 1)
    woo_df['24_mean_TPM'] = woo_df[['SRR2079781_TPM', 'SRR2079795_TPM']].mean(axis = 1)
    woo_df['26_mean_TPM'] = woo_df[['SRR2079781_TPM', 'SRR2079796_TPM']].mean(axis = 1)
    woo_df['28_mean_TPM'] = woo_df[['SRR2079783_TPM', 'SRR2079797_TPM']].mean(axis = 1)
    woo_df['30_mean_TPM'] = woo_df[['SRR2079784_TPM', 'SRR2079798_TPM']].mean(axis = 1)

    if plot_highly_expressed:
        woo_df = woo_df[(woo_df['4_mean_TPM'] > min_mean_TPM) | (woo_df['6_mean_TPM'] > min_mean_TPM) | (woo_df['8_mean_TPM'] > min_mean_TPM) |
        (woo_df['10_mean_TPM'] > min_mean_TPM) | (woo_df['12_mean_TPM'] > min_mean_TPM) | (woo_df['14_mean_TPM'] > min_mean_TPM) |
        (woo_df['16_mean_TPM'] > min_mean_TPM) | (woo_df['18_mean_TPM'] > min_mean_TPM) | (woo_df['20_mean_TPM'] > min_mean_TPM) |
        (woo_df['22_mean_TPM'] > min_mean_TPM) | (woo_df['24_mean_TPM'] > min_mean_TPM) | (woo_df['26_mean_TPM'] > min_mean_TPM) |
        (woo_df['28_mean_TPM'] > min_mean_TPM) | (woo_df['30_mean_TPM'] > min_mean_TPM)]

        if len(woo_df) == 0:
            woo_df = woo_df.append(pd.Series(0, index=woo_df.columns), ignore_index=True)
            woo_df['Name'] = 'All genes'
            at_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    # colours = [(164/255,76/255,64/255,1.0), (133/255,114/255,100/255,1.0), (64/255,70/255,54/255,1.0), (15/2,20/255,12/255,1.0),(74/255,94/255,106/255,1.0)]
    # count = 0
    # for i in arabidopsis_genes:
    #     if i not in at_colour_dict.keys():
    #         at_colour_dict[i] = colours[count]
    #         count += 1

    for i in arabidopsis_genes:
        if i not in at_colour_dict.keys():
            at_colour_dict[i] = (0.5, 0.5, 0.5, 1.0)

    rgb_colours = []
    for i, row in woo_df.iterrows():
        rgb_colours.append(at_colour_dict[row['Name']])
    woo_df['colour'] = rgb_colours

    lines = ["-","--","-.",":"]
    woo_df.sort_values(by=['colour'], inplace=True)
    woo_df.reset_index(inplace=True)
    gene_tp_name = []
    line_styles = []
    count = 1
    for i, row in woo_df.iterrows():
        gene_tp_name.append(orthogroup_location_dict[row['Name']] + ' ' + row['Name'])
        if i == 0:
            line_styles.append(lines[0])
            count = 1
        elif (woo_df.iloc[i]['colour']) == (woo_df.iloc[i-1]['colour']):
            if count == len(lines):
                count = 0
            line_styles.append(lines[count])
            count+=1
        else:
            line_styles.append(lines[0])
            count = 1

    woo_df['gene_tp_name'] = gene_tp_name
    woo_df['line'] = line_styles

    woo_df = woo_df.set_index('gene_tp_name')

    woo_df[['4_mean_TPM', '6_mean_TPM', '8_mean_TPM', '10_mean_TPM', '12_mean_TPM', '14_mean_TPM', '16_mean_TPM', '18_mean_TPM',            '20_mean_TPM', '22_mean_TPM', '24_mean_TPM', '26_mean_TPM', '28_mean_TPM', '30_mean_TPM'           ]].transpose().plot(ax=ax, color=woo_df['colour'], style=list(woo_df['line']), linewidth=1.5)

    N=14
    ax.set_ylabel('TPM')
    ax.set_xticklabels(['4', '6', '8', '10', '12', '14', '16', '18', '20', '22', '24', '26', '28', '30'])
    ax.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))

    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Days after emergence (4th rosette leaf)')



    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# #### Wang et al 2013 data for maize husk/foliar leaf development

# In[10]:


def maize_foliar_husk_wang(ax1, ax2, orthogroup, orthogroup_fasta_file, zm_colour_dict, targetp_dict):
    maize_genes = get_species_genes(orthogroup_fasta_file, 'Zm0', '_', 1)
    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'Zm0' in key:
            orthogroup_location_dict[(key.split('_'))[0]] = targetp_dict[key]

    wang_df = pd.read_csv("data/expression_data/SRP028231_Wang_merged_TPMs_concated_gene_models.csv", delim_whitespace=True)
    wang_df = wang_df.loc[wang_df['Name'].isin(maize_genes)]

    wang_df['F_AM_P1_P2'] = wang_df['SRR942909_TPM']
    wang_df['F_P3_P4'] = wang_df['SRR942915_TPM']
    wang_df['F_P5'] = wang_df['SRR942916_TPM']
    wang_df['F_I'] = wang_df['SRR942911_TPM']
    wang_df['F_E'] = wang_df['SRR942910_TPM']
    wang_df['H_P1_P2'] = wang_df['SRR942912_TPM']
    wang_df['H_P3_P4'] = wang_df['SRR942917_TPM']
    wang_df['H_P5'] = wang_df['SRR942918_TPM']
    wang_df['H_I'] = wang_df['SRR942913_TPM']
    wang_df['H_E'] = wang_df['SRR942914_TPM']

    if plot_highly_expressed:
        wang_df = wang_df[(wang_df['F_AM_P1_P2'] > min_mean_TPM) | (wang_df['F_P3_P4'] > min_mean_TPM) | (wang_df['F_P5'] > min_mean_TPM) |
                   (wang_df['F_I'] > min_mean_TPM) | (wang_df['F_E'] > min_mean_TPM) | (wang_df['H_P1_P2'] > min_mean_TPM) |
                   (wang_df['H_P3_P4'] > min_mean_TPM) | (wang_df['H_P5'] > min_mean_TPM) | (wang_df['H_I'] > min_mean_TPM) |
                   (wang_df['H_E'] > min_mean_TPM)]
        if len(wang_df) == 0:
            wang_df = wang_df.append(pd.Series(0, index=wang_df.columns), ignore_index=True)
            wang_df['Name'] = 'All genes'
            zm_colour_dict['All genes'] = (1, 1, 1, 1)
            orthogroup_location_dict['All genes'] = ''

    for i in maize_genes:
        if i not in zm_colour_dict.keys():
            zm_colour_dict[i] = (0.5, 0.5, 0.5, 1.0)
    rgb_colours = []
    for i, row in wang_df.iterrows():
        rgb_colours.append(zm_colour_dict[row['Name']])
    wang_df['colour'] = rgb_colours

    lines = ["-","--","-.",":"]
    wang_df.sort_values(by=['colour'], inplace=True)
    wang_df.reset_index(inplace=True)
    maize_tp_name = []
    line_styles = []
    count = 1
    for i, row in wang_df.iterrows():
        maize_tp_name.append(orthogroup_location_dict[row['Name']] + ' ' + row['Name'])
        if i == 0:
            line_styles.append(lines[0])
            count = 1
        elif (wang_df.iloc[i]['colour']) == (wang_df.iloc[i-1]['colour']):
            if count == len(lines):
                count = 0
            line_styles.append(lines[count])
            count+=1
        else:
            line_styles.append(lines[0])
            count = 1

    wang_df['maize_tp_name'] = maize_tp_name
    wang_df['line'] = line_styles

    wang_df = wang_df.set_index('maize_tp_name')


    wang_df[['H_P1_P2', 'H_P3_P4', 'H_P5', 'H_I', 'H_E']].transpose().plot(ax=ax1, color=wang_df['colour'], style=list(wang_df['line']), linewidth=1.5)
    wang_df[['F_AM_P1_P2', 'F_P3_P4', 'F_P5', 'F_I', 'F_E']].transpose().plot(ax=ax2, legend=False, color=wang_df['colour'], style=list(wang_df['line']), linewidth=1.5)
    N=5
    ax1.set_ylabel('TPM')
    ax1.set_xticklabels(['P1-P2', 'P3-P4', 'P5', 'I', 'E'])
    ax2.set_xticklabels(['AM-P2', 'P3-P4', 'P5', 'I', 'E'])
    ax1.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))
    ax2.get_xaxis().set_major_locator(mticks.LinearLocator(numticks=N))

    ax1.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax1.set_title('Maize Husk Leaf (Wang et al 2013)')
    ax2.set_title('Foliar Leaf (Wang et al 2013)')

    ax1.set_xlabel('Developmental stage')
    ax2.set_xlabel('Developmental stage')


    # Put a legend to the right of the current axis
    ax1.legend(loc='center', bbox_to_anchor=(0.5, 0.75), fontsize=10)


# #### Gowick et al 2011 for Flaveria C3-C4

# In[11]:


def flaveria_C3_C4_Gowick(ax, orthogroup, orthogroup_fasta_file, targetp_dict, at_colour_dict):
    arabidopsis_genes = get_species_genes(orthogroup_fasta_file, 'AT', '.', 1)

    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'AT' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 1)] = targetp_dict[key]

    gowick_df = pd.read_csv('data/expression_data/Flaveria_Gowik_Supplemental_Dataset_1-1.csv')
    gowick_df = gowick_df.loc[gowick_df['Locus'].isin(arabidopsis_genes)]

    gray = 0.25
    for i in arabidopsis_genes:
        if i not in at_colour_dict.keys():
            at_colour_dict[i] = (gray, gray, gray, 1)
    rgb_colours = []
    for i, row in gowick_df.iterrows():
        rgb_colours.append(at_colour_dict[row['Locus']])
    gowick_df['colour'] = rgb_colours

    gowick_df['F_pringlei_summed_rpm'] = gowick_df['F_pringlei_rpm'].groupby(gowick_df['colour']).transform('sum')
    gowick_df['F_robusta_summed_rpm'] = gowick_df['F_robusta_rpm'].groupby(gowick_df['colour']).transform('sum')
    gowick_df['F_ramosissima_summed_rpm'] = gowick_df['F_ramosissima_rpm'].groupby(gowick_df['colour']).transform('sum')
    gowick_df['F_trinervia_summed_rpm'] = gowick_df['F_trinervia_rpm'].groupby(gowick_df['colour']).transform('sum')
    gowick_df['F_bidentis_summed_rpm'] = gowick_df['F_bidentis_rpm'].groupby(gowick_df['colour']).transform('sum')

    if plot_highly_expressed:
        gowick_df = gowick_df[(gowick_df['F_pringlei_summed_rpm'] > min_mean_TPM) | (gowick_df['F_robusta_summed_rpm'] > min_mean_TPM) |
                    (gowick_df['F_ramosissima_summed_rpm'] > min_mean_TPM) |
                   (gowick_df['F_trinervia_summed_rpm'] > min_mean_TPM) | (gowick_df['F_bidentis_summed_rpm'] > min_mean_TPM)]
        if len(gowick_df) == 0:
            gowick_df = gowick_df.append(pd.Series(0, index=gowick_df.columns), ignore_index=True)
            gowick_df['Locus'] = 'All genes'
            gowick_df['colour'] = [(1, 1, 1, 1)]
            orthogroup_location_dict['All genes'] = ''

    cumval_pri = 0
    cumval_rob = 0
    cumval_ram = 0
    cumval_tri = 0
    cumval_bid = 0

    sorted_pri = gowick_df.set_index('Locus').sort_values('F_pringlei_summed_rpm')['F_pringlei_rpm']
    sorted_rob = gowick_df.set_index('Locus').sort_values('F_robusta_summed_rpm')['F_robusta_rpm']
    sorted_ram = gowick_df.set_index('Locus').sort_values('F_ramosissima_summed_rpm')['F_ramosissima_rpm']
    sorted_tri = gowick_df.set_index('Locus').sort_values('F_trinervia_summed_rpm')['F_trinervia_rpm']
    sorted_bid = gowick_df.set_index('Locus').sort_values('F_bidentis_summed_rpm')['F_bidentis_rpm']

    colours = gowick_df.set_index('Locus')['colour']

    for name, value in sorted_pri.iteritems():
        ax.bar('F. pringeli\nC3', value, bottom=cumval_pri, color=colours[name])
        cumval_pri += value

    for name, value in sorted_rob.iteritems():
        ax.bar('F. robusta\nC3', value, bottom=cumval_rob, color=colours[name])
        cumval_rob += value

    for name, value in sorted_ram.iteritems():
        ax.bar('F. ramosissima\nC3-C4', value, bottom=cumval_ram, color=colours[name])
        cumval_ram += value

    for name, value in sorted_tri.iteritems():
        ax.bar('F. trinervia\nC4', value, bottom=cumval_tri, color=colours[name])
        cumval_tri += value

    for name, value in sorted_bid.iteritems():
        ax.bar('F. bidentis\nC4', value, bottom=cumval_bid, label=(orthogroup_location_dict[name] + ' ' + name), color=colours[name])
        cumval_bid += value


    ax.set_ylabel('RPM')
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# #### Aubry et al 2011 for Gynandropsis gynandra M BS (reads not mapped by me, assumed to be normalised counts

# In[12]:


def gynandropsis_M_BS_Aubry(ax, orthogroup, orthogroup_fasta_file, targetp_dict, at_colour_dict):
    arabidopsis_genes = get_species_genes(orthogroup_fasta_file, 'AT', '.', 1)

    orthogroup_location_dict = {}
    for key in targetp_dict:
        if 'AT' in key:
            orthogroup_location_dict[get_gene_model(key, '.', 1)] = targetp_dict[key]

    aubry_df = pd.read_csv('data/expression_data/gynandropsis_BS_M_Aubry_2014.csv')
    aubry_df = aubry_df.loc[aubry_df['Accession'].isin(arabidopsis_genes)]

    gray = 0.1
    for i in arabidopsis_genes:
        if i not in at_colour_dict.keys():
            at_colour_dict[i] = (gray, gray, gray, 1.0)
            if gray > 0.90:
                gray = 0.1
            gray = gray + 0.1
    rgb_colours = []
    for i, row in aubry_df.iterrows():
        rgb_colours.append(at_colour_dict[row['Accession']])
    aubry_df['colour'] = rgb_colours

    aubry_df['M_summed_mean'] = aubry_df['Mean M'].groupby(aubry_df['colour']).transform('sum')
    aubry_df['BS_summed_mean'] = aubry_df['Mean BS'].groupby(aubry_df['colour']).transform('sum')

    if plot_highly_expressed:
        aubry_df = aubry_df[(aubry_df['M_summed_mean'] > min_mean_TPM) | (aubry_df['BS_summed_mean'] > min_mean_TPM)]

        if len(aubry_df) == 0:
            aubry_df = aubry_df.append(pd.Series(0, index=aubry_df.columns), ignore_index=True)
            aubry_df['Accession'] = 'All genes'
            aubry_df['colour'] = [(1, 1, 1, 1)]
            orthogroup_location_dict['All genes'] = ''

    cumval_M = 0
    cumval_BS = 0

    sorted_M = aubry_df.set_index('Accession').sort_values('M_summed_mean')['Mean M']
    sorted_BS = aubry_df.set_index('Accession').sort_values('BS_summed_mean')['Mean BS']
    colours = aubry_df.set_index('Accession')['colour']

    for name, value in sorted_M.iteritems():
        ax.bar('M', value, bottom=cumval_M, label=(orthogroup_location_dict[name] + ' ' + name), color=colours[name])
        cumval_M += value

    for name, value in sorted_BS.iteritems():
        ax.bar('BS', value, bottom=cumval_BS, color=colours[name])
        cumval_BS += value

    ax.set_ylabel('Normalised read count')
    ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# ### Overall C3 vs C4 expression in mature leaves

# In[13]:


def C3_v_C4_boxplot(ax, orthogroup_fasta_file):
    arabidopsis_genes = get_species_genes(orthogroup_fasta_file, 'AT', '.', 1)

    #arabidopsis_genes = ['AT5G19140', 'AT5G43830']
    species_info_df = pd.read_csv('data/expression_data/species_info.csv')
    C3_code_list = ((species_info_df.loc[species_info_df['Photosynthesis'] == 'C3'])['Mature_leaf_library_ID']).dropna().to_list()
    C4_code_list = ((species_info_df.loc[species_info_df['Photosynthesis'] == 'C4'])['Mature_leaf_library_ID']).dropna().to_list()
    NADP_code_list = ((species_info_df.loc[species_info_df['Type'] == 'NADP-ME'])['Mature_leaf_library_ID']).dropna().to_list()
    NAD_code_list = ((species_info_df.loc[species_info_df['Type'] == 'NAD-ME'])['Mature_leaf_library_ID']).dropna().to_list()

    expression_df = pd.read_csv('data/expression_data/FINAL_results_t10_no_VYNC_new_RBB_sga_assembly_libsize1.csv', skiprows=3)
    expression_df = expression_df.loc[expression_df['Arabidopsis'].isin(arabidopsis_genes)]
    expression_df_summed = expression_df.sum(axis=0, skipna = True)

    C3_expression_list = []
    C4_expression_list =[]
    NADP_expression_list = []
    NAD_expression_list = []
    for code in C3_code_list:
        C3_expression_list.append(expression_df_summed[code])
    for code in C4_code_list:
        C4_expression_list.append(expression_df_summed[code])
    for code in NADP_code_list:
        NADP_expression_list.append(expression_df_summed[code])
    for code in NAD_code_list:
        NAD_expression_list.append(expression_df_summed[code])

    C3_df = pd.DataFrame(columns=['C3'])
    C4_df = pd.DataFrame(columns=['C4'])
    NADP_df = pd.DataFrame(columns=['NADP-ME'])
    NAD_df = pd.DataFrame(columns=['NAD-ME'])

    C3_df['C3'] = C3_expression_list
    C4_df['C4'] = C4_expression_list
    NADP_df['NADP-ME'] = NADP_expression_list
    NAD_df['NAD-ME'] = NAD_expression_list

    C3_df = C3_df.apply(pd.to_numeric, errors='coerce')
    C4_df = C4_df.apply(pd.to_numeric, errors='coerce')
    NADP_df = NADP_df.apply(pd.to_numeric, errors='coerce')
    NAD_df = NAD_df.apply(pd.to_numeric, errors='coerce')

    ax.set_ylabel('TPM')

    result_df = pd.concat([C3_df,C4_df, NADP_df, NAD_df], axis=1)
    result_df.boxplot()


# #### General functions

# In[14]:


def get_species_genes(fasta_file, species_string, delimiter, delimiter_position):
    genes = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if species_string in record.id:
            gene_list = ((record.id).split(delimiter))
            gene = delimiter.join(gene_list[:delimiter_position])
            genes.append(gene)
    return(genes)


def get_rice_orthologues(orthogroup, species):
    for file in glob.glob(f'data/orthologues/Orthologues_Osativa_323_v7.0.protein_primaryTranscriptOnly/*{species}*'):
        df = pd.read_csv(file, sep='\t')
        df = df[df['Orthogroup'] == orthogroup]
        colourmap = plt.get_cmap('Dark2')
        df['colour'] = [colourmap(k) for k in range(len(df))]

        global colour_counter
        colour_counter = len(df)

        return df


# In[16]:


def maize_orthologues(orthogroup, species):
    for file in glob.glob(f'data/orthologues/Orthologues_Zmays_493_RefGen_V4.protein_primaryTranscriptOnly/*{species}*'):
        df = pd.read_csv(file, sep='\t')
        df = df[df['Orthogroup'] == orthogroup]
        return(df)



# In[17]:


def orthologue_colour_dict(orthologue_df):
    os_dict = {}
    zm_dict = {}

    for i, row in orthologue_df.iterrows():
        os_genes = (row['Osativa_323_v7.0.protein_primaryTranscriptOnly']).split(', ')
        for os_gene in os_genes:
            os_dict[os_gene.split('.')[0]] = row['colour']

        zm_genes = (row['Zmays_493_RefGen_V4.protein_primaryTranscriptOnly']).split(', ')
        for zm_gene in zm_genes:
            zm_dict[zm_gene.split('_')[0]] = row['colour']

    return(os_dict, zm_dict)


# In[18]:


def make_targetp_dict(orthogroup_targetP_files):
    targetp_df = pd.read_csv(orthogroup_targetP_files, delimiter = '\t', skiprows=1)
    targetp_dict = pd.Series(targetp_df.Prediction.values, index=targetp_df['# ID']).to_dict()
    return(targetp_dict)


# In[19]:


def get_species_colour_dict(species_string, zm_colour_dict, orthogroup, delimiter, delimiter_position):
    #delimiter of species of interest (not maize)
    maize_species_orthologues = maize_orthologues(orthogroup, species_string)
    species_colour_dict = {}
    global colour_counter
    colourmap = plt.get_cmap('Dark2')
    for i, row in maize_species_orthologues.iterrows():
        zm_genes = row['Zmays_493_RefGen_V4.protein_primaryTranscriptOnly'].split(', ')
        for col in maize_species_orthologues.columns:
            if species_string in col:
                species_genes = (row[f'{col}']).split(', ')

        for zm_gene in zm_genes:
            if zm_gene.split('_')[0] not in zm_colour_dict:
                zm_colour_dict[zm_gene.split('_')[0]] = colourmap(colour_counter)
                colour_counter += 1

        if len(zm_genes) > 1:
            test_colours = []
            for zm_gene in zm_genes:
                colour = zm_colour_dict[zm_gene.split('_')[0]]
                test_colours.append(colour)
            if len(set(test_colours)) == 1:
                for gene in species_genes:
                    species_colour_dict[get_gene_model(gene, delimiter, delimiter_position)] = test_colours[0]
            else:
                for gene in species_genes:
                    species_colour_dict[get_gene_model(gene, delimiter, delimiter_position)] = colourmap(colour_counter)
                colour_counter += 1
        else:
            colour = zm_colour_dict[zm_genes[0].split('_')[0]]
            for gene in species_genes:
                species_colour_dict[get_gene_model(gene, delimiter, delimiter_position)] = colour


    return(species_colour_dict)


# In[20]:


def get_gene_model(gene, delimiter, delimiter_position):
    gene = (gene.split(delimiter))
    gene = delimiter.join(gene[:delimiter_position])
    return(gene)



# ### Make panel figure

# #### Read in orthogroup data

# In[26]:


def panel_fig(orthogroups, orthogroup_fasta_files, orthogroup_targetP_files):
    for i, value in enumerate(orthogroups):
    #i = 0
        global colour_counter
        rice_v_maize_orthologues = get_rice_orthologues(orthogroups[i], 'Zmays')
        os_colour_dict, zm_colour_dict = orthologue_colour_dict(rice_v_maize_orthologues)
        sb_colour_dict = get_species_colour_dict('Sbicolor', zm_colour_dict, orthogroups[i], '.', 2)
        si_colour_dict = get_species_colour_dict('Sitalica', zm_colour_dict, orthogroups[i], '.', 2)
        pv_colour_dict = get_species_colour_dict('Pvirgatum', zm_colour_dict, orthogroups[i], '.', 2)

        at_colour_dict = get_species_colour_dict('Athaliana', zm_colour_dict, orthogroups[i], '.', 1)
        targetp_dict = make_targetp_dict(orthogroup_targetP_files[i])

        plt.rcParams.update({'font.size': 10})
        fig = plt.figure(constrained_layout=False, figsize=(20, 12))
        spec = gridspec.GridSpec(ncols=7, nrows=7, figure=fig)

        f_ax1_1 = fig.add_subplot(spec[0, 0])
        f_ax1_3 = fig.add_subplot(spec[0,2], sharey=f_ax1_1)
        f_ax1_4 = fig.add_subplot(spec[0,3], sharey=f_ax1_1)
        f_ax1_5 = fig.add_subplot(spec[0,4], sharey=f_ax1_1)
        f_ax1_6 = fig.add_subplot(spec[0,5], sharey=f_ax1_1)

        f_ax2_1 = fig.add_subplot(spec[1, 0], sharey=f_ax1_1)
        f_ax2_3 = fig.add_subplot(spec[1, 2], sharey=f_ax1_1)
        f_ax2_5 = fig.add_subplot(spec[1, 4], sharey=f_ax1_1)

        f_ax3_1 = fig.add_subplot(spec[2, 0], sharey=f_ax1_1)
        f_ax3_3 = fig.add_subplot(spec[2, 2], sharey=f_ax1_1)
        f_ax3_5 = fig.add_subplot(spec[2, 4], sharey=f_ax1_1)
        f_ax3_6 = fig.add_subplot(spec[2, 5], sharey=f_ax1_1)

        f_ax4_1 = fig.add_subplot(spec[3,0:2])
        f_ax4_4 = fig.add_subplot(spec[3,3])

        f_ax5_1 = fig.add_subplot(spec[4,0:2])


        maize_M_BS_chang(f_ax1_1, orthogroups[i], orthogroup_fasta_files[i], zm_colour_dict, targetp_dict)
        plt.setp([f_ax1_1], title='Maize (Chang et al 2012)')

        maize_M_BS_tausta(f_ax1_3, f_ax1_4, orthogroups[i], orthogroup_fasta_files[i], zm_colour_dict, targetp_dict)

        maize_M_BS_denton(f_ax1_5, f_ax1_6, orthogroups[i], orthogroup_fasta_files[i], zm_colour_dict, targetp_dict)

        sbicolor_M_BS_Oxford(f_ax2_1, orthogroups[i], orthogroup_fasta_files[i], sb_colour_dict, targetp_dict)
        plt.setp([f_ax2_1], title='Sorghum (Oxford 2016)')

        setaria_M_BS_john(f_ax2_3, orthogroups[i], orthogroup_fasta_files[i], si_colour_dict, targetp_dict)
        plt.setp([f_ax2_3], title='Setaria (John et al 2014)')

        pvirgatum_M_BS_rao(f_ax2_5, orthogroups[i], orthogroup_fasta_files[i], pv_colour_dict, targetp_dict)
        plt.setp([f_ax2_5], title='Panicum (Rao et al 2016)')

        rice_dev_vancampen(f_ax3_1, orthogroups[i], orthogroup_fasta_files[i], os_colour_dict, targetp_dict)
        plt.setp([f_ax3_1], title='Rice leaf development (Van Campen et al 2016)')

        arabidopsis_dev_woo(f_ax3_3, orthogroups[i], orthogroup_fasta_files[i], targetp_dict, at_colour_dict)
        plt.setp([f_ax3_3], title='Arabidopsis leaf development (Woo et al 2016)')

        maize_foliar_husk_wang(f_ax3_5, f_ax3_6, orthogroups[i], orthogroup_fasta_files[i], zm_colour_dict, targetp_dict)

        flaveria_C3_C4_Gowick(f_ax4_1, orthogroups[i], orthogroup_fasta_files[i], targetp_dict, at_colour_dict)
        plt.setp([f_ax4_1], title='Flaveria leaf expression (Gowick et al 2011)')

        gynandropsis_M_BS_Aubry(f_ax4_4, orthogroups[i], orthogroup_fasta_files[i], targetp_dict, at_colour_dict)
        plt.setp([f_ax4_4], title='Gynandropsis gynandra (Aubry et al 2014)')

        C3_v_C4_boxplot(f_ax5_1, orthogroup_fasta_files[i])
        plt.setp([f_ax5_1], title='Expression across selected C3 and C4 dicots (TPMs summed for each arabidopsis gene in orthogroup)')

        for a in fig.axes:
            a.tick_params(
            axis='y',           # changes apply to the x-axis
            which='both',       # both major and minor ticks are affected
            labelleft=True)

        plt.setp(f_ax1_4.get_yticklabels(), visible=False)
        plt.setp(f_ax1_6.get_yticklabels(), visible=False)
        plt.setp(f_ax3_6.get_yticklabels(), visible=False)

        fig.suptitle(f'{orthogroups[i]}', fontsize=24, x=0.15, y=1.88)
        fig.tight_layout()
        subplots_adjust(left=0.125, bottom=-0.5, right=1.1, top=1.8, wspace=0.2, hspace=0.3)

        print(f"Writing panel figure for {orthogroups[i]}")
        plt.savefig(f'figures/{orthogroups[i]}.pdf', bbox_inches='tight')


# #### Phytozome 13 data

# In[27]:
if len(sys.argv) == 1:
    print('Please specify a orthogroup to be analysed (or a file with a list of orthogroups) as a command line variable')
    print('e.g.')
    print('python Panel_figure_C4_gene_expression.py orthogroup_list_example.txt')
    print('')
    exit()

orthogroups = []
if 'OG0' in sys.argv[1]:
    orthogroups.append(sys.argv[1])
else:
    orthogroup_list = sys.argv[1]
    with open(orthogroup_list) as orthogroup_file:
        for line in orthogroup_file:
            orthogroups.append(line.strip())

global plot_highly_expressed
plot_highly_expressed = False
if 'large' in sys.argv[2]:
    plot_highly_expressed = True
    global min_mean_TPM
    min_mean_TPM = 10000


orthogroup_fasta_files = []
orthogroup_targetP_files = []

for group in orthogroups:
    for file in glob.glob(f'data/orthogroup_sequences/*{group}*'):
        orthogroup_fasta_files.append(file)
    for file in glob.glob(f'data/targetp_2_results/*{group}*'):
        orthogroup_targetP_files.append(file)

print('orthogroup_files being used:')
for group, fasta, targetP in zip(orthogroups, orthogroup_fasta_files, orthogroup_targetP_files):
    print(group, fasta, targetP)

panel_fig(orthogroups, orthogroup_fasta_files, orthogroup_targetP_files)
