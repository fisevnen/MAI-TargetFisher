import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import roc_curve, auc
import glob
from scipy import interp
import numpy as np
import os
import random
from matplotlib.ticker import FuncFormatter

def fsave(list_results, file, addn=False):
    if addn == True:
        list_results = [i+'\n' for i in list_results]
    with open(file, 'w') as f:
        f.writelines(list_results)
        f.close()

def fread(file, removen=True, remove1strow=False):
    with open(file, 'r') as f:
        contents = f.readlines()
        f.close()
    if removen == True:
        contents = [i.replace('\n','') for i in contents]
    if remove1strow == True:
        contents.pop(0)
    return contents

def calculate_single(id, method, results):
    with open('BindingDB_human_DATA.csv','r') as f:
        contents = f.readlines()
        f.close()
    targets = []

    for line in contents:
        if id == line.split(',')[0]:
            if float(line.split(',')[2]) < 10000:
                if not ':' in line.split(',')[3][:-1]:
                    targets.append(line.split(',')[3][:-1])
                else:
                    tmp_list = line.split(',')[3][:-1].split(':')
                    for tmp_i in tmp_list:
                        targets.append(tmp_i)
    unique_target_list = list(set(targets))
    unique_target_list = ['P27487', 'P23219', 'P16083', 'O75762', 'P35354', 'P05067', 'Q16678', 'P09884', 'P09917', 'P06276', 'P22303', 'P27338', 'O60341', 'P04798', 'Q13153', 'P21397', 'P11926', 'Q04206', 'P11511', 'P00519', 'P11274', 'P18031', 'P07711', 'P01106', 'P08183', 'P28907', 'P07339', 'P06239', 'Q07817', 'P09960', 'P14679', 'P08684', 'P03372', 'Q9BYT3', 'P17252', 'Q96EB6', 'Q8IXJ6', 'Q9NTG7', 'P00915', 'P00918', 'P07451', 'P22748', 'P35218', 'Q9Y2D0', 'P23280', 'P43166', 'Q16790', 'O43570', 'Q8N1Q1', 'Q9ULX7', 'Q16236', 'P02766']
    if '' in unique_target_list:
        unique_target_list.remove('')
    unique_target_length = len(unique_target_list)
    if method == 'SEA':
        with open('uID-Gene_name.tsv','r') as f:
            uID_Gene_list = f.readlines()
            f.close()
        for line in results:
            if 'HUMAN' not in line:
                results.remove(line)
        count = 0
        x_list = [0]
        y_list = [0]
        appeared_list = []
        Gene_uID_SEA_results = []
        for index in range(len(results)):
            whole_gene_id = results[index].split(',')[1]
            gene_id = whole_gene_id.split('_')[0]
            # print(gene_id)
            for gene_uID in uID_Gene_list:
                if gene_id in gene_uID:
                    Gene_uID_SEA_results.append(gene_uID)
        # print(Gene_uID_SEA_results)
        for index in range(len(Gene_uID_SEA_results)):
            for uID in unique_target_list:
                if uID in Gene_uID_SEA_results[index]:
                    if uID not in appeared_list:
                        appeared_list.append(uID)
                        count += 1
                        y_list.append(count/len(unique_target_list))
                        x_list.append(index/len(Gene_uID_SEA_results))
                    elif uID in appeared_list:
                        pass
        x_list.append(1)
        y_list.append(y_list[-1])
    elif method == 'Swiss':
        count = 0
        x_list = [0]
        y_list = [0]
        appeared_list = []
        for index in range(len(results)):
            for uID in unique_target_list:
                if uID in results[index]:
                    if uID not in appeared_list:
                        appeared_list.append(uID)
                        count += 1
                        y_list.append(count/len(unique_target_list))
                        x_list.append(index/len(results))
                    elif uID in appeared_list:
                        pass
        x_list.append(1)
        y_list.append(y_list[-1])
    elif method == 'MAI_30000_docking':
        if id == '6866':
            pass
        else:
            count = 0
            x_list = [0]
            y_list = [0]
            appeared_list = []
            for index in range(len(results)):
                for uID in unique_target_list:
                    if uID in results[index]:
                        if uID not in appeared_list:
                            appeared_list.append(uID)
                            count += 1
                            y_list.append(count/len(unique_target_list))
                            x_list.append(index/len(results))
                        elif uID in appeared_list:
                            pass
            x_list.append(1)
            y_list.append(y_list[-1])
    elif method == 'MAI_30000_params':
        if id =='6866':
            pass
        else:
            file = f'30000_resveratrol_hybrid_results_sorted'
            with open(file, 'r') as f:
                results = f.readlines()
                f.close()
            count=0
            x_list = [0]
            y_list = [0]
            appeared_list = []
            for index in range(len(results)):
                for uID in unique_target_list:
                    if uID in results[index]:
                        if uID not in appeared_list:
                            appeared_list.append(uID)
                            count += 1
                            y_list.append(count/len(unique_target_list))
                            x_list.append(index/len(results))
                        elif uID in appeared_list:
                            pass
            x_list.append(1)
            y_list.append(y_list[-1])
    return [x_list, y_list]

Swiss_results = fread('SwissTargetPrediction.csv', remove1strow=True)
SEA_results = fread('SEA-resutls.csv', remove1strow=True)
MAI_TF_results = fread('MAI_TargetFisher_results.csv', remove1strow=True)
[swiss_x, swiss_y] = calculate_single('23926', 'Swiss', Swiss_results)
[SEA_x, SEA_y] = calculate_single('23926', 'SEA', SEA_results)
[MAI_x, MAI_y] = calculate_single('23926', 'MAI_30000_docking', MAI_TF_results)
[MAI_h_x, MAI_h_y] = calculate_single('23926', 'MAI_30000_params', None)

from matplotlib.ticker import FuncFormatter
plt.figure(figsize = (6.4,4.8), dpi = 300)
plt.xlim([0,1])
plt.ylim([0,1])
plt.rc('font', family='Helvetica')
def to_percent(temp, position):
  return '%1.0f'%(100*temp) + '%'
plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
plt.gca().xaxis.set_major_formatter(FuncFormatter(to_percent))
plt.plot(MAI_h_x, MAI_h_y, label='MAI-TargetFisher_weighted', linewidth=3,  color='#E47C7C')
plt.plot(MAI_x, MAI_y, label='MAI-TargetFisher_unweighted', linewidth=3, color=(0.22, 0.50, 0.57))
plt.plot(SEA_x, SEA_y, label='SEA', color=(0.96, 0.73, 0.47))
plt.plot(swiss_x, swiss_y, label='SwissTarget', linewidth=3, color=(0.56, 0.77, 0.58))
plt.title(f'Average recall rate of test dataset', fontsize='19')
plt.xlabel('Result', fontsize='16')
plt.ylabel('Identified Tatgets', fontsize='16')
plt.xticks(size='12')
plt.yticks(size='12')
font_legend={'family' : 'Helvetica','style':'normal', 'weight':'bold'}
plt.legend(loc='best', markerfirst=False, shadow=True, fancybox=True, prop=font_legend)
plt.savefig(f'ave_all_4_methods_nolegand.eps')