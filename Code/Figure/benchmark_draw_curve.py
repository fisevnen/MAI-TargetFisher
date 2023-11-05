import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd
from sklearn.metrics import roc_curve, auc
import glob
from scipy import interp
import numpy as np
import os
from matplotlib.ticker import FuncFormatter

all_sampled_list = ['50536679', '50355499', '7533', '50519158', '50355501', '6866', '50536675', '50463479', '50048803', '50028421', '50008735', '31093', '31090', '21079', '16673', '139540', '11639']

def calculate_single(BMID_list, method, out_path):
    x_dict = {} # for drawing total_recall x_dict
    y_dict = {} # for drawing total_recall x_dict
    for id in BMID_list:
        with open('../../BindingDB_human_DATA.csv','r') as f:
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
        if '' in unique_target_list:
            unique_target_list.remove('')
        unique_target_length = len(unique_target_list)
        if method == 'SEA':
            file = f'SEA/SEA_{id}.xls/sea-results.xls'
            with open(file, 'r') as f:
                results = f.readlines()
                f.close()
            results.pop(0)
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
            file = f'swiss-target/{id}.csv'
            with open(file, 'r') as f:
                results = f.readlines()
                f.close()
            results.pop(0)
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
            file = f'../../6new_results_zips/30000_{id}_results/MAI_TargetFisher_results.csv'
            with open(file, 'r') as f:
                results = f.readlines()
                f.close()
            results.pop(0)
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
            file = f'../../30000_{id}_hybrid_results_sorted'
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

        elif method == 'MAI_5000_docking':
            file = f'../{id}_results/MAI_TargetFisher_results.csv'
            with open(file, 'r') as f:
                results = f.readlines()
                f.close()
            results.pop(0)
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

           
        plt.figure(figsize = (6.4,4.8), dpi = 300)
        plt.xlim([0,1])
        plt.ylim([0,1])
        def to_percent(temp, position):
          return '%1.0f'%(100*temp) + '%'
        plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
        plt.gca().xaxis.set_major_formatter(FuncFormatter(to_percent))
        plt.plot(x_list, y_list, label=f'{method}')
        plt.title(f'Recall Curve of BMID:{id}', fontsize='14')
        plt.xlabel('Result', fontsize='12')
        plt.ylabel('Identified Tatgets(%)', fontsize='12')
        plt.legend(loc='upper left', markerfirst=False)
        plt.savefig(f'{out_path}/{method}-{id}.png')
        plt.close()
        x_dict[id] = x_list
        y_dict[id] = y_list

        EF_list = [1, 5, 10, 20, 30, 50, 100]
        for EF in EF_list:
            EF_unique_target_list = unique_target_list
            found_uID_list = []
            EF_count = 0
            for index in range(int(len(results) * 0.01 * EF )):
                for uID in EF_unique_target_list:
                    if uID in results[index]:
                        if not uID in found_uID_list:
                            EF_count += 1
                            found_uID_list.append(uID)
                        break

            try:
                EF_rate = (EF_count / int(len(results) * 0.01 * EF )) / (unique_target_length / len(results))
                with open(f'{out_path}/{method}-{id}_EF.info', 'a') as f:
                    f.write(f'EF{EF}%: {EF_rate}\n')
                    f.write(f'EF{EF}: {EF_count}\n')
                    f.close()
            except:
                pass

            top_list = [10, 20, 30, 50, 100, 500, 1000]
            for top in top_list:
                top_unique_target_list = unique_target_list
                found_uID_list = []
                if top > len(results):
                    pass
                else:
                    top_count = 0
                    for index in range(top):
                        # print(uID)
                        for uID in top_unique_target_list:
                            if uID in results[index]:
                                if not uID in found_uID_list:
                                    top_count += 1
                                    found_uID_list.append(uID)
                                break
                    top_rate = top_count / unique_target_length
                    with open(f'{out_path}/{method}-{id}_top.info', 'a') as f:
                        f.write(f'top{top}%: {top_rate}\n')
                        f.write(f'top{top}: {top_count}\n')
                        f.close()

                x_list.append(0)
                y_list.append(y_list[-1])
    return [x_dict, y_dict]
def process_method_ave_xy(x_dict, y_dict):
    import copy
    all_x_point = []
    all_y_point = []
    for key in x_dict.keys():
        x_points = copy.deepcopy(x_dict[key])
        y_points = copy.deepcopy(y_dict[key])
        all_x_point.append(x_points)
        all_y_point.append(y_points)
    all_x_point = list(set(sum(all_x_point, [])))
    all_x_point = list(set([format(i, '.2f') for i in all_x_point]))
    all_x_point.sort()
    print(len(all_x_point))
    all_x_point = [float(i) for i in all_x_point]
    import copy
    new_all_y_point_dict = {}
    for key in x_dict.keys():
        ori_x_point = copy.deepcopy(x_dict[key])
        ori_y_point = copy.deepcopy(y_dict[key])
        ori_x_point = [format(i, '.2f') for i in ori_x_point]
        ori_x_point = [float(i) for i in ori_x_point]
        set_ori_x_point = list(set([format(i, '.2f') for i in ori_x_point]))
        set_ori_x_point = [float(i) for i in set_ori_x_point]
        ori_y_point = [format(i, '.2f') for i in ori_y_point]
        ori_y_point = [float(i) for i in ori_y_point]
        ori_x_point.sort()
        set_ori_x_point.sort()
        ori_y_point.sort()
        # print(len(ori_x_point), len(set_ori_x_point), len(ori_y_point))
        # print(ori_x_point, set_ori_x_point, ori_y_point)
        remove_count = 0
        for item_index in range(len(set_ori_x_point)):
            if ori_x_point.count(set_ori_x_point[item_index]) > 1:
                item_count = ori_x_point.count(set_ori_x_point[item_index])
                # print(item_count)
                for remove_index in range(item_count-1):
                    ori_x_point.pop(remove_count)
                    ori_y_point.pop(remove_count)
                remove_count += 1
            else:
                remove_count += 1
        tmp = ori_y_point
        for point in all_x_point:
            if point in ori_x_point:
                pass
            else:
                for index in range(len(ori_x_point)-1):
                    if ori_x_point[index] < point < ori_x_point[index+1]:
                        tmp.append(ori_y_point[index] + (ori_y_point[index+1] - ori_y_point[index]) * ((point - ori_x_point[index])/(ori_x_point[index+1] - ori_x_point[index])))
                        break
        if len(tmp) > len(all_x_point):
            print(key)
        new_all_y_point_dict[key] = tmp
    draw_new_all_y_point_dict = {}
    for key in new_all_y_point_dict.keys():
        y_list = new_all_y_point_dict[key]
        y_list.sort()
        draw_new_all_y_point_dict[key] = y_list
    sum_all_y_point = []
    for key in draw_new_all_y_point_dict.keys():
        y_list = draw_new_all_y_point_dict[key]
        y_list.sort()
        sum_all_y_point.append(y_list)
    avg_all_y_point = []
    for i in range(len(all_x_point)):
        tmp = []
        for index in range(len(sum_all_y_point)):
            tmp.append(sum_all_y_point[index][i])
        avg_all_y_point.append(sum(tmp)/len(sum_all_y_point))
    
    return[all_x_point, avg_all_y_point]

# os.mkdir('Swiss_statistics')
[Swiss_x_dict, Swiss_y_dict] = calculate_single(all_sampled_list, 'Swiss', 'Swiss_statistics')
# os.mkdir('SEA_statistics')
[SEA_x_dict, SEA_y_dict] = calculate_single(all_sampled_list, 'SEA', 'SEA_statistics')
# os.mkdir('MAI_30000_statistics')
[MAI_30000_x_dict, MAI_30000_y_dict] = calculate_single(all_sampled_list, 'MAI_30000_docking', 'MAI_30000_statistics')
# os.mkdir('MAI_5000_statistics')
[MAI_5000_x_dict, MAI_5000_y_dict] = calculate_single(all_sampled_list, 'MAI_5000_docking', 'MAI_5000_statistics')
# os.mkdir('MAI_30000_weighted_statistics')
[MAI_30000_weighted_x_dict, MAI_30000_weighted_y_dict] = calculate_single(all_sampled_list, 'MAI_30000_params', 'MAI_30000_weighted_statistics')

[Swiss_x_draw,Swiss_y_draw] =  process_method_ave_xy(Swiss_x_dict, Swiss_y_dict)
[SEA_x_draw,SEA_y_draw] =  process_method_ave_xy(SEA_x_dict, SEA_y_dict)
[MAI_30000_x_draw,MAI_30000_y_draw] =  process_method_ave_xy(MAI_30000_x_dict, MAI_30000_y_dict)
[MAI_5000_x_draw,MAI_5000_y_draw] =  process_method_ave_xy(MAI_5000_x_dict, MAI_5000_y_dict)
[MAI_30000_weighted_x_draw, MAI_30000_weighted_y_draw] = process_method_ave_xy(MAI_30000_weighted_x_dict, MAI_30000_weighted_y_dict)

from matplotlib.ticker import FuncFormatter
plt.figure(figsize = (7,6), dpi = 600)
plt.xlim([0,1])
plt.ylim([0,1])
plt.rc('font', family='Helvetica')
def to_percent(temp, position):
  return '%1.0f'%(100*temp) + '%'
plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
plt.gca().xaxis.set_major_formatter(FuncFormatter(to_percent))
plt.plot(MAI_30000_weighted_x_draw, MAI_30000_weighted_y_draw, linewidth=3, label='MAI-TargetFisher_weighted', color='#E47C7C')
plt.plot(MAI_30000_x_draw, MAI_30000_y_draw, label='MAI-TargetFisher_unweighted', linewidth=3, color=(0.22, 0.50, 0.57))
# plt.plot(MAI_5000_x_draw, MAI_5000_y_draw, label='MAI-TargetFisher_5k', color='black')
plt.plot(SEA_x_draw, SEA_y_draw, label='SEA', linewidth=3, color=(0.96, 0.73, 0.47))
plt.plot(Swiss_x_draw, Swiss_y_draw, label='SwissTarget', linewidth=3, color=(0.56, 0.77, 0.58))
plt.title(f'Average recall rate of test dataset', fontsize='19')
plt.xlabel('Result', fontsize='16')
plt.ylabel('Identified Tatgets', fontsize='16')
plt.xticks(size='12')
plt.yticks(size='12')
font={'family' : 'Helvetica','style':'normal', 'weight':'bold', 'size':14}
plt.legend(loc='lower right', markerfirst=False, shadow=True, fancybox=True, prop=font)
# plt.legend.set_fontsize(12)
plt.savefig(f'results/ave_all_4_methods.eps')
plt.show()
