import os, sys
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

draw_list = []
for pair in [0.18916743, 0.35146541]:
    weight1, weight2 = pair[0], pair[1]
    BMID = 'real_cpd19_10'
    with open('MAI_TargetFisher_results_019_-10.csv', 'r') as f:
        MAI_TF_contents = f.readlines()
        f.close()
    with open(f'{BMID}_ShapeScreen_BDB_results.csv','r') as g:
        shape_contents = g.readlines()
        g.close()
    
    with open(f'../cell_line/HEK293_thd1_uID_gene.tsv', 'r') as f:
        HEK_results = f.readlines()
        f.close()
    
    with open('../FC12P001.result', 'r') as f:
        true_values = f.readlines()
        f.close()
    true_values = [i[:-1] for i in true_values]

    new_contents = pd.DataFrame(data=None, columns=['UniprotID','nDocking_score','similarity', 'hybrid_score', 'true'])

    all_docking_score = []
    for TF_line in MAI_TF_contents:
        if 'uniprotID' not in TF_line:
            all_docking_score.append(float(TF_line.split(',')[2]))
    nomalized_index = max(all_docking_score) - min(all_docking_score)
    higher_score = max(all_docking_score)

    shape_score = {}
    for shape_line in shape_contents:
        shape_line = shape_line[:-1]
        if 'score' not in shape_line:
            uniprotID, score = shape_line.split(',')[-2], float(shape_line.split(',')[-1])
            if uniprotID not in shape_score.keys():
                shape_score[uniprotID] = [score]
            else:
                shape_score[uniprotID].append(score)

    shape_maximum_score = {}
    for key in shape_score.keys():
        shape_maximum_score[key] = max(shape_score[key])

    docking_score = {}
    for TF_line in MAI_TF_contents:
        TF_line = TF_line[:-1]
        if 'score' not in TF_line:
            uniprotID, score = TF_line.split(',')[1], float(TF_line.split(',')[2])
            if uniprotID not in docking_score.keys():
                docking_score[uniprotID] = [score]
            else:
                docking_score[uniprotID].append(score)
    docking_minimum_score = {}
    for key in docking_score:
        docking_minimum_score[key] = min(docking_score[key])

    for uID in docking_minimum_score.keys():
        docking_score = docking_minimum_score[uID]
        normalized_docking_score = (higher_score - docking_score)/nomalized_index
        if uID in shape_maximum_score:
            max_similarity = shape_maximum_score[uID]
        else:
            max_similarity = 0
        if uID in true_values:
            true = 1
        else:
            true = 0
        final_score = round(weight2 * max_similarity + weight1 * normalized_docking_score, 3)
        new_contents.loc[len(new_contents.index)] = [f'{uID}', f'{normalized_docking_score}', f'{max_similarity}', f'{final_score}', f'{true}']

    new_contents = new_contents.sort_values(by='hybrid_score', ascending=False)

    new_contents.to_csv(f'CL2_{BMID}_hybrid_{weight1}_{weight2}_10000_5_FC12_30000_-10.results', index=False)

    with open(f'CL2_real_cpd19_10_hybrid_{weight1}_{weight2}_10000_5_FC12_30000_-10.results', 'r') as f:
        contents = f.readlines()
        MAI_TF_results = [i.replace('\n','') for i in contents]
        f.close()

    MAI_TF_results.pop(0)
    

    with open('019_MS_data_uID', 'r') as f:
        MS_uID = f.readlines()
        f.close()
    
    MS_uID = [i.replace('\n', '') for i in MS_uID]
    MAI_TF_uID_num = 0 #chaji
    for line in MAI_TF_results:
        uID = line.split('.')[1]
        if uID in MS_uID:
            MAI_TF_uID_num += 1

    draw_x = []
    draw_y = []
    hit_num = 0
    for num in range(len(MAI_TF_results)):
        uID  = MAI_TF_results[num].split(',')[0]
        if uID not in MS_uID:
            pass
        else:
            draw_x.append(num/len(MAI_TF_results))
            if uID in true_values:
                hit_num += 1
                draw_y.append(hit_num/len(true_values))
                # print(uID)
            else:
                draw_y.append(hit_num/len(true_values))
    draw_x.append(1)
    draw_y.append(draw_y[-1])
    draw_list.append([draw_x,draw_y])

plt.figure(figsize = (7,6), dpi = 600)
plt.xlim(0,1)
plt.ylim(0,1)
plt.rc('font', family='Helvetica')
def to_percent(temp, position):
    return '%1.0f'%(100*temp) + '%'
plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
plt.gca().xaxis.set_major_formatter(FuncFormatter(to_percent))
plt.plot(draw_list[0][0], draw_list[0][1], label='MAI-TargetFisher_weighted', linewidth=3, color='#E47C7C')
plt.plot(draw_list[1][0], draw_list[1][1], label='MAI-TargetFisher_unweighted', linewidth=3, color=(0.22, 0.50, 0.57))
plt.title(f'Recall rate of compound19', fontsize='19')
plt.xlabel('Result', fontsize='16')
plt.ylabel('Identified Tatgets', fontsize='16')
plt.xticks(size='12')
plt.yticks(size='12')
font={'family' : 'Helvetica','style':'normal', 'weight':'bold', 'size':14}
plt.legend(markerfirst=False, shadow=True, fancybox=True, prop=font)
plt.axhline(draw_y[-1], color='dimgrey',linestyle='dashdot')
plt.savefig(f'CL2_{weight1}_{weight2}_recall_FC12_30000_-10.eps')