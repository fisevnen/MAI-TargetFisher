#导入数据
import pandas as pd # Data analysis
import numpy as np # Scientific computing
import matplotlib.pyplot as plt # Plotting
import matplotlib.colors as colors # Coloring
pd_df = pd.read_csv('results_3columns_hit_weights.csv')
# pd_df = pd.read_csv('results_3columns_hit.csv')
result = pd.DataFrame()
result['x'] = pd_df['logFC']
result['y'] = pd_df['N-log10P']
result['z'] = pd_df['Hit']

x_threshold=1
y_threshold=2

#分组为up, normal, down
result['group'] = 'lightgrey'
result.loc[(result.x > x_threshold)&(result.y > y_threshold),'group'] = 'lightcoral' #x=-+x_threshold直接截断
result.loc[(result.x < -x_threshold)&(result.y > y_threshold),'group'] = 'cornflowerblue' #x=-+x_threshold直接截断
result.loc[result.z == True, 'group'] = 'green'
result.loc[result.y < y_threshold,'group'] = 'lightgrey' #阈值以下点为灰色
# result.loc[abs(result.x) < x_threshold,'group'] = 'lightgrey' #阈值以下点为灰色
result = result.sort_values('z')

print(result.head())

#确定坐标轴显示范围
xmin=-18
xmax=18
ymin=-0.2
ymax=7

#绘制散点图
# fig = plt.figure(figsize=plt.figaspect(7/6)) #确定fig比例（h/w）
fig = plt.figure(figsize=(7, 6), dpi=600)
# plt.tick_params(labelsize=16)
plt.xlim((xmin, xmax)) 
plt.ylim((ymin, ymax))
plt.scatter(result['x'], result['y'], s=6, c=result['group'])

#水平和竖直线
plt.vlines(-x_threshold, ymin, ymax, color='dimgrey',linestyle='dashdot', linewidth=1.5)
plt.vlines(x_threshold, ymin, ymax, color='dimgrey',linestyle='dashdot', linewidth=1.5)
plt.hlines(y_threshold, xmin, xmax, color='dimgrey',linestyle='dashdot', linewidth=1.5)
plt.xticks(range(-16,20,4), fontproperties='Helvetica', size=15)
plt.yticks(range(0,8,2), fontproperties='Helvetica', size=15) 
plt.xlabel('log$_2$(Fold Change)', fontdict={'family':'Helvetica', 'size':16}, fontweight='bold')
plt.ylabel('-log$_{10}$(P-value)', fontdict={'family':'Helvetica', 'size':16}, fontweight='bold')
foo_fig = plt.gcf()
foo_fig.savefig('volcano_FC2_30000_weights.eps', format='eps', dpi=600)
plt.show()
