import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import dim
from bokeh.io import show

cmp = sns.set_palette('Dark2')
sns.set(context = "talk", font_scale=0.9, style = 'white', palette = cmp)

# Plot with holoviews
def boxplot_metrics(data, title, color_by = 'Method', palette = 'Set3', invert_xaxis = True):
    data = data.melt(var_name='Variables', value_name='Metric Value')
    data['Method'] = data['Variables'].apply(lambda x: x.split('-')[0])
    data['Ranking'] = data['Variables'].apply(lambda x: x.split('-')[1])
    data['Metric'] = data['Variables'].apply(lambda x: x.split('-')[2])
    boxwhisker = hv.BoxWhisker(data, ['Metric', 'Ranking', 'Method' ], 
                               'Metric Value', label = title)
    boxwhisker.opts(show_legend = False, width=900, height = 500, 
                    box_fill_color = color_by, cmap=palette, ylim=(0, 1),
                    xrotation=45, invert_xaxis = invert_xaxis, toolbar='above')
    f = hv.render(boxwhisker, backend='bokeh')
    f.toolbar.logo = None
    f.grid.grid_line_color = 'lightgrey'
    f.grid.grid_line_dash = [5, 3]
    return f

def plot_swarm_metrics(data, title, filter_regex, hue = None, hue_name = None, ylim = (0.3, 1),
                      cmap = 'Dark2', legend = True, **kwargs):
    data = data.filter(regex = filter_regex, axis = 1)
    n_cols = len(data.columns)
    data = data.melt(var_name='Var', value_name = 'Values')
    data['Method'] = data['Var'].apply(lambda x: x.split('-')[0])
    if hue is not None:
        data[hue_name] = np.tile(hue, n_cols)
    ax = sns.swarmplot(data = data, x = 'Method', y = 'Values', s = 5,  
                  hue = hue_name, palette = cmap, **kwargs)
    if not legend: ax.legend_.remove()
    plt.title(title)
    plt.grid(linestyle='--', linewidth='0.8')
    plt.ylim(ylim)
    
def plot_swarm_pair(df, title, col = ['DkS', 'DkLEff'], metric='ROC', **kwargs):
    plt.figure(figsize=(17,8))
    plt.subplots_adjust(wspace=0.2)
    metric = '-' + metric
    plt.subplot(1,2,1)
    plot_swarm_metrics(data = df.filter(regex=col[0], axis = 1), 
                  title = F'{title}: {col[0]}{metric}', filter_regex = metric, **kwargs);
    plt.subplot(1,2,2)
    plot_swarm_metrics(data = df.filter(regex=col[1], axis=1), 
                  title = F'{title}: {col[1]}{metric}', filter_regex = metric, **kwargs);