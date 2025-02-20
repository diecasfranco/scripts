#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
 
import numpy as np
import seaborn
import matplotlib.pyplot as plt
 
 
def plot_heatmap(data, xlabels, ylabels, output_filename):
    seaborn.set(color_codes=True)
    plt.figure(1, figsize=(9, 6))
 
    plt.title("Heatmap")
 
    seaborn.set(font_scale=1.4)
    ax = seaborn.heatmap(data, cmap="YlGnBu", cbar_kws={'label': 'Scale'})
 
    ax.set_xticklabels(xlabels, rotation=45)
    ax.set_yticklabels(ylabels, rotation=45)
 
    ax.set(ylabel="X-Label", xlabel="Y-Label")
 
    plt.savefig(output_filename, bbox_inches='tight', dpi=300)
    plt.close()
 
 
# define data
data = np.random.rand(20, 30)
 
# define labels
x_labels = ["x_{}".format(x) for x in range(1, 31)]
y_labels = ["y_{}".format(x) for x in range(1, 21)]
 
# create heat map
plot_heatmap(data, x_labels, y_labels, "heatmap.png")