import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import numpy as np


palette_tab10 = sns.color_palette("tab10", n_colors=10)
palette_dark2 = sns.color_palette("Dark2", n_colors=10)
palette_dark2_sele = palette_dark2[0:1] + palette_dark2[5:6]
palette_css = sns.color_palette(['gold', 'firebrick', 'darkorchid', 'darkcyan', 'olive', 'peru', 'seagreen', 'orangered'])
palette_css2 = sns.color_palette(['deepskyblue', 'deeppink', 'mediumslateblue', 'darkviolet','slategray', 'saddlebrown', 'crimson', 'chocolate', 'teal', 'forestgreen'])

palette = sns.color_palette(palette_tab10+palette_dark2_sele+palette_css+palette_css2)


font_xylabel = {
'family': 'Arial', 
'weight': 'normal', 
'size'  : 12
}
font_xyticks = {
'family': 'Arial', 
'weight': 'normal', 
'size'  : 10
}
font_label = {
'family': 'Arial', 
'weight': 'normal', 
'size'  :10
}

class PlotAgent():
    def __init__(self, nrows=1, ncols=1, figsize_1=[3, 1.8], dpi=200, sharex=True, sharey=True, hspace=0, wspace=0, **kwargs):
        figsize = [figsize_1[0]*ncols, figsize_1[1]*nrows]
        self.fig, self.axs = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex, sharey=sharey, dpi=dpi, **kwargs)
        self.fig.subplots_adjust(hspace=hspace, wspace=wspace)
        plt.rcParams['svg.fonttype'] = 'none'
    def axhlines(self, ys, color='gray', linestyle='--', linewidth=0.8, alpha=0.6, zorder=1):
        for ax in self.axs:
            for y in ys:
                ax.axhline(y=y, color=color, linestyle=linestyle, linewidth=linewidth, alpha=alpha, zorder=zorder)
    def set_xticks(self, x_ticks=None, x_padding=False, fontdict=font_xyticks, rotation=0, padding_scale=1):
        for ax in self.axs.flatten(): 
            if not x_ticks is None:
                ax.set_xticks(x_ticks)
                if x_padding:
                    x_max = max(x_ticks) + (max(x_ticks)-min(x_ticks)) / 18 * padding_scale
                    x_min = min(x_ticks) - (max(x_ticks)-min(x_ticks)) / 18 * padding_scale
                else:
                    x_max = max(x_ticks)
                    x_min = min(x_ticks)
                ax.set_xlim(x_min, x_max)
            xticklabels = [label.get_text() for label in ax.get_xticklabels()]
            if len(ax.get_xticklabels())!=0:
                ax.set_xticklabels(xticklabels, fontdict=fontdict, rotation=rotation)
    def set_yticks(self, y_ticks=None, y_padding=False, fontdict=font_xyticks, padding_scale=1):
        for ax in self.axs.flatten(): 
            if not y_ticks is None:
                ax.set_yticks(y_ticks)
                if y_padding:
                    y_max = max(y_ticks) + (max(y_ticks)-min(y_ticks)) / 18 * padding_scale
                    y_min = min(y_ticks) - (max(y_ticks)-min(y_ticks)) / 18 * padding_scale
                else:
                    y_max = max(y_ticks)
                    y_min = min(y_ticks)
                ax.set_ylim(y_min, y_max)
            yticklabels = [label.get_text() for label in ax.get_yticklabels()]
            if len(ax.get_yticklabels())!=0:
                ax.set_yticklabels(yticklabels, fontdict=fontdict) 
    def set_xlabel_down(self, label, fontdict=font_xylabel):
        self.axs[-1].set_xlabel(label, fontdict=fontdict)
    def set_xlabels_down(self, labels, fontdict=font_xylabel):
        for axs, label in zip(self.axs[-1], labels):
            axs.set_xlabel(label, fontdict=fontdict)
    def set_ylabels(self, tor_agent, labels, fontdict=font_xylabel, position='left', **kwags):
        if position=='left':
            axs = self.axs if len(self.axs.shape)==1 else self.axs[:, 0]
            for ax, (sytem, tor_ss_agent) in zip(axs, tor_agent.items()):
                ax.set_ylabel(getattr(tor_ss_agent, labels), fontdict=fontdict, **kwags)
        elif position=='right':
            axs = self.axs if len(self.axs.shape)==1 else self.axs[:, -1]
            for ax, (sytem, tor_ss_agent) in zip(axs, tor_agent.items()):
                ax.set_ylabel('')
                ax2 = ax.twinx()
                ax2.set_ylabel(getattr(tor_ss_agent, labels), fontdict=fontdict, **kwags)
                ax2.set_yticks([])
    def set_supylabel(self, label, x=0.5, y=0.902, fontdict=font_xylabel, **kwargs):
        self.fig.suptitle(label, fontsize=fontdict['size'], fontweight=fontdict['weight'], fontfamily=fontdict['family'], x=x, y=y, **kwargs)
    def set_legend(self, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=1, frameon=False, prop=font_label, linewidth=None, **kwargs):
        legend = np.array(self.axs).flatten()[0].legend(loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol, frameon=frameon, prop=prop, **kwargs)
        if not legend is None:
            for handle in legend.get_lines():
                handle.set_linewidth(linewidth)
    def savefig(self, filename, transparent=True, folder='pic'):
        self.fig.savefig(f'{folder}/{filename}', transparent=transparent)
        
    
class SinglePlotAgent():
    def __init__(self, nrows=1, ncols=1, figsize_1=[6.0, 3.0], dpi=200, **kwargs):
        figsize = [figsize_1[0]*ncols, figsize_1[1]*nrows]
        self.fig, self.axs = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi, **kwargs)
        plt.rcParams['svg.fonttype'] = 'none'
    def set_xticks(self, x_ticks=None, x_padding=False, fontdict=font_xyticks, rotation=0, padding_scale=1): 
        ax = self.axs
        if not x_ticks is None:
            ax.set_xticks(x_ticks)
            if x_padding:
                x_max = max(x_ticks) + (max(x_ticks)-min(x_ticks)) / 18 * padding_scale
                x_min = min(x_ticks) - (max(x_ticks)-min(x_ticks)) / 18 * padding_scale
            else:
                x_max = max(x_ticks)
                x_min = min(x_ticks)
            ax.set_xlim(x_min, x_max)
        xticklabels = [label.get_text() for label in ax.get_xticklabels()]
        if len(ax.get_xticklabels())!=0:
            ax.set_xticklabels(xticklabels, fontdict=fontdict, rotation=rotation)
    def set_yticks(self, y_ticks=None, y_padding=False, fontdict=font_xyticks, padding_scale=1):
        ax = self.axs
        if not y_ticks is None:
            ax.set_yticks(y_ticks)
            if y_padding:
                y_max = max(y_ticks) + (max(y_ticks)-min(y_ticks)) / 18 * padding_scale
                y_min = min(y_ticks) - (max(y_ticks)-min(y_ticks)) / 18 * padding_scale
            else:
                y_max = max(y_ticks)
                y_min = min(y_ticks)
            ax.set_ylim(y_min, y_max)
        yticklabels = [label.get_text() for label in ax.get_yticklabels()]
        if len(ax.get_yticklabels())!=0:
            ax.set_yticklabels(yticklabels, fontdict=fontdict) 
    def set_xlabel(self, label, fontdict=font_xylabel):
        self.axs.set_xlabel(label, fontdict=fontdict)
    def set_ylabel(self, label, position='left',fontdict=font_xylabel):
        if position=='left':
            self.axs.set_ylabel(label, fontdict=fontdict)
        elif position=='right':
            ax2 = self.axs.twinx()
            ax2.set_ylabel(label, fontdict=fontdict)
            ax2.set_yticks([])
    def set_legend(self, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=1, frameon=False, prop=font_label, linewidth=None, **kwargs):
        legend = self.axs.legend(loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol, frameon=frameon, prop=prop, **kwargs)
        if not legend is None and not linewidth is None:
            for handle in legend.get_lines():
                handle.set_linewidth(linewidth)
    def savefig(self, filename, transparent=True, folder='pic'):
        self.fig.savefig(f'{folder}/{filename}', transparent=transparent)