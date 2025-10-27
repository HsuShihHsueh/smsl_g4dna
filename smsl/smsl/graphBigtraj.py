from .multiWindows import TimeAgent
from .config import ConfAgent
from .plotAgent import font_xylabel, font_xyticks

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ipywidgets
from IPython.display import display


plot_prominent_modes_xylim = {
    'bm' : { 'xlim':[-2, 150], 'ylim':[0.33, 1.03] },
    'hb' : { 'xlim':[-5,  65], 'ylim':[0.37, 1.03] },
    'st' : { 'xlim':[-1,  44], 'ylim':[0.37, 1.03] },
    'rb' : { 'xlim':[-2, 130], 'ylim':[0.33, 1.03] },
    'pr0': { 'xlim':[-2,  60], 'ylim':[0.29, 1.03] },
    'pr1': { 'xlim':[-2,  60], 'ylim':[0.29, 1.03] },
    'cl_st':{'xlim':[-2,  30], 'ylim':[0.37, 1.03] },
}

class BigTrajPlotAgent():
    def plot_eigenvalues(self, is_fullsize_plot=False, strandid=None):
        if is_fullsize_plot:
            plt.figure(figsize=(40,4))
        for st_agent in self.t_agent.values(): 
            if strandid is None:
                w = st_agent.w
            else:
                w = st_agent.w_by_strandid[strandid]
            plt.plot(list(range(len(w))), w, label=st_agent.time_label)
        xticks = self.t_agent.mean.get_mode_xticks(is_fullsize_plot, strandid)
        plt.xticks(xticks)
        plt.xlabel('mode')
        plt.ylabel(r'eigenvalue (kcal/$\rmÅ^2$/mol)')
        plt.title(self.m)
        plt.legend(fontsize=7)
        plt.show()
    
    def plot_prominent_modes(self, strandid=None, alpha_threshold=0.8, sele_mode=None, save_filename=None):
        x = self.t_agent.mean.w
        y = self.mean_r_alpha_array
        fig = plt.figure(figsize=[6.3, 3.95])
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharex=ax0)  # Share the x-axis with ax0
        flierprops = dict(marker="d", markerfacecolor='orange')
        ax0.boxplot(x, whis=1.5, vert=False, flierprops=flierprops, widths=0.7)
        hival = np.quantile(x, 0.75)+1.5*(np.quantile(x, 0.75)-np.quantile(x, 0.25)) # Q3 + 1.5*IQR
        # whishi = min( max(x[x<=hival]), hival ) # whishi = min( max(data_without_outlier), Q3+1.5*IQR )
        ax0.axvline(hival, color="red", linestyle="--")
        ax0.get_xaxis().set_visible(False)
        ax0.get_yaxis().set_visible(False)
        ax1.scatter(x, y, c="dimgray", s=12)
        ax1.set_xlabel(r'$\lambda_{\alpha}^{\rm '+self.m_abbr+r'}\ {\rm(kcal mol/Å^2)}$', fontdict=font_xylabel)
        ax1.set_ylabel(r'$\left< r_{\alpha}^{\rm '+self.m_abbr+r'}\right>$', fontdict=font_xylabel)
        ax1.axvline(hival, color="red", linestyle="--")
        ax1.axhline(alpha_threshold, color="green", linestyle="--")
        ax1.text(0.95, 0.05, self.system, transform=ax1.transAxes, fontsize=20, ha='right', va='bottom')
        x_outlier = x[np.logical_and(x>hival, y>alpha_threshold)]
        y_outlier = y[np.logical_and(x>hival, y>alpha_threshold)]
        ax1.scatter(x_outlier, y_outlier, c="orange",  s=12)
        if not sele_mode is None:
            ax1.scatter(x[sele_mode-1], y[sele_mode-1], c="red", s=12)
        ax1.set_xlim(plot_prominent_modes_xylim[self.m_abbr]['xlim'])
        
        xticklabels = [label.get_text() for label in ax1.get_xticklabels()]
        ax1.set_xticklabels(xticklabels, fontdict=font_xyticks) 
        ax1.set_ylim(plot_prominent_modes_xylim[self.m_abbr]['ylim'])
        yticklabels = [label.get_text() for label in ax1.get_yticklabels()]
        ax1.set_yticklabels(yticklabels, fontdict=font_xyticks) 
        fig.tight_layout()
        plt.rcParams['svg.fonttype'] = 'none'
        if not save_filename is None:
            fig.savefig(save_filename, transparent=True)
        return fig, np.array([ax0, ax1])

    def plot_prominent_modes_and_identify_atom_with_interactive(self, strandid=None, sign=False):
        def plot_prominent_modes_and_identify_atom(sele_mode):
            if strandid != None:
                mode_count_by_strandid = modes_in_strandid[sele_mode-1]
                print(f'mode: {mode_count_by_strandid}')
                self.plot_prominent_modes(sele_mode=sele_mode, strandid=strandid)
                self.t_agent.mean.plot_identify_atom(mode=mode_count_by_strandid, _strandid=strandid, sign=sign)
            else:
                self.plot_prominent_modes(sele_mode=sele_mode)
                self.t_agent.mean.plot_identify_atom(mode=sele_mode, sign=sign)
        if strandid != None:        
            modes_in_strandid = np.array(range(1, len(self.t_agent.mean.eigv_strandid)+1))[self.t_agent.mean.eigv_strandid==strandid]
        ipywidgets.interact(
            plot_prominent_modes_and_identify_atom,
            sele_mode = ipywidgets.IntSlider(min=1, max=self.t_agent.mean.n_node, step=1, value=1)
        )

    def plot_prominent_modes_and_identify_atom_by_atomname_with_interactive(self, strandid=None):
        def plot_prominent_modes_and_identify_atom_by_atomname(sele_mode):
            if strandid != None:
                mode_count_by_strandid = modes_in_strandid[sele_mode-1]
                print(f'mode: {mode_count_by_strandid}')
                self.plot_prominent_modes(sele_mode=sele_mode, strandid=strandid)
                self.t_agent.mean.plot_identify_atom_by_atomname(mode=mode_count_by_strandid, _strandid=strandid)
            else:
                self.plot_prominent_modes(sele_mode=sele_mode)
                self.t_agent.mean.plot_identify_atom_by_atomname(mode=sele_mode)
        if strandid != None:        
            modes_in_strandid = np.array(range(1, len(self.t_agent.mean.eigv_strandid)+1))[self.t_agent.mean.eigv_strandid==strandid]
        ipywidgets.interact(
            plot_prominent_modes_and_identify_atom_by_atomname,
            sele_mode = ipywidgets.IntSlider(min=1, max=self.t_agent.mean.n_node, step=1, value=1)
        )

    def plot_prominent_modes_and_get_identify_pair_with_interactive(self, strandid=None, vlv_threshold=0.1, pe_threshold=0.5, filter=True):
        def plot_prominent_modes_and_get_identify_pair(sele_mode):
            if strandid != None:
                mode_count_by_strandid = modes_in_strandid[sele_mode-1]
                self.plot_prominent_modes(sele_mode=sele_mode, strandid=strandid)
                plt.show()
                identify_pair = self.t_agent.mean.GetIdentifyPair(mode_count_by_strandid, vlv_threshold=vlv_threshold, pe_threshold=pe_threshold, filter=filter)
                display(identify_pair)
                eigenvalue = self.t_agent.mean.get_eigenvalue_by_id(mode_count_by_strandid)
                print(f'origin mode: {mode_count_by_strandid}')
                print(f'λ = {eigenvalue:.1f}')
            else:
                self.plot_prominent_modes(sele_mode=sele_mode)
                plt.show()
                identify_pair = self.t_agent.mean.GetIdentifyPair(sele_mode, vlv_threshold=vlv_threshold, pe_threshold=pe_threshold, filter=filter)
                display(identify_pair)
                eigenvalue = self.t_agent.mean.get_eigenvalue_by_id(sele_mode)
                print(f'λ = {eigenvalue:.1f}')
        if strandid != None:        
            modes_in_strandid = np.array(range(1, len(self.t_agent.mean.eigv_strandid)+1))[self.t_agent.mean.eigv_strandid==strandid]
        ipywidgets.interact(
            plot_prominent_modes_and_get_identify_pair,
            sele_mode = ipywidgets.IntSlider(min=1, max=self.t_agent.mean.n_node, step=1, value=1)
        )
        
class BigTrajGraphAgent(ConfAgent, BigTrajPlotAgent):
    def __init__(self, MAgent):
        ConfAgent.__init__(self)
        self.MAgent = MAgent
        self.t_agent = TimeAgent(self.MAgent, create_mean=True)
        self.n_window = len(self.t_agent)
        self.m        = self.t_agent.mean.m
        self.m_abbr   = self.t_agent.mean.m_abbr
    
    def spectral_decomposition(self):
        for _, st_agent in self.t_agent.items():
            st_agent.spectral_decomposition()
        self.t_agent.mean.spectral_decomposition()
        self.n_eigenvalues = self.t_agent.mean.n_node

    def get_mean_modes_v_mat(self):
        return self.t_agent.mean.v
    
    def get_window_modes_v_mat(self, window_id):
        return list(self.t_agent.values())[window_id].v
    
    def get_r_n_alpha(self):
        mean_modes_v_mat_T = self.get_mean_modes_v_mat().T
        r_n_alpha_mat = np.zeros((self.n_window, self.n_eigenvalues))
        for window_id in range(self.n_window):
            window_modes_v_mat = self.get_window_modes_v_mat(window_id)
            product_mat = np.abs(np.dot(mean_modes_v_mat_T, window_modes_v_mat))
            for eigv_id in range(self.n_eigenvalues):
                r_n_alpha_mat[window_id, eigv_id] = product_mat[eigv_id,:].max()
        return r_n_alpha_mat

    def set_mean_r_alpha_array(self):
        mean_r_alpha_array = np.zeros(self.n_eigenvalues)
        r_n_alpha_mat = self.get_r_n_alpha()
        for eigv_idx in range(self.n_eigenvalues):
            mean_r_alpha_array[eigv_idx] = r_n_alpha_mat[:, eigv_idx].mean()
        self.mean_r_alpha_array = mean_r_alpha_array
