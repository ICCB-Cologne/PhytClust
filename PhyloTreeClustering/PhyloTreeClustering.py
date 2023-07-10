import math
import os
import sys
from copy import deepcopy, copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from PhyloTreeClustering.calculate_pdm import calculate_pdm_from_tree
from PhyloTreeClustering.calinski_harabasz import calinski_harabasz_score
from PhyloTreeClustering.medicc_functions import plot_tree, import_tree

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors)

class PhyloTreeClustering:
    """
    This class creates a dendrogram and clusters the samples on it.

        -   tree_file_path: tree you want to cluster, in new format

        -   output_path: path to save figure produced by algorithm

        -   tree: directly tree to cluster

        -   dist_matrix: panda dataframe, matrix of the distances between all leaves of the tree, index are the names of the leaves 
                The default None means the distance matrix will be calculated by the function calculate_pdm_from_tree.

        -   threshold: manually choose the cutting threshold (not recommended)

        -   calinski_harabasz_score: if True, the algorithm chooses the best threshold based on the clustering with the highest score

        -   n_diff_th: number of repetition of scoring, meaning the number of different threshold cuts will be tested.
                The default None means this number will be calculated according to the size of the tree in init_n_diff_th.
            
        -   CH_matrix: distance matrix that will be used to evaluate the clusering with the CH score
                If no matrix is given, the normal dist_matrix will be used
        
        -   min_size_clus: minimum size of a cluster. 
                If no value is given, it is calculated in init_min_size_clus according to the size of the tree
        
        -   multi_level_clus: if True, make a second level of clusters.
                The default is False, so only one level of clustering will be performed
        
        -   dist_type: distance calculation used to compute the new dendrogram distances.
                The default and recommended distance is 'half_max'. The others are 'mean', 'median' and 'max'

    """

    def __init__(self, tree_file_path= None, output_path = None, tree = None, dist_matrix = None, threshold = None, calinski_harabasz_score = True,
                     n_diff_th = None, CH_matrix = None, min_size_clus = None, multi_level_clus = False, dist_type = 'half_max'):
        self.init_tree( tree_file_path, tree, dist_matrix)
        self.th = self.init_th(threshold)
        self.dendro_tree = deepcopy(self.tree)
        self.dist_matrix_df = self.init_dist_matrix(tree_file_path, dist_matrix)
        self.dist_matrix = self.dist_matrix_df.copy().values
        self.N = len(self.dist_matrix_df)
        self.n_diff_th = self.init_n_diff_th(calinski_harabasz_score , n_diff_th)
        self.name_to_matrix_index = self.init_name_to_index()
        self.name_to_clade = self.init_name_to_clade()
        self.parent_dict_names = self.init_parent()
        self.child_parent_dist = [] # order of tree postorder [child_dist, parent_dist]
        self.min_size_clus = self.init_min_size_clus(min_size_clus)
        self.multi_level_clus = multi_level_clus
        self.CH_matrix = self.init_CH_matrix(CH_matrix, dist_matrix, tree_file_path, threshold)
        self.first_CH_scores = None
        self.dist_type = dist_type
        self.output_path = output_path
        
    def init_tree(self, tree_file_path, tree, dist_matrix):
        '''
        This function initializes the tree depending on where the information comes from.
        '''
        if tree_file_path is None and tree is None:
            raise ValueError('The input data is not sufficient, you must give the tree_file_path ' 
                'or directly the tree ')
        
        if tree is not None : self.tree = deepcopy(tree)
        elif tree_file_path: self.compute_dist_matrix(tree_file_path, only_tree = True) # assigns medicc tree to self.tree

    def init_n_diff_th(self, calinski_harabasz_score, n_diff_th):
        '''
        This function initializes the number of thresholds that will be tested in the th_selection_function.
        '''
        if calinski_harabasz_score :
            if n_diff_th == None or n_diff_th == 0: 
                if self.N > 36 : return round(6*np.sqrt(self.N)) 
                else: return self.N
            elif n_diff_th <= self.N : 
                return n_diff_th
            else: 
                raise ValueError('The number of calinski_harabasz_score repetition ' 
                'must be smaller or equal to the number of internal nodes')
        elif not calinski_harabasz_score: 
            return None

    def init_dist_matrix(self, tree_file_path, dist_matrix):
        '''
        This function initializes the distance matrix, by default the relative distance matrix is calculated.
        '''
        if dist_matrix is None :
            self.pdm_matrix = self.compute_dist_matrix(tree_file_path)
            relative_dist_matrix = copy(self.pdm_matrix)
            leave_mean_dists = self.pdm_matrix.mean()
            relative_dist_matrix = relative_dist_matrix.divide(leave_mean_dists, axis= 'columns')
            relative_dist_matrix = relative_dist_matrix.divide(leave_mean_dists, axis= 'index')
            matrix = relative_dist_matrix
        else: matrix = dist_matrix
        if 'diploid' in matrix.index : 
            matrix.drop(index= 'diploid', columns= 'diploid', inplace=True)
        return matrix

    def init_CH_matrix(self, matrix, dist_matrix, tree_file_path, threshold):
        '''
        This function initializes the distance matrix used to compute the CH score with the cluster labels.
        '''
        if threshold is not None : return None
        if matrix is None: # no CH matrix is given, so we calculate it
            if dist_matrix is not None: 
                self.pdm_matrix = self.compute_dist_matrix(tree_file_path) # verify that it will not be calulated twice
            matrix = self.pdm_matrix
        if 'diploid' in matrix.index : 
            matrix.drop(index= 'diploid', columns= 'diploid', inplace=True)
        return matrix

    def init_name_to_index(self):
        '''
        This function initializes the dictionary that keeps the information of which sample name has which index.
        '''
        self.leave_names = [name for name in list(self.dist_matrix_df.index) if name != 'diploid']
        return {name : index for index, name in enumerate(self.leave_names)}

    def init_name_to_clade(self):
        '''
        This function initializes the dictionary that enables going from the name to the clade in the medicc tree.
        '''
        name_to_clade = {}
        for clade in self.tree.find_clades():
            if clade.name:
                if clade.name in name_to_clade:
                    raise ValueError("Duplicate key: %s" % clade.name)
                name_to_clade[clade.name] = clade
        return name_to_clade

    def init_parent(self):
        '''
        This function initializes the dictionary that enables to find the parent of any sample with its name.
        '''
        parent_dict_names = dict()
        for node in self.tree.get_nonterminals():
            parent_dict_names[node.clades[0].name] = node.name
            parent_dict_names[node.clades[1].name] = node.name
        return parent_dict_names

    def init_min_size_clus(self, min_size_clus):
        '''
        This function initializes the minimun size a cluster must have to be considered as such.
        '''
        if min_size_clus == None:
            min_size_clus = round(np.sqrt(self.N)/2)  
        elif min_size_clus >= self.N:
            raise ValueError("The minimun cluster size entered is bigger than the number of samples")
        return min_size_clus

    def init_th(self, threshold):
        '''
        This function initializes the threshold if one is manually given, if the CH score is not used.
        '''
        if threshold is None:
            return 0
        else: return threshold ## add verification that it is no CH

    def compute_dist_matrix(self, tree_file_path = None, only_tree = False):
        '''
        This function computes the distance matrix from the input files.
        '''
        if tree_file_path:
            self.tree = import_tree(tree_file_path, 'diploid')
        cur_sample_labels = [clade.name for clade in self.tree.get_terminals()]
        if only_tree: return 

        return calculate_pdm_from_tree(self.tree, cur_sample_labels)
    
    def calculate_dendro_dists(self):
        '''
        This function calculates the distances needed to build the dendrogram associated with the tree.
        The default distance chosen is to take the mean of the 50 % biggest distances.
        '''
        self.dendro_dists = dict()
        for leave in self.leave_names : # all the leaves have a dendrogram dist = 0 
            self.dendro_dists[leave] = 0

        self.all_names = []
        self.name_to_index = {}
        n=0

        for clade in self.dendro_tree.find_clades(order='postorder'): # order where children always come before parents

            if clade.name is None:
                self.mrca = [child.name for child in clade.clades if child.name != 'diploid'][0]
                continue

            self.all_names.append(clade.name)
            self.name_to_index[clade.name] = n
            n +=1
            if clade.is_terminal() and clade.name != 'diploid':
                pass
            elif clade.name != 'diploid':
                cur_children = [self.name_to_matrix_index[x.name] for x in clade.get_terminals() if x.name != 'diploid']  # subset of distance matrix          
                flat_subset = copy(self.dist_matrix[np.ix_(cur_children, cur_children)]).flatten() 

                if self.dist_type == 'max': # takes max dist between any terminal leave of the new parent
                    self.dendro_dists[clade.name] = np.max(flat_subset)

                if self.dist_type == 'half_max' or self.dist_type is None: # takes the median of 50 % of the biggest dist 
                    sorted_array = np.sort(flat_subset)[::-1]
                    self.dendro_dists[clade.name] = np.median(sorted_array[:round(0.5*len(sorted_array) +0.5)])

                if self.dist_type == 'mean': # takes the mean of all dist between terminal leaves of the new parent
                    self.dendro_dists[clade.name] = np.mean([value for value in flat_subset if value != 0])
                    if math.isnan(self.dendro_dists[clade.name]): 
                        self.dendro_dists[clade.name] = 0

                if self.dist_type == 'median': # takes the median of all dist between terminal leaves of the new parent
                    self.dendro_dists[clade.name] = np.median([value for value in flat_subset if value != 0]) 
                    if math.isnan(self.dendro_dists[clade.name]): 
                        self.dendro_dists[clade.name] = 0
                
        # fix for half_max, median and mean dist: 
        for clade_name in self.all_names:
            parent = self.parent_dict_names[clade_name]
            if parent is None : continue
            
            if self.dendro_dists[clade_name] > self.dendro_dists[parent]: # if the child's distance is larger than the parent's dist
                self.dendro_dists[parent] = self.dendro_dists[clade_name] 

        self.index_to_name = {index : name for name, index in self.name_to_index.items()}


    def adjust_branch_length_tree(self):
        '''
        This function changes the branch length of the medicc tree to the ones needed for the dendrogram.
        The branch length are not just the distances from the dendrogram because they are relative distances
        to the parent, so they need to be computed.
        '''
        list(self.dendro_tree.find_clades(order='postorder'))[-3].branch_length = 0  # MRCA

        for index, clade_name in enumerate(self.all_names):
            parent = self.parent_dict_names[clade_name] 
            if clade_name != 'diploid' and clade_name != self.mrca: # modify branch length on tree
                list(self.dendro_tree.find_clades(order='postorder'))[index].branch_length = float(self.dendro_dists[parent] - self.dendro_dists[clade_name])
              
        self.tot_branch_length_dict = {name : self.dendro_tree.distance(str(name)) for name in self.all_names if name is not None} # dist from root
        self.tot_branch_length_dict[None] = 0 # parent of diploid and internal 1
        
        for clade_name in self.all_names: # after the changes build child_parent_dist
            parent = self.parent_dict_names[clade_name] 
            self.child_parent_dist.append([self.tot_branch_length_dict[clade_name],
                         self.tot_branch_length_dict[parent]])          
        

    def th_adequate_range(self, nb_tests, sub_tot_branch_length_dict = None):
        '''
        This function computes the range of thresholds that must be evaluated by CH score.
        '''
        if sub_tot_branch_length_dict is None:  # 1st layer clustering
            th_array = sorted(list(self.tot_branch_length_dict.values()))[:nb_tests] 
            self.min_th = np.min(th_array)
            self.max_th = np.max(th_array)
        else :   # 2nd layer clustering
            th_array= sorted(list(sub_tot_branch_length_dict.values()))[:nb_tests+1]
            th_array = [th for th in th_array if th != 0]
        return th_array


    def th_CH_selection(self, subclone = False, sub_tot_branch_length_dict = None, child_parent_dist = None, CH_matrix_sub = None):
        '''
        This function selects the best threshold (cutting point) according to the highest Calinski Harabasz score.
        '''

        if subclone : # clustering of 2nd layer clone, the thresholds and leaves are not the same
            th_range = self.th_adequate_range(self.n_diff_th, sub_tot_branch_length_dict)
            leave_names = self.sub_leave_names
        else: # clustering of 1st layer
            th_range = self.th_adequate_range(self.n_diff_th)
            leave_names =self.leave_names

        if CH_matrix_sub is not None: CH_matrix = CH_matrix_sub
        else: CH_matrix = self.CH_matrix

        self.CH_scores = {}
        self.all_CH_labels = {}
        for i in range(self.n_diff_th):
            dist_th = th_range[i]
            self.cluster(temp_th = dist_th, child_parent_dist= child_parent_dist) # computes color_dict

            labels = [self.color_dict[name][1:] for name in leave_names]

            int_array = [int(label) for label in labels if label != 'lack'] # find max label apart from outliers
            if len(int_array) != 0 : max_label = np.max(int_array)
            else : max_label = 0

            for index, label in enumerate(labels):
                if label == 'lack': # outliers black[1:] --> change cluster so that they are computed as outliers in CH score
                    max_label += 1
                    labels[index] = str(max_label)
            
            if not 1 < len(np.unique(labels)) < self.N: # the score can not be computed 
                score = 0
            else :
                score = calinski_harabasz_score(CH_matrix, labels)
            self.CH_scores[dist_th] = score
            self.all_CH_labels[dist_th] = labels
        self.best_th = max(self.CH_scores.items(), key=lambda x: x[1])[0] # final threshold chosen
        self.best_CH = round(max(self.CH_scores.items(), key=lambda x: x[1])[1],3) # score of the labels chosen
        self.labels_ = self.all_CH_labels[self.best_th] # final labels chosen

        return self.best_th

    def cluster(self, temp_th = None, child_parent_dist = None):
        '''
        This function clusters the samples and computes the color dictionary according to the chosen threshold.
        '''
        if temp_th is not None : self.th = temp_th # if the threshold is given in argument
        self.cluster_clades = [] # index of samples that are clustering clades

        if child_parent_dist is None : # first layer of clustering
            for index, (dist_child, dist_parent) in enumerate(self.child_parent_dist): 
                if dist_child >= self.th and dist_parent < self.th: # sort samples that are the clade starting a cluster
                    self.cluster_clades.append(index)
                n=1
                leave_names = self.leave_names
        else: # 2nd layer of clustering if indicated
            for index, (dist_child, dist_parent) in zip(self.indexes, child_parent_dist): 
                if dist_child >= self.th and dist_parent < self.th:
                    self.cluster_clades.append(index)
                n= self.start_color # to avoid repetition of colors between subclones
                leave_names = self.sub_leave_names
        
        self.color_dict = {}
        self.leaves_added =[]
        self.clus_number = 0
        for clade in self.cluster_clades:
            color = 'C'+ str(n)
            clade_name = self.index_to_name[clade]
            self.color_dict[clade_name] = color
            n +=1
            clades_clusters = list(self.name_to_clade[self.index_to_name[clade]].find_clades())
            if len(clades_clusters) > self.min_size_clus: self.clus_number +=1
            for child in clades_clusters:
                if len(clades_clusters) > self.min_size_clus: # if the cluster is big enough to be considered one
                    self.color_dict[child.name] = color
                    self.leaves_added.append(child.name)
                else:
                    self.color_dict[clade_name] = 'black' # modify color of outliers
                    self.color_dict[child.name] = 'black'
                    self.leaves_added.append(child.name)
            

        # add samples to color_dict that are not in the clusters (outliers):
        for child in leave_names:
            if child not in self.leaves_added:
                self.color_dict[child] = 'black'

    def multi_level_cut(self):
        '''
        This function computes second layer of clones if this option has been chosen. The big clones from the 1st layer 
        go through the same algorithm again, as their own tree, and receive a new cutting threshold each.
        '''
        ## keep first layer of clustering in memory to be able to use it after:
        self.first_CH_scores = copy(self.CH_scores)  
        self.first_best_th = copy(self.best_th)  
        self.first_clus_nb = copy(self.clus_number) 
        self.first_best_CH = copy(self.best_CH)
        self.first_color_dict = self.rearrange_colors(copy(self.color_dict))
        self.second_color_dict = copy(self.color_dict)

        self.multi_level_clus = True
        self.min_size_sub_clus =  round(self.N/6)  # cluster we want to devide must be at least 1/6th of total size
        self.sub_best_th_dict = {}
        self.subclones_clus_nb = 0
        self.start_color = np.max([int(label) for label in self.labels_])

        for clade in self.cluster_clades:
            current_clade = self.name_to_clade[self.index_to_name[clade]]
            clades_clusters = list(current_clade.find_clades())
            if len(clades_clusters) > self.min_size_sub_clus:
                current_leaves = current_clade.get_terminals()
                current_leaves_name = [clade.name for clade in current_leaves]
                sub_tot_branch_length_dict = {clade.name : self.tot_branch_length_dict[clade.name] for clade in clades_clusters}
                parent_name = self.parent_dict_names[self.index_to_name[clade]] 
                sub_tot_branch_length_dict[parent_name] = 0 # parent 
                sub_all_names = []
                for child in current_clade.find_clades(order='postorder'):
                    sub_all_names.append(child.name)
                
                self.indexes =[]
                for index, name in enumerate(self.all_names):
                    if name in sub_all_names:
                        self.indexes.append(index)

                child_parent_dist =[]
                for clade_name in sub_all_names:
                    parent = self.parent_dict_names[clade_name] 
                    child_parent_dist.append([sub_tot_branch_length_dict[clade_name],
                         sub_tot_branch_length_dict[parent]])

                self.n_diff_th = round(0.5*len(clades_clusters)) 
                CH_matrix_sub = self.CH_matrix.loc[current_leaves_name,current_leaves_name] # subset of dist matrix for current clone
                self.sub_leave_names = list(CH_matrix_sub.index)  # modify original not ideal
                self.N = len(self.sub_leave_names) # modify original not ideal
                
                # look for best threshold for this subclone with CH score:
                sub_best_th = self.th_CH_selection(subclone= True, sub_tot_branch_length_dict= sub_tot_branch_length_dict, 
                                            child_parent_dist= child_parent_dist, CH_matrix_sub= CH_matrix_sub)

                self.sub_best_th_dict[self.index_to_name[clade]] = sub_best_th
                self.cluster(temp_th= sub_best_th, child_parent_dist= child_parent_dist) # compute color dict

                # keep track of number of additional clusters made by cutting some clones
                self.subclones_clus_nb += len([color for color in np.unique(list(self.color_dict.values())) if color != 'black']) -1 

                # adjusting the color dictionary by adding the colors for the newly formed clusters 
                for sample_name in sub_all_names:
                    if sample_name in list(self.color_dict.keys()):
                        self.second_color_dict[sample_name] = self.color_dict[sample_name]

                # verify that parents of cluster clades are black, and correct it if not :
                for clade in self.cluster_clades:
                    parent = self.parent_dict_names[self.index_to_name[clade]]
                    if parent not in self.second_color_dict.keys() : continue 
                    if self.second_color_dict[parent] != 'black': 
                        self.second_color_dict[parent] = 'black'
                                
                self.start_color = np.max([int(label) for label in self.labels_]) # update starting color for next subclone

        # compute calinski score of 2nd layer clustering
        self.multi_level_labels = [self.second_color_dict[name][1:] for name in self.leave_names]
        int_array = [int(label) for label in self.multi_level_labels if label != 'lack']
        if len(int_array) != 0 : max_label = np.max(int_array)
        else : max_label =0

        for index, label in enumerate(self.multi_level_labels):
            if label == 'lack': # outliers black[1:] --> change cluster so that they are computed as outliers in CH score
                max_label +=1
                self.multi_level_labels[index] = str(max_label)

        self.multi_level_score = round(calinski_harabasz_score(self.CH_matrix, self.multi_level_labels),3)


    def run(self):
        '''
        This function is the one used to run the whole algorithm.
        '''
        self.calculate_dendro_dists()
        self.adjust_branch_length_tree()
        if self.n_diff_th is not None and self.th is not None: # threshold tuning
            self.th = self.th_CH_selection()
        self.cluster()
        if self.multi_level_clus is True : self.multi_level_cut()

    def rearrange_colors(self, colors_dict):
        '''
        This function rearranges the colors in the color dictionary so that there is no repetition.
        '''
        colors = colors_dict.values()
        modified_colors = {}
        modified_colors['black'] = 'black'
        counter = 0

        for color in colors:
            if color not in modified_colors:
                modified_colors[color] = 'C' + str(counter)
                counter += 1

        result = [modified_colors[color] for color in colors]
        return {name : color for name, color in zip(colors_dict.keys(), result)}

    def label_colors_color_dict(self,label):
        '''
        This function colors the tree leaves according to the clustering provided.
        '''
        if label in self.color_dict:
            return self.color_dict[label]
        elif 'diploid' in label or 'internal' in label:
            return 'black'
        else:
            return 'black'

    def plot_dendro(self, show_multi_level_clus = False, ax = None, annotation = True, **kwargs):
        '''
        This function plots the dendrogram and annotates the cutting points on it, 
        either 1st layer or 2nd layer clustering (show_multi_level_clus = True).
        '''
        if ax is None:
            fig, ax = plt.subplots(figsize = (10,10))

        if not show_multi_level_clus:
            if self.multi_level_clus : 
                self.color_dict = self.first_color_dict
                self.clus_number = self.first_clus_nb
                self.best_th = self.first_best_th
            else: self.color_dict = self.rearrange_colors(self.color_dict)

            plot_tree(self.dendro_tree, label_colors= self.label_colors_color_dict, label_func=lambda x: '', show_branch_lengths=False, 
                title=f'cluster_number = {self.clus_number}\nth = {round(self.th,3)}\nth range = [{round(self.min_th,3)}, {round(self.max_th,3)}]', ax= ax,
                **kwargs)
            
                    
        if show_multi_level_clus:
            if not self.multi_level_clus: 
                raise ValueError('The multi level clustering parameter is set to False so no plotting of it is possible')
            self.color_dict = self.rearrange_colors(self.second_color_dict)
            plot_tree(self.dendro_tree, label_colors= self.label_colors_color_dict, label_func=lambda x: '', show_branch_lengths=False, 
            title=f'cluster_number = {self.first_clus_nb + self.subclones_clus_nb}\nnumber subclone cuts = {len(self.sub_best_th_dict)}', ax =ax,
            **kwargs)

        if annotation :
            ax.axvline(self.best_th, linestyle = '--', color = 'orange', linewidth=2.5)
            ax.axvline(self.min_th, linestyle = '--', color = 'r',linewidth=2.5)
            ax.axvline(self.max_th, linestyle = '--', color = 'r',linewidth=2.5)
            if show_multi_level_clus:
                for clade_name, th in self.sub_best_th_dict.items():
                    ax.axvline(th, linestyle = '--', color = self.first_color_dict[clade_name], linewidth=2.5)
        return ax


    def plot_th_score(self, ax = None, output_path = None ):
        '''
        This function plots the CH score for each threshold tested, the orange dot was the chosen threshold,
        the different colored dots are the thresholds of the 2nd layer clustering.
        '''
        if ax is None:
            fig, ax = plt.subplots(figsize = (5,5))
        elif len(ax) != 1 :
            raise ValueError('The plotting of the score per threshold needs only 1 ax')

        if self.first_CH_scores is not None: 
            self.CH_scores = self.first_CH_scores
            self.best_th = self.first_best_th
        colors = ['black']*len(self.CH_scores)
        keys_list = list(self.CH_scores.keys())
        index_best_th = keys_list.index(self.best_th)
        colors[index_best_th] = 'orange'
        
        if self.multi_level_clus:
            for clade_name, th in self.sub_best_th_dict.items():
                index_th = keys_list.index(th)
                colors[index_th] = self.first_color_dict[clade_name]

        ax.scatter(list(self.CH_scores.keys()), list(self.CH_scores.values()), c= colors, marker='o')
        ax.set_xlabel('Distance threshold')
        ax.set_ylabel('Calinski Harabasz score')

        plt.xticks(rotation=40)
        plt.title('Threshold tuning', size = 15)

        if output_path is not None : self.output_path = output_path

        if self.output_path:
            plt.savefig(os.path.join(self.output_path, 'th_score_plot.png'))
        plt.show()
        return ax

    def plot_summary(self, axs = None, annotation = True, output_path = None, width_scale =1, height_scale = 1 ):
        '''
        This function plots a summary of the algorithm, with the dendrogram next to the medicc tree.
        '''
        plot_height = 5 + height_scale * np.sqrt(10*self.N) * 0.1
        max_leaf_to_root_distances = np.max([np.sum([x.branch_length for x in self.tree.get_path(leaf)])
                            for leaf in self.tree.get_terminals()])
        plot_width = 7 + np.max([0, width_scale * np.log10(max_leaf_to_root_distances / 100) * 5])
       
        if self.multi_level_clus:
            if axs is None: 
                fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(4*min(250, plot_width), min(250, plot_height)))
            elif len(axs) != 4 : 
                raise ValueError('The axs parameter must be of length 4')
            color_dict = self.first_color_dict
        else: 
            if axs is None: 
                fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(2*min(250, plot_width), min(250, plot_height)))
            elif len(axs) != 2 : 
                raise ValueError('The axs parameter must be of length 2')
            self.first_best_CH = self.best_CH
        
        self.plot_dendro(ax = axs[0], annotation= annotation)
        plot_tree(self.tree, label_colors= self.label_colors_color_dict,label_func=lambda x: '',
                show_branch_lengths=False, ax = axs[1],
                title=f'Calinski Harabasz score = {self.first_best_CH}')

        if self.multi_level_clus :        
            self.plot_dendro( ax = axs[2], show_multi_level_clus= True, annotation = annotation)
            plot_tree(self.tree, label_colors= self.label_colors_color_dict,label_func=lambda x: '',
                    show_branch_lengths=False, ax = axs[3],
                    title=f'Calinski Harabasz score = {self.multi_level_score}') 
        #plt.tight_layout()

        if output_path is not None : self.output_path = output_path

        if self.output_path:
            plt.savefig(os.path.join(self.output_path, 'summary_plot.png'))
        #plt.show()
        return axs

    def sample_labels(self, output_path = None, only_terminal = True ):

        if self.multi_level_clus:
            color_dict = self.first_color_dict
        else: 
            color_dict = self.color_dict

        self.labels_= []
        names = []
        for name, label in color_dict.items():
            if only_terminal and 'internal' in name: continue
            if 'C' in label:
                names.append(name)
                self.labels_.append(int(label[1:]))
            else: 
                names.append(name)
                self.labels_.append(-1)

        #sample_labels = [[key, value] for key, value in color_dict.items()]
        sample_labels = np.vstack((names, self.labels_ )).T
        if self.output_path:
            np.savetxt(os.path.join(self.output_path, 'sample_labels.txt'), sample_labels, fmt= '%s', delimiter= ':')

        if self.multi_level_clus:
            color_dict_2 = self.color_dict
            self.labels_2= []
            names_2 = []
            for name, label in color_dict_2.items():
                if only_terminal and 'internal' in name: continue
                if 'C' in label:
                    names_2.append(name)
                    self.labels_2.append(int(label[1:]))
                else: 
                    names_2.append(name)
                    self.labels_2.append(-1)

            sample_labels_2 = np.vstack((names_2, self.labels_2 )).T
            if self.output_path:
                np.savetxt(os.path.join(self.output_path, 'sample_labels_2nd_layer.txt'), sample_labels_2, fmt= '%s', delimiter= ':')

        return sample_labels
