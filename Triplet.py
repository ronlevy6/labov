import diff_treshold
import network
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import sys
import numpy as np

class Triplet():
    ''' A class which represents a tripet in the graph.
    A triplet is a structure that represnets 2 nodes in graph that have a common
    neighbour but they are not.
    There are 5 fields - center node, node1, node2, corr1, coor2'''
    
    def __init__(self, center, node1, node2, corr1, corr2, low_corr):
        #corr1 - the correlation between center and node1
        #corr2 - the correlation between center and node2
        #low_corr - the correlation between node1 and node2 - less than treshold
        self.center = center
        self.node1 = node1
        self.node2 = node2
        self.corr1 = corr1
        self.corr2 = corr2
        self.low_corr = low_corr
        self.score = min(corr1,corr2) / low_corr
    
    def __repr__(self):
        msg = str(self.node1) + " <--("+ str(self.corr1)[:7] +")-- " + str(self.center) + " --("+ str(self.corr2)[:7] +")--> " + str(self.node2)
        
        return msg
    
    def __eq__(self, other):
        if not isinstance(other, Triplet):
            return False
        cen = self.center == other.center 
        corrs = self.corr1 == other.corr1 and self.corr2 == other.corr2 and self.low_corr == other.low_corr
        nodes =  self.node1 == other.node1 and self.node2 == other.node2
        
        return (cen and corrs and nodes)
        
    
    def create_all_triplets(graph, trait_values_path,is_abs = True):
        ''' this method gets a graph and creates all the triplets in it.
        The path of trait_values is needed to get the correlation between
        the unconnected nodes.
        trait_values_file is indicated for the ORIGINAL file!!!'''
            
        diff_treshold.print_log("create_all_triplets - start")
            
        graph_triplets = network.find_all_triplets(graph, "corr")
            
        #create a list with all unnconeted traits
        unconnected_traits = [triplet[1] for triplet in graph_triplets]
        unconnected_traits += [triplet[2] for triplet in graph_triplets]
            
        distinct_unconnected_traits = list(set(unconnected_traits))
            
        triplets_lst = [] # init empty list to contain all triplet objects
            
        diff_treshold.print_log("create_all_triplets - after find_all_triplets, before read file")
            
        df = pd.read_excel(trait_values_path, header = 0, index_col = "FlowJo Subject ID")
            
        diff_treshold.print_log("create_all_triplets - after read file before creating good DF")        
            
        filtered_df = df.loc[distinct_unconnected_traits]
            
        corr_df = diff_treshold.create_corr_between_all(filtered_df)
            
        diff_treshold.print_log("create_all_triplets - DF is good, before loop")        
        
        for i in range(0 , len(graph_triplets)):
            if i % 800 == 0:
                diff_treshold.print_log("create_all_triplets - in loop, i = " + str(i))
            
            (center,node1,node2,corr1,corr2) = graph_triplets[i]
            unconnected_corr = corr_df.loc[node1,node2]
            if is_abs:
                unconnected_corr = abs(unconnected_corr)
            curr_triplet = Triplet(center,node1, node2, corr1, corr2, unconnected_corr)
            triplets_lst.append(curr_triplet)
            
        diff_treshold.print_log("create_all_triplets - after loop - finished") 
        return triplets_lst
        
    def get_graph_Triplets_heatmap(graph, trait_values_path, heatmap_save_path, is_save):
        ''' this method gets a graph (usually, subgraph - a conneceted component)
        find its triplets and creates an heatmap whit X axis as Triplet centers,
        Y axis as all nodes in graph. The value in eah place in the matrix is the 
        Triplet's score, or 0 when the nodes don't have Triplet'''
        
        diff_treshold.print_log("Triplet.get_graph_Triplets_heatmap - start, before create triplets")
        triplets_lst = []
        triplets_lst = Triplet.create_all_triplets(graph, trait_values_path)       
        if len(triplets_lst) == 0:
            #no trios..
            print("no trios here..")
            return triplets_lst,None,None
        diff_treshold.print_log("Triplet.get_graph_Triplets_heatmap - after create triplets")
        
        #get all centers and nodes to create a proper DF
        centers = [triplet.center for triplet in triplets_lst]
        centers_no_dups = network.remove_dups(centers)
        all_nodes = list(graph.nodes())
        #create a dataframe of zeros with triplets centers as indexs and nodes as columns
        df = pd.DataFrame(0.0,index = centers_no_dups, columns = all_nodes)
        
        diff_treshold.print_log("Triplet.get_graph_Triplets_heatmap - df ready, fill it")
                
        for triplet in triplets_lst:
            #check the score exists and there is a correlation result between the nodes
            curr_score = triplet.score
            if curr_score != curr_score:
                curr_score = -1
            
            df_node1_score = df.loc[triplet.center, triplet.node1]
            if df_node1_score != 0:
                #more than one triplet with these cells, do average
                new_score = (curr_score + df_node1_score) / 2
            else:
                new_score = curr_score
            df = df.set_value(triplet.center, triplet.node1, new_score)

            #same for node2
            df_node2_score = df.loc[triplet.center, triplet.node2]
            if df_node2_score != 0:
                new_score = (curr_score + df_node2_score) / 2
            else:
                new_score = curr_score
            
            df = df.set_value(triplet.center, triplet.node2, new_score)
        
        diff_treshold.print_log("Triplet.get_graph_Triplets_heatmap - after updating DF, create heatmap")
        
        
        #Way number 1 to do it, didn't work beee
        #now df is filled, create an heatmap from it
        plt.new_figure_manager(0)
        heatmap = sns.heatmap(df)
        heatmap.set_title("heatmap of scores between triplets in the graph")
        heatmap_fig = heatmap.get_figure()
        heatmap_fig.show()
        
        '''
        #wWay number 2 - lets hope
        plt.imshow(df)
        plt.colorbar()
        plt.show()     
        heatmap_fig = plt.figure()
        '''
        if is_save:
            #way number 1
            heatmap_fig.savefig(heatmap_save_path, format = "png", dpi = 1200)
            
        #return for number1
        plt.close(0)
        return triplets_lst,heatmap_fig,df


def list_of_triplets_to_df(lst_of_triplets,graph):
    ''' this method gets a list of triplets and returns a dataframe with those triplets.
    It also gets the graph in order to add the markers of current the nodes'''
    lst = []
    for triplet in lst_of_triplets:
        #make sure the triplet is valid - no nan as score
        if triplet.score == triplet.score:
            center_node = triplet.center
            node1 = triplet.node1            
            node2 = triplet.node2
            center_markers = graph.node[center_node]["markers_lst"]
            node1_markers = graph.node[node1]["markers_lst"]
            node2_markers = graph.node[node2]["markers_lst"]
            #make it pretier
            center_markers = ', '.join(center_markers)
            node1_markers = ', '.join(node1_markers)
            node2_markers = ', '.join(node2_markers)
            
            new_val = (center_node,node1,node2,triplet.corr1,triplet.corr2,triplet.low_corr,triplet.score,center_markers,node1_markers,node2_markers)
            lst.append(new_val)
    
    res_df = pd.DataFrame(data = lst, columns = ["center","node1","node2","corr1","corr2","low corr","triplet score","center markers","node1 markers","node2 markers"])
    res_df.sort_values(["triplet score","center"],inplace = True,ascending=False)
    
    return res_df
    
    


def shuffle_triplets(triplets_df_path, trait_vals_path, ret_df_path, num_of_iterations = 1000, is_return = False):    
    ''' this method check the of a triplet is robust or not by shuffling its data 100 times and createing 100 different triplets
    of same cells, each time saving the score of the triplet'''
    
    diff_treshold.print_log("Triplet.shuffle_triplets - start, before read files")
    
    trait_vals_df = pd.read_csv(trait_vals_path, header=0,index_col="FlowJo Subject ID")
    
    df_no_zero =  trait_vals_df.replace(0,np.nan)
    
    
    triplet_lst = []
    #this list will contain the result of the permutaions and original data   
    data_lst = []
    
    triplets_df = pd.read_excel(triplets_df_path, header = 0)
    needed_df = triplets_df[["center","node1","node2","corr1","corr2","low corr", "triplet score"]]
    
    diff_treshold.print_log("Triplet.shuffle_triplets - after read file,before putting real data")
    for i in range(0, len(needed_df)):
        var = tuple(needed_df.iloc[i])
        var = (0,0) + var
        var = var[1:]
        data_lst.append(var)
        triplet = (var[1],var[2],var[3])
        triplet_lst.append(triplet)
            
    diff_treshold.print_log("Triplet.shuffle_triplets - before loop")
    
    for i in range(0, num_of_iterations):
        if i % 10 == 0:
            diff_treshold.print_log("Triplet.shuffle_triplets - in loop, i = " + str(i))
        shuf_df = shuffle_df(df_no_zero, "row")
        corr_df = shuf_df.corr(method="pearson",min_periods = 90) #changed from 100!!
        for c,n1,n2 in triplet_lst:
            corr1 = abs(corr_df[c][n1])
            corr2 = abs(corr_df[c][n2])
            low_corr = abs(corr_df[n1][n2])
            score = min(corr1,corr2) / low_corr
            perm_var = (i+1,c,n1,n2,corr1,corr2,low_corr,score)
            data_lst.append(perm_var)
    
    diff_treshold.print_log("Triplet.shuffle_triplets - after loop")
    ret_df = pd.DataFrame(data = data_lst,columns = ["perm num","center","node1","node2","corr1","corr2","low corr","triplet score"])
    
    ret_df.to_csv(ret_df_path,index = False)
    
    if is_return:
        return ret_df
            

def shuffle_df(df, axis):     
    ''' this method shuffles a df by row or by col'''
    n = 1
    df = df.copy()
    
    if axis == "row":
        axis_num = 0
    elif axis == "col":
        axis_num = 1        
    for _ in range(n):
        df = df.apply(np.random.permutation, axis=axis_num)
        
    return df


def check_triplets_significance(triplets_shuffled_data_path, pct_to_check ,result_lst_path, is_return = False):
    ''' this method gets as input the output of shuffle_triplets and checks 
    for each triplet (perm num == 0) if the score is bigger than the score in pct_to_check location'''
    
    diff_treshold.print_log("Triplet.check_triplets_significance - start, before read file")
    triplets_shuffled_data = pd.read_csv(triplets_shuffled_data_path,header = 0, index_col = ["center","node1","node2"])
    
    result_lst = []
    all_triplets = list(set(triplets_shuffled_data.index))
    
    
    diff_treshold.print_log("Triplet.check_triplets_significance - after read file, before loop")    
    
    for i in range(0, len(all_triplets)):
        if i % 15 == 0:
            diff_treshold.print_log("Triplet.check_triplets_significance - in loop, i = " + str(i))                
            
        center,node1,node2 = triplet = all_triplets[i]
        triplet_data = triplets_shuffled_data.loc[triplet]
        #get the triplet real score
        orig_score = triplet_data[triplet_data["perm num"] == 0]
        orig_score = float(orig_score["triplet score"])
        #get all the randomized scores
        triplet_data = triplet_data[triplet_data["perm num"] != 0]     
        rand_scores = triplet_data["triplet score"]
        rand_scores_not_nan = rand_scores.dropna()
        rand_scores_good = list(rand_scores_not_nan)
        #check if the triplet is significant
        rand_scores_good.sort()
        n = len(rand_scores_good)
        needed_pos = int(n * pct_to_check)
        if orig_score >= rand_scores_good[needed_pos]:
            is_sig = True
        else:
            is_sig = False
        #real_pos = binary_search(rand_scores,orig_score,0,n)
        #real_pct = real_pos / n
        real_pct = get_real_pct(rand_scores_good,orig_score)
        var = (center, node1,node2, orig_score,is_sig,real_pct)
        result_lst.append(var)
        
    df = pd.DataFrame(data = result_lst, columns = ["center","node1","node2","triplet score", "is significant","real pct"])
    
    df.to_csv(result_lst_path, index = False)
    
    if is_return:
        return df
    

def binary_search(lst,val,left,right):
    ''' this method will return the position of the triplet score real pct'''
    mid = (left + right)//2
    
    if left >= right:
        #we just passed the value
        return mid
        
    if lst[mid] == val:
        return mid
    else:
        if lst[mid] > val:
            return binary_search(lst,val,left,mid - 1)
        else:
            return binary_search(lst,val,mid+1,right)

def get_real_pct(lst,val):
    
    lst.sort()
    n = len(lst)
    i = 0
    while i < n and lst[i] < val:
        i += 1
    
    pct = i / n
    return pct
            