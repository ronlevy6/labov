import diff_treshold
import pandas as pd
import os
import re
import numpy as np
import networkx as nx
import Triplet

regex = 'P\d\d_P\d\d_treshold.*corr.csv'
PATTERN = re.compile(regex) #check if need to add re.DOTALL

SUB_PANELS_LST = ["P11","P12","P13","P14","P15","P16","P21","P22","P23","P24","P25",
                  "P26", "P31","P32","P41","P42","P51","P61","P71","P72"]


def get_file_suffix(filename):
    ''' This method gets a filename and returns its suffix'''
    
    diff_treshold.print_log("get_file_suffix start")
    
    rev = filename[::-1] # reverse the file
    pos = rev.find(".") #last "." in original filename
    rev_ret = rev[:pos]
    ret = rev_ret[::-1]
    
    diff_treshold.print_log("get_file_suffix end")    
    
    return ret


def create_sub_panel_id_connections(treshold_dir,sub_panel_id,index_column,filter_col):
    ''' this method gets a directory, and sub - panel ID in which there are the traits that their
    correlation is higher than treshold. It connects all of them to single DF
    and creates a dictionary'''
    
    diff_treshold.print_log("create_sub_panel_id_connections start")

    is_first = True
    did_something = False
    #check which pattern to take
    if sub_panel_id != None:
        pattern = sub_panel_id + '_P\d\d_treshold.*corr.csv'
    else:
        #no sub pannel id given - run on all data
        pattern = PATTERN
    
    for filename in os.listdir(treshold_dir):
        if re.match(pattern,filename) != None:
            #correct file
            full_path = treshold_dir + filename
            diff_treshold.print_log("create_sub_panel_id_connections - before read file")
            #check if should be with index or not!!!!!!!!!!!!!
            if is_first:
                df = pd.read_csv(full_path,header = 0, index_col = index_column)
                is_first = False
            else:
                tmp = pd.read_csv(full_path,header = 0, index_col = index_column)
                df = df.append(tmp, ignore_index = False)
            diff_treshold.print_log("create_sub_panel_id_connections - after read file")
            did_something = True
    #now the df is full with all the data
    diff_treshold.print_log("create_sub_panel_id_connections - created united DF, before mega dict")
    
    if not (did_something):
        #no data for specific sub_panel_id
        return
        
    d = {} #init empty dict
    df.sort_values(by = "corr", ascending = False, inplace = True)
    idx_lst = remove_dups(df.index)
    for idx in idx_lst:
        vals = df.loc[idx][filter_col] #get trait2 values for specific trait
        
        if type(vals) == str:
            vals = [vals] #one variable list
        else:
            vals = list(vals)
        d[idx] = vals
    
    diff_treshold.print_log("create_sub_panel_id_connections - after dict - return it")           
    
    return d
            
            
def create_network(treshold_dir):
    ''' This method gets a dir of tresholded files and creates the connections 
    for each sub panel id and returns a united dictionary'''
    diff_treshold.print_log("create_network - start")           
    final_dict = {}
    for sub_panel_id in SUB_PANELS_LST:
        diff_treshold.print_log("create_network - now for sub panel id " + sub_panel_id)           
        d = create_sub_panel_id_connections(treshold_dir,sub_panel_id,"trait1","trait2")
        diff_treshold.print_log("create_network - after first create connections, before combine")           
        final_dict = combine_dicts(final_dict,d)
        diff_treshold.print_log("create_network - after first combine")
        
        ## second run on the same sub panel id with different traits as index
        diff_treshold.print_log("create_network - second run for sub panel id " + sub_panel_id)
        d = create_sub_panel_id_connections(treshold_dir,sub_panel_id,"trait2","trait1")
        diff_treshold.print_log("create_network - after second create connections, before combine")
        final_dict = combine_dicts(final_dict,d)
        diff_treshold.print_log("create_network - after second combine")
    
    diff_treshold.print_log("create_network - finish, before return")           
    return final_dict


    
def combine_dicts(dict1, dict2):
    ''' This method gets 2 dictionaries and combines them'''
    
    if dict1 == None:
        return dict2
    if dict2 == None:
        return dict1
    
    diff_treshold.print_log("combine dicts - start")
    combined = {}
    keys = list(dict1.keys()) + list(dict2.keys())
    # remove duplicates using the set datatype    
    keys_non_dup = remove_dups(keys) 
    #gather values for each key and remove duplicates
    for key in keys_non_dup:
        vals = dict1.get(key, []) + dict2.get(key, [])
        vals_non_dup = remove_dups(vals)
        combined[key] = vals_non_dup
        
    diff_treshold.print_log("combine dicts - after gathering all data")
    return combined
        
def remove_dups(lst):
    ''' this method removes duplicates from a list '''
    t = type(lst) #lst can be string or list
    n = len(lst)
    if t == str:
        ret = [lst] # a list with one var        
       
    else:
        if (n == 1):
            ret = lst
        else:
            ret = list(set(lst))

    
    return ret
    
def get_network_hubs(network):
    ''' This method gets a dictionary network and returns a list of nodes that have 
    more than avg + 2sd neighbours'''
    
    diff_treshold.print_log("get network hubs start - before get_num_of_neis_per_key")
    (key_num_of_neis,lst_num_of_neis) = get_num_of_neis_per_key(network)
    
    avg_of_nei = sum(lst_num_of_neis) / len (lst_num_of_neis)
    sd = np.std(lst_num_of_neis)    
    treshold = avg_of_nei + 2 * sd
    diff_treshold.print_log("get network hubs - after get_num_of_neis_per_key and calc of statsitacal data")
    
    hub_lst = []
    i = 0
    #key_num_of_neis is ordered
    while i < len(key_num_of_neis) and key_num_of_neis[i][1] > treshold:
        hub_lst.append(key_num_of_neis[i])
        i += 1
    
    diff_treshold.print_log("get network hubs - after appending hub list - return it")
    return hub_lst

def get_hubs(graph):
    ''' this method gets a graph and returns list of nodes with amount of neighbours 
    thats higher than avg + 2 * sd'''

    diff_treshold.print_log("get hubs start - before nodes loop")
    nodes_nei_lst = []  
    only_nei_num = []
    
    
    for node in graph.nodes():
        deg = graph.degree(node)
        nodes_nei_lst.append((node, deg))
        only_nei_num.append(deg)
        
    
    diff_treshold.print_log("get hubs after nodes loop")    
    
    # calculate treshold
    nei_avg = sum(only_nei_num) / len(only_nei_num)
    sd = np.std(only_nei_num)
    treshold = nei_avg + sd * 2
    nodes_nei_lst.sort(key = lambda x: x[1], reverse = True)
    
    hubs_lst = []
    i = 0
    while i < len(nodes_nei_lst) and nodes_nei_lst[i][1] > treshold:
        hubs_lst.append(nodes_nei_lst[i])
        i += 1
    
    return hubs_lst
        
    

def get_num_of_neis_per_key(network):
    ''' This method gets a network and returns a list with each var a tuple
    in the way - (key, degree). It also returns a list with num of neighbours only'''
    ret_lst = []
    for key in network.keys():
        nei_num = len(network[key])
        val = (key,nei_num)
        ret_lst.append(val)    
    
    ret_lst.sort(key = lambda x: x[1], reverse = True)
    num_of_neis_only = [x[1] for x in ret_lst]    
    return (ret_lst,num_of_neis_only)


def create_graph_from_network(network):
    ''' This method gets a network and creates a graph from it'''
    
    diff_treshold.print_log("create graph - start")
    
    g = nx.Graph()
    for key in network.keys():
        g.add_node(key, panel_id = diff_treshold.get_panel_id_substr(key))
    #g.add_nodes_from(network.keys())
    
    diff_treshold.print_log("create graph - before adding edges")
    #create list of verteces
    neigbours_lst = []    
    for edge in network.keys():
        for val in network[edge]:
            neigbours_lst.append((edge, val))
        
    g.add_edges_from(neigbours_lst)
    
    diff_treshold.print_log("create graph - after adding edges, return graph")
    return g

def add_markers_for_traits(graph, trait_analysis_path):
    ''' this method adds the markers of each trait (node) in the graph'''
    
    diff_treshold.print_log("add_markers_for_traits start before read file")
    
    trait_analysis_df = pd.read_excel(trait_analysis_path, header = 0, index_col = "Trait ID")
    
    diff_treshold.print_log("add_markers_for_traits after read file before loop")
    # run on all traits, for each one, check if has a node in the graph. 
    
    for trait in trait_analysis_df.index:
        if trait in graph:
            #trait has a node, take only the markers that are "+" and add as attributes
            trait_data = trait_analysis_df.loc[trait]
            filtered_data = trait_data[trait_data == "+"]
            markers = list(filtered_data.index)
            markers.sort()
            graph.add_node(trait, markers_lst = markers)
    
    diff_treshold.print_log("add_markers_for_traits finished")

def add_corr_to_graph(graph, united_pass_path):
    ''' this method gets the graph and a path to a file contains all 
    pairs of traits in the graph and add the correlation as data of the edge'''
    
    diff_treshold.print_log("add_corr_to_graph start")    
    df = pd.read_csv(united_pass_path, header = 0, index_col = ["trait1","trait2"])
    diff_treshold.print_log("add_corr_to_graph after read file, before loop on edges")    
    
    for idx_tup in df.index:
        u = idx_tup[0]
        v = idx_tup[1]
        if graph.has_edge(u, v):
            attr = df.loc[idx_tup]["corr"]
            graph.add_edge(u, v, corr = attr)
        else:
            #sanity check
            print("something is tahat!!!")
    
    diff_treshold.print_log("add_corr_to_graph after loop - finish")    
    

def find_all_triplets(graph):    
    ''' this method runs on the graph and finds triplets of nodes x,i,j such that
    the edges x-i and x-j exist but the edge i-j doesn't'''
    
    
    diff_treshold.print_log("find all triplets start")
    triplets_lst = []
    
    for node_x in graph.nodes():
        #diff_treshold.print_log("find all triplets new node - " + str(node_x))
        nei_lst = graph.neighbors(node_x)
        for i in range(0 , len(nei_lst)):
            node_i = nei_lst[i]
            panel_i = diff_treshold.get_panel_id_substr(node_i)
            node_i_neis = graph.neighbors(node_i)
            corr_i = graph.get_edge_data(node_x,node_i)["corr"]
            for j in range(i + 1, len(nei_lst)):
                if nei_lst[j] not in node_i_neis and panel_i != diff_treshold.get_panel_id_substr(nei_lst[j]):
                    #node_x and node_i are neighbours, also node_x and node_j
                    #but node_i and node_j not and they are from different panels
                    corr_j = graph.get_edge_data(node_x,nei_lst[j])["corr"]
                    var = (node_x,node_i,nei_lst[j],corr_i,corr_j)                    
                    triplets_lst.append(var)
                    #print("enter var")
    
    diff_treshold.print_log("find all triplets after loop return list")        
    return triplets_lst