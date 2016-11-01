import diff_treshold
import pandas as pd
import os
import re
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#import Triplet

regex = 'P\d\d_P\d\d_treshold.*corr.csv'
PATTERN = re.compile(regex) #check if need to add re.DOTALL

SUB_PANELS_LST = ["P11","P12","P13","P14","P15","P16","P21","P22","P23","P24","P25",
                  "P26", "P31","P32","P41","P42","P51","P61","P71","P72"]
                  
MAX_P_VAL_TRESHOLD_95 = 0.000002         

MAX_P_VAL_TRESHOLD_09 = 0.00000025

MAX_P_VAL_ALL_DATA = 10**6 # the max treshold is 10^5, this will have all in it

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


def create_p_val_network_one_side(p_val_file, max_p_val,index_col,filter_col):
    ''' this method gets an outfile R script file and creates a network according to it 
    using index col as edge start and filter col as edge end'''
    
    diff_treshold.print_log("create_p_val_network_one_side - start")    

    network_dict = {}   
    
    df = pd.read_csv(p_val_file,header = 0,index_col = index_col)
    
    good_p_vals = df[df["needed_p_val"] < max_p_val]
    
    good_p_vals.sort_values(by = "needed_p_val", ascending = True, inplace = True)
    
    diff_treshold.print_log("create_p_val_network_one_side - after read file and filter DF")    
    
    idx_lst = remove_dups(good_p_vals.index)
    for idx in idx_lst:
        vals = good_p_vals.loc[idx][filter_col] #get trait2 values for specific trait
        
        if type(vals) == str:
            vals = [vals] #one variable list
        else:
            vals = list(vals)
        network_dict[idx] = vals
    
    diff_treshold.print_log("create_p_val_network_one_side - after dict - return it")           
    
    return network_dict
    
def create_p_val_network(p_val_file, max_p_val):
    ''' creates a full network'''
    
    diff_treshold.print_log("create_p_val_network - start")           
    final_dict = {}
    
    d = create_p_val_network_one_side(p_val_file, max_p_val, "trait1","trait2")
    final_dict = combine_dicts(final_dict,d)
    
    d = create_p_val_network_one_side(p_val_file, max_p_val, "trait2","trait1")
    final_dict = combine_dicts(final_dict,d)
    
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
            plus_filtered_data = trait_data[trait_data == "+"]
            plus_markers = list(plus_filtered_data.index)
            minus_filtered_data = trait_data[trait_data == "-"]
            minus_markers = list(minus_filtered_data.index)
            markers = []
            for i in range(len(plus_markers)):
                var = plus_markers[i]
                insert_var = var + "+"
                markers.append(insert_var)
            
            for j in range(len(minus_markers)):
                var = minus_markers[j]
                insert_var = var + "-"
                markers.append(insert_var)
                
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
        attr = df.loc[idx_tup]["corr"]
        if graph.has_edge(u, v):
            #attr = df.loc[idx_tup][col_to_add]
            graph.add_edge(u, v, corr = attr)
    
    diff_treshold.print_log("add_corr_to_graph after loop - finish")    
    

def find_all_triplets(graph,edge_attr_name):    
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
            corr_i = graph.get_edge_data(node_x,node_i)[edge_attr_name]
            for j in range(i + 1, len(nei_lst)):
                if nei_lst[j] not in node_i_neis and panel_i != diff_treshold.get_panel_id_substr(nei_lst[j]):
                    #node_x and node_i are neighbours, also node_x and node_j
                    #but node_i and node_j not and they are from different panels
                    corr_j = graph.get_edge_data(node_x,nei_lst[j])[edge_attr_name]
                    var = (node_x,node_i,nei_lst[j],corr_i,corr_j)                    
                    triplets_lst.append(var)
                    
                    #print("enter var")
    
    diff_treshold.print_log("find all triplets after loop return list")        
    return triplets_lst
    


def create_full_graph(p_val_file, max_p_val, united_corr_pass_path, trait_analysis_path):
    ''' this method runs all the needed stuff to get a graph'''
    
    diff_treshold.print_log("create_full_graph - before create network")    
    
    network = create_p_val_network(p_val_file,max_p_val)
    
    diff_treshold.print_log("create_full_graph - before create graph")        
    
    graph = create_graph_from_network(network)
    
    diff_treshold.print_log("create_full_graph - before add corr")    
    
    add_corr_to_graph(graph,united_corr_pass_path)
    
    diff_treshold.print_log("create_full_graph - before add markers")        
    
    add_markers_for_traits(graph,trait_analysis_path)
    
    return graph       



def find_trios_in_graph_same_panel_also(graph,trait_vals_path):
    ''' gets all trios and adds here also the nodes that are from same panel'''

    diff_treshold.print_log("find_trios_in_graph_same_panel_also - start")
    trios_lst = []
    
    df = pd.read_excel(trait_vals_path, header = 0, index_col = "FlowJo Subject ID")

    corr_df = diff_treshold.create_corr_between_all(df)    
    
    for node_x in graph.nodes():
        
        nei_lst = graph.neighbors(node_x)
        for i in range(0 , len(nei_lst)):
            node_i = nei_lst[i]
            node_i_neis = graph.neighbors(node_i)
            corr_i = graph.get_edge_data(node_x,node_i)["corr"]
            for j in range(i + 1, len(nei_lst)):
                if nei_lst[j] not in node_i_neis :
                    #node_x and node_i are neighbours, also node_x and node_j
                    #but node_i and node_j not and they are from different panels
                    corr_j = graph.get_edge_data(node_x,nei_lst[j])["corr"]
                    unconnected_corr = corr_df.loc[node_i,nei_lst[j]]
                    var = (node_x,node_i,nei_lst[j],corr_i,corr_j,unconnected_corr)                    
                    trios_lst.append(var)
                    #print("enter var")
    
    diff_treshold.print_log("find_trios_in_graph_same_panel_also - after loop return list")        
    return trios_lst    

def analyze_trios(trios_lst, treshold):
    
    same_panel_high_corr = 0
    same_panel_low_corr = 0
    diff_panel_high_corr = 0
    diff_panel_low_corr = 0
    
    for trio in trios_lst:
        node1 = trio[1]
        node2 = trio[2]
        unconnected_corr = trio[5]
        diff_panel = False
        low_corr = unconnected_corr < treshold
        if diff_treshold.get_panel_id_substr(node1) != diff_treshold.get_panel_id_substr(node2):
            diff_panel = True
            
        if diff_panel:
            if low_corr:
                diff_panel_low_corr += 1
            else:
                diff_panel_high_corr += 1
        else:
            if low_corr:
                same_panel_low_corr += 1
            else:
                same_panel_high_corr += 1
    
    return (same_panel_high_corr,same_panel_low_corr,diff_panel_high_corr,diff_panel_low_corr)
        

def merge_graph_nodes(graph,trait_val_path, min_unconneced_treshold, min_common_markers):
    ''' this method unites nodes in the graph when the following conditions are met:
    1 - they have at least min_common_markers (3) common markers and above 66% common markers
    2 - they have the same neighbours. OR - for each uncommon neighbour the correlation
    between the non-connected traits is above min_unconneced_treshold (0.85).
    When merging 2 nodes the following attributes will be changed:
        markers_lst - combine both traits markers
        correlation - for each edge the new correlation will be the average of correlations
        merged_with - list with previous nodes before merge'''
    
    diff_treshold.print_log("merge_graph_nodes - start, before read file")
    
    traits_val_df = pd.read_csv(trait_val_path, header = 0, index_col = "FlowJo Subject ID")
    
    df_no_zero =  traits_val_df.replace(0,np.nan)
    
    to_corr_df = df_no_zero.astype('float32')
    
    corr_df = to_corr_df.corr(method="pearson",min_periods = 100)
    
    diff_treshold.print_log("merge_graph_nodes - after read file and calc corr, before loop")
    
    merged_nodes_lst = [] # a list of tuples to contain all merged nodes
    
    all_nodes = list(graph.nodes())
    n = len(all_nodes)
    
    for i in range(0, n):
        if i % 100 == 0:
            diff_treshold.print_log("merge_graph_nodes - in loop, i = " + str(i))            
        curr_node = all_nodes[i]
        if (curr_node in graph):
            #curr node wasn't merged
            j = i + 1
            while j < n:
                #second node wasn't merged as well
                if (all_nodes[j] in graph):
                    if can_merge_nodes(graph,curr_node,all_nodes[j], corr_df, min_unconneced_treshold, min_common_markers):
                        merged_tup = (curr_node, all_nodes[j])
                        merged_nodes_lst.append(merged_tup)
                        merge_2_nodes(graph, curr_node, all_nodes[j], corr_df)
                        #remove when doing real merging
                        j = n #exit loop and start merging the next node
                j += 1
                        
    return merged_nodes_lst
    

def can_merge_nodes(graph, node1, node2, corr_df, min_unconneced_treshold, min_common_markers):
    ''' this method gets a graph, 2 nodes and traits val DF and returns True if the
    nodes can be merged according to conditions described in unite_graph_nodes method'''
    
    # check neighberood        
    node1_neis = graph.neighbors(node1)    
    node2_neis = graph.neighbors(node2)
    
    node1_only_neis = list(set(node1_neis) - set(node2_neis))    
    node2_only_neis = list(set(node2_neis) - set(node1_neis))
    
    for nei in node1_only_neis:
        if corr_df.loc[node2,nei] < min_unconneced_treshold:
            return False
    
    for nei in node2_only_neis:
        if corr_df.loc[node1, nei] < min_unconneced_treshold:
            return False
    
    #if got here - passed the neighbours condition, now check the markers
    node1_markers = graph.node[node1]["markers_lst"]
    node2_markers = graph.node[node2]["markers_lst"]
    
    common_markers = set(node1_markers).intersection(node2_markers)
    num_of_common = len(common_markers)
    #check for enough common markers
    if num_of_common < min_common_markers:
        return False
    
    #check for common markers pct
    if num_of_common / len(node1_markers) < 0.66 or num_of_common / len(node2_markers) < 0.66:
        return False
    
    #if got here - passed all checks and nodes can be merged
    return True
    
    
def merge_2_nodes(graph, node1, node2, corr_df):
    ''' this method will merge node1 with node2 and make the changes written in merge_graph_nodes'''
    
    #unite markers_lst
    full_markers_lst = list(graph.node[node1]["markers_lst"]) + list(graph.node[node2]["markers_lst"])
    distinct_markers_lst = list(set(full_markers_lst))
    
    #get all prevoius nodes
    node1_prev = graph.node[node1].get("merged_with",[])
    if type(node1_prev) == str:
        node1_prev = [node1_prev]
    
    node2_prev = graph.node[node2].get("merged_with",[])
    if type(node2_prev) == str:
        node2_prev = [node2_prev]
    
    all_prevs = node1_prev + node2_prev + [node2] #add node2 to the list
    distinct_prevs = list(set(all_prevs)) 
    
    #fix all correlations between nodes
    node1_neis = graph.neighbors(node1)    
    node2_neis = graph.neighbors(node2)
    all_neis = node1_neis + node2_neis
    distinct_neis = list(set(all_neis))
    
    #print("all_neis = ")
    #print(all_neis)
    
    #print("distinct_neis = " )
    #print(distinct_neis)
    
    new_edges_lst = []
    for nei_node in distinct_neis:
        corr1 = corr_df.loc[node1,nei_node]
        corr2 = corr_df.loc[node2,nei_node]
        if corr1 == corr1 and corr2 == corr2:
            #both are real vals
            new_corr = (corr1 + corr2) / 2
        elif corr1 == corr1:
            #only corr1 is good
            new_corr = corr1
        else:
            #only corr2 is valid
            new_corr = corr2
        new_edge = ( nei_node, new_corr)
        new_edges_lst.append(new_edge)
    
    #update node fields
    graph.add_node(node1, markers_lst = distinct_markers_lst, merged_with = distinct_prevs)
    
    #create and update new edges
    for edge_to_add in new_edges_lst:
        graph.add_edge(node1, edge_to_add[0], corr = edge_to_add[1])
    
    #remove the merged node
    graph.remove_node(node2)
    
    #return new_edges_lst


def get_graph_all_data(combined_file, trait_anal_path, n):
    ''' this method creates a graph and runs the following checks on it:
    p - val - 0.01 / n
    #nodes
    #edges
    #trios
    #connected components.
    The tests are ran twice - first time for all data, the second one is
    only for the connections below the p- val.
    combined file - united file with p val and with correlation'''
    
    p_val = 0.01 / n  
    
    full_graph = create_full_graph(combined_file, MAX_P_VAL_ALL_DATA, combined_file, trait_anal_path)
    full_graph_dict = fetch_graph_data(full_graph)
    
    p_val_graph = create_full_graph(combined_file, p_val, combined_file, trait_anal_path)
    p_val_graph_dict = fetch_graph_data(p_val_graph)
    
    #fetch graphs data
    
    ret = (p_val, full_graph, full_graph_dict, p_val_graph, p_val_graph_dict)    
    return ret
    
def fetch_graph_data(graph):
    
    nodes_num = len(graph.nodes())
    edges_num = len(graph.edges())
    
    trios = find_all_triplets(graph, "corr")
    trios_num = len(trios)
    
    graph_connected_components = list(nx.connected_components(graph))
    connected_components_num = len(graph_connected_components)        
    
    ret_dict = {}
    ret_dict["num_of_nodes"] = nodes_num
    ret_dict["num_of_edges"] = edges_num
    ret_dict["num_of_trios"] = trios_num
    ret_dict["num_of_connected_components"] = connected_components_num
    
    return ret_dict

def try_to_print():
    ''' main for printing the graph'''
    united_corr_pass_path = p_val_file = r'C:\Ron\אוניברסיטה\תל אביב\שנה ב\פרויקט\twins data\changed_data\666_Ron\tresh_098\tresh_98_super_duper_file.csv'
    max_p_val = 0.01 / 794
    trait_analysis_path = r'C:\Ron\אוניברסיטה\תל אביב\שנה ב\פרויקט\twins data\changed_data\4_Trait_Analysis\4. Trait Analysis_good.xlsx'
    graph = create_full_graph(p_val_file,max_p_val, united_corr_pass_path, trait_analysis_path)

    pos=nx.spring_layout(graph)
    nx.draw(graph,pos,node_color='#A0CBE2',edge_color='#BB0000',width=2,edge_cmap=plt.cm.Blues,with_labels=True)
    
    #option number1
    plt.savefig("graph.png", dpi=500, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1)
    
    #option number 2
    plt.savefig("graph2.png",dpi = 1000)
    
    #option 3
    plt.savefig("graph2.pdf")
    
    #different figure type
    plt.figure(3,figsize=(12,12)) 
    nx.draw(graph,pos,node_size=60,font_size=8) 
    plt.savefig("graph3.png",dpi = 1000)
    
    
    return graph
    
    
    
    
    