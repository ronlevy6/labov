import diff_treshold
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt 
import network
import networkx as nx


def load_df(united_file_path):
    ''' this method loads the df of all data with correlation above treshold'''
    df = pd.read_csv(united_file_path,header = 0, index_col = ["trait1","trait2"])
    return df


def filter_df_according_to_panels(united_df, panel1, panel2):
    ''' This method gets the DF of above treshold and returns a DF with traits only from 
    panel1 and panel2'''
    
    # condition for needed indecis
    good_idx = [x[0].startswith(panel1) and x[1].startswith(panel2) 
    for x in united_df.index]
    
    filtered_df = united_df[good_idx]
    
    return filtered_df

def corr_dist_between_pannels(filtered_df,panel1,panel2,is_show):
    ''' this method calculates the distribution of the correlation between 2 panels'''
    
    corr_dist = np.array(filtered_df["corr"])
    plt.xlabel("correlation")
    plt.ylabel("frequency")
    plt.title("frequency for correlation between panels " + panel1 + " and " + panel2)
    plt.hist(corr_dist)
    fig_name = panel1 + "_" + panel2 + "_corr_dist.png"
    if is_show:
        plt.savefig(fig_name)
        plt.show()
    

def calc_corr_dist_per_all_panels(united_df):
    ''' this method calculates the distributaion of correaltions between all panels'''
    
    n = len(diff_treshold.PANELS_LST)
    for i in range(0, n):
        panel1 = diff_treshold.PANELS_LST[i]
        for j in range(i + 1, n):
            #no need to calc between same panels
            panel2 = diff_treshold.PANELS_LST[j]
            filtered_df = filter_df_according_to_panels(united_df,panel1,panel2)        
            corr_dist_between_pannels(filtered_df,panel1,panel2,is_show = False)
                                      
                                      
def pie_chart_neigbours_panel(united_df,is_pie,is_save):
    ''' this method gets a united_df and returns the distribution of panels of neigbours'''
        
    n = len(diff_treshold.PANELS_LST)
    nei_mat = np.zeros((n,n)) #matrix with cell for each two panels combination
    
    for i in range(0, n):
        panel1 = diff_treshold.PANELS_LST[i]
        for j in range(i + 1, n):
            #no need to calc between same panels
            panel2 = diff_treshold.PANELS_LST[j]
            #get data relevant only for these 2 panels
            filtered_df = filter_df_according_to_panels(united_df,panel1,panel2)
            nei_mat[i][j] = len(filtered_df)
            nei_mat[j][i] = len(filtered_df) # symetry
    
    #print and save results
    labels = diff_treshold.PANELS_LST
    colors = ['gold','yellow','blue','red','green','brown']
    for i in range(0, n):
        #set the plot
        plt.xlabel("panel ID")
        plt.ylabel("amount of neigbours")
        plt.title("histogram of neigbours of traits from panel P" + str(i+1))
        if is_pie:
            plt.pie(nei_mat[i], colors = colors,labels = labels,autopct='%1.1f%%')
        else:
            #bar plot
            plt.bar(np.arange(n), nei_mat[i],align='center')
            plt.xticks(np.arange(n), labels)
            
        if is_save:
            fig_name = "P" + str(i+1) + "_neigbours_distribution_full.png"
            plt.savefig(fig_name)
        plt.show()
        
def degree_distribution(graph,panel_id,is_save):
    ''' this method gets a graph as input and prints for each panel the distribution
    of neighbours number'''
    
    num_of_nei_lst = []
    #get only nodes of needed panel id
    panel_id_nodes = [n for n,attrdict in graph.node.items() 
    if attrdict ['panel_id'] == panel_id]
        
    for node in panel_id_nodes:
        num_of_nei_lst.append(graph.degree(node))
        
    plt.title("amount of neighbours histogram for panel id " + panel_id)
    plt.ylabel("amount of nodes")
    plt.xlabel("amount of neighbours")
    plt.hist(num_of_nei_lst)
    fig_name = panel_id + "_amount_of_neigbours_dist.png"
    if is_save:
        plt.savefig(fig_name)
    plt.show()
        
def all_pannels_degree_distribution(graph,is_save):
    
    for panel_id in diff_treshold.PANELS_LST:
        degree_distribution(graph,panel_id,is_save)


def calc_hub_and_neighbours_markers(graph, is_save):
    ''' this method gets a graph, uses get hubs to find the center nodes 
    and calculates the markeres of the neighbours of those nodes'''    
    
    diff_treshold.print_log("calc_hub_and_neighbours_markers start")
    
    hubs_lst = network.get_hubs(graph)
    
    diff_treshold.print_log("calc_hub_and_neighbours_markers after get hubs before loop")
    for hub in hubs_lst:
        hub_name = hub[0]
        markers = [] # init empty list to be populated with all markers
        for neighbour in graph.neighbors(hub_name):
            markers += graph.node[neighbour]["markers_lst"]
        
        # create a dictionary from the list when key is marker and value is count
        diff_treshold.print_log("calc_hub_and_neighbours_markers in loop after list creation")
        d = create_dict_from_list(markers)
        label = "histogram of markers of the neighbours of the hub " + hub_name
        xlabel = "markers"
        ylabel = "count"
        file_name = "markers_of_neigbours_of_hub_" + hub_name.replace(":","_") + ".png"
        histogram_dict(d, file_name, label , xlabel, ylabel ,'green', is_save)
    
    diff_treshold.print_log("calc_hub_and_neighbours_markers finished")
        

def calc_markers_per_connected_components(graph, is_save):
    ''' this methdo gets a graph and creates an plotbar for each connected component
    with the count of each marker in the nodes in the connected component'''
    
    diff_treshold.print_log("calc_markers_per_connected_components start")
    
    connected_components = list(nx.connected_components(graph))
    for i in range(0, len(connected_components)):
        group = connected_components[i]
        markers = []
        for node_group in group:
            markers += graph.node[node_group]["markers_lst"]
            
        #create a dictionary from the list and plot it
        diff_treshold.print_log("calc_markers_per_connected_components before plotting")
        d = create_dict_from_list(markers)
        label = "markers per connected component with " + str(len(group)) + " cells"
        xlabel = "markers"
        ylabel = "count"
        file_name = "markers_for_connected_component_number_" + str(i) + ".png"
        histogram_dict(d, file_name, label, xlabel, ylabel, 'yellow', is_save)
    

          
def create_dict_from_list(lst):
    ''' this method gets a list and returns a dictionary with the count of each var'''
    ret_dict = {}
    for val in lst:
        ret_dict[val] = ret_dict.get(val, 0) + 1
    
    return ret_dict
    

def histogram_dict(d, filename , label, xlabel, ylabel ,color, is_save):
    ''' this method gets a dictionary in the form of key:count and print histogram of it'''
    
    
    axis_font = { 'size':'5'}    
    plt.bar(range(len(d)), d.values(), align='center',color = color)
    plt.xticks(range(len(d)), d.keys(), **axis_font)
    plt.title(label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    
    if is_save:
        plt.savefig(filename,dpi = 420)
    
    plt.show()
    

def main():
    # get parameters from user
    args = sys.argv
    united_data_path = args[1]
    treshold_dir = args[2]
    
    multi_idx_df = load_df(united_data_path)
    calc_corr_dist_per_all_panels(multi_idx_df)
    
    pie_chart_neigbours_panel(multi_idx_df,is_pie=False,is_save=False)
    
    above_treshold_network = network.create_network(treshold_dir)
    graph = network.create_graph_from_network(above_treshold_network)
    
    
    
    

'''                                      
if __name__ == "__main__":
    main()                                      
'''