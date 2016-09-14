import diff_treshold
import network
import pandas as pd

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
        
    
    def create_all_triplets(graph, trait_values_path):
        ''' this method gets a graph and creates all the triplets in it.
        The path of trait_values is needed to get the correlation between
        the unconnected nodes'''
            
        diff_treshold.print_log("create_all_triplets - start")
            
        graph_triplets = network.find_all_triplets(graph)
            
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
            triplet = Triplet(center,node1, node2, corr1, corr2, unconnected_corr)
            triplets_lst.append(triplet)
            
        diff_treshold.print_log("create_all_triplets - after loop - finished") 
        return triplets_lst