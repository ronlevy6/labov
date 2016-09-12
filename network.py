import diff_treshold
import pandas as pd
import os
import re

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


def create_sub_panel_id_connections(treshold_dir,sub_panel_id):
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
                df = pd.read_csv(full_path,header = 0, index_col = "trait1")
                is_first = False
            else:
                tmp = pd.read_csv(full_path,header = 0, index_col = "trait1")
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
        vals = df.loc[idx]["trait2"] #get trait2 values for specific trait
        
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
        d = create_sub_panel_id_connections(treshold_dir,sub_panel_id)
        diff_treshold.print_log("create_network - after create connections, before combine")           
        final_dict = combine_dicts(final_dict,d)
        diff_treshold.print_log("create_network - after combine")           
    
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