import pandas as pd
import diff_treshold
import random
import sys
import datetime as dt
import numpy as np
import os

RAND_INPUT_DIR = r'/groups/igv/ronlevy/data/random_data/broken_by_panels//'

RAND_OUTPUT_DIR = r'/groups/igv/ronlevy/data/random_data/rand_corr_results/with_nans//'


def randomize_trait_vals_file(trait_vals_path,output_path,is_return):
    ''' this method gets the path to traits values file and shuffels the data of 
    certain panel id'''
    
    rand_print_log("randomize_trait_vals_file - start")
    
    df = pd.read_excel(trait_vals_path, header = 0, index_col = "FlowJo Subject ID")
    
    rand_print_log("randomize_trait_vals_file - after read file before loop")
    
    for i in range(0, len(diff_treshold.PANELS_LST)):
        panel_id = diff_treshold.PANELS_LST[i]
        rand_print_log("randomize_trait_vals_file - panel id is " + panel_id)
        copy = df.copy()
        shuffled_panel_id_df = randomize_panel_id(copy,panel_id)
        
        if i == 0:
            ret_df = shuffled_panel_id_df
        else:
            #not the first panel id, so append
            ret_df = ret_df.append(shuffled_panel_id_df)
    
    ret_df.to_excel(output_path,header= True,index= True)    
    if is_return:
        return ret_df   
    
    
    
def randomize_panel_id(df,panel_id):
    ''' this method does the shuffleing'''
    
    needed_df = df[df["Panel ID"] == panel_id].copy()
    
    rand_print_log("randomize_panel_id - after get panel_id data")
    
    lst_to_shuffle = list(needed_df.columns)
    
    lst_to_shuffle.remove("Panel ID")
    
    random.shuffle(lst_to_shuffle)
    
    rand_print_log("randomize_panel_id - after shuffle list, now shuffle data")
    
    n= len(lst_to_shuffle)    
    
        
    j = 0
    while j < n - 1:
        id1 = lst_to_shuffle[j]
        id2 = lst_to_shuffle[j + 1]
        
        tmp1 = df.loc[:,id1]        
        tmp2 = df.loc[:,id2]        
        
        
        needed_df.loc[:,id1] = tmp2.loc[:,]
        needed_df.loc[:,id2] = tmp1.loc[:,]
        
        j += 2
    
    return needed_df

def rand_create_corr_between_all(df):
    ''' this method gets a FD to a united vals and returns dataframe(!!!!) with 
    all the correlations between all traits'''
    
    rand_print_log("create coor - start - before create needed dataframe")
    #create the dataframe in the needed way
    df_no_zero =  df.replace(0,np.nan)
    trans = df_no_zero.transpose()
    trans = trans.drop("Panel ID")
    trans_float = trans.astype('float32')
    
    
    # run the corr function
    rand_print_log("create coor - before coor")
    corr_res = trans_float.corr(method="pearson",min_periods = 1)
    
    rand_print_log("create coor - before return")
    
    return corr_res


def rand_create_corr_between_two_files(path1, path2):
    ''' Because of Memory Error when running the correlation on all data once this method
    gets 2 paths, combine them to create new DF and sends it into create_corr_between_all'''
    
    rand_print_log("create_corr_between_two_files - start")
    if (path1 == path2):
        # same file, send it directly
        rand_print_log("create_corr_between_two_files - same file")
        df = pd.read_excel(path1, 0, index_col = "FlowJo Subject ID")
    else:
        #create single file
        rand_print_log("create_corr_between_two_files - different files")
        df1 = pd.read_excel(path1, 0, index_col = "FlowJo Subject ID")
        df2 = pd.read_excel(path2, 0, index_col = "FlowJo Subject ID")
        #combine DataFrames
        rand_print_log("create_corr_between_two_files - combine DataFrames")
        df = df1.append(df2, ignore_index = False)
    #now there is a DF, send it's path to create_corr_between_all
    corr_df = rand_create_corr_between_all(df)
        
    rand_print_log("create_corr_between_two_files - before return")
    return corr_df





def rand_create_correlation_map(path1, path2, is_return,to_overwrite):
    ''' This function does calls function described in create_jobs_functions'''
    
    rand_print_log("create_correlation_map - start before create_corr_between_two_files")
    
    corr_df = rand_create_corr_between_two_files(path1, path2)
    tresh_df = "same panel, no tresholded corr DF" 

    panel1 = diff_treshold.get_panel_from_path(path1)
    panel2 = diff_treshold.get_panel_from_path(path2)
    output_path = RAND_OUTPUT_DIR + panel1.upper() + "_" + panel2.upper() + '_corr.csv'
    tresh_output_path = RAND_OUTPUT_DIR + panel1.upper() + "_" + panel2.upper() + '_treshold_corr.csv'
    rand_print_log("create_correlation_map - after calculation of panels and path")
    if (panel1 != panel2):
        #different panels, create tresholded correlation map
        rand_print_log("create_correlation_map - different panels, treshold map")
        if ( not(os.path.isfile(tresh_output_path)) or to_overwrite == True):
            #OK overwrite existing file (if there is one..)
            tresh_df = diff_treshold.create_tresholded_corr_df(corr_df, treshlod=0.6,
                                             compare_only_diff_panels=True,
                                             output_file=tresh_output_path,
                                             print_indicator=800,is_return=True)
        #write tresholded correlation map
    
    #write regular correlation map, the treshold data is written in create_tresholded_corr_df
    rand_print_log("create_correlation_map - before write correlation data")
    if ( not(os.path.isfile(output_path)) or to_overwrite == True):
        #OK overwrite existing file (if there is one..)
        corr_df.to_csv(output_path, header = True, index = True)
    
    if is_return:
        return(corr_df,tresh_df)
    

def rand_get_all_original_nans(corr_res_df,output_path):
    ''' this method gets the path of the dir with all corr data, and writes a csv
    with all pairs of traits that thier correlation wasn't calculated because
    there were not enough observations'''
    
    nans_lst = []
    idx_lst = list(corr_res_df.index)
    n = len(idx_lst)
    for i in range(0,n):
        trait1 = idx_lst[i]
        tmp_lst = list(corr_res_df[trait1].index[corr_res_df[trait1].apply(np.isnan)])
        for trait2 in tmp_lst:
            val = (trait1, trait2)
            nans_lst.append(val)
    
    nans_df = pd.DataFrame(nans_lst)
    
    nans_df.to_csv(output_path,index = False)
            
def rand_break_vals_file_to_files_by_panel(vals_path, res_dir):
    ''' This method gets a path to a SINGELES PANELED vals file and creates 
    different files for each panel ID'''
    
    rand_print_log("break vals file - start")    
       
    
    # read files by Panel ID as index
    vals_df = pd.read_excel(vals_path,0,index_col = "Panel ID")
    
    rand_print_log("break vals file - after file read")    
    
    list_of_panel_id = list(vals_df.index.unique())
    
    rand_print_log("break vals file - before loop")    
    for panel_id in list_of_panel_id:
        rand_print_log("break vals file - in loop for panel id = " + str(panel_id)) 
        #get only the data of specific panel ID
        tmp_df = vals_df.loc[panel_id]
        tmp_df_idx = tmp_df.reset_index()
        #set back to original index and write file
        tmp_df_idx_final = tmp_df_idx.set_index('FlowJo Subject ID')
        new_path = res_dir + "2. Trait Values_SINGLES_PANELED_RANDOMIZED_" + panel_id + ".xlsx"
        tmp_df_idx_final.to_excel(new_path)
        del tmp_df
        del tmp_df_idx
        del tmp_df_idx_final            

def rand_print_log(msg):
    ''' prints massage and current time to screen'''
    print(msg + " " + str(dt.datetime.now()))   
    
def put_nan_in_random_files_according_to_original(random_file,nan_file,rand_file_nans_output):
    ''' this method cputs nan in a correlation in the random data if there is one
    in the original corr file'''
    
    rand_print_log("put_nan_in_random_files_according_to_original - start, read files")
    
    rand_corr_df = pd.read_csv(random_file, header = 0, index_col = "FlowJo Subject ID")
    nans_df = pd.read_csv(nan_file, header = 0, index_col = ["0","1"])
    traits_to_nan = list(nans_df.index) #the indexes are the nan traits
    for idx_tup in traits_to_nan:
        trait1 = idx_tup[0]
        trait2 = idx_tup[1]
        rand_corr_df = rand_corr_df.set_value(trait1, trait2,np.nan) 
        rand_corr_df = rand_corr_df.set_value(trait2, trait1,np.nan) #symetry
    
    rand_corr_df.to_csv(rand_file_nans_output)


def check_random_corr_file(df):
        
    a95 = 0
    a90 = 0
    a85 = 0
    a80 = 0
    a75 = 0
    a70 = 0
    a65 = 0
    a60 = 0
    a_nans = 0   
    a_nans_orig = 0
 
    all_idx = list(df.index)
    n = len(all_idx)
    lst60 = []
    
    for i in range (0, n):
        trait1 = all_idx[i]
        if i % 500 == 0:
            print(i)
        for j in range(i, n):
            trait2 = all_idx[j]
            val = df[trait1][trait2]
            if trait1[:2] != trait2[:2]:
                if val > 0.6:
                    a60 += 1
                    lst60.append((trait1,trait2,df[trait1][trait2]))
                if val > 0.65:
                    a65+= 1
                if val > 0.7:
                    a70 += 1
                if val > 0.75:
                    a75 += 1
                if val > 0.8:
                    a80 += 1
                if val > 0.85:
                    a85 += 1
                if val > 0.9:
                    a90 += 1
                if val > 0.95:
                    a95 += 1
            if  val == np.nan:
                a_nans_orig += 1
            if val == -666:
                a_nans += 1
                
    d ={}
    d[60]  = a60
    d[65]  = a65
    d[70]  = a70
    d[75]  = a75
    d[80]  = a80
    d[85]  = a85
    d[90]  = a90
    d[95]  = a95
    d["nans_orig"] = a_nans_orig
    d["nans"] = a_nans
    
    return d,lst60
    
    
####################mains

def main_for_tresholded_data():
    # get parameters from user
    args = sys.argv
    corr_df_path = args[1]
    treshold = args[2]
    output_dir = args[3]
    panel1 = args[4]
    panel2 = args[5]
    
    is_return = False  
    compare_only_diff_panels = True
    print_indicator = 800    
    
    tresh = float(treshold)
    
    tresh_output_path = output_dir + panel1.upper() + "_" + panel2.upper() + '_treshold_'+ treshold.replace(".","_") + '_corr.csv'
    
    corr_df_path = RAND_OUTPUT_DIR + corr_df_path

    #the : is added here so the real panel ID will be get and not sub panel id
    if (panel1[:2] != panel2[:2]):
        rand_print_log("MAIN - tresh is needed, before read DF")
        corr_df = pd.read_csv(corr_df_path,header = 0, index_col = "FlowJo Subject ID")
        rand_print_log("MAIN - after read DF, before tresh")
        
        diff_treshold.create_tresholded_corr_df(corr_df, tresh, compare_only_diff_panels ,tresh_output_path, print_indicator,is_return)
    

def main_put_nan_in_rand_corr_results():
    #directories
    random_corr_dir = '/groups/igv/ronlevy/data/random_data/rand_corr_results/'
    nans_dir = '/groups/igv/ronlevy/data/corr_results/all_nans/'
    output_dir = '/groups/igv/ronlevy/data/random_data/rand_corr_results/with_nans/'
    
    #read data    
    args = sys.argv
    random_corr_res_path = args[1]
    nan_file_path = args[2]
    
    #combine it
    full_rand_path = random_corr_dir + random_corr_res_path
    full_nans_path = nans_dir + nan_file_path
    full_output_path = output_dir + random_corr_res_path #same filename, different dir
    
    
    put_nan_in_random_files_according_to_original(full_rand_path, full_nans_path, full_output_path)

    
def main_randomize_panels():
    '''
    #read_files    
    args = sys.argv
    traits_val_path = args[1]
    panel_id = args[2]
    output_path = args[3]
    '''
    
    input_path = '/groups/igv/ronlevy/data/5_Trait_Values/2. Trait Values_SINGLES_PANELED.xlsx'
    output_path = '/groups/igv/ronlevy/data/5_Trait_Values/2. Trait Values_SINGLES_PANELED_RANDOMIZED.xlsx'
    
    randomize_trait_vals_file(input_path, output_path,False)
    
def main_get_idx_of_nans():
    #read args

    args = sys.argv
    input_path = args[1]
    output_path = args[2]
    
    corr_df = pd.read_csv(input_path,header = 0, index_col = "FlowJo Subject ID")    
    
    rand_get_all_original_nans(corr_df, output_path)
    
def main_correlation_map_of_random():
    
    #read_files 
    args = sys.argv
    filename1 = args[1]
    filename2 = args[2]
    is_return = args[3]
    to_overwrite = args[4]
    
    path1 = RAND_INPUT_DIR + filename1
    path2 = RAND_INPUT_DIR + filename2
    
    rand_create_correlation_map(path1, path2, is_return,to_overwrite)

    
'''
if __name__ == "__main__":
    main()
'''        
    
        
    
    