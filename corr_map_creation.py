import pandas as pd
import datetime as dt
import os
import sys

TWINS_DATA = r'/groups/igv/ronlevy/data//'

DEST_DIR = r'/groups/igv/ronlevy/data//'

TMP_DIR = r'/groups/igv/ronlevy/data/tmp//'

PANELS_LST = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']

pd.set_option('precision',4)

OUTPUT_DIR = r'/groups/igv/ronlevy/data/corr_results//'

INPUT_DIR = r'/groups/igv/ronlevy/data/5_Trait_Values/By_Panels/broken//'

# Consts for creationg correlation maps
DEFAULT_TRESHOLD = 0.2
COMPARE_DIFF_PANELS_ONLY = True
PRINT_INDICATOR = 400


################# irrelevant part #########################

def get_idx_of_samples(demographic_path):
    ''' This method runs on the demograpic xlsx file and returns numpy nd-array
    of indexes represnts the sample ID's of chosen samples while 
    making sure to not take two sisters'''
    
    full_path = TWINS_DATA + demographic_path
    df = pd.read_excel(full_path,0)  # read first sheet only
    
    grouped = df.groupby("Unique Family ID").min() # get min sample ID per Family ID
    grouped = grouped["FlowJo Sample ID"] # Sample ID is the only needed data
    ret_val = grouped.values
    
    #ret_val = list(grouped.values)
    
    return ret_val


def create_panel_dictionary(mmc1_path):
    ''' This method gets mmc1.xlsx path and returns a dictionary which built as:
    Key -> 1..7 represents each panel.
    Value -> list of markers for specefic panel according to the mmc1 file'''
    
    hard_coded_data_to_remove = ["Reagent","Viability"]
    # init empty dictionary
    ret_dict = {}
       
    full_path = TWINS_DATA + mmc1_path
    df = pd.read_excel(full_path, skiprows = 2) # to get panels as indexes
    #run on all 7 panels
    for i in range(1,8):
        tmp_lst = []
        key = 'P' + str(i)
        for val in df[key]:
            if (val == val):
                # way to check the value is not NaN
                tmp_lst.append(val)
        
        # remove unneeded data
        tmp_lst.pop(0) # remove tmp value
        
        for check_remove in hard_coded_data_to_remove:
            if check_remove in tmp_lst:
                tmp_lst.remove(check_remove)        
            
        ret_dict[i] = tmp_lst
    
    return ret_dict
            
        
def get_panel_number_by_dict(trait_analysis_path,mmc1_path):
    ''' This method creates new trait_analysis file with extra column - panel number.
    It does so according to the columns in the end of the original trait_analysis file
    and mmc1 file'''
    panel_lst = [] #Empty list to be populated with panel Id for each trait
    panel_dict = create_panel_dictionary(mmc1_path)
    panel_numbers = list(panel_dict.keys()) # get number of panels
    full_path = TWINS_DATA + trait_analysis_path
    
    # Read file and use trait id as index column
    df_idx = pd.read_excel(full_path,0,index_col = "Trait ID") 
    
    for trait_id in df_idx.index:
        panel_id = get_panel_id_dict(trait_id,df_idx.loc[trait_id],panel_numbers,panel_dict)
        panel_lst.append(panel_id)
        
    df_idx.insert(0,"Panel ID",panel_lst,True) #Add panel id column to file
    
    new_path = DEST_DIR + trait_analysis_path
    df_idx.to_excel(new_path, index = True) #check field vals




def get_panel_id_dict(trait_id, single_trait_data,panel_numbers,panel_dict ):
    ''' This method gets trait_id and it's data. It also gets all panels exists 
    and the dictionary of the pannels and returns the panel ID for specific trait ID'''
    
    copy_panel_numbers = panel_numbers[::] #copy of the list
    reduced_data = single_trait_data[13:] #remove unneeded data 
    is_non_null = False
        
    
    indexes = reduced_data.index
    # check over all optional markers
    for i in range(0, len(reduced_data)):
        if reduced_data[i] == reduced_data[i]:
            is_non_null = True
            #current field is not nan
            marker = indexes[i] #get marker name
            for j in copy_panel_numbers:
                if marker not in panel_dict[j]:
                    copy_panel_numbers.remove(j)
    if (is_non_null == True):
        ret_val = str(copy_panel_numbers).strip('[]') # change panels into string
    else:
        ret_val = "All Markers are Null"
    return ret_val














################# only this is good #########################



def add_panel_id_column_to_file(path,header,mark):
    ''' This method gets path to trait analysis or values file and creates a new one
    with a new column - panel ID. Does it according to first letters of trait ID.
    For analysis file the header is Trait ID. For values file its FlowJo Subject ID'''
    
    panel_lst = [] #Empty list to be populated with panel Id for each trait
    full_path = TWINS_DATA + path
    
    print_log("add panel id - before file read") 
    
    # Read file and use trait id as index column
    df_idx = pd.read_excel(full_path,0,index_col = header )
    for trait in df_idx.index:
        panel_id = get_panel_id_substr(trait) # get substring of panel ID only
        
        panel_lst.append(panel_id)
    
    print_log("add panel id - after list creation")     
    
    df_idx.insert(0,"Panel ID",panel_lst,True) #Add panel id column to file
    
    print_log("add panel id - before new path and write") 
    
    # create new file and write it
    #new_path = DEST_DIR + path
    new_path = create_new_path(DEST_DIR,path,mark)
    df_idx.to_excel(new_path, index = True) #check field vals
    
    print_log("add panel id - after file write") 
        
    
def get_panel_id_substr(trait_id):
    ''' This method gets trait id and returns the pannel it belongs to according to 
    the id. 
    DOESNT WORK FOR MFI AND LIN'''
     # 2 separators - space and :
    pos_dot = trait_id.find(":")
    pos_space = pos = trait_id.find(" ")
    if pos_space == -1:
        pos = pos_dot
    elif pos_dot == -1:
        pos = pos_space
    else:
        # trait ID has both separators
        pos = min(pos_space,pos_dot)
    
    ret_val = trait_id[:pos]
    return ret_val
    

def create_new_demographic_file(demographic_path):
    ''' This method runs on the demograpic xlsx file writes a new one containing only
    one sample from each family'''
    
    full_path = TWINS_DATA + demographic_path
    df = pd.read_excel(full_path,0)  # read first sheet only
    
    grouped = df.groupby("Unique Family ID").min() # get min sample ID per Family ID
    
    # write new data into file
    new_path = create_new_path(DEST_DIR ,demographic_path,"CHANGED")
    grouped.to_excel(new_path)
    
    

def create_new_path(directory, path,mark):
    ''' This method get as input a path to a file and creates a new one '''
        
    upper_mark = str(mark).upper()
    
    res_path = directory
    
    path_rev = path[::-1] # reverse path
    add_place = len(path) - path_rev.find(".") - 1 # position of last "." in orirginal path
    
    new_path = path[0:add_place] + "_" + upper_mark + path[add_place:]
    
    res_path = res_path + new_path
    
    return res_path

def remove_family_members(values_path,changed_demographic_path):
    ''' This method updated the values file according to the new demographic file'''

    print_log("remove family members - before load")
    
    # Load the data
    changed_demo_full_path = DEST_DIR + changed_demographic_path
    demographic_df = pd.read_excel(changed_demo_full_path,0)
    
    values_full_path = TWINS_DATA + values_path
    values_df = pd.read_excel(values_full_path,0)
    
    print_log("remove family members - after load")
    
    # Extracat sample IDs from both files and sort
    sample_id_lst = list(demographic_df["FlowJo Sample ID"] ) 
    sample_id_lst.sort()
    
    vals_sample_data = list(values_df.columns)
    vals_sample_data.remove("FlowJo Subject ID") #remove string
    vals_sample_data.sort()

    print_log("remove family members - after sort")
    
    # Delete columns of samples which not appear in demographic file
    i = 0 # index on sample_id_lst
    j = 0 # index on vals_sample_data
    sample_ids_to_remove = []
    while (i < len(sample_id_lst) and j < len(vals_sample_data)):
        if (sample_id_lst[i] == vals_sample_data[j]):
            i += 1
            j += 1
        elif (sample_id_lst[i] > vals_sample_data[j]):
            # vals_sample_data[j] doesnt exist in new demographic file
            # so it will be removed from values file
            sample_ids_to_remove.append(vals_sample_data[j])
            j += 1
        elif (sample_id_lst[i] < vals_sample_data[j]):
            i += 1
    
    print_log("before delete")
    # remove marked IDs
    sample_ids_to_remove += vals_sample_data[j:] #rest of samples IDs
    values_df.drop(sample_ids_to_remove,axis = 1, inplace = True)
    
    ''' this is a column by column delete
    for val_sample_id in sample_ids_to_remove:
        values_df.drop(val_sample_id, axis = 1, inplace = True)           
    '''
    new_path = create_new_path(DEST_DIR,values_path, "singeles")
    
    print_log("remove family members - before write to file")
    
    values_df.to_excel(new_path, index = False)
    
    print_log("remove family members - after write to file")
    
    
def print_log(msg):
    ''' prints massage and current time to screen'''
    print(msg + " " + str(dt.datetime.now()))
    
    
def comapre_values_file(before_df, after_df):
    dif_lst = []
    for sample_id in after_df.columns:
        cnt = 0
        for trait_id in after_df.index:
            af_val = after_df[sample_id][trait_id]
            bef_val = before_df[sample_id][trait_id]
            if (af_val == af_val or bef_val == bef_val):
                #one of the items is not NULL
                if af_val != bef_val:
                    #print("for sample " + str(sample_id) + " and trait " + str(trait_id)
                    #+ " the before value is " + str(bef_val) + " and the after val is " + str(af_val))
                    new_var  = (sample_id, trait_id, bef_val,af_val)
                    dif_lst.append(new_var)
                    cnt += 1
        if cnt > 0:
            print("for sample id " + str(sample_id) + " there are " + str(cnt) + " diffs")
            
    return dif_lst
                
                
def run_function_on_united_traits_file(vals_file_path, header ,func, is_numerical, treshold, outuput_file_path):
    ''' This method gets a singles paneled values file and runs the function func
    on its traits. The field is_numberical is here to change type of values of
    the data from float64 to float.
    The result is stored in a list in the following way:
    (trait1,trait2,func_result)'''
    
    total = 0 
    cnt = 0
    
    print_log("Start of run_func - before read file")
    full_path = DEST_DIR + vals_file_path
    result_lst = []
    
    vals_df = pd.read_excel(full_path,0,index_col = header)
    
    print_log("run_func - after read file")
    
    idx_lst = list(vals_df.index)
    n = len(idx_lst)
    print_log("Run_func - start of mega loop")
    for i in range(0, n):
        if i % 1000 == 0:
            print_log("run func - mega loop i = " + str(i))
        trait1_name = idx_lst[i]
        trait1 = vals_df.loc[trait1_name]
        panel1 = trait1["Panel ID"]
        trait1 = trait1.drop("Panel ID")        
        if is_numerical == True:
            trait1 = trait1.astype(float)
        #func_to_run = getattr(type(trait1),func)  
        for j in range(i + 1, n):
            total += 1
            trait2_name = idx_lst[j]
            trait2 = vals_df.loc[trait2_name]
            panel2 = trait2["Panel ID"]
            if panel1 != panel2:
                # different panels - run function
                trait2 = trait2.drop("Panel ID")
                if is_numerical == True:
                    trait2 = trait2.astype(float)
                result = run_function(trait1, trait2, func, treshold,trait1_name,trait2_name)
                if result != None:
                    result_lst.append(result)
                    cnt += 1
                #result = func_to_run(trait1,trait2)
                #if result > treshold:
                    #result_var = (trait1_name, trait2_name, result)
                    #result_lst.append(result_var)
                    #cnt += 1
    
    print_log("run_func - after loop, create dataframe")
    print_log("run_func summary - total checks " + str(total) 
    + " passed treshlod " + str(cnt))
    
    df = pd.DataFrame(data = result_lst, columns = ['trait1', 'trait2', func])
    full_output_path = DEST_DIR + outuput_file_path
    
    print_log("run_func - write file")
    df.to_excel(full_output_path, index = False)
    
    print_log("run_func - after write file")
    

def break_vals_file_to_files_by_panel(vals_path, res_dir):
    ''' This method gets a path to a SINGELES PANELED vals file and creates 
    different files for each panel ID'''
    
    print_log("break vals file - start")    
    full_path = DEST_DIR + vals_path    
    
    # read files by Panel ID as index
    vals_df = pd.read_excel(full_path,0,index_col = "Panel ID")
    
    print_log("break vals file - after file read")    
    
    list_of_panel_id = list(vals_df.index.unique())
    
    print_log("break vals file - before loop")    
    for panel_id in list_of_panel_id:
        print_log("break vals file - in loop for panel id = " + str(panel_id)) 
        #get only the data of specific panel ID
        tmp_df = vals_df.loc[panel_id]
        tmp_df_idx = tmp_df.reset_index()
        #set back to original index and write file
        tmp_df_idx_final = tmp_df_idx.set_index('FlowJo Subject ID')
        new_path = create_new_path(DEST_DIR,vals_path,panel_id)
        tmp_df_idx_final.to_excel(new_path)
        del tmp_df
        del tmp_df_idx
        del tmp_df_idx_final



    
def run_function(x, y, func, treshold, x_name, y_name):
    ''' This method gets 2 variables and names and runs the attribute function func 
    on them. if the result is higher than treshold than returns the value, else None'''
    func_to_run = getattr(type(x),func) 
    result = func_to_run(x,y)
    if result > treshold:
        result_var = (x_name, y_name, result)
        return result_var
    else:
        return None
    


def run_function_on_traits(files_dir, panel1, panel2, func, is_numerical, treshold):
    ''' Gets the directory of values files separated by panel id in SINGELES_PANELED format 
    and runs an attribute function func on each trait of traits of panel1 and panel2'''
    
    total = 0
    cnt = 0
    
    is_same_panel = (panel1 == panel2)    
    
    print_log("Start of run_function_on_traits - before read file")
    result_lst = []
    
    for file in os.listdir(files_dir):
        if file.endswith(panel1 + ".xlsx"):
            path1 = files_dir + file
            file_name = file
        if file.endswith(panel2 + ".xlsx"):
            path2 = files_dir + file
    #print(path1)
    #print("") #\n
    #print(path2)    
    #print("") #\n
    
    df1 = pd.read_excel(path1,0,index_col = "FlowJo Subject ID")
    df1.drop('Panel ID', axis=1, inplace=True)
    if is_same_panel:
        df2 = pd.read_excel(path2,0,index_col = "FlowJo Subject ID")
        df2.drop('Panel ID', axis=1, inplace=True)        
    else:
        df2 = df1 #same panel - same file
        
    print_log("run_function_on_traits - after read file")
    
    idx1_lst = list(df1.index)
    n = len(idx1_lst)
    
    idx2_lst = list(df2.index)
    m = len(idx2_lst)    
    
    print_log("run_function_on_traits - before super loop")
    for i in range(0, n):
        if (i % 200 == 0):
            print_log("run_function_on_traits - in super loop, i = " + str(i))
        trait1_name = idx1_lst[i]
        trait1 = df1.loc[trait1_name]
        if is_numerical == True:
            trait1 = trait1.astype(float)
        # get starting index for calculation in the second file
        if is_same_panel:
            k = i
        else:
            k = 0
        for j in range(k, m):
            total += 1
            trait2_name = idx2_lst[j]
            trait2 = df2.loc[trait2_name]
            if is_numerical == True:
                trait2 = trait2.astype(float)
            result = run_function(trait1, trait2, func, treshold,trait1_name,trait2_name)
            if result != None:
                result_lst.append(result)
                cnt += 1
    
    print_log("run_function_on_traits - after loop, create dataframe")
    print_log("run_function_on_traits summary - total checks " + str(total) 
    + " passed treshlod " + str(cnt))        
    
    df = pd.DataFrame(data = result_lst, columns = ['trait1', 'trait2', func])
    
    full_output_path = create_new_path(files_dir, file_name, 
                                       "_" + panel2 + "_" + func)
    #print(full_output_path)
    print_log("run_function_on_traits - write file")
    df.to_excel(full_output_path, index = False)
    
    print_log("run_function_on_traits - after write file")
        
    

def create_corr_between_all(df):
    ''' this method gets a FD to a united vals and returns dataframe(!!!!) with 
    all the correlations between all traits'''
    
    print_log("create coor - start - before create needed dataframe")
    #create the dataframe in the needed way
    trans = df.transpose()
    trans = trans.drop("Panel ID")
    trans_float = trans.astype('float32')
    
    # run the corr function
    print_log("create coor - before coor")
    corr_res = trans_float.corr()
    
    print_log("create coor - before return")
    
    return corr_res
    
    
def create_tresholded_corr_df(corr_df, treshlod, compare_only_diff_panels ,output_file, print_indicator,is_return):
    ''' This method gets dataframe contains the correlations of all traits 
    and creates a new one with those from different panels 
    with higher correlation than treshold'''
    
    print_log("start of create_tresholded_corr_df - before loop")
    corr_res_lst = []    
    traits = corr_df.index    
    
    
    for i in range(0, len(traits)):
        trait1 = traits[i]
        if i % print_indicator == 0:
            print_log("treshold loop - i is " + str(i))
        for j in range(i, len(traits)):
            trait2 = traits[j]
            is_same_panel = get_panel_id_substr(trait1) == get_panel_id_substr(trait2)
            if ( not(is_same_panel and compare_only_diff_panels) and 
            corr_df[trait1][trait2] > treshlod):
                var = (trait1,trait2,corr_df[trait1][trait2])
                corr_res_lst.append(var)
    
    print_log("create_tresholded_corr_df - after loop, before create final df")
            
    res_df = pd.DataFrame(data = corr_res_lst, columns = ["trait1","trait2", "corr"])
    
    print_log("create_tresholded_corr_df - after final df before write file")
    
    res_df.to_csv(output_file,index = False) #changed to CSV to be limitless
    
    print_log("create_tresholded_corr_df - after write file - finished")
    
    if is_return:
        return res_df


def create_corr_between_two_files(path1, path2):
    ''' Because of Memory Error when running the correlation on all data once this method
    gets 2 paths, combine them to create new DF and sends it into create_corr_between_all'''
    
    print_log("create_corr_between_two_files - start")
    if (path1 == path2):
        # same file, send it directly
        print_log("create_corr_between_two_files - same file")
        df = create_corr_between_all(path1)
    else:
        #create single file
        print_log("create_corr_between_two_files - different files")
        df1 = pd.read_excel(path1, 0, index_col = "FlowJo Subject ID")
        df2 = pd.read_excel(path2, 0, index_col = "FlowJo Subject ID")
        #combine DataFrames
        print_log("create_corr_between_two_files - combine DataFrames")
        df = df1.append(df2, ignore_index = False)
    #now there is a DF, send it's path to create_corr_between_all
    corr_df = create_corr_between_all(df)
        
    print_log("create_corr_between_two_files - before return")
    return corr_df
    

def break_files_into_parts(path, new_amount_of_rows):
    ''' This method takes a file and breaks it into new files, 
    each with new_amount_of_rows rows'''
    
    print_log("break_files_into_parts - start, before read file")    
    
    df = pd.read_excel(path, 0, index_col = "FlowJo Subject ID")
    print_log("break_files_into_parts - after read file, before loop")
    num_of_rows = len(df.index)
    
    i = 0
    cycle = 1
    while i < num_of_rows:
        upper_limit = min(i + new_amount_of_rows, num_of_rows) 
        tmp = df[i:upper_limit]
        rev = path[::-1]
        new_path = path[:(len(path) - rev.find(".")) - 1] + "_" + str(cycle) + path[(len(path) - rev.find(".")) - 1:]
        #path till .suffix + _CYCLE_NUM + suffix
        print_log("break_files_into_parts - before write cycle number " + str(cycle))
        tmp.to_excel(new_path,index = True)
        print_log("break_files_into_parts - after write cycle number " + str(cycle))
        cycle += 1
        i = upper_limit
    
        


def create_correlation_map(path1, path2, is_return,to_overwrite):
    ''' This function does calls function described in create_jobs_functions'''
    
    print_log("create_correlation_map - start before create_corr_between_two_files")
    
    corr_df = create_corr_between_two_files(path1, path2)
    
    panel1 = get_panel_from_path(path1)
    panel2 = get_panel_from_path(path2)
    output_path = OUTPUT_DIR + panel1.upper() + "_" + panel2.upper() + '_corr.csv'
    tresh_output_path = OUTPUT_DIR + panel1.upper() + "_" + panel2.upper() + '_treshold_corr.csv'
    print_log("create_correlation_map - after calculation of panels and path")
    if (panel1 != panel2):
        #different panels, create tresholded correlation map
        print_log("create_correlation_map - different panels, treshold map")
        if ( not(os.path.isfile(tresh_output_path)) or to_overwrite == True):
            #OK overwrite existing file (if there is one..)
            tresh_df = create_tresholded_corr_df(corr_df, treshlod=DEFAULT_TRESHOLD,
                                             compare_only_diff_panels=COMPARE_DIFF_PANELS_ONLY,
                                             output_file=tresh_output_path,
                                             print_indicator=PRINT_INDICATOR,is_return=True)
        #write tresholded correlation map
    
    #write regular correlation map, the treshold data is written in create_tresholded_corr_df
    print_log("create_correlation_map - before write correlation data")
    if ( not(os.path.isfile(output_path)) or to_overwrite == True):
        #OK overwrite existing file (if there is one..)
        corr_df.to_csv(output_path, header = True, index = True)
    
    if is_return:
        return(corr_df,tresh_df)
    

def create_jobs_functions(to_run,to_print,is_return,to_overwrite):
    ''' This function runs on all broken files in INPUT_DIR, calculates the correlation 
    between all traits in each two files using create_corr_between_two_files. 
    Then, for the result, in case the panels are different, it calculates
    also the tresholded correlatin using create_tresholded_corr_df. 
    Both of files are written into OUTPUT_DIR. 
    to_run and to_print are 2 booleans to indicate wheter to run or print (or both)
    the commands'''
    
    print_log("create_jobs_functions - start")
    input_files = os.listdir(INPUT_DIR)
    n = len(input_files) #amount of files
    for i in range(0, n):
        print_log("create_jobs_functions - i = " + str(i) + " n = " + str(n))
        path1 = INPUT_DIR + input_files[i]
        for j in range(i, n):
            path2 = INPUT_DIR + input_files[j]
            if to_print:
                print("\n")
                print("create_correlation_map(r'" + path1 + "',r'" + 
                path2 + "'," + str(is_return) + "," + str(to_overwrite) +")")
                print("\n")
            if to_run:
                create_correlation_map(path1, path2,is_return,to_overwrite)
            
    
    
def get_panel_from_path(singeles_paneled_broken_filename):
    ''' This method gets a filename of a singles broken filename (INPUT_DIR files)
    and returns the panel of the data inside in file concataneted to file number'''
    
    rev_path = singeles_paneled_broken_filename[::-1]
    #find the last "_" and the one before it which wraps the panel
    rev_start_pos = rev_path.find(".")
    rev_end_pos = rev_path.find("_") # first accurence in reverse
    rev_end_pos = rev_path.find("_",rev_end_pos + 1) # second accurence in reverse
    rev_panel = rev_path[rev_start_pos + 1 : rev_end_pos]
    ret_val = rev_panel[::-1] #reverse back
    ret_val = ret_val.replace("_","") # remove "_" 
    return ret_val


def main():
    # get parameters from user
    args = sys.argv
    filename1 = args[1]
    filename2 = args[2]
    is_return = args[3]
    to_overwrite = args[4]
    
    path1 = INPUT_DIR + filename1
    path2 = INPUT_DIR + filename2
    
    create_correlation_map(path1, path2, is_return,to_overwrite)
    


if __name__ == "__main__":
    main()
    