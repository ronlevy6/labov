library("lmtest") 
library("nlme")
library("stringr")

#read data
united_traits = read.delim('C:/Ron/tmp/r_files_tresh_08/united_08.txt', header = TRUE)
#traits_with_age = read.csv('C:/Ron/tmp/r_files_tresh_08/traits_file_for_R_tresh_08.csv')
#add -1 to NULLS
traits_with_age[is.na(traits_with_age)] <- -1

trait1_vec = as.vector(united_traits$trait1)
trait2_vec = as.vector(united_traits$trait2)
n = length(trait1_vec)

for (i in c(j:n)){
  trait1 = trait1_vec[i]
  trait2 = trait2_vec[i]
  #fix data to fit column name
  trait1 = str_replace_all(trait1,":",".")
  trait1 = str_replace_all(trait1," ",".")
  
  trait2 = str_replace_all(trait2,":",".")
  trait2 = str_replace_all(trait2," ",".")
  
  #create subset for linear model
  subset_for_linear = as.data.frame(cbind(traits_with_age$Age,traits_with_age[trait1]))
  colnames(subset_for_linear) = c("Age", "trait1") 
  
  linear_model = lm(trait1~Age,data=subset_for_linear)
  
  #create subset for mixed model
  subset_for_mixed = as.data.frame(cbind(traits_with_age$Age,traits_with_age[trait1],traits_with_age[trait2]))
  colnames(subset_for_mixed) = c("Age", "trait1","trait2") 
  mixed_model = lme(trait1~Age,data=subset_for_mixed, random=~1|trait2, control=lmeControl(returnObject=TRUE))
  
  
  test_result = lrtest(linear_model,mixed_model)
  
  p_val = test_result$`Pr(>Chisq)`
  needed_p_val = p_val[2]
  var_to_keep = cbind(trait1,trait2,needed_p_val)
  if (i == 1){
    ret_val = var_to_keep
  }
  else{
    ret_val = rbind(ret_val, var_to_keep)
  }
}
ret_df = as.data.frame(ret_val)
write.csv(x = ret_df,file = 'c:/ron/tmp/r_files_tresh_08/r_output_tresh_08.csv')

