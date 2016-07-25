import os, sys

def check_nan(param):
        #Check there is nan which is not a valid value.
        #It has nan, change to 99999 value. Otherwise, return param value
        if str(param).find("nan") >= 0:
                return "9999999"
        else:
                return param
        
def main():
    obj_1_file_name = "Gyrate_0.fit"
    obj_2_file_name = "Hydrophilic_0.fit"
    obj_3_file_name = "Potential_0.fit"

    obj_5_file_name = "Hydrophobic_0.fit"
    obj_6_file_name = "H_Bond_Main_0.fit"
    obj_7_file_name = "Electrostatic_0.fit"
    obj_8_file_name = "Van_der_Waals_0.fit"
    obj_9_file_name = "GBSA_Solvatation_0.fit"

    path_fit = sys.argv[1] #path of files
    number_individuals = int(sys.argv[2]) 
    number_generation = int(sys.argv[3])
    number_objectives = int(sys.argv[4])
    sufix_output_file = sys.argv[5]

    path_file_obj_1 = os.path.join(path_fit,obj_1_file_name)
    path_file_obj_2 = os.path.join(path_fit,obj_2_file_name)
    path_file_obj_3 = os.path.join(path_fit,obj_3_file_name)

    path_file_obj_5 = os.path.join(path_fit,obj_5_file_name)
    path_file_obj_6 = os.path.join(path_fit,obj_6_file_name)
    path_file_obj_7 = os.path.join(path_fit,obj_7_file_name)
    path_file_obj_8 = os.path.join(path_fit,obj_8_file_name)
    path_file_obj_9 = os.path.join(path_fit,obj_9_file_name)

    #Loading file Objective 1
    dic_obj_1 = {}
    g = 0
    ind = 1    
    f_1 = open(path_file_obj_1,"r")    
    for line in f_1:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_1[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_1.close()

    #Loading file Objective 2
    dic_obj_2 = {}
    g = 0
    ind = 1    
    f_2 = open(path_file_obj_2,"r")    
    for line in f_2:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_2[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_2.close()

    #Loading file Objective 3
    dic_obj_3 = {}
    g = 0
    ind = 1    
    f_3 = open(path_file_obj_3,"r")    
    for line in f_3:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_3[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_3.close()

    #Objective 4 was removed. It was based on Stride program

    #Loading file Objective 5
    dic_obj_5 = {}
    g = 0
    ind = 1    
    f_5 = open(path_file_obj_5,"r")    
    for line in f_5:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_5[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_5.close()

    #Loading file Objective 6
    dic_obj_6 = {}
    g = 0
    ind = 1    
    f_6 = open(path_file_obj_6,"r")    
    for line in f_6:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_6[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_6.close()

    #Loading file Objective 7
    dic_obj_7 = {}
    g = 0
    ind = 1    
    f_7 = open(path_file_obj_7,"r")    
    for line in f_7:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_7[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_7.close()

    #Loading file Objective 8
    dic_obj_8 = {}
    g = 0
    ind = 1    
    f_8 = open(path_file_obj_8,"r")    
    for line in f_8:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_8[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_8.close()

    #Loading file Objective 9
    dic_obj_9 = {}
    g = 0
    ind = 1    
    f_9 = open(path_file_obj_9,"r")    
    for line in f_9:
        if line.find("#") < 0:
            key = str(str(g) + '_'+ str(ind))
            dic_obj_9[key] = check_nan(line.split()[1])
            ind = ind + 1
    f_9.close()


    #Saving files for each generation
    #Here all objectives will be join into one column for each generation    
    path_file_out = os.path.join(path_fit,sufix_output_file+str(g)+'_'+str(number_individuals)+'.txt')
    f_out = open(path_file_out,"w")
    f_out.write('n= '+str(number_individuals)+'\n')
    f_out.write('k= '+str(number_objectives)+'\n')
    ind = 1
    while ind < number_individuals:          
       key = str(str(g) + '_'+ str(ind))
       line_value = str(ind) + ' '+ dic_obj_1[key] + ' '+ dic_obj_2[key]+' '+ dic_obj_3[key]+' '+ dic_obj_5[key]+' '+ dic_obj_6[key]+' '+ dic_obj_7[key]+' '+ dic_obj_8[key]+' '+ dic_obj_9[key] +'\n'
       f_out.write(line_value) 
       ind = ind + 1

    f_out.write('EOF') 
    f_out.close()


main()
