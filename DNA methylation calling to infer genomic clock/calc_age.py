#!/usr/bin/env python3

# just libaries
import sys
import pandas as pd 
from copy import deepcopy

def calc_age():
    # get the file path
    file_path_bismark = sys.argv[1]
    file_path_weights = sys.argv[2]

    # using pandas to read in the file and
    # create columns for methylation ratio and methylation percentage
    bisOrig = pd.read_csv(file_path_bismark, sep='\t', header=None)
    bisOut = deepcopy(bisOrig)
    bisOut = bisOut.drop(columns=[0,1])

    header = ['CpG_loc','Meth_Percentage','methC', 'unmethC']
    bisOut.columns = header[:len(bisOut.columns)]
    bisOut = bisOut.set_index("CpG_loc")

    bisOut['Meth_Ratio'] = bisOut['Meth_Percentage']/100
    bisOut['Total_Reads'] = bisOut['methC'] + bisOut['unmethC']

    # uncomment this to see the table
    #display(bisOut.head())

    # loading parameters weights
    rDNAloci_original = pd.read_csv(file_path_weights, header=None)
    rDNAloci_original.columns = ['CpG_loc', 'weight']
    rDNAloci = deepcopy(rDNAloci_original)

    # uncomment this to see weights
    #print(rDNAloci.head())

    # adjust CpG loci
    for i in range(len(rDNAloci)):
        rDNAloci.loc[i,'CpG_loc'] = int(rDNAloci['CpG_loc'][i].split(',')[1])
        if rDNAloci['CpG_loc'][i] < 0:
            rDNAloci.loc[i,'CpG_loc'] += 501
        else: 
            rDNAloci.loc[i,'CpG_loc'] += 500

    # uncomment below to display table
    #print(rDNAloci.head())


    # merge weights and methylation tables
    clockdtfr = pd.merge(rDNAloci, bisOut, on='CpG_loc', how='inner')

    # uncomment blew to display table
    #print(clockdtfr.head())

    # cacluate the ages
    intercept= -2.223915739
    prod=0
    for i in range(len(clockdtfr)):
        prod+=(clockdtfr['weight'][i]*clockdtfr['Meth_Ratio'][i])
    rDNAage=2**(intercept+prod)

    print(file_path_bismark, "has age:", rDNAage)

    # write it out to a file
    with open("./output_age.txt", "w") as  file:
        file.write(f"{file_path_bismark} has rDNA age: {rDNAage}")

if __name__ == "__main__":
    calc_age()











