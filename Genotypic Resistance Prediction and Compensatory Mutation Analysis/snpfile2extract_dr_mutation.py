#!/usr/local/bin/python3
import sys
### prepair resistance source database
## add index ID （2726196_C/T）
dr_list_single = [] # single resistance mutation database
dr_list_multi = [] # multiple resistance mutations database
with open(sys.argv[1]) as f1:
    next(f1)
    for line in f1:
        line=line.strip().split('\t')
        if line[2] == '-': # drop useless data
            pass
        else:
            # split the dr_list into single and multiple mutation list
            if len(line[2].split('/')) == 1:
                line.append(line[2]+"_"+line[4])
                dr_list_single.append(line)
            else:
                line[2] = line[2].split('/')
                num_muti_mutation = len(line[2])
                line[4] = [line[4].split('/')[2*i]+"/"+line[4].split('/')[2*i+1] for i in range(num_muti_mutation)]
                line.append([line[2][i]+"_"+line[4][i] for i in range(num_muti_mutation)])
                dr_list_multi.append(line)

## add gene name （ahpC P2S）
for line in dr_list_single: # add gene name into single mutation list
    if len(line)==6+1:
        name_row = line[1] + " " + line[4].split('/')[0]+line[3]+line[4].split('/')[1]
    elif len(line)==10+1:
        if line[7] == '-' or line[5] == '-':
            name_row = line[1] + " " + line[4].split('/')[0]+line[3]+line[4].split('/')[1]
        else:
            name_row = line[1] + " " + line[7].split('/')[0]+line[5]+line[7].split('/')[1]
    line.append(name_row)
for line in dr_list_multi: # add gene name into multiple mutations list
    name_row = line[1] + " " + line[7].split('/')[0]+line[5]+line[7].split('/')[1]
    line.append(name_row)
    
### prepair SNP file
## add index ID
SNP_list = {}
with open(sys.argv[2]) as f2:
    for line in f2:
        line=line.strip().split('\t')  
        line.append(line[0]+"_"+line[1]+"/"+line[2]) # add index ID
        SNP_list[line[-1]] = line
        
### extract resistance mutation based on source database
## generate index ID list
dr_list_single_ID = [line[-2] for line in dr_list_single] # add single resistance infor
for SNP_ID in SNP_list.keys():
    ID_matches = [i for i,x in enumerate(dr_list_single_ID) if x==SNP_list[SNP_ID][-1]]
    if len(ID_matches) != 0:
        drug_matched = []
        for ID_match in ID_matches:
            name_matched = dr_list_single[ID_match][-1]
            drug_matched.append(dr_list_single[ID_match][0])
        SNP_list[SNP_ID].append('single')
        SNP_list[SNP_ID].append(name_matched)
        SNP_list[SNP_ID].append("/".join(drug_matched))

for line in dr_list_multi:
    flag_all_detected = True
    for line_ID in line[-2]:
        if line_ID not in SNP_list.keys():
            flag_all_detected = False
    if flag_all_detected is True:
        for line_ID in line[-2]:
            SNP_list[line_ID].append('multiple')
            SNP_list[line_ID].append(line[-1])
            SNP_list[line_ID].append(line[0])

for lines in SNP_list.values():
    if len(lines) != 4:
       	#print(sys.argv[1].strip('.snp'), end = '\t') 
        for number in [0,1,2,4,5,6]:
            print(lines[number], end = '\t')
        print('')
