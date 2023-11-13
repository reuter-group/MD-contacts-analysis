#!/usr/bin/env python

### Hbond post-processing script
# this script will extract inforamtion from the output file of hbond analysis

import numpy as np 
import sys


## The lipid types used in your hbond analysis
LIPIDS = ['POPC', 'POPE', 'POPI', 'POPS', 'CER1', 'CHL1'] 


def exclude_lipids(hbond_file, exc_lipids, total_frame, output_lipid_contact=True):
    '''
    Report hbond occupancy, contacts information but excluding the specified lipids
    :: exc_lipids: ["RESNAME-ID", ...]
    :: total_frame: the number of frames in DCD file  used for the analysis
    :: output_lipid_contact: set False if you do not need that information
    '''
    import sys 
    sys.path.insert(1, '../')
    from hbond_candi import ATOMS_DICT

    with open(hbond_file) as f:    
        contents = f.readlines()        

    new_file_lines=[contents[0]]
    for line in contents[1:]:
        hbond = line.strip().split()
        l_index = 3 if hbond[3] in LIPIDS else 7                
        lipid = hbond[l_index] + "-" + hbond[l_index+1]
        if lipid not in exc_lipids: 
            new_file_lines.append(line)

    # writ to a new hbond file
    new_hbond_file = hbond_file[:-4]+"_new.txt"
    with open(new_hbond_file, 'w') as new_file:
        new_file.writelines(new_file_lines) 

    # output information for each amino acid
    line_count = len(new_file_lines)
    for i, line in enumerate(new_file_lines):    
        if i == 0: # init, 1st line is useless
            pho_dict, head_dict, gly_dict = [{k: [] for k in LIPIDS} for i in range(3)]
            continue 

        hbond = line.strip().split()
        if i == 1:
            res = hbond[0] 
        
        if res != hbond[0]:
            report_occupancy(res, pho_dict, gly_dict, head_dict, total_frame)
            report_contact(res, pho_dict, gly_dict, head_dict, line_count-1)

            if output_lipid_contact:
                report_lipid_percent([pho_dict, gly_dict, head_dict], 
                                    ["PHO", "GLY", "HEAD"])
            res = hbond[0]
            pho_dict, head_dict, gly_dict = [{k: [] for k in LIPIDS} for i in range(3)]

        l_index = 3 if hbond[3] in LIPIDS else 7                
        lipid = hbond[l_index] + "-" + hbond[l_index+1]
        lipid_name = hbond[l_index]
        lipid_atom = hbond[l_index+2]

        frame = int(hbond[-1])                    
        if lipid_atom in ATOMS_DICT[lipid_name]['pho'].split():
            pho_dict[lipid_name].append(frame)
        elif lipid_atom in ATOMS_DICT[lipid_name]['gly'].split():
            gly_dict[lipid_name].append(frame)
        elif lipid_atom in ATOMS_DICT[lipid_name]['head'].split():
            head_dict[lipid_name].append(frame)
        else:
            print("WARNING: found contact from parts other than PHO, GLY and HEAD")
    
        if i == line_count-1: # last line 
            report_occupancy(res, pho_dict, gly_dict, head_dict, total_frame)  
            report_contact(res, pho_dict, gly_dict, head_dict, line_count-1)
            if output_lipid_contact:
                report_lipid_percent([pho_dict, gly_dict, head_dict], 
                                    ["PHO", "GLY", "HEAD"])

    return new_hbond_file

def report_occupancy(res, pho_dict, gly_dict, head_dic, total_frames):  
    pho_num = len(set(sum(pho_dict.values(), []))) 
    pho_occ = pho_num / total_frames * 100  

    gly_num = len(set(sum(gly_dict.values(), []))) 
    gly_occ = gly_num / total_frames * 100

    head_num = len(set(sum(head_dic.values(), []))) 
    head_occ = head_num / total_frames * 100
    
    res_num = len(set(sum(list(pho_dict.values()) + list(gly_dict.values()) + list(head_dic.values()), []))) 
    res_occ = res_num / total_frames * 100
    print("\n{} Occupancy count {}, {}/{} = {:.2f}% ; PHO {:.2f}%, GLY {:.2f}%, HEAD {:.2f}%".format(
        res, res_num, res_num, total_frames, res_occ, pho_occ, gly_occ, head_occ))


def report_contact(res, pho_dict, gly_dict, head_dic, total_contacts):
    pho_con = len(sum(pho_dict.values(), []))
    pho_per = pho_con / total_contacts * 100  

    gly_con = len(sum(gly_dict.values(), []))
    gly_per = gly_con / total_contacts * 100

    head_con = len(sum(head_dic.values(), []))
    head_per = head_con / total_contacts * 100
    
    res_con = pho_con + gly_con + head_con
    res_per = res_con / total_contacts * 100

    print("{} Contacts count {}, {}/{} = {:.2f}% ; PHO {:.2f}%, GLY {:.2f}%, HEAD {:.2f}%".format(
        res, res_con, res_con, total_contacts, res_per, pho_per, gly_per, head_per))


def report_lipid_percent(group_dicts, group_names):
    for group_dict, group_name in zip(group_dicts, group_names):
        all_indices = sum(group_dict.values(), [])
        occ_count = len(set(all_indices))
        con_count = len(all_indices)

        if con_count != 0:
            print("  - {}: \t occupancy count {} \t;\t contact count {} ".format(
                group_name, occ_count, con_count))

            for lipid, indices in group_dict.items():
                output = ""
                print("    -- {}: occupancy\t{}\t{:.2f}%\t;\tcontacts\t{}\t{:.2f}%".format(
                    lipid, len(set(indices)), len(set(indices))/occ_count * 100,
                           len(indices), len(indices)/con_count * 100))                
                    

def report_statistics(hbond_file, total_frame):
    """
    Output some statistics information from the given hbond file
    also write the contact number of each lipid part to "frame_contacts.txt" 
    """
    import sys 
    sys.path.insert(1, '../')
    from hbond_candi import ATOMS_DICT

    with open(hbond_file) as f:    
        contents = f.readlines()

    pho_dict, head_dict, gly_dict = [{k: [] for k in LIPIDS} for i in range(3)]
    frame_info = np.zeros((total_frame, 3)) # [[pho_count, gly_count, head_count],...]

    for line in contents[1:]:
        hbond = line.strip().split()
        l_index = 3 if hbond[3] in LIPIDS else 7                
        lipid_name = hbond[l_index]
        lipid_atom = hbond[l_index+2]

        frame = int(hbond[-1])  
        if lipid_atom in ATOMS_DICT[lipid_name]['pho'].split():
            pho_dict[lipid_name].append(frame)
            frame_info[frame][0] += 1
        elif lipid_atom in ATOMS_DICT[lipid_name]['gly'].split():
            gly_dict[lipid_name].append(frame)
            frame_info[frame][1] += 1
        elif lipid_atom in ATOMS_DICT[lipid_name]['head'].split():
            head_dict[lipid_name].append(frame)
            frame_info[frame][2] += 1
        else:
            print("WARNING: contacts from parts other than PHO, GLY and HEAD")
    
    print("\n\n **** TOTAL lines of contact: {}".format( len(contents)-1))
    print("--------------------------------------------------------------")
    print("\t" + "\t".join(LIPIDS)+"\tTotal")
    for d, name in [(pho_dict, "PHO"), (gly_dict, "GLY"), (head_dict, "HEAD")]:
        output = ""
        total = 0
        for lipid, frames in d.items():
            output += "\t" + str(len(frames))
            total += len(frames)
        print("{}{}\t{}".format(name, output, total))
    print("--------------------------------------------------------------")

    print("Average number of contacts per frame: {}/{}={:.2f}".format(
        len(contents)-1, total_frame, (len(contents)-1)/(total_frame)))

    indices = np.arange(total_frame).reshape(total_frame, 1)
    frame_output = np.concatenate((indices, frame_info), axis=1)
    np.savetxt("frame_contacts.txt", frame_output, fmt='%d', 
                    header="index,PHO,GLY,HEAD")


if __name__ == "__main__":
    ## Set these two arguments 
    hbond_file = "out_hbonds.txt"
    total_frame = 22861

    # This line for Feature 1, 2
    # will output hbonds without the specified lipids to a new file
    new_hbond_file = exclude_lipids(hbond_file,                                    
                                    ["POPS-150", "POPC-212"], # set excluded lipids here
                                    total_frame, 
                                    output_lipid_contact=True # set this False to show less info
                                    )
    
    # This line will give you information including Feature 3
    # contacts in each frame will be output to file "frame_contacts.txt"
    report_statistics(new_hbond_file,  # this uses the newly generated hbond file
                        total_frame)

    # shortcut: use command-line tool to exclude specific lipid
    #$ grep -v " POPC 158 " out_hbonds.txt > out_hbonds_new.txt
    