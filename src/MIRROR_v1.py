import pandas as pd
from collections import defaultdict as dd
from subprocess import Popen,PIPE
import sys
import os
import argparse
from time import strftime, localtime
from multiprocessing import Pool
import multiprocessing as mp
import RNA
import numpy as np
from Bio import SeqIO
from itertools import groupby,count
from collections import defaultdict as dd
from typing import Literal, get_args
import itertools
from copy import deepcopy
from pathlib import Path
_bases = {"A", "T", "C", "G", "N", "*", " ", "U", "a", "t", "c", "g", "n", "u"}
_modes = Literal["raw", "A-A", "A-C", "CA-AC",'GA-GC',"G-G","GAP-GCP"]
SCRIPT_DIR = Path(__file__).resolve().parent
SUBSTRATES_DIR = SCRIPT_DIR.parent / 'substrates'
RNAHYBRID_BIN = SCRIPT_DIR / 'RNAhybrid-2.1.2' / 'bin' / 'bin' / 'RNAhybrid'


class Hybrid_str():
    """Parse transformed RNAhybrid-output(Editing site was marked as * in es_pair or es_unpair)

    Args:
        ecs_unpair: 5' unpaired Editing complementary sequence 3'
        ecs_pair  : 5' paired Editing complementary sequence   3'
        es_pair   : 3' paired Editing sequence                 5'
        es_unpair : 3' unpaired Editing sequence               5'

    Returns:
        a Hybrid_str object to get base-pair information."""

    def __init__(self, ecs_unpair: str, ecs_pair: str, es_pair: str, es_unpair: str):

        for input_str in [ecs_unpair, ecs_pair, es_pair, es_unpair]:
            if isinstance(
                    input_str,
                    str) == False or set(input_str).issubset(_bases) == False:
                raise ValueError(
                    "Invalid input! Input str should only contain A,C,T,G,N,*,' '(space),"
                )

        if len(ecs_unpair) == len(ecs_pair) == len(es_pair) == len(es_unpair):
            pass
        else:
            raise ValueError("Input should be same length!")
        ecs_unpair = ecs_unpair.upper().replace("T","U")
        ecs_pair = ecs_pair.upper().replace("T","U")
        es_pair = es_pair.upper().replace("T","U")
        es_unpair = es_unpair.upper().replace("T","U")
        self.ecs_unpair = ecs_unpair
        self.ecs_pair = ecs_pair
        self.es_pair = es_pair
        self.es_unpair = es_unpair
        bp_dict = dd(dict)

        ###bp annotation
        for i in range(len(ecs_pair)):
            if ecs_pair[i] != " ":
                if (es_pair[i] == "U" and ecs_pair[i] == "G") or (es_pair[i] == "G" and ecs_pair[i] == "U"):
                    bp_dict[i]["info"] = "wobble_paired"
                else:
                    bp_dict[i]["info"] = "paired"
                bp_dict[i]["ecs_base"] = ecs_pair[i]
                bp_dict[i]["es_base"] = es_pair[i]
            elif ecs_pair[i] == " ":
                bp_dict[i]["info"] = "unpaired"
                bp_dict[i]["ecs_base"] = ecs_unpair[i]
                bp_dict[i]["es_base"] = es_unpair[i]
        unpaired_locs = [key for key,value in bp_dict.items() if value["info"] == "unpaired"]
        grouped_locs = [list(g) for _,g in groupby(unpaired_locs,lambda x,c=count(): x-next(c))]
        for g in grouped_locs:
            if len(g) == 1:
                loc = g[0]
                if bp_dict[loc]["es_base"]  != " " and bp_dict[loc]["ecs_base"] == " ":
                    bp_dict[loc]["info"] = "es_bulge"
                elif bp_dict[loc]["es_base"] == " " and bp_dict[loc]["ecs_base"] != " ":
                    bp_dict[loc]["info"] = "ecs_bulge"
                else:
                    bp_dict[loc]["info"] = "loop"
            elif len(g) > 1:
                es_bases = [bp_dict[i]["es_base"] for i in g]
                ecs_bases = [bp_dict[i]["ecs_base"] for i in g]
                for loc in g:
                    if es_bases.count(" ") == len(g):
                        bp_dict[loc]["info"] = "ecs_bulge"
                    elif ecs_bases.count(" ") == len(g):
                        bp_dict[loc]["info"] = "es_bulge"
                    else:
                        bp_dict[loc]["info"] = "loop"
        self.bp_dict = bp_dict
                

    def __iter__(self):
        return iter(
            [self.ecs_unpair, self.ecs_pair, self.es_pair, self.es_unpair])

    def __eq__(self,other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __len__(self):
        return len(self.ecs_pair)

    def show(self, out=False,reverse=False):
        """show all base pair information of Hybrid_str"""
        if reverse == False:
            out_s1 = "5'-" + self.ecs_unpair + "-3'"
            out_s2 = "5'-" + self.ecs_pair + "-3'"
            out_s3 = "3'-" + self.es_pair + "-5'"
            out_s4 = "3'-" + self.es_unpair + "-5'"
        elif reverse == True:
            out_s1 = "3'-" + self.ecs_unpair[::-1] + "-5'"
            out_s2 = "3'-" + self.ecs_pair[::-1] + "-5'"
            out_s3 = "5'-" + self.es_pair[::-1] + "-3'"
            out_s4 = "5'-" + self.es_unpair[::-1] + "-3'"   
        if out == True:
            return "\n".join(
                [out_s1, out_s2, out_s3, out_s4])
        elif out == False:
            print("\n".join([out_s1, out_s2, out_s3, out_s4]))

    def get_edit_site_index(self):
        """return 0_based coordinate of editing site(marked as * in es_pair or es_unpair"""
        es_idx = 0
        try:
            try:
                es_idx = self.es_pair.index("*")
            except:
                es_idx = self.es_unpair.index("*")
        except ValueError:
            raise ValueError("Editing site was not marked as * in input!")
        return es_idx

    def get_edit_site_str(self):
        """return base pair information of editing site"""
        es_idx = self.get_edit_site_index()
        return self.get_base_pair_by_loc(es_idx)

    def offset_to_loc(self, five_offset=0, three_offset=0):
        """return 0 based coordinate of input base offset to editing sites
        
        Args:
            five_offset : es base offset upstream to editing sites.
            three_offset: es base offset downstream to editing sites.

        Returns:
            0_based coordinate of offset on the input str!
            (three_offset_loc,five_offset_loc)"""
        es_idx = self.get_edit_site_index()
        right_idx = self.get_edit_site_index()
        left_idx = self.get_edit_site_index()
        right_n = 0
        left_n = 0

        if isinstance(five_offset, int) == False or isinstance(
                three_offset, int) == False:
            raise TypeError("Input should be int!")

        if three_offset > 0:  #downstream to editing site!
            for i in range(es_idx - 1, -1, -1):
                if self.es_pair[i] != " " or self.es_unpair[i] != " ":
                    left_n += 1
                if left_n == three_offset:
                    left_idx = i
                    break
            if left_n != three_offset:
                raise ValueError("Input three_offset is beyond the hybrid_str!")

        if five_offset > 0:  #upstream to editing site!
            for i in range(es_idx + 1, self.__len__()):
                if self.es_pair[i] != " " or self.es_unpair[i] != " ":
                    right_n += 1
                if right_n == five_offset:
                    right_idx = i
                    break
            if right_n != five_offset:
                raise ValueError("Input five_offset is beyond the hybrid_str!")

        return [left_idx, right_idx]

    def get_base_pair_by_loc(self, loc_0: int):
        """return base pair informatin based on the input 0_based coordinate"""
        ecs_base = " "
        es_base = " "

        try:
            es_base = self.bp_dict[loc_0]["es_base"]
            ecs_base = self.bp_dict[loc_0]["ecs_base"]
            info = self.bp_dict[loc_0]["info"]
        except:
            raise KeyError("Input loc_0 is beyond the hybrid_str!")

        return "-".join([es_base, ecs_base, info])

    def get_base_pair_by_three_offset(self, three_offset: int):
        """get base pair information by 1_based base five_offset to editing site

        Args:
            three_offset: base span downstream to editing sites.

        Returns:
            Base pair information"""
        left_idx, _ = self.offset_to_loc(three_offset=three_offset)
        return self.get_base_pair_by_loc(left_idx)

    def get_base_pair_by_five_offset(self, five_offset: int):
        """get base pair information by 1_based base five_offset to editing site

        Args:
            five_offset: base span upstream to editing sites.

        Returns:
            Base pair information"""
        _, right_idx = self.offset_to_loc(five_offset=five_offset)
        return self.get_base_pair_by_loc(right_idx)

    def truncated_by_es_base_offset(self,
                                    five_offset=0,
                                    three_offset=0,
                                    mode: _modes = "raw"):
        """return a truncated Hybrid_str object based on input base's 1_based five_offset or 1_based three_offset
        
        Args:
            five_offset : base offset upstream to editing site,
            three_offset: base offset downstream to editing site,

        Returns:
            a truncated Hybrid_str object."""
        assert mode in get_args(
            _modes
        ), f"'{mode}' is invalid - valid options are {get_args(_modes)}"

        if isinstance(five_offset, int) == False or isinstance(
                three_offset, int) == False:
            raise ValueError("Input offset should be int!")
        else:
            if five_offset < 0 or three_offset < 0:
                raise ValueError("Input offset should >= 0!")

        # print(self.offset_to_loc(five_offset=five_offset,three_offset=three_offset))
        left_idx, right_idx = self.offset_to_loc(five_offset=five_offset,
                                                 three_offset=three_offset)
        es_idx = self.get_edit_site_index()
        if mode == "raw":
            new_s1 = self.ecs_unpair[left_idx:right_idx + 1]
            new_s2 = self.ecs_pair[left_idx:right_idx + 1]
            new_s3 = self.es_pair[left_idx:right_idx + 1]
            new_s4 = self.es_unpair[left_idx:right_idx + 1]
        elif mode == "A-A":
            new_s1 = self.ecs_unpair[left_idx:es_idx] + "A" + self.ecs_unpair[es_idx + 1:right_idx + 1]
            new_s2 = self.ecs_pair[left_idx:es_idx] + " " + self.ecs_pair[es_idx + 1:right_idx + 1]
            new_s3 = self.es_pair[left_idx:es_idx] + " " + self.es_pair[es_idx + 1:right_idx + 1]
            new_s4 = self.es_unpair[left_idx:es_idx] + "*" + self.es_unpair[es_idx + 1:right_idx + 1]
        elif mode == "A-C":
            new_s1 = self.ecs_unpair[left_idx:es_idx] + "C" + self.ecs_unpair[es_idx + 1:right_idx + 1]
            new_s2 = self.ecs_pair[left_idx:es_idx] + " " + self.ecs_pair[es_idx + 1:right_idx + 1]
            new_s3 = self.es_pair[left_idx:es_idx] + " " + self.es_pair[es_idx + 1:right_idx + 1]
            new_s4 = self.es_unpair[left_idx:es_idx] + "*" + self.es_unpair[es_idx + 1:right_idx + 1]
        elif mode == "G-G": 
            if  five_offset < 1:
                raise ValueError(
                    "five offest shoud > 1 for G-G mode")
            else:
                if five_offset == 1:
                    new_s1 = self.ecs_unpair[left_idx:es_idx + 1] + "G" 
                    new_s2 = self.ecs_pair[left_idx:es_idx + 1] + " " 
                    new_s3 = self.es_pair[left_idx:es_idx + 1] + " " 
                    new_s4 = self.es_unpair[left_idx:es_idx + 1] + "G" 
                else:
                    three_n = 0
                    for i in range(es_idx+2, len(self.es_unpair)):
                        if self.es_unpair[i] != " " or self.es_pair[i] != " ":
                            three_n += 1
                        if three_n == five_offset - 1:
                            three_idx = i
                            break
                    if three_idx != None:
                        new_s1 = self.ecs_unpair[left_idx:es_idx + 1] + "G" + self.ecs_unpair[es_idx+2:three_idx+1]
                        new_s2 = self.ecs_pair[left_idx:es_idx + 1] + " " + self.ecs_pair[es_idx+2:three_idx+1]
                        new_s3 = self.es_pair[left_idx:es_idx + 1] + " " + self.es_pair[es_idx+2:three_idx+1]
                        new_s4 = self.es_unpair[left_idx:es_idx + 1] + "G" + self.es_unpair[es_idx+2:three_idx+1]
                    else:
                        raise ValueError("The current hs do not support.")
        elif mode == "GA-GC":
            if  five_offset < 1:
                raise ValueError(
                    "five offest shoud > 1 for GA-GC mode")
            else:
                if five_offset == 1:
                    new_s1 = self.ecs_unpair[left_idx:es_idx] + "CG"  
                    new_s2 = self.ecs_pair[left_idx:es_idx] + "  " 
                    new_s3 = self.es_pair[left_idx:es_idx] + "  " 
                    new_s4 = self.es_unpair[left_idx:es_idx] + "*G" 
                else:
                    three_n = 0
                    for i in range(es_idx+2, len(self.es_unpair)):
                        if self.es_unpair[i] != " " or self.es_pair[i] != " ":
                            three_n += 1
                        if three_n == five_offset - 1:
                            three_idx = i
                            break
                    if three_idx != None:
                        new_s1 = self.ecs_unpair[left_idx:es_idx] + "CG" + self.ecs_unpair[es_idx+2:three_idx+1]
                        new_s2 = self.ecs_pair[left_idx:es_idx] + "  " + self.ecs_pair[es_idx+2:three_idx+1]
                        new_s3 = self.es_pair[left_idx:es_idx] + "  " + self.es_pair[es_idx+2:three_idx+1]
                        new_s4 = self.es_unpair[left_idx:es_idx] + "*G" + self.es_unpair[es_idx+2:three_idx+1]
                    else:
                        raise ValueError("The current hs do not support.")
        elif mode == "GAP-GCP":
            if three_offset < 1 or five_offset < 1:
                raise ValueError(
                    "Three offset value should be >= 2 and five offset value should be >= 1for GAP-GCP mode!")
            else:
                if five_offset == 1:
                    if three_offset == 1:
                        new_s1 =  " CG"  
                        new_s2 =  "N  " 
                        new_s3 =  "N  " 
                        new_s4 =  " *G" 
                    else:
                        five_n = 0
                        for i in range(es_idx-2, -1, -1):
                            if self.es_unpair[i] != " " or self.es_pair[i] != " ":
                                five_n += 1
                            if five_n == three_offset - 1:
                                five_idx = i
                                break
                        if five_idx != None:
                            new_s1 = self.ecs_unpair[five_idx:es_idx - 1] + " CG" 
                            new_s2 = self.ecs_pair[five_idx:es_idx - 1] + "N  " 
                            new_s3 = self.es_pair[five_idx:es_idx -1] + "N  " 
                            new_s4 = self.es_unpair[five_idx:es_idx - 1] + " *G" 
                        else:
                            raise ValueError("The current hs do not support.")
                else:
                    # five_offset > 1
                    three_n = 0
                    for i in range(es_idx+2, len(self.es_unpair)):
                        if self.es_unpair[i] != " " or self.es_pair[i] != " ":
                            three_n += 1
                        if three_n == five_offset - 1:
                            three_idx = i
                            break

                    if three_offset == 1:
                        if three_idx != None:
                            new_s1 =  " CG" + self.ecs_unpair[es_idx+2:three_idx+1]
                            new_s2 =  "N  " + self.ecs_pair[es_idx+2:three_idx+1]
                            new_s3 =  "N  " + self.es_pair[es_idx+2:three_idx+1]
                            new_s4 =  " *G" + self.es_unpair[es_idx+2:three_idx+1]
                        else:
                            raise ValueError("The current hs do not support.")
                    else:
                        five_n = 0
                        for i in range(es_idx-2, -1, -1):
                            if self.es_unpair[i] != " " or self.es_pair[i] != " ":
                                five_n += 1
                            if five_n == three_offset - 1:
                                five_idx = i
                                break
                        if five_idx != None and three_idx != None:
                            new_s1 = self.ecs_unpair[five_idx:es_idx - 1] + " CG" + self.ecs_unpair[es_idx+2:three_idx+1]
                            new_s2 = self.ecs_pair[five_idx:es_idx - 1] + "N  " + self.ecs_pair[es_idx+2:three_idx+1]
                            new_s3 = self.es_pair[five_idx:es_idx -1] + "N  " + self.es_pair[es_idx+2:three_idx+1]
                            new_s4 = self.es_unpair[five_idx:es_idx - 1] + " *G" + self.es_unpair[es_idx+2:three_idx+1]
                        else:
                            raise ValueError("The current hs do not support.")

        new_out = Hybrid_str(new_s1, new_s2, new_s3, new_s4)
        return new_out

    def es_paired_base_count(self):
        """return paired es base count"""
        pair_count = 0
        for _, value in self.bp_dict.items():
            if (value["info"] == "paired" or value["info"] == "wobble_paired") and value["es_base"] != " ":
                pair_count += 1
        return pair_count

    def es_base_count(self):
        """return es base count"""
        n = 0
        for _, value in self.bp_dict.items():
            if value["es_base"] != " ":
                n += 1
        return n

    def ecs_base_count(self):
        """return ecs base count"""
        n = 0
        for _, value in self.bp_dict.items():
            if value["ecs_base"] != " ":
                n += 1
        return n

    def es_seq(self, show_space=False):
        """return 3' > 5' es sequence"""
        out = ""
        for i in range(len(self.es_pair)):
            if self.es_pair[i] != " ":
                out += self.es_pair[i]
            elif self.es_unpair[i] != " ":
                out += self.es_unpair[i]
            else:
                if show_space == False:
                    pass
                else:
                    out += " "
        return out

    def ecs_seq(self, show_space=False, ignore_es_base=True):
        """return 5' > 3' ecs sequence"""
        
        out = ""
        if ignore_es_base == True:
            for i in range(len(self.ecs_pair)):
                if self.ecs_pair[i] != " ":
                    out += self.ecs_pair[i]
                elif self.ecs_unpair[i] != " ":
                    out += self.ecs_unpair[i]
                else:
                    if show_space == False:
                        pass
                    else:
                        out += " "
        elif ignore_es_base == False:
            for _,value in self.bp_dict.items():
                if value["es_base"] != " ":
                    ecs_base = value["ecs_base"]
                    if show_space == True:
                        out += ecs_base
                    elif show_space == False:
                        if ecs_base == " ":
                            pass
                        else:
                            out += ecs_base
                else:
                    pass
        return out
    def loc_to_offset(self,loc):
        """return 1 based offset of input 0_based coordinate on the hybrid_str"""
        es_idx = self.get_edit_site_index()
        five_offset = 0
        three_offset = 0
        if isinstance(loc, int) == False:
            raise TypeError("Input should be int!")
        if loc < 0 or loc >= self.__len__():
            raise ValueError("Input loc is beyond the hybrid_str!")
        if loc < es_idx: # downstream to the editing site
            for i in range(es_idx - 1, -1, -1):
                if self.es_pair[i] != " " or self.es_unpair[i] != " ":
                    three_offset += 1
                if i == loc:
                    break
            return ["three",three_offset]
        elif loc > es_idx: #upstream to the editing site
            for i in range(es_idx + 1, self.__len__()):
                if self.es_pair[i] != " " or self.es_unpair[i] != " ":
                    five_offset += 1
                if i == loc:
                    break
            return ["five",five_offset]
        else:
            return ["site",0]

def reverse_complement(s):
    return str.translate(s,str.maketrans("AUCGNaucgn","UAGCNuagcn"))[::-1]

def reverse_complement_wobble(s):
    return str.translate(s,str.maketrans("UGug","GUgu"))[::-1]

def get_mimic_region(target_seq, ref_bp_df):
    '''return the mimic region of ASO 5\' AUCG......AUCG 3\''''
    base_list = list(target_seq)
    s1,s2,s3,s4 = "","","",""
    for idx, row in ref_bp_df.iterrows():
        if row["es_base"] == " " and row["ecs_base"] != " ":
            base_list.insert(idx, " ")
    if ref_bp_df.shape[0] != len(base_list):
        print(ref_bp_df)
        print("".join(base_list))
        ref_bp_df.to_csv('./test.csv',index=0)

    # print(ref_bp_df.shape[0])
    # print(len(base_list))
    for idx, row in ref_bp_df.iterrows():
        raw_base = base_list[idx]
        ecs_base = row["ecs_base"]
        es_base = row["es_base"]
        if row["info"] == "loop":
            s4 += raw_base
            s2 += " "
            s3 += " "
            if row["ecs_base"] != " " and row["es_base"] != " ":
                if raw_base == es_base or (raw_base == "A" and es_base == "*"):
                    s1 += ecs_base
                elif raw_base == ecs_base:
                    s1 += es_base
                else:
                    if ecs_base == es_base:
                        s1 += raw_base
                    else:
                        if raw_base == "U":
                            s1 += "C"
                        elif raw_base == "G":
                            s1 += "A"
                        elif raw_base == "A":
                            s1 += "G"
                        elif raw_base == "C":
                            if ecs_base != "G":
                                s1 += ecs_base
                            else:
                                s1 += "A"
            else:
                s1 += ecs_base
                
        elif row["info"] == "es_bulge":
            s1 += " "
            s2 += " "
            s3 += " "
            s4 += raw_base
        elif row["info"] == "ecs_bulge":
            s1 += ecs_base
            s2 += " "
            s3 += " "
            s4 += raw_base
        elif row["info"] == "paired" or row["info"] == "wobble_paired":
            s1 += " "
            s3 += raw_base
            s4 += " "
            if raw_base == es_base or (raw_base == "A" and es_base == "*"):
                s2 += ecs_base
            elif raw_base == ecs_base:
                s2 += es_base
            else:
                if raw_base == "U":
                    s2 += "A"
                elif raw_base == "A":
                    s2 += "U"
                elif raw_base == "G":
                    s2 += "C"
                elif raw_base == "C":
                    s2 += "G"
    return (s1,s2,s3,s4)

def mark_editing_sites(target_pair,target_unpair,pos_1=None):
    n = 0
    es_idx = None
    s1 = list(target_pair[::-1])
    s2 = list(target_unpair[::-1])
    for i,value in enumerate(s1):
        if s2[i] != " " or value != " ":
            n += 1
        if n == pos_1:
            es_idx = i
            break
    if s1[es_idx] != " ":
        s1[es_idx] = "*"
    elif s2[es_idx] != " ":
        s2[es_idx] = "*"
    return ("".join(s1[::-1]),"".join(s2[::-1]))

def add_arm(es_five_arm,es_three_arm,s1,s2,s3,s4):
    out_s1,out_s2,out_s3,out_s4 = "","","",""
    for i in range(len(es_three_arm)):
        out_s1 += " "
        out_s2 += reverse_complement(es_three_arm[::-1][i])
        out_s3 += es_three_arm[::-1][i]
        out_s4 += " "
    for i in range(len(s1)):
        out_s1 += s1[i]
        out_s2 += s2[i]
        out_s3 += s3[i]
        out_s4 += s4[i]
    for i in range(len(es_five_arm)):
        out_s1 += " "
        out_s2 += reverse_complement(es_five_arm[::-1][i])
        out_s3 += es_five_arm[::-1][i]
        out_s4 += " "
    return (out_s1,out_s2,out_s3,out_s4)

def get_RNAhybrid_res(ASO_seq=None,target_seq=None,five_arm=None,five_mimic=None):
    tmp_cmd =  "{hybrid} -m 1000 -n 1000 -s 3utr_human {seq1} {seq2}".format(hybrid=options.hybrid,seq1=ASO_seq,seq2=target_seq)
    tmp_process = Popen(tmp_cmd,shell=True,stdin=PIPE,stdout=PIPE,stderr=PIPE)
    tmp_lines = [i.decode(encoding="utf-8") for i in tmp_process.stdout.read().splitlines()]
    hybrid_aso_unpair,hybrid_aso_pair,hybrid_target_pair,hybrid_target_unpair = tmp_lines[-6][10:-3],tmp_lines[-5][10:-3],tmp_lines[-4][10:-3],tmp_lines[-3][10:-3]
    hybrid_mfe = tmp_lines[5].split(": ")[1]
    
    hybrid_target_pair,hybrid_target_unpair = mark_editing_sites(hybrid_target_pair,hybrid_target_unpair,pos_1=five_arm+five_mimic+1)
    tmp_hs = Hybrid_str(hybrid_aso_unpair,hybrid_aso_pair,hybrid_target_pair,hybrid_target_unpair)
    return (tmp_hs,hybrid_mfe)

def write_details(df_out):
    for _,row in df_out.iterrows():
        sys.stderr.write(
            "===================================================================\n")
        sys.stderr.write("Name: %s\n" % row["site"])
        sys.stderr.write("Target seq(5'>3'): %s\n" % row["target_seq"])
        sys.stderr.write("ASO seq(5'>3')   : %s\n" % row["ASO_seq"])
        sys.stderr.write("ASO type: %s\n" % row["ASO_type"])
        sys.stderr.write("Mimic   : %s\n" % row["mimic"])
        sys.stderr.write("Reference structure:\n")
        sys.stderr.write("Ecs: %s \n" % row["ecs"])
        sys.stderr.write("Es : %s \n" % row["es"])
        sys.stderr.write(
            "5\'%s3\'\n5\'%s3\'\n3\'%s5\'\n3\'%s5\'\n\n" % (row["ref_s1"], row["ref_s2"], row["ref_s3"], row["ref_s4"]))
        sys.stderr.write("Output ASO structure:\n")
        sys.stderr.write("5\'%s3\'\n5\'%s3\'\n3\'%s5\'\n3\'%s5\'\n" % (
            row["hybrid_aso_unpair"],row["hybrid_aso_pair"],row["hybrid_target_pair"],row["hybrid_target_unpair"]))
        sys.stderr.write(
            "===================================================================\n\n") 
## low throughput
def get_control_aso(target_seq,fa=0,fm=0,GAN=False):
    if GAN == False:
        target_bases = list(target_seq)
        out_bases = [reverse_complement(i) for i in target_bases]
        out_bases[fa+fm] = 'C'
    else:
        target_bases = list(target_seq)
        out_bases = [reverse_complement(i) for i in target_bases]
        out_bases[fa+fm] = 'C'
        out_bases[fa+fm-1] = 'G'

    return ''.join(out_bases)[::-1]

def get_stru_anno(hs,wobble=True):
    out_info = []
    bp_dict = hs.bp_dict.copy()
    bp_df = pd.DataFrame.from_dict(bp_dict,orient='index')
    unpaired_locs = bp_df[bp_df["info"].isin(["loop","es_bulge","ecs_bulge"])].index.to_list()
    grouped_locs = [list(g) for _,g in groupby(unpaired_locs,lambda x,c=count(): x-next(c))]

    # for loop,bulge
    if len(grouped_locs) >0:
        for g in grouped_locs:
            if len(g) > 1:
                stru = bp_df.loc[g[0]]["info"]
                if stru == 'ecs_bulge':
                    if g[-1] < hs.get_edit_site_index():
                        offset = hs.loc_to_offset(g[-1] + 1)
                        if offset[1] >0:
                            offset = '+' + str(offset[1])
                        elif offset[1] == 0:
                            offset = '0'
                    elif g[0] > hs.get_edit_site_index():
                        offset = hs.loc_to_offset(g[0] - 1)
                        if offset[1] >0:
                            offset = '-' + str(offset[1])
                        elif offset[1] == 0:
                            offset = '0'
                    g_info = offset + '|' + stru + '|' + str(len(g))
                else:
                    g = [i for i in g if bp_df.loc[i]["es_base"] != ' ']
                    if len(g) > 1:
                        offset_3 = hs.loc_to_offset(g[0])
                        offset_5 = hs.loc_to_offset(g[-1])
                        if offset_3[1] > 0:
                            if offset_3[0] == 'three':
                                offset_3 = '+' + str(offset_3[1])
                            elif offset_3[0] == 'five':
                                offset_3 = '-' + str(offset_3[1])
                        elif offset_3[1]== 0:
                            offset_3 = '0'
                        if offset_5[1] > 0:
                            if offset_5[0] == 'three':
                                offset_5 = '+' + str(offset_5[1])
                            elif offset_5[0] == 'five':
                                offset_5 = '-' + str(offset_5[1])
                        elif offset_5[1] ==0:
                            offset_5 = '0'
                        es_base_count = len(g) - bp_df.loc[g[0]:g[-1]]["es_base"].to_list().count(' ')
                        g_info = offset_5 + ',' + offset_3 + '|' + stru  + '|' + str(es_base_count)
                    elif len(g) == 1:
                        offset = hs.loc_to_offset(g[0])
                        if offset[0] == 'five':
                            offset = '-' + str(offset[1])
                        elif offset[0] == 'three':
                            offset = '+' + str(offset[1])
                        elif offset[0] == 'site':
                            offset = '0'
                        g_info = offset + '|' + stru 
                out_info.append(g_info)

            elif len(g) == 1:
                stru = bp_df.loc[g[0]]["info"]
                if stru != "ecs_bulge":
                    offset = hs.loc_to_offset(g[0])
                    if offset[0] == 'five':
                        offset = '-' + str(offset[1])
                    elif offset[0] == 'three':
                        offset = '+' + str(offset[1])
                    elif offset[0] == 'site':
                        offset = '0'
                    g_info = offset + '|' + stru
                elif stru == "ecs_bulge":
                    if g[0] > hs.get_edit_site_index():
                        offset = hs.loc_to_offset(g[0]-1)
                        if offset[0] != 'site':
                            offset = '|' + str(offset[1])
                        else:
                            offset = '0'
                    elif g[0] < hs.get_edit_site_index():
                        offset = hs.loc_to_offset(g[0]+1)
                        if offset[0] != 'site':
                            offset = '+' + str(offset[1])
                        else:
                            offset = '0'
                    g_info = offset + '|' + stru
                out_info.append(g_info)
    
    # for wobble pair
    if wobble == True:
        for i,row in bp_df.iterrows():
            if row["info"] == 'wobble_paired':
                offset = hs.loc_to_offset(i)
                if offset[0] == 'five':
                    offset = '-' + str(offset[1])
                elif offset[0] == 'three':
                    offset = '+' + str(offset[1])
                elif offset[0] == 'site':
                    offset = '0'
                g_info = offset + '|' + "wobble_paired"
                out_info.append(g_info)

    return out_info

def add_IA_name(row):
    return '@'.join([row["es"],row["ecs"]])

def add_stru_stats(row):
    hs = Hybrid_str(row["hybrid_aso_unpair"],row["hybrid_aso_pair"],row["hybrid_target_pair"],row["hybrid_target_unpair"]).truncated_by_es_base_offset(five_mimic,three_mimic)

    # all features
    stru_all = get_stru_anno(hs,wobble=True)

    # unpair features
    stru_unpair = get_stru_anno(hs,wobble=False)

    # wobble_count 
    wobble_n = 0
    for stru in stru_all:
        if 'wobble_paired' in stru:
            wobble_n += 1

    # unpaired baes in the mimic region of target
    unpaired_es_base = 0
    for _,bp_row in pd.DataFrame.from_dict(hs.bp_dict,orient='index').iterrows():
        if bp_row['es_base'] not in ['*',' '] and bp_row['info'] not in ['paired','wobble_paired']:
            unpaired_es_base += 1
    # print(wobble_n)
    return (stru_all,len(stru_unpair),wobble_n,unpaired_es_base)

####################################################################################
def worker(row):
    global out,site_motif
    ecs_unpair = row["ecs_unpair"]
    ecs_pair = row["ecs_pair"]
    es_pair = row["es_pair"]
    es_unpair = row["es_unpair"]
    modes = None
    # get reference structure based on the five_mimic and three_mimic
    if site_motif.startswith("U") or site_motif.startswith("A"):
        modes = ["raw","A-C"]
    elif site_motif.startswith("G"):
        modes = ["raw","GAP-GCP","G-G","GA-GC"]
    else:
        modes = ["raw","A-C"] # add CA-AC if necessary
    
    for mode in modes:
        hs = None
        # step1
        try:
            hs = Hybrid_str(ecs_unpair,ecs_pair,es_pair,es_unpair)
            hs = hs.truncated_by_es_base_offset(five_offset=five_mimic,three_offset=three_mimic,mode=mode)
        except:
            continue
        if five_mimic + five_arm > loc or len(fa) - 1 < loc + three_arm + three_mimic or hs == None:
            pass
        else:
            # editing sequence(es) five_arm and three_arm region,5'>3'
            es_five_arm = fa[loc-five_mimic-five_arm:loc-five_mimic]
            es_three_arm = fa[loc+three_mimic + 1:loc+three_mimic+three_arm+1]

            # ASO target seq 5'>3'
            target_seq = fa[loc-five_mimic - five_arm:loc+three_mimic+three_arm+1]
            # ref-structure,ecs:5'>3',es:3'>5'
            ref_bp_df = pd.DataFrame.from_dict(hs.bp_dict,orient="index")

            # step2
            # ASO-mimic target seq 3'>5'
            mimic_target_seq = fa[loc-five_mimic:loc + three_mimic+1][::-1]

            # ASO-mimic
            new_s1,new_s2,new_s3,new_s4 = get_mimic_region(mimic_target_seq, ref_bp_df)
            new_s3,new_s4 = mark_editing_sites(new_s3,new_s4,pos_1=five_mimic + 1)

            # step3
            # add paired arm 
            new_s1,new_s2,new_s3,new_s4 = add_arm(es_five_arm,es_three_arm,new_s1,new_s2,new_s3,new_s4)

            # expected structure
            expected_mimic_hs = Hybrid_str(new_s1,new_s2,new_s3,new_s4)
            expected_ASO_seq = expected_mimic_hs.ecs_seq(show_space=False)

            # output reference structrure
            row["ref_s1"],row["ref_s2"],row["ref_s3"],row["ref_s4"] = hs

            # step4
            # RNAhybrid output
            hybrid_hs,hybrid_mfe = get_RNAhybrid_res(expected_ASO_seq,target_seq,five_arm,five_mimic)
            is_mimic = True
            if hybrid_hs.truncated_by_es_base_offset(five_mimic,three_mimic) != expected_mimic_hs.truncated_by_es_base_offset(five_mimic,three_mimic):
                is_mimic = False
                raw_dict = deepcopy(expected_mimic_hs.bp_dict)
                hybrid_dict = deepcopy(hybrid_hs.bp_dict)

                feature_locs = [k for k,v in raw_dict.items() if v["info"] != "paired"]
                conflict_locs = []
                if len(raw_dict) == len(hybrid_dict):
                    for i in feature_locs:
                        if hybrid_dict[i]["info"] == raw_dict[i]["info"] or (hybrid_dict[i]["info"] == "wobble_paired" and raw_dict[i]["info"] == "paired") \
                            or (hybrid_dict[i]["info"] == "paired" and raw_dict[i]["info"] == "wobble_paired"):
                            pass
                        else:
                            conflict_locs.append(i)
                else:
                    conflict_locs = [i for i in feature_locs]

                base_combinations = []
                for i in conflict_locs:
                    combination = ""
                    # here only for A-C,A-A modes
                    # raw mismatches were always in the first rank.
                    if raw_dict[i]["es_base"] != "*":
                        if raw_dict[i]["info"] == "wobble_paired":
                            combination += reverse_complement_wobble(raw_dict[i]["es_base"])
                            combination += reverse_complement(raw_dict[i]["es_base"])
                        elif raw_dict[i]["info"] == "es_bulge":
                            combination += " "
                        elif raw_dict[i]["info"] == "ecs_bulge":
                            combination += raw_dict[i]["ecs_base"]
                            combination += "".join(list({"A","U","C","G"} - {raw_dict[i]["ecs_base"]}))
                        elif raw_dict[i]["info"] == "loop":
                            if raw_dict[i]["ecs_base"]  == " ":
                                combination += " "
                            else:
                                combination += raw_dict[i]["ecs_base"]
                                combination += "".join(list({"A","U","C","G"} - {reverse_complement(raw_dict[i]["es_base"])} - {raw_dict[i]["ecs_base"]}))
                    else:
                        combination += raw_dict[i]["ecs_base"]

                    base_combinations.append(combination)
                
                for base_combination in itertools.islice(itertools.product(*base_combinations),options.n):
                    s1,s2,s3,s4 = "","","",""
                    for i,value in enumerate(base_combination):
                        changed_loc = conflict_locs[i]
                        raw_dict[changed_loc]["ecs_base"] = value

                    for _,value in raw_dict.items():
                        info = value["info"]
                        if info == "paired" or info == "wobble_paired":
                            s1 += " "
                            s2 += value["ecs_base"]
                            s3 += value["es_base"]
                            s4 += " "
                        else:
                            s1 += value["ecs_base"]
                            s2 += " "
                            s3 += " "
                            s4 += value["es_base"]

                    new_hs = Hybrid_str(s1,s2,s3,s4)
                    new_ASO_seq = new_hs.ecs_seq(show_space=False)
                    new_hybrid_hs,new_hybrid_mfe = get_RNAhybrid_res(new_ASO_seq,target_seq,five_arm,five_mimic)
                    if new_hybrid_hs.truncated_by_es_base_offset(five_mimic,three_mimic) == new_hs.truncated_by_es_base_offset(five_mimic,three_mimic):
                        expected_mimic_hs = Hybrid_str(*new_hs)
                        hybrid_hs,hybrid_mfe = new_hybrid_hs,new_hybrid_mfe
                        expected_ASO_seq = new_ASO_seq
                        is_mimic = True
                        break
            else:
                pass

            row["ASO_unpair"],row["ASO_pair"],row["target_pair"],row["target_unpair"] = expected_mimic_hs
            row["editing_site_str"] = expected_mimic_hs.get_edit_site_str()
            row["hybrid_aso_unpair"],row["hybrid_aso_pair"],row["hybrid_target_pair"],row["hybrid_target_unpair"] = hybrid_hs
            row["editing_site_str_hybrid"] = hybrid_hs.get_edit_site_str()
            row["hybrid_mfe"] = hybrid_mfe
            row["target_seq"] = target_seq
            row["target_mfe"] = round(RNA.fold(row["target_seq"])[1], 5)
            row["ASO_seq"] = expected_ASO_seq
            row["ASO_seq_mfe"] = round(RNA.fold(row["ASO_seq"])[1],5)
            row["ASO_type"] = mode
            row["ASO_length"] = len(row["ASO_seq"])
            row["mimic"] = is_mimic
            out.append(row.to_frame().T)

if __name__ == "__main__":
    des = """
            Design ASO sequence based on the input GTEX Alu editing sites dsRNA structure,
            ASO(ECS):  5' |----paired----|----------------------|N|----------------------|----paired-----| 3'
            target(ES):3' |--three_arm---|--three_mimic_region--|*|---five_mimic_region--|---five_arm----| 5'
            * is the editing site,N can be A,T,C,G or None.
            
            Steps:
                1. Get truncated RNAhybrid double strand structures from substrates based on the input fm,tm parameters.         
                   truncated-structure example:
                            5'  A               CC 3'
                            5'UA CUGGGAUUACAGGCA  C3'
                            3'AU GACCCUG*UGUCCGU  G5'
                            3'  C               AC 5'
                2. Simulate ASO based on the feature types of truncated structures, introducing structure to ASO-target at same locations.
                    ASO-target example:
                            5'  C               CC 3'
                            5'CC GGGGGGUUAAGCAGU  G3'
                            3'GG CCCCCCG*UUCGUCA  C5'
                            3'  U               AC 5'
                   
                3. Add five_arm(complementary sequence) and three_arm to the output ASO.
                    ASO-target example:
                  5'            C               CC           3'
                  5'AUGACCUUGGCC GGGGGGUUAAGCAGU  GUGGUGCAGGA3'
                  3'UACUGGAACCGG CCCCCCG*UUCGUCA  CACCACGUCCU5'
                  3'            U               AC           5'

                4. Run RNAhybrid based on the ASO and target sequence. If the structure is different to the truncated structure,
                   then simulate new ASOs by introduing all possible bases at the location of features. Max simulating times
                   could be set by -n, default is 1000.
                   RNAbybrid output:
                  5'            C               CC           3'
                  5'AUGACCUUGGCC GGGGGGUUAAGCAGU  GUGGUGCAGGA3'
                  3'UACUGGAACCGG CCCCCCG*UUCGUCA  CACCACGUCCU5'
                  3'            U               AC           5'

                5. Output the final ASO sequence and its structure.
        """

    parser = argparse.ArgumentParser(description=des, formatter_class=argparse.RawTextHelpFormatter)
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-t", dest="target", metavar="target",help="ASO target,fasta format,one seq only", required=True)
    group_required.add_argument("-l", dest="loc", metavar="loc",help="1_based editing_site coordinate on the input target fasta", required=True, type=int)
    group_required.add_argument("-o", dest="out", metavar="output", help="output file name prefix", required=True)
    group_required.add_argument("-fm", dest="five_mimic", metavar="five_mimic_length",help="five_mimic region length of the editing site.", type=int, required=True)
    group_required.add_argument("-tm", dest="three_mimic", metavar="three_mimic_length", type=int,help="three_mimic region length of the editing site", required=True)

    group_optional = parser.add_argument_group("Optional")
    group_optional.add_argument('--hybrid',dest='hybrid',metavar='RNAhybird',type=str, required=False, help=f'Path to RNAhybrid binary (default: {RNAHYBRID_BIN})',default=RNAHYBRID_BIN)
    group_optional.add_argument('--substrates',dest='substrates',metavar='substrates',type=str,required=False,help=f'Path to substrates directory (default: {SUBSTRATES_DIR})',default=SUBSTRATES_DIR)
    group_optional.add_argument("-m", dest="motif", metavar="motif",help="input motif structure to mimic,default is the site's orginal motif", required=False)
    group_optional.add_argument("-fa", dest="five_arm", metavar="five_arm_length",help="length of the ES 5' binding region,default=25,if you do not want it,just set it be 0.", type=int, default=25, required=False)
    group_optional.add_argument("-ta", dest="three_arm", metavar="three_arm_length",help="length of the ES 3' binding region,default=25,if you do not want it,just set it be 0.", type=int, default=25, required=False)
    group_optional.add_argument("-p", dest="process", metavar="process",default=1, type=int, help="processor number,default=1", required=False)
    group_optional.add_argument("-n", dest="n", help="Max n times for simulating ASOs for RNAhybrid, default=1000",type=int,default=1000,required=False) 
    group_optional.add_argument("-L", dest="low_throughput",help="Filter for low_throughput ASO synthesis",default=False,required=False,action="store_true")
    options = parser.parse_args()

    parent_pid = os.getpid()
    # reporter
    fa_dict = {record.id: str(record.seq).upper().replace("T","U")
               for record in SeqIO.parse(options.target, "fasta")}
    if len(fa_dict) != 1:
        print("One target seq only.")
        sys.exit("1")

    # raw parameters
    key = list(fa_dict.keys())[0]
    fa = fa_dict[key]
    loc_1 = options.loc

    loc = loc_1 - 1
    five_arm = options.five_arm
    three_arm = options.three_arm
    five_mimic = options.five_mimic
    three_mimic = options.three_mimic
    p = options.process
    total_length = five_arm + five_mimic + 1 + three_mimic + three_arm
    # check mimic parameters
    if three_mimic == 0 and five_mimic == 0:
        print("Mimic length should be bigger than 0!")
        sys.exit(1)

    available_motifs = [i+"A"+j for i in ["A", "U", "C", "G"]
                        for j in ["A", "U", "C", "G"]]
    if not options.motif:
        site_motif = fa[loc-1:loc+2].upper().replace("T", "U")
    elif options.motif:
        site_motif = options.motif.upper().replace("T", "U")

    # format input
    if site_motif not in available_motifs:
        print("Input editing site motif is not avaiable! Check the input editing site location on the target fasta!")
        sys.exit(1)
    df = pd.read_csv("{substrates}/{site_motif}.csv".format(substrates=options.substrates,site_motif=site_motif))
    df["level"] = df["site"].apply(lambda x: float(x.split("|")[-1].split(":")[1]))
    rows = [row for _, row in df.iterrows()]

    sys.stderr.write("[%s] CMD: %s, pid: %d \n" % (strftime("%Y-%m-%d %H:%M:%S", localtime()), " ".join(sys.argv), parent_pid))


    out_csv_name = options.out + "_len{length}_fa{fa}_ta{ta}_fm{fm}_tm{tm}.csv".format(
        length=total_length, fa=five_arm, ta=three_arm, fm=five_mimic, tm=three_mimic)
    out_csv_log = options.out + "_len{length}_fa{fa}_ta{ta}_fm{fm}_tm{tm}.log".format(
        length=total_length, fa=five_arm, ta=three_arm, fm=five_mimic, tm=three_mimic)
    manager = mp.Manager()
    out = manager.list()
    sys.stderr = open(out_csv_log, "w")
    with Pool(p) as pool:
        pool.map(worker, rows)
        pool.close()
        pool.join()
    df_out = pd.concat(out, ignore_index=True)
    df_out = df_out.sort_values(by="level",ascending=False)

    # if options.low_throughput:
    if not options.low_throughput:
        # write out log
        write_details(df_out)
        # df_out = df_out.drop_duplicates(subset="ASO_seq",keep="first",ignore_index=True)
        df_out.to_csv(out_csv_name, index=0)
        sys.stderr.write("[%s] Completed!\n" % strftime("%Y-%m-%d %H:%M:%S", localtime()))

    if options.low_throughput:
        df_out = df_out[df_out['mimic'] == True].copy(deep=True)
        df_out["IA_name"] = df_out.apply(add_IA_name,axis=1)
        df_out = df_out.sort_values(by=["top","IA_name","ASO_type"],ascending=[True,True,False])
        df_out = df_out.drop_duplicates(subset=["ASO_seq"],keep='first')
        target_seq = df_out["target_seq"].values[0]
        control_aso = None
        if target_seq[five_arm+five_mimic-1] != 'G':
            # force A-C
            control_aso = get_control_aso(target_seq,fa=five_arm,fm=five_mimic,GAN=False)
            df_out = df_out[df_out["editing_site_str_hybrid"] == '*-C-loop'].copy(deep=True)
        else:
            # force G-G,G-A
            control_aso = get_control_aso(target_seq,fa=five_arm,fm=five_mimic,GAN=True)
            keeped_idx = []
            for i,row in df_out.iterrows():
                hs = Hybrid_str(row["hybrid_aso_unpair"],row["hybrid_aso_pair"],row["hybrid_target_pair"],row["hybrid_target_unpair"])
                if hs.get_base_pair_by_five_offset(five_offset=1) in ["G-G-loop","G-A-loop"]:
                    keeped_idx.append(i)
            df_out = df_out.loc[keeped_idx].copy(deep=True)
        df_out = df_out[df_out["ASO_seq"] != control_aso].copy(deep=True)
        df_out = df_out.reset_index(drop=True)
        df_out[['All_feature','Loop_bulge_count','wobble_count','unpaired_es_base_count']] = df_out.apply(add_stru_stats,axis=1,result_type='expand')
        ## short data
        df_out = df_out[(df_out["Loop_bulge_count"] >= 2) & (df_out['Loop_bulge_count'] <= 4)].copy(deep=True)
        

        # down sampling,each bin select 5 gRNAs
        df_output = None
        down_sampling_dfs = []
        df_fully_edited = df_out[df_out['level'] == 1].copy(deep=True)
        df_fully_edited = df_fully_edited.sort_values(by=["unpaired_es_base_count"],ascending=[True])
        df_fully_edited = df_fully_edited[df_fully_edited["wobble_count"] <= 1].head(5).copy(deep=True)
        bins = np.arange(1,0,-0.05)[::-1]
        labels = [f"{round(bins[i], 2)} <= level < {round(bins[i+1], 2)}" for i in range(len(bins)-1)]
        
        df_out = df_out[df_out['level'] < 1].copy(deep=True)
        df_out['bin'] = pd.cut(df_out['level'], bins=bins, right=False,labels=labels)
        for label in labels:
            bin_df = df_out[df_out['bin'] == label].sort_values(by='unpaired_es_base_count',ascending=True)
            bin_df = bin_df[bin_df["wobble_count"] <= 1].head(5).copy(deep=True)
            down_sampling_dfs.append(bin_df)

        down_sampling_dfs = down_sampling_dfs[::-1]
        down_sampling_dfs.insert(0,df_fully_edited)
        df_output = pd.concat(down_sampling_dfs,ignore_index=True)

        with open(options.out + '_low_throughput_structure.txt','w') as f:
            for i,row in df_output.iterrows():
                hs = Hybrid_str(row["hybrid_aso_unpair"],row["hybrid_aso_pair"],row["hybrid_target_pair"],row["hybrid_target_unpair"])
                
                stru_anno = get_stru_anno(hs,wobble=True)
                f.write('====================================================================================\n')
                f.write("          Top: " + str(i+1) + '\n')
                f.write("      ASO seq: " + row["ASO_seq"] + '\n')
                f.write("   target seq: " + row["target_seq"] + '\n')
                f.write("    structure: " +  '  '.join(stru_anno) + '\n')
                f.write(hs.show(out=True) + '\n')
                f.write('====================================================================================\n')
                f.write('\n')

        output_columns = ['target_seq','ASO_seq','All_feature']
        control_row = pd.DataFrame([[target_seq, control_aso, 'control']], columns=output_columns)
        df_output = df_output[output_columns]
        df_final = pd.concat([control_row,df_output],ignore_index=True)
        df_final.to_csv(options.out + "_low_throughput.csv",index=0)
        sys.stderr.write("[%s] Completed!\n" % strftime("%Y-%m-%d %H:%M:%S", localtime()))
