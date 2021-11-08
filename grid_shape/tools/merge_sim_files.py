import os
import json
import glob
import numpy as np
import argparse

'''
Merge simulation outputs to reduce number of files
Usage:
    python3 merge_sim_files.py sim_path merge_config_file

Parameters:
    sim_path: path of the main folder containing the Stshp_XmaxLibrary_* directories
    merge_config_file: .json parameter files describing primaries and layout geometry and step size

Output:
    a bunch of .json files of the form Primary_geom_step.json
'''



def parse_args():
    ap = argparse.ArgumentParser(description="aggregate simulation files")

    ap.add_argument(
        "sim_path",
        type=str,
        help="Folder containing the simulation",
    )
    ap.add_argument(
        "merge_config_file",
        type=str,
        help="layoutvconfig file ",
    )
    ap.add_argument(
        "output_path",
        type=str,
        help="output folder",
    )   
    args = ap.parse_args()
    return args


def isfile(filepath):
    if os.path.isfile(filepath):
        return filepath
    else:
        return None


def make_combined_layout(layouts):
    combined_layout_list = []
    for k, v in layouts.items(): 
        for i, step in enumerate(v):
            d = {}
            d["name"] = k+'_'+step
            d["id"] = k+'_%d'%i
            combined_layout_list.append(d)

    return combined_layout_list


def merge_sims(path, output_path, primary_list, combined_layout_list):

    for primary in primary_list:
        for cl in combined_layout_list:
            globals()["%s_%s"%(primary, cl["id"])] = {}
        
    for primary in primary_list:
        l = glob.glob(os.path.join(path  ,'*_%s_*'%primary))
        for cl in combined_layout_list:
            print("Processing %s"%("%s_%s"%(primary, cl["name"])))
            for dirr in l:
                bn = os.path.basename(dirr)
            
                ef_showerinfo = isfile(
                    os.path.join(
                        dirr,
                        bn +".Interpolated."+cl["name"]+"_efield.P2Pdat.showerinfo"
                    )
                )
                ef = isfile(
                    os.path.join(
                        dirr,
                        bn +".Interpolated."+cl["name"]+"_efield.P2Pdat"
                    )
                )
                vo_showerinfo = isfile(
                    os.path.join(
                        dirr,
                        bn +".Interpolated."+cl["name"]+"_voltage.P2Pdat.showerinfo"
                    )
                )
                vo = isfile(
                    os.path.join(
                        dirr,
                        bn +".Interpolated."+cl["name"]+"_voltage.P2Pdat"
                    )
                )
                d = {}
                if ef_showerinfo:
                    d["ef_showerinfo"] = open(ef_showerinfo).readlines()
                if vo_showerinfo:
                    d["vo_showerinfo"] = open(vo_showerinfo).readlines()
                if ef:
                    d["ef"] = np.loadtxt(ef).tolist()
                if vo:
                    d["vo"] = np.loadtxt(vo).tolist()

                globals()["%s_%s"%(primary, cl["id"])][bn] = d    

            output_file = os.path.join(
                output_path,
                "%s.json"%("%s_%s"%(primary, cl["name"]))
            )
            with open(output_file, "w") as f:
                json.dump(globals()["%s_%s"%(primary, cl["id"])], f, indent=4)       
        

if __name__ == "__main__":

    args = parse_args()

    with open(args.merge_config_file) as f:
        config = json.load(f)
 
    combined_layout_list = make_combined_layout(config["layouts"])
    
    os.makedirs(args.output_path, exist_ok=True)
    merge_sims(
        args.sim_path,
        args.output_path,
        config["primaries"],
        combined_layout_list
    )

