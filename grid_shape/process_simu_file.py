import os
import json
import glob
import numpy as np


path = "/Users/benoitl/Documents/GRAND/P2PDataNew/P2PdataNew/"


layout_list = ["rect", "hexhex"]

step_hexhex_list = ["125",]# "250", "500", "750", "1000"]
step_rect_list = ["137.5",]# "275", "550", "825", "1100"]
primary_list = ["Proton",]# "Gamma", "Iron"]


combined_layout_list = [] # will be a list of dictionaries
for i, step in enumerate(step_hexhex_list):
    d = {}
    d["name"] = 'hexhex'+'_'+step
    d["id"] = 'hexhex'+'_%d'%i
    combined_layout_list.append(d)
  
for i, step in enumerate(step_rect_list):
    d = {}  
    d["name"] = 'rect'+'_'+step
    d["id"] = 'rect'+'_%d'%i
    combined_layout_list.append(d)



def isfile(filepath):
    if os.path.isfile(filepath):
        return filepath
    else:
        return None



for primary in primary_list:
    for cl in combined_layout_list:
        globals()["%s_%s"%(primary, cl["id"])] = {}
    


for primary in primary_list:
    l = glob.glob(os.path.join(path  ,'*_%s_*'%primary))
    for cl in combined_layout_list:
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
    
        with open("./%s.json"%("%s_%s"%(primary, cl["name"])), "w") as f:
            json.dump(globals()["%s_%s"%(primary, cl["id"])], f, indent=4)       
    