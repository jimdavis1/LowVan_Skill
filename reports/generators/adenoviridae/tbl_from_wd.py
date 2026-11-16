#!/usr/bin/env python3
"""Recover <genome>.feature.tbl from a run_gto_eval work dir (annotator stdout), adding
the splice/special rows exactly as run_gto_eval --tbl-dir does.  usage: tbl_from_wd.py <anno_dir> <tbl_dir>"""
import sys,os,glob,shutil
sys.path.insert(0,'/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/LowVan_Skill/skill/scripts')
from run_gto_eval import _append_special_rows
a,t=sys.argv[1:3]; os.makedirs(t,exist_ok=True); n=0
for wd in sorted(glob.glob(os.path.join(a,'_wd','*'))):
    b=os.path.basename(wd); src=os.path.join(wd,b+'.stdout.txt'); ann=os.path.join(a,b+'.ann.gto')
    if os.path.exists(src) and os.path.getsize(src)>0 and os.path.exists(ann):
        dst=os.path.join(t,b+'.feature.tbl'); shutil.copyfile(src,dst); _append_special_rows(dst,ann); n+=1
print(a,'->',t,n,'tables')
