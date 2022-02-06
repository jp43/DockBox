import os
import sys
import pandas as pd
from glob import glob
import shutil
import argparse

parser = argparse.ArgumentParser(description="Computes hit rates and enrichment factors from output of extract_dbx_best_poses (each folder generated with extract_dbx_best_poses should contain a file called best_poses.csv)")

parser.add_argument('-n',
    dest='nactives',
    type=int,
    required=True,
    help='Number of active compounds in the set (in best_poses.csv files, all the active compounds should be provided first, followed by decoys")')

nhits = 100 # number of hits that should be considered for hit rates

# update parsers with arguments
args = parser.parse_args()

dirs = []
for dir in glob("*"):
    if os.path.isfile(dir+"/best_poses.csv"):
        dirs.append(dir)

topdir = "top_hits"
shutil.rmtree(topdir, ignore_errors=True)
os.mkdir(topdir)

info = {"method": [], "EF": [], "hit-rate": []}

for dir in dirs:
    df = pd.read_csv(dir+"/best_poses.csv")
    df['status'] = "decoy"
    df.iloc[:args.nactives, df.columns.get_loc('status')] = "active"

    if 'consensus' in df.columns:
        df_groupby = df.groupby(['status'])[['consensus']].sum()
        tp = int(df_groupby.ix[0]['consensus']) # True Positives
        fn = args.nactives - tp # False Negatives
        fp = int(df_groupby.ix[1]['consensus']) # False Negatives
        tn = len(df)-args.nactives - fp # True Negatives

        nctot = tp + fn + fp + tn
        nc = tp + fp
        ratio = tp*1./fp
        ef = tp*1./(tp+fn)*nctot*1./nc

        df = df[df['consensus']]
    else:
        ratio = 100.
        ef = 1.

    if dir.startswith("docking"):
        column = ["score"]

    elif dir.startswith("rescoring"):
        column = [dir.split("_")[-1]]

    elif dir.startswith("cd"):
        column = []
        for prgm in dir.split("_")[1:]:
            column.append("score_"+prgm)

    elif dir.startswith("sbcd"):
        column =[]
        for prgm in dir.split("_")[1:]:
            column.append(prgm)

    for cc in column:
        df_top_hits = df.sort_values(by=cc).head(nhits)
        if dir.startswith(("cd", "sbcd")):
            ccs = cc.split('_')
            method = dir + "_scored_with_" + ccs[-1]
        else:
            method = dir
        csvfile = topdir + "/" + method + ".csv"
        df_top_hits[['ligID', cc]].to_csv(csvfile, index=False, float_format="%.3f")
        info["method"].append(method)
        info["EF"].append(ef)
        info["hit-rate"].append(len(df_top_hits[df_top_hits['status']=='active']))

df_info = pd.DataFrame(info)
df_info = df_info.sort_values('hit-rate', ascending=False)
df_info[["method", "hit-rate", "EF"]].to_csv(topdir+"/ranking.csv",index=False, float_format="%.3f")
