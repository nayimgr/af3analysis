import pandas as pd
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob
import string
import argparse

# Parse input
parser = argparse.ArgumentParser(description='Script for plotting pae from all models from AF3 server')
parser.add_argument('--dir', type=str, help='directory containing all json files', default="." )
args = parser.parse_args()
directory_path = args.dir

def get_full_data(path):
    job_parameters_file = glob.glob(path+"/*job_request.json")[0]
    with open(job_parameters_file, 'r') as p:
        parameters = json.load(p)[0]
    name = parameters['name']
    sequences = parameters['sequences']
    protein_lengths = []
    for i in sequences:
        if "proteinChain" in i:
            for n in range(int(i['proteinChain']['count'])):
                protein_lengths.append(len(i['proteinChain']['sequence']))
        elif "dnaChain" in i:
            for n in range(int(i['dnaChain']['count'])):
                protein_lengths.append(len(i['dnaChain']['sequence']))
        elif "rnaChain" in i:
            for n in range(int(i['rnaChain']['count'])):
                protein_lengths.append(len(i['rnaChain']['sequence']))
        elif "ion" in i:
            for n in range(int(i['ion']['count'])):
                protein_lengths.append(len(i['ion']['sequence']))
        elif "ligand" in i:
            for n in range(int(i['ligand']['count'])):
                protein_lengths.append(len(i['ligand']['sequence']))
        else: print("There is a molecule I cannot recognize in "+str(i))
    file_list = sorted(glob.glob(path+"/*full_data_?.json"))
    all_preds = []
    for file in file_list:
        with open(file, 'r') as f:
            all_preds.append(json.load(f))
    return [all_preds,protein_lengths,name]

def plot_pae(data):
    name = data[2]+" model "
    proteins_lengths = data[1]
    num_models = len(data[0])
    total_length = sum(proteins_lengths)
    n_model = "0"
    fig, axes = plt.subplots(1, num_models, figsize=(8*num_models, 6))
    for i, model in enumerate(data[0]):
        pae_matrix = np.array(model['pae'])
        sns.heatmap(pae_matrix, cmap='viridis', ax=axes[i], square=True, cbar=True, xticklabels=False, yticklabels=False, cbar_kws={"shrink":0.2})
        axes[i].set_title(name+n_model, y=1.05)
        axes[i].set_aspect('equal', adjustable='box')
        previous_l = 0
        palette = plt.cm.get_cmap("viridis")
        palette = [palette(a / (len(proteins_lengths) - 1)) for a in range(len(proteins_lengths))] # Adapting to discrete palette
        chain_ids = list(string.ascii_uppercase[:len(proteins_lengths)])
        c = 0
        for l in proteins_lengths:
            rect_x = patches.Rectangle((previous_l, 0), l, -0.03*total_length, linewidth=1, edgecolor='none', facecolor=palette[c], clip_on=False)
            rect_y = patches.Rectangle((total_length, previous_l), 0.03*total_length, l, linewidth=1, edgecolor='none',facecolor=palette[c], clip_on=False)
            axes[i].add_patch(rect_x)
            axes[i].add_patch(rect_y)
            axes[i].text(previous_l + l / 2, -0.015*total_length, chain_ids[c], ha='center', va='center', color='black')
            axes[i].text(1.015*total_length, previous_l + l / 2, chain_ids[c], ha='center', va='center', color='black')
            previous_l += l
            c += 1
        n_model = str(int(n_model)+1)
    plt.savefig(data[2]+"_pae.png", dpi=300.0)


data = get_full_data(directory_path)
plot_pae(data)

