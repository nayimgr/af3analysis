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
    chain_lengths = []
    PTM_locations = []

    for sequence in sequences:
        if "proteinChain" in sequence:
            for n in range(int(sequence['proteinChain']['count'])):
                additional_atoms = 0
                if "modifications" in sequence["proteinChain"]:
                    for ptm in sequence["proteinChain"]["modifications"]:
                        if ptm["ptmType"] == "CCD_SEP": additional_atoms += 10
                        elif ptm["ptmType"] == "CCD_TPO": additional_atoms += 11
                        elif ptm["ptmType"] == "CCD_PTR": additional_atoms += 16
                        elif ptm["ptmType"] == "CCD_NEP": additional_atoms += 14
                        elif ptm["ptmType"] == "CCD_HIP": additional_atoms += 14
                        elif ptm["ptmType"] == "CCD_ALY": additional_atoms += 12
                        elif ptm["ptmType"] == "CCD_MLY": additional_atoms += 11
                        elif ptm["ptmType"] == "CCD_M3L": additional_atoms += 12
                        elif ptm["ptmType"] == "CCD_MLZ": additional_atoms += 10
                        elif ptm["ptmType"] == "CCD_2MR": additional_atoms += 13
                        elif ptm["ptmType"] == "CCD_AGM": additional_atoms += 12
                        elif ptm["ptmType"] == "CCD_MCS": additional_atoms += 12
                        elif ptm["ptmType"] == "CCD_HYP": additional_atoms += 8
                        elif ptm["ptmType"] == "CCD_HY3": additional_atoms += 8
                        elif ptm["ptmType"] == "CCD_LYZ": additional_atoms += 10
                        elif ptm["ptmType"] == "CCD_AHB": additional_atoms += 9
                        elif ptm["ptmType"] == "CCD_P1L": additional_atoms += 22
                        elif ptm["ptmType"] == "CCD_SNN": additional_atoms += 7
                        elif ptm["ptmType"] == "CCD_SNC": additional_atoms += 8
                        elif ptm["ptmType"] == "CCD_TRF": additional_atoms += 16
                        elif ptm["ptmType"] == "CCD_KCR": additional_atoms += 14
                        elif ptm["ptmType"] == "CCD_CIR": additional_atoms += 11
                        elif ptm["ptmType"] == "CCD_YHA": additional_atoms += 11

                chain_lengths.append(len(sequence['proteinChain']['sequence'])+additional_atoms)
        elif "dnaSequence" in sequence:
            for n in range(int(sequence['dnaSequence']['count'])):
                chain_lengths.append(len(sequence['dnaSequence']['sequence']))
        elif "rnaSequence" in sequence:
            for n in range(int(sequence['rnaSequence']['count'])):
                chain_lengths.append(len(sequence['rnaSequence']['sequence']))
        elif "ion" in sequence:
            continue
            # for n in range(int(sequence['ion']['count'])):
            #     chain_lengths.append(1)
        elif "ligand" in sequence:
            # Number of atoms in each ligand
            for n in range(int(sequence['ligand']['count'])):
                if sequence['ligand']['ligand'] == "CCD_ADP": chain_lengths.append(27)
                elif sequence['ligand']['ligand'] == "CCD_ATP": chain_lengths.append(31)
                elif sequence['ligand']['ligand'] == "CCD_AMP": chain_lengths.append(23)
                elif sequence['ligand']['ligand'] == "CCD_GTP": chain_lengths.append(32)
                elif sequence['ligand']['ligand'] == "CCD_GDP": chain_lengths.append(28)
                elif sequence['ligand']['ligand'] == "CCD_FAD": chain_lengths.append(53)
                elif sequence['ligand']['ligand'] == "CCD_NAD": chain_lengths.append(44)
                elif sequence['ligand']['ligand'] == "CCD_NAP": chain_lengths.append(48)
                elif sequence['ligand']['ligand'] == "CCD_NDP": chain_lengths.append(48)
                elif sequence['ligand']['ligand'] == "CCD_HEM": chain_lengths.append(43)
                elif sequence['ligand']['ligand'] == "CCD_HEC": chain_lengths.append(43)
                elif sequence['ligand']['ligand'] == "CCD_PLM": chain_lengths.append(18)
                elif sequence['ligand']['ligand'] == "CCD_OLA": chain_lengths.append(20)
                elif sequence['ligand']['ligand'] == "CCD_MYR": chain_lengths.append(16)
                elif sequence['ligand']['ligand'] == "CCD_CIT": chain_lengths.append(13)
                elif sequence['ligand']['ligand'] == "CCD_CLA": chain_lengths.append(65)
                elif sequence['ligand']['ligand'] == "CCD_CHL": chain_lengths.append(66)
                elif sequence['ligand']['ligand'] == "CCD_BCL": chain_lengths.append(66)
                elif sequence['ligand']['ligand'] == "CCD_BCB": chain_lengths.append(66)
                else: print("Could not recognize ligand "+str(sequence))
        else: print("There is a molecule I cannot recognize: "+str(sequence))

    file_list = sorted(glob.glob(path+"/*full_data_?.json"))
    all_preds = []

    for file in file_list:
        with open(file, 'r') as f:
            all_preds.append(json.load(f))

    return [all_preds,chain_lengths,name]

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
        chain_number = 0

        for length in proteins_lengths:
            rect_x = patches.Rectangle((previous_l, 0), length, -0.03*total_length, linewidth=1, edgecolor='none', facecolor=palette[chain_number], clip_on=False)
            rect_y = patches.Rectangle((total_length, previous_l), 0.03*total_length, length, linewidth=1, edgecolor='none',facecolor=palette[chain_number], clip_on=False)
            axes[i].add_patch(rect_x)
            axes[i].add_patch(rect_y)
            axes[i].text(previous_l + length / 2, -0.015*total_length, chain_ids[chain_number], ha='center', va='center', color='black')
            axes[i].text(1.015*total_length, previous_l + length / 2, chain_ids[chain_number], ha='center', va='center', color='black')
            previous_l += length
            chain_number += 1

        n_model = str(int(n_model)+1)

    plt.savefig(data[2]+"_pae.png", dpi=300.0)


data = get_full_data(directory_path)
plot_pae(data)

