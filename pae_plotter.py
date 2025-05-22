import json
import numpy as np
import matplotlib.pyplot as plt
import glob
import string
import argparse

# Parse input
parser = argparse.ArgumentParser(description='Script for plotting pae from all models from AF3 server')
parser.add_argument('--dir', type=str, help='directory containing all json files', default="." )
parser.add_argument('--palette', type=str, help='color scheme for the plot', default="bwr")
args = parser.parse_args()
directory_path = args.dir
color_scheme = args.palette


# How many additional atoms to add per PTM or ligand
protPTMs = {"CCD_SEP":10,
            "CCD_TPO":11,
            "CCD_PTR":16,
            "CCD_NEP":14,
            "CCD_HIP":14,
            "CCD_ALY":12,
            "CCD_MLY":11,
            "CCD_M3L":12,
            "CCD_MLZ":10,
            "CCD_2MR":13,
            "CCD_AGM":12,
            "CCD_MCS":12,
            "CCD_HYP":8,
            "CCD_HY3":8,
            "CCD_LYZ":10,
            "CCD_AHB":9,
            "CCD_P1L":22,
            "CCD_SNN":7,
            "CCD_SNC":8,
            "CCD_TRF":16,
            "CCD_KCR":14,
            "CCD_CIR":11,
            "CCD_YHA":11}

# dnaPTMs = { "CCD_5CM":
#             "CCD_C34":
#             "CCD_5HC":
#             "CCD_6OG":
#             "CCD_6MA":
#             "CCD_1CC":
#             "CCD_8OG":
#             "CCD_5FC":
#             "CCD_3DR":}

ligands_dict = {"CCD_ADP":27,
                "CCD_ATP":31,
                "CCD_AMP":23,
                "CCD_GTP":32,
                "CCD_GDP":28,
                "CCD_FAD":53,
                "CCD_NAD":44,
                "CCD_NAP":48,
                "CCD_NDP":48,
                "CCD_HEM":43,
                "CCD_HEC":43,
                "CCD_PLM":18,
                "CCD_OLA":20,
                "CCD_MYR":16,
                "CCD_CIT":13,
                "CCD_CLA":65,
                "CCD_CHL":66,
                "CCD_BCL":66,
                "CCD_BCB":66}

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
                        additional_atoms += protPTMs[ptm["ptmType"]]

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
                try: chain_lengths.append(ligands_dict[sequence['ligand']['ligand']])
                except: print("Could not recognize ligand "+str(sequence))
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
    n = 0

    fig, axes = plt.subplots(nrows=1, ncols=num_models, figsize=(6*num_models, 5))

    for i, model in enumerate(data[0]):
        pae_matrix = np.array(model['pae'])
        im = axes[i].imshow(pae_matrix, cmap=color_scheme, vmin=0, vmax=30)
        axes[i].set_title(name+n_model, y=1.05)
        previous_l = 0
        palette = plt.cm.get_cmap(color_scheme)
        palette = [palette(a / (len(proteins_lengths) - 1 if len(proteins_lengths) > 1 else 1)) for a in range(len(proteins_lengths))]  # Fixing division by zero
        chain_ids = list(string.ascii_uppercase[:len(proteins_lengths)])
        chain_number = 0

        for length in proteins_lengths[:-1]: # trimming list to avoid painting the last line
            axes[i].axvline(previous_l+length, color="black", linewidth=1.5)
            axes[i].axhline(previous_l+length, color="black", linewidth=1.5)
            axes[i].text(previous_l + length / 2, -0.015*total_length, chain_ids[chain_number], ha='center', va='center', color='black')
            axes[i].text(1.015*total_length, previous_l + length / 2, chain_ids[chain_number], ha='center', va='center', color='black')
            previous_l += length
            chain_number += 1

        n_model = str(int(n_model)+1)
    
    cax = fig.add_axes([axes[-1].get_position().x1+0.012, axes[-1].get_position().y0, 0.01, axes[-1].get_position().height])
    fig.colorbar(im, cax=cax)
    plt.subplots_adjust(wspace=0.06)
    plt.savefig(data[2]+"_pae.svg", dpi=300.0)
    plt.savefig(data[2]+"_pae.png", dpi=300.0)

def plot_plddt(data, window_size=8):
    
    name = data[2]+" model "
    proteins_lengths = data[1]
    num_models = len(data[0])
    total_length = sum(proteins_lengths)
    n_model = "0"
    chains, chain_lengths = np.unique(data[0][0]["atom_chain_ids"], return_counts=True)

    palette = plt.cm.get_cmap(color_scheme)
    palette = [palette(a / (len(proteins_lengths) - 1 if len(proteins_lengths) > 1 else 1)) for a in range(len(proteins_lengths))]  # Fixing division by zero
    prev_chain_length = 0

    fig, axes = plt.subplots(figsize=(24, 6))

    for i, model in enumerate(data[0]):
        plddt_array = np.array(model["atom_plddts"])
        smoothed_plddt = np.convolve(plddt_array, np.ones(window_size)/window_size, mode='valid')
        label = name+n_model
        axes.plot(smoothed_plddt, label=label, color=palette[i % len(palette)])  # Using modulo to avoid index errors
        axes.set_xticks([])
        axes.set_ylabel("plddt") 
        

        n_model = str(int(n_model)+1)
        prev_chain_length += chain_lengths

    length = 0
    for i in range(len(chain_lengths)):
        # axes.text(length + chain_lengths[i] / 2, np.argmin(plddt_array)*0.9, chains[i], ha='center', va='center', color='black')
        length += chain_lengths[i]
        axes.axvline(x=length, color='black', linestyle='--', linewidth=2)

    
    plt.savefig(data[2]+"_plddt.png", dpi=300.0)
    plt.savefig(data[2]+"_plddt.svg", dpi=300.0)
        
print("starting")
data = get_full_data(directory_path)
print("data loaded")
plot_pae(data)
print("pae")
plot_plddt(data)
print("plddt")
