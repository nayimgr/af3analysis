import pandas as pd
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob


# Replace 'file_path.json' with the actual path to your JSON file
directory_path = '.'
# file_path = 'fold_cst_telojunction_job_request.json'

def get_full_data(path):
    job_parameters_file = glob.glob(path+"/*job_request.json")[0]
    with open(job_parameters_file, 'r') as p:
        parameters = json.load(p)[0]
    name = parameters['name']
    sequences = parameters['sequences']
    print(sequences[0]['proteinChain'])
    protein_lengths = []
    for i in sequences:
        try: protein_lengths.append(len(i['proteinChain']['sequence'])*i['proteinChain']['count'])
        except:
            try: protein_lengths.append(len(i['dnaChain']['sequence'])*i['dnaChain']['count'])
            except: 
                try: protein_lengths.append(len(i['ion']['sequence'])*i['ion']['count'])
                except:
                    try: protein_lengths.append(len(i['rnaChain']['sequence'])*i['rnaChain']['count'])
                    except:
                        print("wtf")
                        exit()
    file_list = sorted(glob.glob(path+"/*full_data_?.json"))
    all_preds = []
    for file in file_list:
        with open(file, 'r') as f:
            all_preds.append(json.load(f))
    return [all_preds,protein_lengths]

all_data = get_full_data(directory_path)[0]
proteins_lengths =  get_full_data(directory_path)[1]
total_length = sum(proteins_lengths)

name = "pae_"
n_model = "0"
num_models = len(all_data)
fig, axes = plt.subplots(1, num_models+1, figsize=((8*num_models)+1, 6))

color_bar_bool = False
for i, model in enumerate(all_data):
    pae_matrix = np.array(model['pae'])
    pred_len = pae_matrix.shape[0]
    if i == len(all_data)-1: color_bar_bool = True
    sns.heatmap(pae_matrix, cmap='viridis', ax=axes[i], square=True, cbar=color_bar_bool, cbar_ax=axes[-1], xticklabels=False, yticklabels=False, cbar_kws={"shrink":0.2})
    axes[i].set_title(name+n_model)
    axes[i].set_aspect('equal', adjustable='box')
    previous_l = 0
    for l in proteins_lengths:
        rect = patches.Rectangle((previous_l, 0), l, -40, linewidth=1, edgecolor='r', facecolor='none', clip_on=False)
        axes[i].add_patch(rect)
        previous_l += l
    n_model = str(int(n_model)+1)  

axes[-1].set_aspect(1.5)

plt.tight_layout()

# for model in get_full_data(directory_path):
#     pae_matrix = np.array(model['pae'])
#     sns.heatmap(pae_matrix, cmap='viridis')
#     ax = plt.gca()
#     rect = patches.Rectangle((0, 0), 1000, -40, linewidth=1, edgecolor='r', facecolor='none', clip_on=False)
#     ax.add_patch(rect)
#     plt.title(name+n_model)
#     plt.savefig(name+n_model+'.png')
#     plt.clf()
#     n_model = str(int(n_model)+1)  



# # Show the plot
plt.show()