import pandas as pd
import numpy as np
import seaborn as sns
from strcture_alignment import Alignment
import scipy.cluster.hierarchy as sch
from collections import defaultdict
import os
from PIL import Image
import io
from Bio import SeqIO
from seqfold import fold, dg, dg_cache, dot_bracket
import RNA  # https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/api_python.html
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
import forgi  # install Cython before forgi !!!!!
import multiprocessing as mp
import time
from multiprocessing import Pool
from joblib import Parallel, delayed


# mfold
# http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi?PAGE=3&ID=Lgs50derfi
# https://alphafoldserver.com/
# https://github.com/LinearFold/LinearFold
# https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/MaxExpect/MaxExpect.html
# https://viennarna.github.io/forgi/graph_tutorial.htmlpip


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)

# pip install ViennaRNA
# pip install seqfold


def get_relative_cdna_position(strand, pos, five_start, five_end, three_start, three_end):
    if strand == '-':
        pass
        # three_start, three_end = three_end, three_start
        # five_start, five_end = five_end, five_start

    three_len = three_end - three_start
    cds_len = three_start - five_end
    five_len = five_end - five_start
    if three_len == 0 or five_len == 0 or cds_len == 0:
        return np.nan
    if pos > three_start:
        return 2 + (pos - three_start) / three_len
    elif pos > five_end:
        return 1 + (pos - five_end) / cds_len
    else:
        return 0 + (pos - five_start) / five_len


def get_seq_from_fasta(target_transcript_id="ENST00000000233",
                       fasta_file="Homo_sapiens.GRCh38.cdna.all.fa/Homo_sapiens.GRCh38.cdna.all.fa"):
    # Path to your downloaded FASTA file
    # The transcript ID you're searching for
    # target_transcript_id = "ENST00000000233"  # Replace with your actual transcript ID

    # Search for the transcript
    found_transcript = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.description.split()[0].split('.')[0] == target_transcript_id:
            found_transcript = record
            break

    # Print the result
    if found_transcript:
        return found_transcript.seq
    else:
        return ""


def get_mfe_structur_from_seqfold(sequence):
    print("get_mfe_structur_from_seqfold")
    # sequence = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
    sequence = sequence.replace("T", "U")
    # just returns minimum free energy
    mfe = dg(sequence, temp=37.0)  # -13.4
    # print(mfe)
    # `fold_object` returns a list of `seqfold.Struct` from the minimum free energy structure
    structs = fold_object(sequence)
    # List[Struct]: A list of structures. Stacks, bulges, hairpins, etc.
    # print(structs)

    mfe = sum(s.e for s in structs)
    print(mfe)  # -13.4, same as dg()

    """for struct_elem in structs:
        break
        print(struct_elem.ij)  # prints the i, j, ddg, and description of each structure
        print(struct_elem.desc)
        # print(struct_elem.fmt)
        print(struct_elem.e)"""

    # `dot_bracket` returns a dot_bracket representation of the folding
    mfe_structure = dot_bracket(seq=sequence, structs=structs)
    print(mfe_structure)  # ((((((((.((((......))))..((((.......)))).))))))))

    # `dg_cache` returns a 2D array where each (i,j) combination returns the MFE from i to j inclusive
    cache = dg_cache(sequence)
    cache = np.asarray(cache)
    cache[np.where(np.isinf(cache))] = 0
    temp = np.where(cache > 0)
    # x = np.asarray(list(zip(list(temp[0]), list(temp[1]))))
    # print(x)

    G = nx.Graph()
    # G = nx.from_numpy_array(cache)
    for idx, i in enumerate(sequence):
        G.add_node(str(idx + 1))

    e_ = []
    for struct_elem in structs:
        i, j = struct_elem.ij[0][0], struct_elem.ij[0][1]
        # print(struct_elem.ij)  # prints the i, j, ddg, and description of each structure
        # print(struct_elem.desc)
        # print(struct_elem.e)

        if "STACK" in struct_elem.desc or "BIFURCATION" in struct_elem.desc:
            e_.append(abs(struct_elem.e))
            G.add_edge( str(i + 1), str(j + 1), weight=abs(struct_elem.e))
        if "HAIRPIN" in struct_elem.desc:
            G.add_edge(str(i + 1), str(j + 1), weight=abs(struct_elem.e) / 2)
        else:
            G.add_edge(str(i + 1), str(j + 1), weight=abs(struct_elem.e))
        # print(struct_elem.ij)  # prints the i, j, ddg, and description of each structure
        # print(struct_elem.desc)
        # print(struct_elem.fmt)
        # print(struct_elem.e)
    for idx, b in enumerate(sequence[:-1]):
        G.add_edge(str(idx + 1), str(idx + 2), weight=np.mean(e_))

    node_pos = nx.kamada_kawai_layout(G)  # Layout for better visualization

    # print(node_pos)

    # color_array = [color_map[b] for b in sequence]
    color_array = ["lightgray"] * len(sequence)
    for idx in range(left_window - 2, left_window + 3):
        color_array[idx] = color_map[sequence[idx]]

    plt.subplot2grid((1, 2), (0, 0))
    # print(len(list(G.nodes)))
    # print(len(color_array))
    nx.draw(G, node_pos,
            with_labels=True,
            node_size=70,
            node_color=color_array,
            font_size=6,
            font_weight='regular')
    plt.title('seqfold', pad=0)


    # cache = [[0 for x in line if np.isinf(x)] for line in cache]
    # print(cache)
    # print(np.max(cache))
    # sns.heatmap(cache)
    # plt.show()

    return mfe_structure


def get_mfe_structur_from_vianna_rna(sequence_):
    # print("get_mfe_structur_from_vianna_rna")

    # sequence = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
    sequence_ = sequence_.replace("T", "U")

    # Fold the RNA sequence
    """
    The differences in the 2D structures produced by RNA.fold_object and 
    RNA.fold_compound in the ViennaRNA Package arise from their 
    underlying mechanisms and intended use cases.
    
    RNA.fold_object: This function is designed to take a single RNA sequence as input and directly compute its minimum 
    free energy (MFE) structure. It simplifies the process for users who need a quick prediction for a single RNA strand
    
    (pred_structure, mfe) = RNA.fold_object(sequence)
    mfe = round(mfe, 2)
    
    RNA.fold_compound: This function allows for more complex inputs, including single sequences and sequence alignments. 
    It is designed to handle multiple sequences and can compute MFE for these alignments, potentially leading to 
    different structural predictions due to the consideration of interactions among multiple RNA strands.

    """

    fc = RNA.fold_compound(sequence_)
    (ss, ensemble_energy) = fc.pf()
    (centroid_struct, dist) = fc.centroid()
    mfe_centroid = fc.eval_structure(centroid_struct)
    mfe_centroid = round(mfe_centroid, 4)
    return mfe_centroid, centroid_struct, ss, fc, ensemble_energy


def define_rna_element(dot_bracket_seq):
    """
    In forgi, the RNA secondary structure is divided into the following elements:

    f0: The unpaired nucleotides at the 5' end of the molecule.
    s0, s1, s2: Stems, which are regions of contiguous canonical Watson-Crick base-paired nucleotides.
    h0, h1: Hairpin loops, which are single-stranded regions enclosed by a stem.
    m0: A multiloop segment, which is a single-stranded region between two stems.
    i0: An interior loop, which is a bulged out region or a loop with unpaired bases on either strand, flanked by stems on either side.
    t0: The unpaired nucleotides at the 3' end of the molecule.
    """
    rna_structure = forgi.load_rna(dot_bracket_seq)
    return "".join([rna_structure[0].get_elem(elem) for elem in range(1, len(dot_bracket_seq) + 1)])


def plot_rna_structure_4_seq(centroid_struct):
    # fc from get_mfe_structur_from_vianna_rna()
    fc = fold_object

    # Get the centroid structure from get_mfe_structur_from_vianna_rna()
    """
    The centroid structure in RNA secondary structure prediction refers to a 
    representative structure derived from a Boltzmann-weighted ensemble of possible RNA conformations. 
    This concept is particularly useful in improving the accuracy of RNA structure predictions compared 
    to traditional methods that focus solely on minimum free energy (MFE) structures.
    """

    # Compute partition function and base pair probabilities
    # need ensemble_energy for fc.centroid
    # ensemble_energy = ensembl_es[structures_index - 1]

    # Convert structure to pair table format
    pt = RNA.ptable(centroid_struct)

    # Evaluate the energy of the structure in pair table format
    # Function to get the energy of a specific base pair
    def get_base_pair_energy(fc, pt, i, j):
        # Temporarily break the base pair
        pt[i] = 0
        pt[j] = 0
        # Evaluate the energy of the structure without the base pair
        energy_without_bp = fc.eval_structure_pt(pt)
        # Restore the base pair
        pt[i] = j
        pt[j] = i
        # The energy contribution of the base pair is the difference
        energy_with_bp = fc.eval_structure_pt(pt)
        return energy_with_bp - energy_without_bp

    structure_plot = nx.Graph()
    for idx_ in range(len(centroid_struct)):
        structure_plot.add_node(str(idx_ + start))

    # print("\nEnergy contributions of base pairs:")
    e_ = []
    found_bounds = False
    for i in range(1, len(pt)):
        if pt[i] > 0:
            energy_bp = get_base_pair_energy(fc, pt, i, pt[i])
            if energy_bp == 0:  # replace 0 with 1 otherwise 2D structure is compromised
                energy_bp = 1
            e_.append(abs(energy_bp))
            # print(f"Base pair ({i},{pt[i]}): {energy_bp:.2f} kcal/mol")
            found_bounds = True
            edge = str(i + start), str(pt[i] + start)
            structure_plot.add_edge(str(i + start - 1), str(pt[i] + start - 1),
                                    # weight=abs(energy_bp)
                                    weight=200  # between bases in stack
                                    )

    # print(e_)
    # print(np.mean(e_))
    for idx_ in range(len(centroid_struct) - 1):
        structure_plot.add_edge(str(idx_ + start), str(idx_ + start + 1),
                                # weight=np.mean(e_) * 0.55,
                                weight=100  # between bases in sequence
                                )  # default 0.55

    base_pos = nx.kamada_kawai_layout(structure_plot)  # Layout for better visualization

    color_array = ["lightgray"] * (len(sequence) + 2)
    color_array = color_array[:len(centroid_struct)]
    node_list = list(structure_plot.nodes())
    if len(node_list) > len(centroid_struct):
        structure_plot.remove_node(node_list[-1])


    # structure_plot.remove_nodes_from(list(nx.isolates(structure_plot)))
    for idx in range(pos - start - 2, pos - start + 3):
        color_array[idx] = color_map[sequence[idx]]

    if found_bounds or True:
        # plt.figure(figsize=(7, 7))
        # ViannaRNA
        # plt.title(f'{t_id} - {kmer}',
        #           pad=-5, fontsize=14, loc="right")

        # plt.text(x=0.1, y=1.1, s=f"Structure: {structures_index}\nWindow: {start}-{end - 1}\nMFE = {mfe} kcal/mol",
        #          fontsize=10, transform=ax1.transAxes, va="center")

        new_labels = []
        node_sizes = []

        # node_list = node_list[:-1]
        for idx, n in enumerate(node_list):
            if show_bases:
                new_labels.append(sequence[idx])
                node_sizes.append(50)
            elif pos - 2 <= int(n) <= pos + 2:
                # new_labels.append(n)
                # new_labels.append(sequence[col_idx])
                new_labels.append("")
                node_sizes.append(50)
            elif int(n) % 10 == 0:
                new_labels.append(n)
                node_sizes.append(30)
            elif idx == 0:
                new_labels.append(n)
                node_sizes.append(30)
            else:
                new_labels.append("")
                node_sizes.append(5)

        # new_labels[-1] = node_list[-1]
        """nx.draw_networkx_nodes(structure_plot,
                               base_pos,
                               node_color=color_array,
                               node_size=node_sizes,
                               )"""
        while len(color_array) > len(node_sizes):
            node_sizes.append(10)

        while len(color_array) < len(node_sizes):
            node_sizes = node_sizes[:-1]

        nx.draw_networkx(structure_plot,
                         base_pos,
                         labels=dict(zip(node_list, new_labels)),
                         with_labels=True,
                         node_color=color_array,
                         font_size=7,
                         node_size=node_sizes,
                         font_weight='regular')

        plt.axis('off')
        # nx.draw_networkx(structure_plot, pos=base_pos, labels=dict(zip(node_list, new_labels)))
        # nx.relabel_nodes(structure_plot, dict(zip(node_list, new_labels)), copy=False)
    return found_bounds


def plot_rna_structures(centroid_struct):
    # fc from get_mfe_structur_from_vianna_rna()
    fc = folds[structures_index - 1]

    # Get the centroid structure from get_mfe_structur_from_vianna_rna()
    """
    The centroid structure in RNA secondary structure prediction refers to a 
    representative structure derived from a Boltzmann-weighted ensemble of possible RNA conformations. 
    This concept is particularly useful in improving the accuracy of RNA structure predictions compared 
    to traditional methods that focus solely on minimum free energy (MFE) structures.
    """

    # Compute partition function and base pair probabilities
    # need ensemble_energy for fc.centroid
    # ensemble_energy = ensembl_es[structures_index - 1]

    # Convert structure to pair table format
    pt = RNA.ptable(centroid_struct)

    # Evaluate the energy of the structure in pair table format
    # Function to get the energy of a specific base pair
    def get_base_pair_energy(fc, pt, i, j):
        # Temporarily break the base pair
        pt[i] = 0
        pt[j] = 0
        # Evaluate the energy of the structure without the base pair
        energy_without_bp = fc.eval_structure_pt(pt)
        # Restore the base pair
        pt[i] = j
        pt[j] = i
        # The energy contribution of the base pair is the difference
        energy_with_bp = fc.eval_structure_pt(pt)
        return energy_with_bp - energy_without_bp

    structure_plot = nx.Graph()
    for idx_ in range(len(centroid_struct)):
        structure_plot.add_node(str(idx_ + start))

    # print("\nEnergy contributions of base pairs:")
    e_ = []
    found_bounds = False
    for i in range(1, len(pt)):
        if pt[i] > 0:
            energy_bp = get_base_pair_energy(fc, pt, i, pt[i])
            if energy_bp == 0:  # replace 0 with 1 otherwise 2D structure is compromised
                energy_bp = 1
            e_.append(abs(energy_bp))
            # print(f"Base pair ({i},{pt[i]}): {energy_bp:.2f} kcal/mol")
            found_bounds = True
            edge = str(i + start), str(pt[i] + start)
            structure_plot.add_edge(str(i + start - 1), str(pt[i] + start - 1),
                                    # weight=abs(energy_bp)
                                    weight=200  # between bases in stack
                                    )

    # print(e_)
    # print(np.mean(e_))
    for idx_ in range(len(centroid_struct) - 1):
        structure_plot.add_edge(str(idx_ + start), str(idx_ + start + 1),
                                # weight=np.mean(e_) * 0.55,
                                weight=100  # between bases in sequence
                                )  # default 0.55

    base_pos = nx.kamada_kawai_layout(structure_plot)  # Layout for better visualization

    color_array = ["lightgray"] * (len(sequence) + 2)
    color_array = color_array[:len(centroid_struct)]
    node_list = list(structure_plot.nodes())
    if len(node_list) > len(centroid_struct):
        structure_plot.remove_node(node_list[-1])


    # structure_plot.remove_nodes_from(list(nx.isolates(structure_plot)))
    for idx in range(pos - start - 2, pos - start + 3):
        color_array[idx] = color_map[sequence[idx]]

    if found_bounds or True:
        # plt.figure(figsize=(7, 7))
        # ViannaRNA
        plt.title(f'{t_id} - {kmer}',
                  pad=-5, fontsize=14, loc="right")

        plt.text(x=0.1, y=1.1, s=f"Structure: {structures_index}\nWindow: {start}-{end - 1}\nMFE = {mfe} kcal/mol",
                 fontsize=10, transform=ax1.transAxes, va="center")

        new_labels = []
        node_sizes = []

        # node_list = node_list[:-1]
        for idx, n in enumerate(node_list):
            if show_bases:
                new_labels.append(sequence[idx])
                node_sizes.append(50)
            elif pos - 2 <= int(n) <= pos + 2:
                # new_labels.append(n)
                # new_labels.append(sequence[col_idx])
                new_labels.append("")
                node_sizes.append(50)
            elif int(n) % 10 == 0:
                new_labels.append(n)
                node_sizes.append(30)
            elif idx == 0:
                new_labels.append(n)
                node_sizes.append(30)
            else:
                new_labels.append("")
                node_sizes.append(5)

        # new_labels[-1] = node_list[-1]
        """nx.draw_networkx_nodes(structure_plot,
                               base_pos,
                               node_color=color_array,
                               node_size=node_sizes,
                               )"""
        while len(color_array) > len(node_sizes):
            node_sizes.append(10)

        while len(color_array) < len(node_sizes):
            node_sizes = node_sizes[:-1]

        nx.draw_networkx(structure_plot,
                         base_pos,
                         labels=dict(zip(node_list, new_labels)),
                         with_labels=True,
                         node_color=color_array,
                         font_size=7,
                         node_size=node_sizes,
                         font_weight='regular')

        plt.axis('off')
        # nx.draw_networkx(structure_plot, pos=base_pos, labels=dict(zip(node_list, new_labels)))
        # nx.relabel_nodes(structure_plot, dict(zip(node_list, new_labels)), copy=False)
    return found_bounds


def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))
    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [int(den[label][i]) for i in l]
        if len(i_l) > 3:
            cluster_classes[c] = i_l
    del_c = []
    for k in cluster_classes.keys():
        if len(cluster_classes[k]) < 8 and max(cluster_classes[k]) - min(cluster_classes[k]) > 8:
            del_c.append(k)
    for c_ in del_c:
        cluster_classes.pop(c_)
    if "C0" in cluster_classes.keys():
        cluster_classes.pop("C0")

    for k in cluster_classes.keys():
        cluster_classes[k] = sorted(cluster_classes[k])
    return cluster_classes


def plot_heatmap_of_alignment():
    # cluster_fig = plt.figure(2)

    # plt.title('Alignment score')
    labels = [" "] * len(score_matrix)
    labels[structures_index - 1] = str(structures_index)
    ax = sns.heatmap(score_matrix, cbar_kws={'ticks': [0, 0.5, 1], 'pad': 0.03},
                     yticklabels=labels,
                     xticklabels=labels)

    cbar = ax.collections[0].colorbar
    cbar.ax.yaxis.set_tick_params(labelsize=7)
    ax.yaxis.set_tick_params(labelsize=9)
    ax.xaxis.set_tick_params(labelsize=9)

    cbar.set_label('Similarity score', rotation=270, labelpad=15, fontsize=9)
    # Visualize the clustered matrix
    plt.hlines(y=structures_index - 0.5, xmin=0, xmax=structures_index - 0.5, linestyle='solid',
               color="steelblue", linewidth=1)
    plt.vlines(x=structures_index - 0.5, ymin=structures_index - 0.5, ymax=len(score_matrix), linestyle='solid',
               color="steelblue", linewidth=1)
    cluster_mfes = dict()
    for name, c_ in cluster.items():
        min_ = min(c_)
        max_ = max(c_) + 1
        center = (max_ + min_) / 2
        rect = plt.Rectangle(xy=(min_, min_), width=max_ - min_, height=max_ - min_,
                             fill=False, edgecolor='steelblue', linewidth=1.5)
        ax.add_patch(rect)
        avg_mfe = round(np.mean(energies_[min_: max_]), 2)
        cluster_mfes[name] = avg_mfe
        """plt.text(x=center, y=center, s=round(avg_mfe, 1), ha='center', va='center',
                 fontsize=5, color="black", fontweight='bold',
                 bbox=dict(facecolor='white', edgecolor='none',
                           boxstyle='round', alpha=0.6))"""


def plot_mfe_profile():

    pos_energies_ = [-x for x in energies_]
    plt.fill_between(x=range(1, len(energies_) + 1),
                     y1=pos_energies_,
                     alpha=0.5, color='steelblue')
    plt.vlines(x=structures_index, ymin=0,
               ymax=pos_energies_[structures_index - 1],
               color='steelblue', linewidth=3, zorder=3)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)

    plt.text(x=structures_index, y=pos_energies_[structures_index - 1] * 1.02,
             s=energies_[structures_index - 1], ha='center', va='bottom', fontsize=8.5, zorder=3,
             bbox=props)
    c_name = 1
    for name, c_ in cluster.items():
        c = [int(x) for x in c_]
        min_ = min(c) + 1
        max_ = max(c) + 1
        center = (max_ + min_) / 2
        if c_name != 1:
            plt.axvline(x=min_, color='gray', linestyle='--')
            plt.axvline(x=max_, color='gray', linestyle='--')
        else:
            plt.axvline(x=max_, color='gray', linestyle='--')
        plt.text(x=center, y=max(pos_energies_) * 0.04, s=str(c_name), ha='center', va='bottom', fontsize=6.5)
        c_name += 1
    # plt.axvline(x=max_, color='gray', linestyle='--')
    plt.xlim(1, len(energies_))
    plt.ylim(0, max(pos_energies_) * 1.12)
    plt.ylabel("-MFE (kcal/mol)")
    plt.xlabel("RNA structure index")


def add_params_for_window():
    for i in range(window - crop, crop, -steps):
        left_window = i
        if pos - left_window < 0:
            continue
        start = pos - left_window
        end = start + window
        seq_temp = seq[start: end].replace("T", "U")
        # print("get_mfe_structur_from_vianna_rna for", i)
        e, s, ps, fold_object, ensembl_e = get_mfe_structur_from_vianna_rna(seq_temp)
        # print(s)
        energies_.append(e)
        structures_.append(s)
        pseudostructures.append(ps)
        starts.append(start)
        ends.append(end)
        seq_window.append(seq_temp)
        folds.append(fold_object)
        ensembl_es.append(ensembl_e)


def calculate_score(i, j, structures_, energies_):
    if abs(energies_[i] - energies_[j]) < 5 and (energies_[i] < -10 and energies_[j] < -10):
        score = Alignment(structures_[i], structures_[j]).needleman_wunsch()
    else:
        score = 0
    return score


class Align:

    def __init__(self, rna_structures : list, energies : list):

        self.dim = len(rna_structures)
        self.score_matrix = np.ones((self.dim, self.dim))
        self.rna_structures = rna_structures
        self.energies = energies


    def compute_score_matrix(self):

        # Create a list of arguments for each pair of structures
        args_list = []
        # Prepare the tasks
        for i in range(self.dim - 1):
            for j in range(i + 1, self.dim):
                args_list.append((i, j, self.rna_structures, self.energies))

        # Use multiprocessing to compute scores
        print("CPUs: ", mp.cpu_count())
        with mp.Pool(mp.cpu_count()) as pool:
            scores = pool.starmap(calculate_score, args_list)

        # Fill the score matrix with results
        idx = 0
        for i in range(self.dim - 1):
            for j in range(i + 1, self.dim):
                score = scores[idx]
                self.score_matrix[i, j] = score
                self.score_matrix[j, i] = score
                idx += 1

        return self.score_matrix


# Assuming Alignment is defined elsewhere
def compute_alignment(i, j):
    # return Alignment(structures_[i][crop:-crop], structures_[j][crop:-crop]).needleman_wunsch()
    compare_window = 40
    start_temp_i = kmer_pos_in_window[i] - compare_window  # default 40
    start_temp_j = kmer_pos_in_window[j] - compare_window

    if start_temp_i < 0:
        start_temp_i = 0
    if start_temp_j < 0:
        start_temp_j = 0

    return Alignment(structures_[i][start_temp_i:kmer_pos_in_window[i] + compare_window],
                     structures_[j][start_temp_j:kmer_pos_in_window[j] + compare_window]).needleman_wunsch()


if __name__ == '__main__':

    #####
    """
    test
     e, s, pvals, fold_object, ensembl_e = get_mfe_structur_from_vianna_rna(seq_temp)
    """

    color_map = {"A": "#109648",
                 "C": "#255C99",
                 "G": "#F7B32B",
                 "U": "#D62839",
                 }

    ### test section  ###

    # test_seq = sequence = "TATTTCTAGCTCCACAAGTGTGTGGCCCCGCCCGAGCCCCTGCCCACGCCCTTGGAGCCTTCCACCGGCACTCATGACGGCCTGCCTGCAAACCTGCTGGTGGGGCAGACCCGAAAATCCAGCGTGCACCCCGCCGGAGGAAGGTCCCAT"
    # e, centroid_struct_, ps, fold_object, ensembl_e = get_mfe_structur_from_vianna_rna(test_seq)
    # print(e, centroid_struct_, ps, fold_object, ensembl_e)
    # print(test_seq[53:93])
    # print(centroid_struct_[53:93])
    # start = 0
    # pos = 74
    # kmer = "ATGA"
    # show_bases = True
    # sequence = sequence.replace("T", "U")
    # plot_rna_structure_4_seq(centroid_struct_)
    # plt.show()


    xpore_file = "../Nanopore/direct_RNA_seq/Xpore/DATA/"\
                 "xpore_cdna.all_group_compare_2025/"\
                 "diffmod_s3-s4/"\
                 "majority_direction_kmer_diffmod_padj_modrate_025_GENOME_cDNA_UTR.tsv"
    xpore_data = pd.read_table(xpore_file)

    # filter
    ens_id_filer = ["ENST00000417122",
                    # "ENST00000617924", "ENST00000225964"
                    ]
    pos_filter = [3092,
                  # 5043, 5517, 5561
                  ]
    kmer_filter = ["GUUUC",
                   ]

    # xpore_data = xpore_data[xpore_data["len"] < 300]
    #  seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
    # seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCCGTGACCTCCCTGGGTGA"  #  + A not working anymore


    """
    "ENST00000216037" in t_id and kmer == "GUUUU" and pos == 1629 > use for example of switch !!!! 
    """

    xpore_data = xpore_data.dropna(subset="relative_cdna_pos")
    print(xpore_data.info())

    plot_structures = True
    show_bases = False
    save_svg = True
    overwrite = False
    take_random_pos = False
    close_gaps = False

    params = {150: (0.85, 20, 4),
              170: (0.85, 20, 4),
              200: (.8, 30, 6),
              300: (0.8, 40, 8)}  # (cluster_thresh, energy diff, min_num_cluster)


    for window in [150]:  # default 150
        print(window)
        cluster_thresh = params[window][0]
        steps = 2  # default 2
        crop = 10  # default 10

        if take_random_pos:
            folder_out = f"RNAfolds/{window}bp_window_RANDOM_CDS_3UTR"
        else:
            folder_out = f"RNAfolds/{window}bp_window_run3_corrected"

        os.mkdir(folder_out) if not os.path.exists(folder_out) else None

        # plot_extra
        # tbl = pd.read_csv('RNAfolds/100bp_window/seq_params_best_5_combined.csv')
        # tbl = tbl[(tbl["rel_pos"] >= .8) & (tbl.kmer_in_structure == "sss")]
        # print(tbl)
        num_kmers = len(xpore_data)
        time_needed = []
        first_plot = True
        for seq, pos, kmer, t_id, cdna_pos, five_utr, three_utr in zip(xpore_data.cdna_seq,
                                                                       xpore_data.position,
                                                                       xpore_data.kmer,
                                                                       xpore_data.id,
                                                                       xpore_data.relative_cdna_pos,
                                                                       xpore_data.cDNA_end_fiveUTR,
                                                                       xpore_data.cDNA_start_threeUTR
                                                                       ):

            print(kmer, pos, t_id)
            pos_temp = pos
            kmer_temp = kmer
            if take_random_pos:
                if 1 > cdna_pos > 0:
                    if int(five_utr) > 0:
                        pos = np.random.randint(0, int(five_utr))
                elif 2 > cdna_pos >= 1:
                    pos = np.random.randint(int(five_utr), int(three_utr))
                elif cdna_pos >= 2:
                    if three_utr < len(seq) - 10:
                        pos = np.random.randint(int(three_utr), int(len(seq) - 10))
                else:
                    continue
                print("OLD_NEW_POS: ",  pos_temp, pos)
                kmer = seq[pos - 2: pos + 3].replace("T", "U")
                print(seq[pos], kmer)
                # xpore_data.loc[
                #     (xpore_data["id"] == t_id) &
                #     (xpore_data["kmer"] == kmer) &
                #     (xpore_data["position"] == pos_temp),
                #     "Random_kmer"] = new_kmer
                #
                # xpore_data.loc[
                #     (xpore_data["id"] == t_id) &
                #     (xpore_data["kmer"] == kmer) &
                #     (xpore_data["position"] == pos_temp),
                #     "Random_pos"] = pos

            start_time_main = time.time()
            num_kmers -= 1
            # GGAAC 1283 ENST00000323851.13
            # if "ENST00000226225" not in t_id or kmer != "ACUGG":
            # print(t_id, kmer)
            # ENST00000225964 > COL1A1
            plot_structures = False
            if not (t_id in ens_id_filer and kmer in kmer_filter and pos in pos_filter):
                continue
            if t_id in ens_id_filer and kmer in kmer_filter and pos in pos_filter:
                plot_structures = True
            if kmer != "GUUUU":
                pass

            print(num_kmers, "\n")
            t_id_name = t_id.replace(".", "_")
            # os.startfile(f"RNAfolds/{window}bp_window/{t_id_name}")
            # quit()
            if not os.path.exists(f"{folder_out}/{t_id_name}"):
                os.mkdir(f"{folder_out}/{t_id_name}")
            if not os.path.exists(f"{folder_out}/{t_id_name}/{kmer}_{pos}"):
                os.mkdir(f"{folder_out}/{t_id_name}/{kmer}_{pos}")

            if (os.path.exists(f"{folder_out}/{t_id_name}/{kmer}_{pos}/seq_params_best_10.csv")
                    and not overwrite and not plot_structures):
                continue

            output_main = f"{folder_out}/{t_id_name}/{kmer}_{pos}/seq_params.csv"
            if (os.path.exists(f"{folder_out}/{t_id_name}/{kmer}_{pos}/seq_params_best_10.csv")
                    and not overwrite and not plot_structures):
                seq_params = pd.read_csv(output_main)
            else:
                # fetch_paramsmfe_mean
                start_time = time.time()
                energies_, structures_, pseudostructures, starts, ends, seq_window = [], [], [], [], [], []
                folds, ensembl_es = [], []

                # MAIN FUNCTION
                add_params_for_window()

                print("Structure prediction finished in: ", time.time() - start_time, "s")

                if len(structures_) == 0:
                    continue
                seq_params = pd.DataFrame({"seq_window": seq_window,
                                           "structure": structures_,
                                           "pseudo_structure": pseudostructures,
                                           "mfe": energies_,
                                           "start": starts,
                                           "end": ends})
                seq_params["kmer"] = kmer
                seq_params["pos"] = pos
                if take_random_pos:
                    seq_params["kmer_org"] = kmer_temp
                    seq_params["pos_org"] = pos_temp

                seq_params["kmer_pos_in_window"] = seq_params["pos"] - seq_params["start"]
                # a = Align(rna_structures=structures_, energies=energies_)  # works worse!
                # score_matrix = a.compute_score_matrix()

            start_time = time.time()
            structures_ = seq_params["structure"].to_list()
            kmer_pos_in_window = seq_params["kmer_pos_in_window"].to_list()

            energies_ = seq_params["mfe"].to_list()
            score_matrix = np.zeros((len(structures_), len(structures_)))
            avg_e = np.mean(energies_)

            # Pre-filter structures based on energy criteria
            filtered_indices = [
                i for i in range(len(energies_))
                if energies_[i] < avg_e * 0.25  # default 0.5
            ]

            energy_diff = params[window][1]  # default 10 !!! don't calculate alignment to accelerate analysis

            filtered_pairs = [
                (i, j) for idx, i in enumerate(filtered_indices)
                for j in filtered_indices[idx + 1:]
                if abs(energies_[i] - energies_[j]) < energy_diff
                   and (structures_[i][steps * (j - i):] != structures_[j][:-steps * (j - i)])
            ]

            equal_pairs = [
                (i, j) for i in range(len(structures_) - 1)
                for j in range(i, len(structures_))
                if structures_[i][steps * (j - i):] == structures_[j][:-steps * (j - i)]
            ]
            # print(equal_pairs)

            # Compute alignment scores in parallel
            scores = Parallel(n_jobs=-1)(delayed(compute_alignment)(i, j) for i, j in filtered_pairs)  # works best!
            # Fill the score matrix
            for (i, j), score in zip(filtered_pairs, scores):
                score_matrix[i, j] = score
                score_matrix[j, i] = score

            for i in range(len(structures_)):
                score_matrix[i, i] = 1

            for (i, j) in equal_pairs:
                # print(structures_[i][steps:])
                # print(structures_[j][:-steps])
                score_matrix[i, j] = 1
                score_matrix[j, i] = 1

            """for i in range(len(structures_) - 1):
                break
                for j in range(i + 1, len(structures_)):
                    if abs(energies_[i] - energies_[j]) < 3 and (energies_[i] < -10 and energies_[j] < -10):  # to speed up code
                        score = Alignment(structures_[i], structures_[j]).needleman_wunsch()
                    else:
                        score = 0
                    score_matrix[i, j] = score
                    score_matrix[j, i] = score"""

            print("Alignment finished in: ", time.time() - start_time, "s")
            start_time = time.time()
            # print(score_matrix)
            # Perform hierarchical clustering
            if len(score_matrix) < 2:
                continue
            # score_matrix[score_matrix < 0.4] = 0
            pairwise_distances = sch.distance.pdist(score_matrix)
            linkage_matrix = sch.linkage(pairwise_distances, method='single')  # default 'single'
            # print(linkage_matrix)

            clusters = sch.fcluster(linkage_matrix, t=cluster_thresh,  # default 0.85
                                    criterion='distance'
                                    )
            # print(clusters)
            cluster = dict()
            structure_ids = np.arange(0, len(structures_))
            cluster_name = 0
            for c in sorted(np.unique(clusters)):
                c_elements = list(np.where(clusters == c)[0])
                if len(c_elements) < params[window][2]:
                    continue
                cluster_name += 1
                # cluster[f"C{c}"] = list(structure_ids[c_elements])
                cluster[f"C{cluster_name}"] = list(structure_ids[c_elements])
            cluster = dict(sorted(cluster.items(), key=lambda x: min(x[1]), reverse=False))
            # print(cluster)

            # Create a dendrogram and get the order of rows/columns
            # dendro = sch.dendrogram(linkage_matrix, no_plot=False, color_threshold=0.65,)  # default 0.65  0.65*max(linkage_matrix[:, 2])
            # cluster = get_cluster_classes(dendro)

            seq_params["cluster"] = "."

            # get most largest and most stable structure ensembl
            cluster_name = 0
            for c, s in cluster.items():
                cluster_name += 1
                for elm in s:
                    seq_params.loc[seq_params.index == int(elm), "cluster"] = f"C{cluster_name}"

            # close_gaps

            if close_gaps:
                for c in cluster.keys():
                    c_temp = seq_params[seq_params["cluster"] == c]
                    c_min, c_max = c_temp.start.min(), c_temp.start.max()
                    seq_params.loc[seq_params.start.between(c_min, c_max), "cluster"] = c

            # seq_params = seq_params[seq_params["cluster"] != "."]
            # print(seq_params)

            # MFE is negative !!!
            mean_e_per_cluster = seq_params.groupby("cluster")["mfe"].quantile(0.1).reset_index().sort_values(
                by="mfe", ascending=True).rename(columns={"mfe": "mfe_mean"})

            std_e_per_cluster = seq_params.groupby("cluster")["mfe"].std().reset_index().sort_values(
                by="mfe", ascending=True).rename(columns={"mfe": "mfe_std"})
            counts = seq_params.groupby("cluster")["mfe"].count().reset_index().sort_values(
                by="mfe", ascending=True).rename(columns={"mfe": "mfe_counts"})
            stats_df = pd.merge(mean_e_per_cluster, std_e_per_cluster, on="cluster")
            stats_df = pd.merge(stats_df, counts, on="cluster")
            # print(mean_e_per_cluster)
            # print(std_e_per_cluster)
            # print(counts)
            stats_df = stats_df[stats_df["cluster"] != "."]
            # stats_df = stats_df[stats_df["mfe_mean"] <= stats_df["mfe_mean"].min() * 0.7]
            # stats_df = stats_df[stats_df.mfe_std < 4]
            # stats_df = stats_df.nsmallest(3, "mfe_mean")
            stats_df = stats_df.nlargest(3, "mfe_counts")
            stats_df.mfe_mean = abs(stats_df.mfe_mean)
            # stats_df.sort_values(by="mfe_mean", ascending=False, inplace=True)
            stats_df["mfe_mean_score"] = stats_df["mfe_mean"] / stats_df["mfe_mean"].max()
            # stats_df.sort_values(by="mfe_std", ascending=False, inplace=True)
            # stats_df["mfe_std_score"] = range(1, len(stats_df) + 1)
            stats_df["mfe_std"] = stats_df["mfe_std"] + 1
            stats_df["mfe_std_score"] = stats_df["mfe_std"].max() / stats_df["mfe_std"]
            stats_df["mfe_count_score"] = stats_df["mfe_counts"] / stats_df["mfe_counts"].max()
            # stats_df.sort_values(by="mfe_counts", ascending=True, inplace=True)
            # stats_df["mfe_count_score"] = range(1, len(stats_df) + 1)
            # stats_df["mfe_count_score"] = stats_df["mfe_count_score"] * 3
            stats_df["all_score"] = stats_df["mfe_mean_score"] + stats_df["mfe_count_score"]

            stats_df = stats_df.sort_values(by=["all_score", "mfe_mean"], ascending=False)
            print(stats_df)
            if stats_df.empty:
                continue
            mfe_cluster = stats_df.iloc[0, 0]
            print(mfe_cluster)
            # print(mfe_cluster)
            seq_params["is_canonical_cluster"] = False
            seq_params.loc[seq_params.cluster == mfe_cluster, "is_canonical_cluster"] = True
            best_10 = seq_params.loc[seq_params.is_canonical_cluster].nsmallest(10, 'mfe')
            best_10["forgi_structure"] = best_10["structure"].apply(define_rna_element)
            seq_params.loc[seq_params.is_canonical_cluster, "forgi_structure"] = seq_params.loc[seq_params.is_canonical_cluster, "structure"].apply(define_rna_element)

            seq_params.loc[seq_params.is_canonical_cluster, "forgi_structure"] = (
                seq_params.loc[seq_params.is_canonical_cluster, "forgi_structure"].apply(
                    lambda x: "".join(i for i in x if i not in [str(x) for x in range(0, 10)])))

            best_10["forgi_structure"] = best_10["forgi_structure"].apply(
                lambda x: "".join(i for i in x if i not in [str(x) for x in range(0, 10)]))

            best_10["rel_pos"] = best_10["pos"] - best_10["start"]
            best_10["kmer_in_structure"] = best_10.apply(lambda x: x["forgi_structure"][x["rel_pos"] - 2:
                                                                                        x["rel_pos"] + 3], axis=1)

            best_10["kmer2"] = best_10.apply(lambda x: x["seq_window"][x["rel_pos"] - 2:
                                                                       x["rel_pos"] + 3], axis=1)

            if not best_10[best_10["kmer"] != best_10["kmer2"]].empty:
                print("CAUTION !!!!! WRONG Position !!! ")

            seq_params["rel_pos"] = seq_params["pos"] - seq_params["start"]
            seq_params.loc[seq_params.is_canonical_cluster, "kmer_in_structure"] = (
                seq_params.loc[seq_params.is_canonical_cluster].apply(lambda x: x["forgi_structure"][x["rel_pos"] - 2:
                                                                                                     x["rel_pos"] + 3],
                                                                      axis=1))
            # print(best_10)
            # print(best_10)

            print("Clustering finished in: ", time.time() - start_time, "s")
            for col in seq_params.columns:
                if "Unnamed" in col:
                    seq_params.drop(col, axis=1, inplace=True)

            for col in best_10.columns:
                if "Unnamed" in col:
                    best_10.drop(col, axis=1, inplace=True)

            # print(best_10)
            seq_params.to_csv(f"{folder_out}/{t_id_name}/{kmer}_{pos}/seq_params.csv", index=False)
            best_10.to_csv(f"{folder_out}/{t_id_name}/{kmer}_{pos}/seq_params_best_10.csv", index=False)
            print(found_forgi := best_10.kmer_in_structure.unique())

            if plot_structures:
                print("plot_structures")
                structures_index = 0
                out_files = []

                for cs, sequence, mfe, start, end, is_cluster in zip(seq_params.structure, seq_params.seq_window,
                                                         seq_params.mfe, seq_params.start, seq_params.end,
                                                         seq_params.is_canonical_cluster):
                    if not first_plot:
                        break
                    structures_index += 1

                    if not is_cluster:
                        pass
                        # continue

                    if structures_index < 40:
                        pass
                    # print(structures_index)
                    # print(cs)
                    if cs not in best_10["structure"].values:
                        pass
                    plt.figure(figsize=(9, 9))
                    print(structures_index, "/", len(seq_params))
                    plt.subplots_adjust(left=0, right=0.95, top=0.9, bottom=0.25, wspace=0.45)
                    # ax1 = plt.subplot2grid((12, 5), (0, 0), colspan=3, rowspan=12)
                    ax1 = plt.subplot2grid((6, 10), (0, 0), colspan=10, rowspan=4)
                    plot_rna_structures(cs)
                    # ax2 = plt.subplot2grid((12, 5), (1, 3), colspan=2, rowspan=5)
                    ax2 = plt.subplot2grid((6, 10), (4, 6), colspan=3, rowspan=2)
                    plot_heatmap_of_alignment()
                    # ax3 = plt.subplot2grid((12, 5), (7, 3), colspan=2, rowspan=4)
                    ax3 = plt.subplot2grid((6, 10), (4, 1), colspan=4, rowspan=2)
                    plot_mfe_profile()
                    out_file = f"{folder_out}/{t_id_name}/{kmer}_{pos}/{kmer}_{pos}_{t_id}_{structures_index}.png"
                    out_files.append(out_file)
                    plt.savefig(out_file, dpi=300)
                    if save_svg:
                        plt.savefig(out_file.replace(".png", ".svg"))
                    # plt.show()
                    plt.close()

                # quit()
                # create GIF

                if first_plot:
                    frames = [Image.open(png_file) for png_file in out_files]
                    frames[0].save(f"{folder_out}/{t_id_name}/{kmer}_{pos}/{kmer}_{pos}_{t_id}.gif",
                                   save_all=True, append_images=frames[1:], loop=3, duration=750)
                    # os.startfile(os.getcwd() + f"/RNAfolds/{window}bp_window/{t_id_name}/{kmer}_{pos}/")
                first_plot = False
                os.startfile(os.getcwd() + f"/{folder_out}/{t_id_name}")
                # quit()

            time_needed.append(time.time() - start_time_main)
            print("\nestimated residual time:", (np.mean(time_needed) * num_kmers / 60) / 60, "h")
            # quit()
            # os.startfile(os.getcwd() + f"/RNAfolds/{window}bp_window/{t_id_name}")

    # xpore_data.to_csv(f"RNAfolds/xpore_C3_C4_RANDOM_positions_in_5_CDS_3.tsv", sep="\t",
    #                   index=False)
