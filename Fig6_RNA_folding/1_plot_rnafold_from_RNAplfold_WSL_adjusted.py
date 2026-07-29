# pip install networkx matplotlib
import math
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors, colormaps
import subprocess
import time
import forgi  # install Cython before forgi !!!!!
import pandas as pd


class Alignment:
    def __init__(self, seq1: str, seq2: str,
                 match: int = 1,
                 wobble: float = 0.5,
                 gap_open: float = -2,
                 gap_extend: float = -0.5,
                 gap_penalty: float = -2,
                 mismatch_penalty: float = -1,
                 only_score: bool = False,
                 allow_wobble: bool = False,
                 allow_gaps: bool = False):


        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
        self.wobble = wobble  # e.g. +0.5 or +1.0
        self.gap_open = gap_open  # e.g. -2.0 or -3.0
        self.gap_extend = gap_extend  # e.g. -0.5 or -1.0
        self.only_score = only_score
        self.allow_wobble = allow_wobble # True/False
        self.allow_gaps = allow_gaps # True/False

    def needleman_wunsch(self):
        # still reverse complement as input needed !!!
        complement = {"A": "U",
                      "U": "A",
                      "G": "C",
                      "C": "G"}
        # reverse complement for direct match
        self.seq2 = "".join([complement[b] for b in self.seq2][::-1])

        n, m = len(self.seq1), len(self.seq2)

        # DP score matrix
        score = [[0] * (m + 1) for _ in range(n + 1)]
        # Optional traceback matrix: 'D' diag, 'U' up, 'L' left
        tb = [[''] * (m + 1) for _ in range(n + 1)]

        # Initialize first row/column
        for i in range(1, n + 1):
            score[i][0] = i * self.gap_penalty
            tb[i][0] = 'U'
        for j in range(1, m + 1):
            score[0][j] = j * self.gap_penalty
            tb[0][j] = 'L'

        # Fill
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                sdiag = score[i - 1][j - 1] + (self.match
                                               if (self.seq1[i - 1] == self.seq2[j - 1])

                                               else self.mismatch_penalty)

                if self.allow_wobble and self.seq1[i - 1] != self.seq2[j - 1]:
                    sdiag = score[i - 1][j - 1] + (self.wobble + self.mismatch_penalty
                                                   if (((self.seq1[i - 1] == "U") & (self.seq2[j - 1] == "G")) |
                                                       ((self.seq1[i - 1] == "G") & (self.seq2[j - 1] == "U")))

                                                   else self.mismatch_penalty)
                sup = score[i - 1][j] + self.gap_penalty    # gap_penalty in seq2
                sleft = score[i][j - 1] + self.gap_penalty    # gap_penalty in seq1
                best = max(sdiag, sup, sleft)
                score[i][j] = best
                # tie-breaking order: diag > up > left
                if best == sdiag:
                    tb[i][j] = 'D'
                elif best == sup:
                    tb[i][j] = 'U'
                else:
                    tb[i][j] = 'L'

        alignment_score = score[n][m]

        if self.only_score:
            return alignment_score, None, None

        # Traceback
        align1, align2 = [], []
        i, j = n, m
        while i > 0 or j > 0:
            move = tb[i][j] if i >= 0 and j >= 0 else ''
            if i > 0 and j > 0 and move == 'D':
                align1.append(self.seq1[i - 1])
                align2.append(self.seq2[j - 1])
                i -= 1; j -= 1
            elif i > 0 and (j == 0 or move == 'U'):
                align1.append(self.seq1[i - 1])
                align2.append('-')
                i -= 1
            else:  # j > 0 and (i == 0 or move == 'L')
                align1.append('-')
                align2.append(self.seq2[j - 1])
                j -= 1

        align1 = ''.join(reversed(align1))
        align2 = ''.join(reversed(align2))

        # Optional identity metrics
        matches = sum(a == b and a != '-' for a, b in zip(align1, align2))
        aln_len = len(align1)
        gaps = sum(a == '-' or b == '-' for a, b in zip(align1, align2))
        pid = matches / aln_len if aln_len else 0.0

        return {
            "score": alignment_score,
            "align1": align1,
            "align2": align2,
            "percent_identity": pid,
            "alignment_length": aln_len,
            "matches": matches,
            "gaps": gaps,
        }

    def smith_waterman(self):

        # still reverse complement as input needed !!!
        complement = {"A": "U",
                      "U": "A",
                      "G": "C",
                      "C": "G"}
        # reverse complement for direct match
        self.seq2 = "".join([complement[b] for b in self.seq2][::-1])

        n, m = len(self.seq1), len(self.seq2)
        score = [[0] * (m + 1) for _ in range(n + 1)]
        tb = [[None] * (m + 1) for _ in range(n + 1)]  # traceback

        max_score = 0
        max_pos = (0, 0)

        # Fill the score matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag = score[i - 1][j - 1] + (self.match if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty)
                up = score[i - 1][j] + self.gap_penalty
                left = score[i][j - 1] + self.gap_penalty
                score[i][j] = max(0, diag, up, left)
                if score[i][j] == 0:
                    tb[i][j] = None
                elif score[i][j] == diag:
                    tb[i][j] = 'D'
                elif score[i][j] == up:
                    tb[i][j] = 'U'
                else:
                    tb[i][j] = 'L'
                if score[i][j] > max_score:
                    max_score = score[i][j]
                    max_pos = (i, j)

        # Traceback from max_pos
        i, j = max_pos
        align1, align2 = [], []
        while i > 0 and j > 0 and tb[i][j]:
            if tb[i][j] == 'D':
                align1.append(self.seq1[i - 1])
                align2.append(self.seq2[j - 1])
                i -= 1
                j -= 1
            elif tb[i][j] == 'U':
                align1.append(self.seq1[i - 1])
                align2.append('-')
                i -= 1
            elif tb[i][j] == 'L':
                align1.append('-')
                align2.append(self.seq2[j - 1])
                j -= 1
            if score[i][j] == 0:
                break

        align1 = ''.join(reversed(align1))
        align2 = ''.join(reversed(align2))

        return {
            "score": max_score,
            "align1": align1,
            "align2": align2,
            "start_in_seq1": i,
            "end_in_seq1": max_pos[0],
        }

    def smith_waterman_affine(self):
        """
        Local alignment (Smith–Waterman).
        - self.allow_gaps=False → gapless local alignment (diagonal only)
        - self.allow_gaps=True  → affine gaps (gap_open + gap_extend)
        - self.allow_wobble=True/False toggles G:U scoring
        """

        def is_wc(a, b):
            return (a == 'A' and b == 'U') or (a == 'U' and b == 'A') or (a == 'G' and b == 'C') or (
                        a == 'C' and b == 'G')

        def is_gu(a, b):
            return (a == 'G' and b == 'U') or (a == 'U' and b == 'G')

        s1, s2 = self.seq1, self.seq2
        n, m = len(s1), len(s2)

        # helper: substitution score
        def sub_score(a, b):
            if is_wc(a, b):
                return float(self.match)
            if is_gu(a, b):
                return float(self.wobble) if getattr(self, "allow_wobble", True) else float(self.mismatch_penalty)
            return float(self.mismatch_penalty)

        # ========== GAPLESS (no gaps allowed) ==========
        if not getattr(self, "allow_gaps", True):
            S = [[0.0] * (m + 1) for _ in range(n + 1)]  # score
            CONT = [[False] * (m + 1) for _ in range(n + 1)]  # continue-diagonal flag

            best = 0.0
            bi = bj = 0  # best-end indices (ints)

            for i in range(1, n + 1):
                ai = s1[i - 1]
                for j in range(1, m + 1):
                    bjc = s2[j - 1]
                    diag = S[i - 1][j - 1] + sub_score(ai, bjc)
                    if diag > 0.0:
                        S[i][j] = diag
                        CONT[i][j] = True
                        if diag > best:
                            best = diag
                            bi, bj = i, j
                    else:
                        S[i][j] = 0.0
                        CONT[i][j] = False

            # Traceback using separate variables to avoid collisions
            ti, tj = int(bi), int(bj)
            a1, a2 = [], []
            while ti > 0 and tj > 0 and CONT[ti][tj]:
                a1.append(s1[ti - 1])
                a2.append(s2[tj - 1])
                ti -= 1
                tj -= 1

            a1 = ''.join(reversed(a1))
            a2 = ''.join(reversed(a2))
            return {
                "score": best,
                "align1": a1,
                "align2": a2,
                "start_in_seq1": ti,
                "start_in_seq2": tj,
                "end_in_seq1": bi,
                "end_in_seq2": bj,
            }

        # ========== AFFINE GAPS (allow_gaps=True) ==========
        go = float(self.gap_open)
        ge = float(self.gap_extend)

        M = [[0.0] * (m + 1) for _ in range(n + 1)]
        X = [[0.0] * (m + 1) for _ in range(n + 1)]
        Y = [[0.0] * (m + 1) for _ in range(n + 1)]
        PM = [[None] * (m + 1) for _ in range(n + 1)]
        PX = [[None] * (m + 1) for _ in range(n + 1)]
        PY = [[None] * (m + 1) for _ in range(n + 1)]

        best_score = 0.0
        best_matrix = 'M'
        best_pos = (0, 0)

        for i in range(1, n + 1):
            ai = s1[i - 1]
            for j in range(1, m + 1):
                bjc = s2[j - 1]
                sub = sub_score(ai, bjc)

                x_open = M[i - 1][j] + go + ge
                x_ext = X[i - 1][j] + ge
                X[i][j] = max(0.0, x_open, x_ext)
                PX[i][j] = None if X[i][j] == 0.0 else ('M' if X[i][j] == x_open else 'X')

                y_open = M[i][j - 1] + go + ge
                y_ext = Y[i][j - 1] + ge
                Y[i][j] = max(0.0, y_open, y_ext)
                PY[i][j] = None if Y[i][j] == 0.0 else ('M' if Y[i][j] == y_open else 'Y')

                m_from_M = M[i - 1][j - 1] + sub
                m_from_X = X[i - 1][j - 1] + sub
                m_from_Y = Y[i - 1][j - 1] + sub
                M[i][j] = max(0.0, m_from_M, m_from_X, m_from_Y)
                PM[i][j] = None if M[i][j] == 0.0 else (
                    'M' if M[i][j] == m_from_M else ('X' if M[i][j] == m_from_X else 'Y'))

                if M[i][j] > best_score:
                    best_score, best_matrix, best_pos = M[i][j], 'M', (i, j)
                if X[i][j] > best_score:
                    best_score, best_matrix, best_pos = X[i][j], 'X', (i, j)
                if Y[i][j] > best_score:
                    best_score, best_matrix, best_pos = Y[i][j], 'Y', (i, j)

        ti, tj = best_pos
        matrix = best_matrix
        a1, a2 = [], []

        while ti > 0 and tj > 0:
            if matrix == 'M':
                if M[ti][tj] == 0.0 or PM[ti][tj] is None:
                    break
                a1.append(s1[ti - 1]);
                a2.append(s2[tj - 1])
                matrix = PM[ti][tj]
                ti -= 1;
                tj -= 1
            elif matrix == 'X':
                if X[ti][tj] == 0.0 or PX[ti][tj] is None:
                    break
                a1.append(s1[ti - 1]);
                a2.append('-')
                matrix = PX[ti][tj]
                ti -= 1
            else:  # 'Y'
                if Y[ti][tj] == 0.0 or PY[ti][tj] is None:
                    break
                a1.append('-');
                a2.append(s2[tj - 1])
                matrix = PY[ti][tj]
                tj -= 1

        a1 = ''.join(reversed(a1))
        a2 = ''.join(reversed(a2))
        return {
            "score": best_score,
            "align1": a1,
            "align2": a2,
            "start_in_seq1": ti,
            "start_in_seq2": tj,
            "end_in_seq1": best_pos[0],
            "end_in_seq2": best_pos[1],
        }


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            if line.startswith(">"): continue
            seq.append(line.strip().upper())
    return "".join(seq)


def load_plp_to_pairs(plp_path, n, threshold=0.35):
    """Read i j p lines and keep (i-1, j-1, p) with p >= threshold."""
    pairs = []
    all_pairs = []
    with open(plp_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("/ubox"):
                continue
            if s.startswith("%"):
                continue

            if not s or s[0] in "#%":  # comments
                continue
            i, j, p = s.split()[:3]
            i, j, p = int(i)-1, int(j)-1, float(p)
            all_pairs.append((i, j, p))
            if 0 <= i < j < n and p >= threshold:
                pairs.append((i, j, p))
    # highest probability first
    pairs.sort(key=lambda x: x[2], reverse=True)
    return pairs, all_pairs


def choose_non_crossing_pairs(candidates, n, allow_pseudoknots=False):
    """Greedy ThreshKnot-like selection. Returns list of (i,j,p)."""
    chosen = []
    used = [False]*n

    def crosses(a,b,c,d):
        return (a < c < b < d) or (c < a < d < b)

    for i, j, p in [(i,j,p) for (i,j,p) in candidates]:
        if used[i] or used[j]:
            continue
        if not allow_pseudoknots:
            # forbid crossing with any chosen pair
            conflict = False
            for x,y,_ in chosen:
                if crosses(i,j,x,y):
                    conflict = True
                    break
            if conflict:
                continue
        # accept this pair
        chosen.append((i,j,p))
        used[i] = used[j] = True
    # sort by i for nicer downstream handling
    chosen.sort(key=lambda t: t[0])
    return chosen

def dotbracket_from_pairs(pairs, n):
    s = ["." for _ in range(n)]
    for i,j,_ in pairs:
        s[i] = "("
        s[j] = ")"
    return "".join(s)


def build_graph(seq, pairs):
    """Build NetworkX graph with backbone and base-pair edges."""
    n = len(seq)
    G = nx.Graph()
    # nodes
    for i, base in enumerate(seq):
        G.add_node(i, base=base, pos=i+1)  # 1-based label for aesthetics
    # backbone edges
    for i in range(n-1):
        G.add_edge(i, i+1, kind="backbone", weight=1, prob=0)
    # base-pair edges (prob attribute)
    for i, j, p in pairs:
        G.add_edge(i, j, kind="pair", weight=4, prob=p)

    return G


def draw_arc_diagram(G, seq, pairs, figsize=(14, 4)):
    """
    Linear layout (x = position) with arcs for base pairs above the line.
    Backbone edges are straight segments; base pairs are curved.
    """
    n = len(seq)
    # positions along x-axis
    pos = {i: (i, 0) for i in range(n)}
    labels = {i: f"{seq[i]}{i+1}" for i in range(len(seq))}

    # separate edge lists
    backbone = [(u,v) for u,v,d in G.edges(data=True) if d.get("kind")=="backbone"]
    bpedges  = [(u,v,d) for u,v,d in G.edges(data=True) if d.get("kind")=="pair"]

    fig, ax = plt.subplots(figsize=figsize)

    # draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=100, ax=ax)
    # draw backbone as simple straight lines
    nx.draw_networkx_edges(G, pos, edgelist=backbone, width=1, ax=ax)

    # draw base pairs as arcs with curvature related to span
    for u, v, d in bpedges:
        span = v - u
        # map span to a reasonable curvature (rad): larger span -> higher arc
        rad = min(0.6, 0.06 + 0.0025*span)
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)], ax=ax, width=1.5,
            connectionstyle=f"arc3,rad={rad}",
        )

    # optional: labels (comment out if too busy)
    # nx.draw_networkx_labels(G, pos, {i: seq[i] for i in range(n)}, font_size=6, ax=ax)

    ax.set_xlim(-2, n+1)
    ax.set_ylim(-0.1, 1.5)   # room for arcs
    ax.set_xticks(range(0, n, max(1, n//20)))
    ax.set_yticks([])
    ax.set_xlabel("Position")
    ax.set_title("RNA secondary structure (NetworkX arc diagram)")
    plt.tight_layout()
    plt.show()


def find_helices_strict(pairs_ijp):
    P = {(i,j): p for (i,j,p) in pairs_ijp if i < j}
    seen = set()
    helices = []

    for (i,j) in sorted(P.keys()):
        if (i,j) in seen:
            continue
        # walk backwards to the helix start
        a,b = i,j
        while (a-1, b+1) in P:
            a -= 1
            b += 1
        # extend forward collecting the stack
        stack = []
        x,y = a,b
        while (x,y) in P:
            stack.append((x,y,P[(x,y)]))
            seen.add((x,y))
            x += 1
            y -= 1
        helices.append(stack)
    return helices


def filter_helices(helices, min_len=3, min_mean_p=0.35):
    """
    Keep only helices with at least min_len stacked pairs and mean(p) >= min_mean_p.
    Returns filtered helices and a flat list of all pairs kept.
    """
    kept = []
    kept_pairs = []
    for h in helices:
        if len(h) < min_len:
            continue
        mean_p = sum(p for _,_,p in h) / len(h)
        if mean_p >= min_mean_p:
            kept.append(h)
            kept_pairs.extend(h)
    return kept, kept_pairs


def pairs_to_matrix(pairs, n):
    P = np.zeros((n, n), dtype=float)
    for i, j, p in pairs:
        P[i, j] = P[j, i] = p
    return P


def draw_arc_diagram_colored(G, seq, P, pairs,
                             color_mode="unpaired",   # "unpaired" or "paired"
                             u_values=None,           # 1D array from read_lunp(...), if color_mode="unpaired"
                             u_len=1,
                             cmap_name="viridis",
                             figsize=(14,4)):
    """
    Draw an arc diagram; color nodes by either:
      - unpaired probability (u_len), or
      - paired probability Σ_j P_ij
    Also color arcs by their base-pair probability.
    """
    n = len(seq)
    pos = {i: (i, 0) for i in range(n)}
    labels = {i: f"{seq[i]}{i+1}" for i in range(n)}
    backbone = [(u, v) for u, v, d in G.edges(data=True) if d.get("kind") == "backbone"]
    bpedges = [(u, v, d) for u, v, d in G.edges(data=True) if d.get("kind") == "pair"]

    # node values
    if color_mode == "unpaired":
        if u_values is None or len(u_values) != n:
            raise ValueError("Provide u_values = read_lunp(lunp_path, u=...) matching sequence length.")
        node_vals = u_values
        cbar_label = f"Unpaired probability (u={u_len})"
    elif color_mode == "paired":
        node_vals = P.sum(axis=1)          # Σ_j P_ij
        node_vals = np.clip(node_vals, 0, 1)
        cbar_label = "Paired probability (Σ Pij)"
    else:
        raise ValueError("color_mode must be 'unpaired' or 'paired'.")

    norm = colors.Normalize(vmin=0.0, vmax=1.0)
    cmap = colormaps[cmap_name]
    # node_colors = [cmap(norm(v)) for v in node_vals]

    fig, ax = plt.subplots(figsize=figsize)

    # nodes and backbone
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, ax=ax, linewidths=0.5)
    nx.draw_networkx_edges(G, pos, edgelist=backbone, width=1, ax=ax, arrows=False)

    # arcs for base pairs (width + color by pair prob)
    for u, v, d in bpedges:
        span = v - u
        rad = min(1, 0.06 + 0.0025 * span)
        prob = d["prob"]
        width = 0.5 + 1.0 * prob
        edge_color_ = [cmap(norm(prob))]
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)], ax=ax, width=width,
            edge_color=edge_color_,
            connectionstyle=f"arc3,rad={-rad *4}", arrows=True
        )

    # labels (optional)
    nx.draw_networkx_labels(G, pos, labels=dict(zip(G.nodes(), new_labels)), font_size=6, ax=ax)

    ax.set_xlim(-2, n + 1)
    ax.set_ylim(-0.1, 1.6)
    ax.set_yticks([])
    ax.set_xlabel("Position")
    ax.set_title(f"RNA secondary structure\nNode color = {cbar_label}")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(cbar_label)
    plt.tight_layout()
    plt.ylim(-0.05, 0.4)
    ax.set_frame_on(False)
    plt.show()


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
    temp = "".join([rna_structure[0].get_elem(elem) for elem in range(1, len(dot_bracket_seq) + 1)])
    exclude = [str(x) for x in range(0, 10)]
    return "".join([i for i in temp if i not in exclude])


def rc(seq_temp):
    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G",
                  "-": "-"}
    seq_out = "".join([complement[x] for x in seq_temp][::-1])
    return seq_out


def get_complement(seq_temp):
    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G",
                  "-": "-"}
    seq_out = "".join([complement[x] for x in seq_temp])
    return seq_out


def find_mirna_binding():


    # reverse complement for direct match
    seed_rc = rc(seed)
    mirna_seq_rc = rc(mirna_seq)
    print("mirna_seq: ", mirna_seq)
    print("mirna_seq_rc:",  mirna_seq_rc)
    print("seed: ", seed)
    print("seed_rc: ", seed_rc)
    print("seed_r: ", seed[::-1])
    # local alignment of miRNA seed to target
    aligner = Alignment(seq, mirna_seq, match=1, mismatch_penalty=-2, gap_penalty=-2,
                        )
    smith_waterman_result = aligner.smith_waterman()

    mirna_seq_c = mirna_seq_rc[::-1]
    start_mirna = mirna_seq_c.find(smith_waterman_result["align2"][::-1])
    end_mirna = start_mirna + len(smith_waterman_result["align2"])

    # mirna_pairs = list(zip(np.arange(smith_waterman_result["start_in_seq1"], smith_waterman_result["end_in_seq1"]),
    #                        ["m" + str(x) for x in np.arange(start_mirna, end_mirna)]))
    # print(mirna_pairs)
    print(smith_waterman_result)

    # find direct match for seed without wobble base
    # mirna_binding_start = seq.find(seed_rc)
    # mirna_binding_end = mirna_binding_start + len(seed_rc)
    mirna_binding_start = smith_waterman_result["start_in_seq1"]
    mirna_binding_end = smith_waterman_result["end_in_seq1"]

    print(mirna_binding_start, mirna_binding_end)
    print(seq[mirna_binding_start:mirna_binding_end])

    upstream_window = 20
    # check binding upstream match für residual miRNA binding in 35 nt window
    target_seq2 = seq[mirna_binding_start - upstream_window:mirna_binding_end]

    # get non_seed_mirna
    residual_mirna = mirna_seq[8:]
    print("residual_mirna: ", residual_mirna)
    print("target_seq2", target_seq2)

    aligner2 = Alignment(target_seq2, mirna_seq, match=1, wobble=1,
                         mismatch_penalty=-1.5, gap_penalty=-2,
                         gap_open=0, gap_extend=-0.5)
    aligner2.allow_wobble = True
    aligner2.allow_gap = True

    needleman_wunsch_result = aligner2.needleman_wunsch()

    print(needleman_wunsch_result["align1"])
    print(get_complement(needleman_wunsch_result["align2"]))
    # print(rc(needleman_wunsch_result["align2"]))

    target_idx = seq.find(target_seq2) - 1
    mirna_idx = len(mirna_seq)
    mirna_pairs = []
    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G"}
    for idx, (elem1, elem2) in enumerate(zip(needleman_wunsch_result["align2"], needleman_wunsch_result["align1"])):
        if elem2 != "-":
            target_idx += 1
        if elem2 == "-":
            mirna_idx -= 1
        if elem1 != "-" and elem2 != "-":
            mirna_idx -= 1
            print(elem2, complement[elem1], is_wc(complement[elem1], elem2), is_gu(complement[elem1], elem2))
            if is_wc(complement[elem1], elem2) or is_gu(complement[elem1], elem2):
                mirna_pairs.append((target_idx, "m" + str(mirna_idx)))

    print(mirna_pairs)
    return mirna_pairs


# miRNA binding edges
def is_wc(a, b):
    return ((a == 'A' and b == 'U') or
            (a == 'U' and b == 'A') or
            (a == 'G' and b == 'C') or
            (a == 'C' and b == 'G'))


def is_gu(a, b):
    return ((a == 'G' and b == 'U') or (a == 'U' and b == 'G'))


if __name__ == '__main__':

    # run RNAplfold in WSL
    """
    Windows > command / cmd
    WSL 
    conda activate rnafold 
    cd rnafold
    
    > store fasta in .fa-file 
     cat > seq.fa << 'EOF'
    >myseq
    GUGCCGUAACAACUGUGGUCAAUCCGAAGUAUGAGGGAAAAUGAGUACUGCCCGUGCAAAUCCCACAACACUGAAUGCAAAGUAGCAAUUUCCAUAGUCACAGUUAGGUAGCUUUAGGGCAAUAUUGCCAUGGUUUUACUCAUGUGCAGGUUUUGAAAAUGUACAAUAUGUAUAAUUUUUAAAAUGUUUUAUUAUUUUGAAAAUAAUGUUGUAAUUCAUGCCAGGGACUGACAAAAGACUUGAGACAGGAUGGUUACUCUUGUCAGCUAAGGUCACAUUGUGCCUUUUUGACCUUUUCUUCCUGGACUAUUGAAAUCAAGCUUAUUGGAUUAAGUGAUAUUUCUAUAGCGAUUGAAAGGGCAAUAGUUAAAGUAAUGAGCAUGAUGAGAGUUUCUGUUAAUCAUGUAUUAAAACUGAUUUUUAGCUUUACAAAUAUGUCAGUUUGCAGUUAUGCAGAAUCCAAAGUAAAUGUCCUGCUAGCUAGUUAAGGAUUGUUUUAA
    EOF
    
    > run RNAplfold 
     RNAplfold -W 200 -L 150 -u 10 --id-prefix myseq < seq.fa
     
    > extract base pairs
    awk '/ubox/{print $1, $2, ($3)^2}' myseq_0001_dp.pvals > myseq_0001_basepairs.txt
    
    > run script
    """
    """
    or 
    > run first in wsl without conda
    sudo add-apt-repository multiverse   # enable the multiverse repo
    sudo apt update
    sudo apt install vienna-rna
    
     > run directly in python
    """

    color_map = {"A": "#109648",  # green
                 "C": "#255C99",  # blue
                 "G": "#F7B32B",  # yellow
                 "U": "#D62839"  # red
                 }

    with open("FASTA/FASTA_ITB1_204.fa") as fa:
        fasta = ""
        for line in fa.readlines():
            if not line.startswith(">"):
                fasta = line.strip("\t")
                break

    positions = [2859,
                 3092,
                 ]

    pos1 = 2859
    pos2 = 3092
    start, end = min(positions) - 70, max(positions) + 50
    print(seq := fasta[start:end])
    rel_pos = [x - start - 1for x in positions]

    w_window = 150  # default 200   −−winsize=size
    l_span = 100  # default 150  # −L, −−span=size
    u_unpaired_len = 31  # default 31

    # load miRNA
    mirnas = pd.read_table(
        r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\miRNA_Expression\TargetScan_miR_Family_Info\miR_Family_Info.txt")
    name = "hsa-miR-493-5p"
    seed, mirna_seq = mirnas.loc[mirnas["MiRBase ID"] == name, ["Seed+m8", "Mature sequence"]].to_numpy()[0]
    mirna_nodes = find_mirna_binding()

    # 1) Set your paths
    fasta_path = r"C:\Users\tobia\rnafold\seq.fa"           # the FASTA you ran RNAplfold on
    plp_path = r"C:\Users\tobia\rnafold\myseq_0001_basepairs.txt"   # produced by RNAplfold (id-prefix may vary)
    with open(fasta_path, "w") as fh:
        fh.write(">myseq\n")
        fh.write(seq)
    # 2) Load data
    subprocess.run(["wsl", "bash", "-c", "echo Writing Fasta"])
    out_path = "/mnt/c/Users/tobia/rnafold"
    cmd = (f"cd {out_path} && RNAplfold -W {w_window} -L {l_span} -u {u_unpaired_len} --id-prefix myseq < /mnt/c/Users/tobia/rnafold/seq.fa")

    subprocess.run(["wsl", "bash", "-c", "echo Running RNAplfold"])
    subprocess.run(["wsl", "bash", "-c", cmd])

    time.sleep(2)

    with open(r"C:\Users\tobia\rnafold/myseq_0001_dp.ps") as f:
        file = f.read()
        s = "%start of base pair probability data"
        ppl = file[file.find(s):file.find("showpage\nend")]
        ppl = ppl[len(s) + 1:]
        with open(plp_path, "w") as bp:
            bp.write(ppl)

    seq = read_fasta(fasta_path)
    n = len(seq)

    # Balanced default (most use cases): 0.35–0.40 → good sensitivity with tolerable noise.
    # High-confidence, publication figures: 0.45–0.55 → fewer, more stable stems.

    pairs_candidates, pairs_all = load_plp_to_pairs(plp_path, n, threshold=0.5)
    print(pairs_candidates)

    # 2) build and filter helices
    helices = find_helices_strict(pairs_candidates)
    helices_kept, pairs_kept = filter_helices(helices, min_len=3, min_mean_p=0.35)

    # 3) Choose a consistent set of pairs (no pseudoknots)
    pairs = choose_non_crossing_pairs(pairs_kept, n, allow_pseudoknots=False)
    print(pairs)

    # 4) Build dot-bracket (optional, for saving/reporting)
    ss = dotbracket_from_pairs(pairs, n)
    print(ss)

    rna_elements = define_rna_element(ss)
    for pos in rel_pos:
        print(rna_elements[pos - 2: pos + 3])

    G = build_graph(seq, pairs)

    # add miRNA nodes
    for i in range(len(mirna_seq)):
        G.add_node("m" + str(i), base=mirna_seq[i],
                   pos=i + 1000)
    # print(mirna_nodes)
    # for pair in mirna_nodes:
    #     mirna_node = pair[1]
    #     target_node = pair[0]
    #     G.add_node(mirna_node, base=mirna_seq[int(mirna_node.replace("m", ""))],
    #                pos=target_node + 1000)  # 1-based label for aesthetics

    # miRNA backbone edges
    for i in range(len(mirna_seq) - 1):
        G.add_edge("m" + str(i), "m" + str(i + 1), kind="mirna_backbone", weight=1, prob=0)

    for i in range(0, len(mirna_nodes)):
        n1, n2 = mirna_nodes[i][1], mirna_nodes[i][0]
        b1, b2 = mirna_seq[int(n1[1:])], seq[n2]
        print(b1, b2)
        if is_wc(b1, b2):
            G.add_edge(n1, n2, kind="mirna_target_pair", weight=2, prob=1)
            print("is_WC")
        elif is_gu(b1, b2):
            G.add_edge(n1, n2, kind="mirna_target_pair", weight=6, prob=1)
        else:
            G.add_edge(n1, n2, kind="mirna_target_pair", weight=8, prob=1)
        # base-pair edges (prob attribute)
    edges = list(G.edges())

    print(G.edges())
    print(G.nodes())

    n = len(list(G.nodes()))
    # add color
    cmap_name = "Reds"
    cmap = colormaps[cmap_name]
    norm = colors.Normalize(vmin=0, vmax=1)
    node_colors = ["gray"] * n

    edge_color = ["gray"] * len(edges)
    node_sizes = [5] * n
    new_labels = [""] * n
    nodes = list(G.nodes())
    # new_labels = list(G.nodes())
    base_pos = np.arange(start, end)
    base_pos_show = [""] * n


    def highlight_pos(pos_temp, only_dmr=True, label_pos=False):
        if only_dmr:
            for idx in range(pos_temp - 2, pos_temp + 3):
                node_colors[idx] = color_map[seq[idx]]
                node_sizes[idx] = 20
                new_labels[idx] = seq[idx]
        else:
            for (rna_idx, mirna_idx) in mirna_nodes:
                node_colors[rna_idx] = color_map[seq[rna_idx]]
                node_sizes[rna_idx] = 20
                new_labels[rna_idx] = seq[rna_idx]

            # label miRNA
            for idx in range(1, len(mirna_seq) + 1):
                node_colors[-idx] = color_map[mirna_seq[int(nodes[-idx][1:])]]
                new_labels[- idx] = mirna_seq[int(nodes[-idx][1:])]
                node_sizes[- idx] = 20
        # label DMR base
        node_sizes[pos_temp] = 80

        if label_pos:
            # label pos
            for idx, _ in enumerate(base_pos):
                if int(idx) % 10 == 0 and not base_pos[idx]  in  np.arange(pos_temp - 2, pos_temp + 3):
                    new_labels[idx] = str(base_pos[idx])

    def color_edges():
        edge_width = [1] * len(edges)
        for idx, edge in enumerate(edges):
            if "m" in str(edge[1]):
                edge_color[idx] = "lightgray"
                continue
            if edge[0] == edge[1] + 1:
                continue
            for n1, n2, p in pairs:
                if edge[0] == n1 and edge[1] == n2:
                    edge_color[idx] = cmap(norm(p))
                    edge_width[idx] = p * 2
                    break


    highlight_pos(rel_pos[0], only_dmr=False)
    highlight_pos(rel_pos[1], only_dmr=True)

    color_edges()

    base_pos = nx.kamada_kawai_layout(G)  # Layout for better visualization
    # pos_spring = nx.spring_layout(G, pos=base_pos, fixed=list(helix_ends), iterations=200, seed=1)
    # base_pos = nx.kamada_kawai_layout(G, pos=pos_spring)  # Layout for better visualization
    # 5) Build graph and draw

    plt.figure(figsize=(22, 18))
    nx.draw_networkx(G,
                     base_pos,
                     labels=dict(zip(G.nodes(), new_labels)),
                     with_labels=True,
                     node_color=node_colors,
                     edge_color=edge_color,
                     font_size=4,
                     font_color="white",
                     node_size=node_sizes,
                     font_weight='regular')

    plt.axis('off')
    plt.tight_layout()
    plt.show()

    quit()
    P = pairs_to_matrix(pairs_all, n)
    draw_arc_diagram_colored(G, seq, P, pairs,
                             color_mode="paired",
                             cmap_name=cmap_name)
