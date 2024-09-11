import pandas as pd
import glob
from sklearn.neighbors import kneighbors_graph
import networkx as nx
from collections import defaultdict
import numpy as np
import plotly.express as px
from scipy.spatial import distance_matrix
import math
from matplotlib import pyplot as plt
import scipy.sparse as sp
import sys
import os

### INPUTS AND PARAMETERS
cells = pd.read_csv("/path/to/data/synthetic_mrf_neighborhoods_v1.csv")
intermediate_path = "/path/to/intermediates/SyntheticV1/"
IMAGE_NAME = sys.argv[1]
N = 2 # Recommended value N = 2
NUCLEUS_PENALTY = 5

# Column labels
im_label = "Image"
x_label = "x"
y_label = "y"
cell_type_label = "CellType"
neighborhood_label = "Neighborhood"


# Threshold for kNN graph filtering by spatial proximity (50um recommended)
threshold = 50

# Number of n-orbit samples to represent each neighborhood
sample_size = 10000


# List of all cell types (alphabetical)
uniqueCellTypes = list(sorted(set(cells[cell_type_label])))

# Create directories for intermediates if they don't already exist
if not os.path.exists(intermediate_path+"norbits/"):
    os.makedirs(intermediate_path+"norbits/")
if not os.path.exists(intermediate_path+"neighborhood_vectors/"):
    os.makedirs(intermediate_path+"neighborhood_vectors/")

def create_graph(image, threshold=50):
    # Create kNN graph from spatial coordinates
    A = kneighbors_graph(image[[x_label,y_label]].values, 6, mode='connectivity', include_self=False)
    print("KNN graph created.")
    D = distance_matrix(image[[x_label,y_label]].values, image[[x_label,y_label]].values)
    D[D <= threshold] = 1
    D[D > threshold] = 0
    D = sp.csr_matrix(D)
    A = A.multiply(D.astype(int))
    A.eliminate_zeros()
    G = nx.from_scipy_sparse_array(A)
    print("NX Graph created.")
    attr = {}
    for index, row in image.iterrows():
        attr[index] = row[cell_type_label]
    nx.set_node_attributes(G,attr,"Cell Type")
    attr2 = {}
    for index, row in image.iterrows():
        attr2[index] = row[neighborhood_label]
    nx.set_node_attributes(G,attr2,"Neighborhood")
    print("Graph attributes added.")
    return G

def filter_by_attribute(G, attr, value):
    return [x for x,y in G.nodes(data=True) if y[attr]==value]

def norbit(G,node, n, filter_by_neighborhood = False):
    # Construct the norbit from a specific node
    nodes = defaultdict(set)
    nodes[0] = set([node])
    all_nodes = set([node])
    for i in range(1,n+1,1):
        for j in nodes[i-1]:
            neighbors = set(G.neighbors(j))
            if filter_by_neighborhood:
                neighbors = neighbors.intersection(filter_by_attribute(G,"Neighborhood",G.nodes[node]["Neighborhood"]))
            nodes[i].update(neighbors.difference(all_nodes))
        all_nodes.update(nodes[i])
    return nodes

def norbit_to_attr(G, norbit, attr):
    # Convert an norbit to show its attribute rather than identity
    norbit_alt = defaultdict(list)
    for i in norbit.keys():
        for j in norbit[i]:
            norbit_alt[i].append(G.nodes[j][attr])
    return norbit_alt

def norbit_comp(G, norbit, uniqueCellTypes, mode = "raw"):
    # Convert an norbit to show its composition
    norbit1 = norbit_to_attr(G, norbit, "Cell Type")
    norbit_alt = defaultdict(list)
    for i in norbit1.keys():
        norbit_alt[i] = [0]*len(uniqueCellTypes)
        for j in range(len(uniqueCellTypes)):
            norbit_alt[i][j] = norbit1[i].count(uniqueCellTypes[j])
            if mode == "set" and norbit_alt[i][j] > 0:
                already = False
                for k in range(1,i):
                    if norbit_alt[k][j] > 0:
                        already = True
                if not already:
                    norbit_alt[i][j] = 1
                else:
                    norbit_alt[i][j] = 0
        if mode == "prop":
            norbit_alt[i] = np.asarray(norbit_alt[i])/sum(norbit_alt[i])
    return norbit_alt

def count_norbits(G, n = N):
    # Count the n-orbit abundance in a given graph
    counts = defaultdict(int)
    for celltype in uniqueCellTypes:
        nodes_of_type = filter_by_attribute(G, "Cell Type", celltype)
        for node in nodes_of_type:
            nodeorbit = norbit(G,node,2, filter_by_neighborhood = True)
            nodeorbit_comp = norbit_comp(G, nodeorbit, uniqueCellTypes, mode = "set")
            counts[norbit_encode(nodeorbit_comp, n)] += 1
    countsDF = pd.DataFrame({"Code": counts.keys(), "Count": counts.values()})
    countsDF = countsDF.sort_values(by="Count",ascending=False).reset_index(drop=True)
    return countsDF

def filter_by_attribute(G, attr, value):
    return [x for x,y in G.nodes(data=True) if y[attr]==value]

def norbit(G,node, n, filter_by_neighborhood = False):
    # Construct the norbit from a specific node
    nodes = defaultdict(set)
    nodes[0] = set([node])
    all_nodes = set([node])
    for i in range(1,n+1,1):
        for j in nodes[i-1]:
            neighbors = set(G.neighbors(j))
            if filter_by_neighborhood:
                neighbors = neighbors.intersection(filter_by_attribute(G,"Neighborhood",G.nodes[node]["Neighborhood"]))
            nodes[i].update(neighbors.difference(all_nodes))
        all_nodes.update(nodes[i])
    return nodes

def norbit_to_attr(G, norbit, attr):
    # Convert an norbit to show its attribute rather than identity
    norbit_alt = defaultdict(list)
    for i in norbit.keys():
        for j in norbit[i]:
            norbit_alt[i].append(G.nodes[j][attr])
    return norbit_alt

def norbit_comp(G, norbit, uniqueCellTypes, mode = "raw"):
    # Convert an norbit to show its composition
    norbit1 = norbit_to_attr(G, norbit, "Cell Type")
    norbit_alt = defaultdict(list)
    for i in norbit1.keys():
        norbit_alt[i] = [0]*len(uniqueCellTypes)
        for j in range(len(uniqueCellTypes)):
            norbit_alt[i][j] = norbit1[i].count(uniqueCellTypes[j])
            if mode == "set" and norbit_alt[i][j] > 0:
                already = False
                for k in range(1,i):
                    if norbit_alt[k][j] > 0:
                        already = True
                if not already:
                    norbit_alt[i][j] = 1
                else:
                    norbit_alt[i][j] = 0
        if mode == "prop":
            norbit_alt[i] = np.asarray(norbit_alt[i])/sum(norbit_alt[i])
    return norbit_alt

def norbit_encode(norbit_comp, n = N):
    # Encode norbit as string
    flattened = np.ndarray.flatten(np.asarray([np.array(x) for x in norbit_comp.values()]))
    string = ""
    for i in flattened:
        string += str(i)
    while len(string) < len(uniqueCellTypes) * (n+1):
        string += "0"
    return string

def norbit_decode(norbit_encoded):
    # Decode an norbit string to a dictionary
    decoded = {}
    for i in range(0,len(norbit_encoded),len(uniqueCellTypes)):
        j = i/len(uniqueCellTypes)
        string = norbit_encoded[i:i+len(uniqueCellTypes)]
        decoded[int(i/len(uniqueCellTypes))] = np.where(np.asarray(list(string)) == "1")
        decoded[int(i/len(uniqueCellTypes))] = [uniqueCellTypes[x] for x in decoded[int(i/len(uniqueCellTypes))][0]]
    return decoded

def orbit_code(code, n=N, inclusion_penalty = 0):
    # Orbit encoding of n-orbit
    code_l = np.asarray([*code]).astype(np.uint8)
    code_orbits = [0]*(len(uniqueCellTypes))
    digit = n+inclusion_penalty
    for r in range(len(uniqueCellTypes),(n+1)*len(uniqueCellTypes), len(uniqueCellTypes)):
        code_r = code_l[r:r+len(uniqueCellTypes)]
        code_orbits += digit*code_r
        digit -= 1
    return code_orbits

def explore_graph(image, name = "", n = N):
    # Get n-orbits for all neighborhoods in a graph
    G = create_graph(image,threshold=threshold)
    graphDFs = []
    for neighborhood in set(nx.get_node_attributes(G, "Neighborhood").values()):
        print(IMAGE_NAME, neighborhood)
        G_n = G.subgraph(filter_by_attribute(G, "Neighborhood", neighborhood))
        df = count_norbits(G_n, n)
        df["image"] = [name]*len(df)
        df["neighborhood"] = [str(neighborhood)]*len(df)
        graphDFs.append(df)
    return pd.concat(graphDFs)

# Enumerate all n-orbits
image_names = sorted(set(cells[im_label]))
allNorbits = []
for name in [IMAGE_NAME]:
    cells_im = cells[cells[im_label]==name].reset_index()
    print(name, len(cells_im))
    graphDF = explore_graph(cells_im, name=name)
    allNorbits.append(graphDF)
df = pd.concat(allNorbits)

df["Code Orbit"] = [orbit_code(code) for code in df["Code"]]
df["Code Orbit Sum"] = [sum(codeOrbit) for codeOrbit in list(df["Code Orbit"])]
df = df[df["Code Orbit Sum"] > 0]
df["Nucleus"] = [uniqueCellTypes[code[:len(uniqueCellTypes)].find("1")] for code in list(df["Code"])]
df['percentage'] = df.groupby(['image', 'neighborhood'])['Count'].transform(lambda z: z / z.sum() * 100)
df['image'] = df['image'].astype(str)
df['neighborhood'] = df['neighborhood'].astype(str)
df["image_neighborhood"] = df["image"] + "_" + df["neighborhood"]
df.to_csv(intermediate_path + "norbits/neighborhood_norbits_"+IMAGE_NAME + ".csv")

def norbit_code(code, n = N, nucleus_penalty = 10):
    # Full vector encoding of n-orbit
    code
    code_l = np.asarray([*code]).astype(np.uint8)
    code_norbit = [0]*(len(uniqueCellTypes))
    digit = n
    for r in range(len(uniqueCellTypes),(n+1)*len(uniqueCellTypes), len(uniqueCellTypes)):
        code_r = code_l[r:r+len(uniqueCellTypes)]
        code_norbit += digit*code_r
        digit -= 1
    code_norbit = np.append(code_l[:len(uniqueCellTypes)]*nucleus_penalty, code_norbit)
    return code_norbit

# Store sampled n-orbit vectors
neighborhood_list = sorted(set(df["image_neighborhood"]))

for i in neighborhood_list:
    print(i)
    image1 = "_".join(i.split("_")[:-1])
    n1 = i.split("_")[-1]
    df1 = df[(df["image"]==image1) & (df["neighborhood"] == n1)]
    df1 = df1.sample(sample_size,weights = "percentage",replace=True)
    df1["Code NOrbit"] = [norbit_code(code, nucleus_penalty=NUCLEUS_PENALTY) for code in df1["Code"]]
    df1_vectors = np.array([row["Code NOrbit"] for index, row in df1.iterrows()])
    df1_vectorsDF = pd.DataFrame(df1_vectors)
    df1_vectorsDF.to_csv(intermediate_path+"neighborhood_vectors/"+i+".csv")
