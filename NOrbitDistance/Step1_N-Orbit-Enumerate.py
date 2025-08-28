import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import networkx as nx
from scipy.spatial import KDTree
import numpy as np
from collections import defaultdict
import sys
import os

### INPUTS AND PARAMETERS
input_file_path = "/path/to/data/synthetic_mrf_neighborhoods_v1.csv"
intermediate_path = "/path/to/intermediates/SyntheticV1/"
MODE = "Neighborhood"

# Maximum radius for distance measurement between cells
radius = 100 # default 100
# Maximum radius for connected components in Instance Mode
instanceRadius = 50  # default 50
# Minimum neighborhood size (number of cells) to be recorded
minSize = 40  # default 40
# Nucleus penalty p (higher means more weight toward cell type composition)
nucleusPenalty = 1 # default 1
# Number of n-orbit samples to represent each neighborhood
sample_size = 1000 #default 1000

# Default 10 microns
distMin = 10

# Column labels
im_label = "Image"
x_label = "x"
y_label = "y_"
cell_type_label = "CellType"
neighborhood_label = "Neighborhood"
### ---

# List of all cell types (alphabetical)
cells = pd.read_csv(input_file_path)
uniqueCellTypes = list(sorted(set(cells[cell_type_label])))
numCellTypes = len(uniqueCellTypes)

# Create directories for intermediates if they don't already exist
if not os.path.exists(intermediate_path+"sampled_vectors/"):
	os.makedirs(intermediate_path+"sampled_vectors/")

cell_type_indices = {}
for i in range(len(uniqueCellTypes)):
    cell_type_indices[uniqueCellTypes[i]]=i
    
if MODE == "Image":
    cells["unit"] = cells[im_label]
else:
    cells["unit"] = cells[im_label].astype(str) + "_" + cells[neighborhood_label].astype(str)
cells["unit"] = [i.replace(" ","-") for i in cells["unit"]]


def generate_spatial_graph(df, radius=radius):
    '''Generate spatial graph for sample or neighborhood whose metadata is in df'''
    coords = df[[x_label,y_label]].to_numpy()
    info = df[[x_label,y_label,cell_type_label]].to_numpy()
    # Radius is maximum distance (in μm) between cells for neighbors
    tree = KDTree(coords)
    pairs = tree.query_pairs(radius)
    
    # Create graph
    G = nx.Graph()
    i=0
    for cell in info:
        x = cell[0]
        y = cell[1]
        ct = cell[2]
        G.add_node(i, pos=(x, y), CellType=ct)
        i+=1
    
    for i, j in pairs:
        dist = np.linalg.norm(np.array(G.nodes[i]["pos"])-np.array(G.nodes[j]["pos"]))
        G.add_edge(i, j, weight = (distMin/max(distMin,dist)))
    return G

def generate_instances(df,unit,radius=radius):
    coords = df[[x_label,y_label]].to_numpy()
    info = df[[x_label,y_label,cell_type_label]].to_numpy()

    # Radius is maximum distance (in μm) between cells for neighbors
    tree = KDTree(coords)
    pairs = tree.query_pairs(radius)
    
    # Create graph
    G = nx.Graph()
    i=0
    for cell in info:
        x = cell[0]
        y = cell[1]
        ct = cell[2]
        G.add_node(i, pos=(x, y), CellType=ct)
        i+=1
    
    for i, j in pairs:
        G.add_edge(i, j)

    ccs = nx.connected_components(G)
    df = {x_label:[],y_label:[],cell_type_label:[], "unit":[]}
    
    i=1
    for cc in ccs:
        if len(cc) > minSize:
            for node in cc:
                df[x_label].append(G.nodes[node]["pos"][0])
                df[y_label].append(G.nodes[node]["pos"][1])
                df[cell_type_label].append(G.nodes[node]["CellType"])
                df["unit"].append(unit+"-"+str(i))
            i+=1
    
    print(str(i),"instances created.")
    return pd.DataFrame(df)

def generate_vectors(G):
    ''' Generate list of N-Orbit vectors from graph '''
    vectors = []
    for cell in G.nodes:
#         if cell % 1000 == 0:
#             print(cell,"/",len(G.nodes))
        neighbors = G.neighbors(cell)
        vector = np.zeros((numCellTypes*2))

        nucCellType = G.nodes[cell]["CellType"]
        vector[cell_type_indices[nucCellType]] = nucleusPenalty

        for n in neighbors:
            nbrCellType = G.nodes[n]["CellType"]
            nbrIndex = numCellTypes + cell_type_indices[nbrCellType]
            vector[nbrIndex] = np.max([vector[nbrIndex], G.edges[cell,n]["weight"], G.edges[n,cell]["weight"]])
        vectors.append(vector)
    return np.array(vectors)

def process_unit(df, unit):
    G = generate_spatial_graph(df, radius=radius)
    vectors = pd.DataFrame(generate_vectors(G))
    if len(vectors) > minSize:
        sampled_vectors = vectors.sample(sample_size,replace=True)
        sampled_vectors.to_csv(intermediate_path+"sampled_vectors/"+unit+".csv")
        vectors["unit"] = unit
        return vectors

def main(cells, instances = False, instanceRadius = instanceRadius):
    vectors = []
    if instances:
        instanceDF = []
    for unit in sorted(set(cells["unit"])):
        print(unit)
        df = cells[cells["unit"]==unit]
        if instances:
            df = generate_instances(df,unit,radius=instanceRadius)
            for instance in sorted(set(df["unit"])):
                print(instance)
                vectors.append(process_unit(df[df["unit"]==instance],instance))
            instanceDF.append(df)
        else:
            unit_vectors = process_unit(df,unit)
            if unit_vectors is not None:
                vectors.append(unit_vectors)
    vectors = pd.concat(vectors)
    vectors.to_csv(intermediate_path+"norbits.csv")
    if instances:
        pd.concat(instanceDF).to_csv(intermediate_path+"instance_metadata.csv")

import time
start_time = time.time()
main(cells,instances=(MODE=="Instance"))
elapsed_time = time.time()-start_time
elapsed_hours = np.floor(elapsed_time/3600)
elapsed_minutes = np.floor((elapsed_time-elapsed_hours*3600)/60)
elapsed_seconds = elapsed_time-elapsed_hours*3600-elapsed_minutes*60

with open(intermediate_path+"Step1_Runtime.txt", "w") as file:
    file.write(str(elapsed_hours)+" hours, " + str(elapsed_minutes) + " minutes, " + str(elapsed_seconds) + " seconds.")
