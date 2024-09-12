# **Distance-Based Neighborhood Exploration with N-Orbits**
![header]("https://github.com/barbara-xiong/N-Orbit/blob/main/images/NOrbitLogo-03.png")

## **Contents**

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Citation](#citation)

## **Overview**

![N-Orbit Schematic]("https://github.com/barbara-xiong/N-Orbit/blob/main/images/Fig1_NOrbitSchematic.png")

Tissue cellular neighborhoods (TCNs) are spatially contiguous regions of homogeneous and distinct cell type composition. Formation of TCNs may be indicative of various cell types coordinating to form functional niches. Studying spatial relationships on the level of TCNs may allow for identification of localized patterns that may otherwise be diluted and overlooked when analyzing in bulk.

Several neighborhood detection methods have been developed in recent years including our lab's own method, [CytoCommunity](https://github.com/huBioinfo/CytoCommunity). However, downstream methods to quantify TCN changes across conditions (e.g. time, clinical subtypes) beyond the level of cell type enrichment are limited. To achieve this goal, we introduce a distance-based approach, centered around a novel N-Orbit formalism of a neighborhood.

The N-Orbit structure, derived from that of the n-hop, encodes the level of proximity of each possible cell type to each cell in a TCN. The focus on proximity, rather than exact subgraph connections, allows for increased computational efficiency and greater instances of each structure. N-Orbits can be formulated as vectors, allowing for calculation of pairwise distance between two N-Orbit structures using their Manhattan distance. Distances between neighborhoods can then be computed from the minimum total distance between their representative set of N-Orbits. Pairwise neighborhood distances can be compiled into an overall neighborhood distance matrix that can be used for clustering, projection, visualization, and other methods for analyzing neighborhoods on a global scale.

## **Installation**

### **Hardware requirement**

CPU: i7

Memory: 16G or more

Storage: 10GB or more

### **Software requirement**

Conda version: 4.13.0

Python version: 3.10.4

Clone this repository and cd into it as below

```bash
git clone https://github.com/barbara-xiong/N-Orbit.git
cd N-Orbit
```

**For Linux**

**Preparing the virtual environment**

```bash
conda env create -f environment_norbit.yml
conda activate CytoCommunity
```

## **Usage**

This package requires single-cell resolution spatial omics data, where cell types have already been annotated.
By default, the N-Orbit distance pipeline is intended for neighborhood-level use. However, small adjustments may be made to calculate distances on the sample level and are italicized below.

**Prepare input data**

To prepare the input data, compile a CSV table where each row represents a cell and columns are provided for **(1) x coordinates, (2) y coordinates, (3) cell type labels, (4) sample/image labels, and (5) neighborhood labels**. Column names are flexible and specified in Step1 below. Synthetic data examples are provided under the directory examples/synthetic_mrf_neighborhoods_v1.csv and examples/synthetic_mrf_neighborhoods_v2.csv.

*For calculating sample-level distances, use a column of a constant values (e.g. all zeros) for the neighborhood label column.*

### **N-Orbit Distance**

Scripts for calculating N-Orbit distance are in the NOrbitDistance folder.

**Run the following steps in the Windows Powershell or Linux Bash shell:**

**Step 1a. Enumerate all N-Orbits in each sample and bootstrap representative vectors.**

This step generates a folder "norbits" and counts the number of occurrences of all N-Orbits for each neighborhood in a given sample. In addition, this step generates a folder "neighborhood_vectors" and bootstraps the specified number of N-Orbits, storing them in vector form, for each neighborhood of each sample. This step needs to be run for each image in the dataset, with the image name passed as a command line argument.

```bash
conda activate NOrbit
cd NOrbitDistance
python Step1a_N-Orbit-Enumerate.py Image1     # replace Image1 with your image name
```

**Hyperparameters**

* input_file_path: The path to your input dataset CSV file.

* intermediate_path: The path to the folder where your intermediate files will be stored. This path should be different for each dataset / analysis performed to avoid file conflicts and ensure each is performed independently.

* im_label: The name of your image/sample label column

* x_label: The name of your x-coordinate column

* y_label: The name of your y-coordinate column

* neighborhood_label: The name of your neighborhood label column

* N: Your desired N-Orbit depth *N* (recommended 2)

* NUCLEUS_PENALTY: Your desired nucleus change penalty *p* (recommended 10)

* K: The number of neighbors in the k-NearestNeighbors. This should be 6 in most cases, but for extremely dense or 3D data this value may be higher.

* threshold: The maximum distance for two cells to be considered neighbors after k-NearestNeighbors graph construction. This should be the equivalent of 50 micrometers in whichever unit is being used for the x and y coordinates.

* sample_size: The number of boostrapped N-Orbits used to represent each neighborhood. (recommended at least 1000)

*For calculating sample-level distances, specify the neighborhood_label to the column of constant values, as mentioned in Preparing Inputs.*

**Step 1b: Compile N-Orbit enumerations into a single CSV**

This step generates a compiled CSV called "neighborhood_norbits.csv" of all N-Orbit enumeration CSVs in the "norbits" folder.

```bash
python Step1b_Compile-N-Orbits.py
```

**Hyperparameters**

* intermediate_path: The path to the folder where your intermediate files are stored, as in Step1a.

**Step 2a: Compute pairwise neighborhood distances**

This step creates a "neighborhood_dists" folder, and computes a distance matrix between representative N-Orbit vectors of each neighborhood pair before calculating the minimum total distance via cost matrix optimization. This step may be parallelized in chunks by image/sample or by neighborhood, or run serially on the entire dataset **("All Mode")**. When parallelizing in neighborhood chunks **("Neighborhood Mode")**, distances are calculated from the specified neighborhood to all other neighborhoods lexicographically before it. For image/sample chunks **("Image Mode")**, the same is done for all neighborhoods in the specific image/sample.

For **"Neighborhood Mode"**, this step needs to be run for every neighborhood of every image in the dataset. Both the image and neighborhood labels need to be passed as a single command line argument. For example, for Neighborhood1 of Image1, the command would be as follows. 

```bash
python Step2a_Neighborhood-Distances.py Image1_Neighborhood1
```

For **"Image Mode"**, this step needs to be run for every image in the dataset. The image is passed as a command line argument as follows.

```bash
python Step2a_Neighborhood-Distances.py Image1
```

For **"All Mode"**, this step only needs to be run once, without passing any command line arguments.

```bash
python Step2a_Neighborhood-Distances.py
```

**Hyperparameters**

* UNIT_MODE: The parallelization mode as just described. This can take on the values "Image", "Neighborhood", or "All".

* intermediate_path: The path to the folder where your intermediate files are stored, as earlier.

*For calculating sample-level distances, use Image Mode, substituting the constant value for the neighborhood label, e.g. Image1_0.*

**Step 2b: Compile individual distance calculation runs into a single distance matrix**

This step creates a CSV of the compiled neighborhood distance matrix from the individual CSVs in the neighborhood_dists folder. This step needs to be run regardless of UNIT_MODE used in Step2a.

```bash
python Step2b_Neighborhood-Distances.py
```

**Hyperparameters**

* input_file_path: The file path to the original input file, as in Step1a.

* intermediate_path: The folder path where intermediates are stored, as earlier.

* im_label: The name of the image label column, as in Step1a.

* neighborhood_label: The name of the neighborhood label column, as in Step1a.

### **N-Orbit Enrichment**

Given a neighborhood, or neighborhood cluster of interest, enriched N-Orbit may be identified using a bootstrap-permutation test. Scripts for this procedure are located in the NOrbitEnrichment folder.

**Run the following steps in the Windows Powershell or Linux Bash shell:**

**Setup**

The following should compile the Cython code and create permutation_test.c and permutation_test.html folders.

```bash
cd ../NOrbitEnrichment
python setup.py build_ext --inplace
```

**Preparing inputs**

To prepare inputs for the neighborhoods of interest, compile the corresponding CSVs in the neighborhood_vectors folder into a single CSV with an additional column annotating which image and neighborhood each vector initially came from. An example is provided in "examples/NR_cluster.csv".

**Calculating N-Orbit Enrichment**

Bootstrap-permutation tests may be run as multiple jobs to allow for parallelization. The number of tests should be kept consistent for all runs. Each run (or batch of tests) is designated through a command line argument as follows.

```bash
python n-orbit-enrichment.py 1  # example for Run 1
```

**Hyperparameters**

* vectors_path: The path to your input file storing the neighborhood vectors of interest.

* trials_path: The folder path where results will be stored for each set of trials for this neighborhood cluster.

* num_permutations: The number of bootstrap-permutation tests for each run.

* numCellTypes: The number of unique cell types in your dataset. It should be half the length of each N-Orbit vector.

**Compiling bootstrap permutation tests across multiple runs**

After completing all runs, the following command will compile records from trials_path into final enrichment results.

```bash
python compile.py
```

**Hyperparameters**

* input_file_path: The path to your original input files for calculating N-Orbit distance, as in Step1a.

* trials_path: The folder path where results are stored for each set of trials for this neighborhood cluster, as earlier.

* cell_type_label: The name of the column for cell type labels, as in Step1a.

* output_path: The file path where the final enrichment results will be stored.

## Maintainers

Barbara Xiong (barbara.xiong@pennmedicine.upenn.edu)

Yuxuan Hu (huyuxuan@xidian.edu.cn)

Kai Tan (tank1@chop.edu)

## Citation
TBD
