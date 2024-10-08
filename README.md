# Evolutionary-Arms-Race
The repository contains scripts and processed data for **"Landscape of Evolutionary arms race between Transposable elements and KRAB-ZFP family"**.

## 1. Dependencies
We have tested our script in the following environments.
```
OS: macOS 12.5.1
CPU: Intel Xeon W
conda: 23.1.0
```
You can reconstruct the virtual environment following these command lines:
```
git clone https://github.com/hmdlab/Evolutionary-Arms-Race.git
cd Evolutionary-Arms-Race
conda env create -f dependencies_all.yml
conda activate arms_race_all
```

if you can't reconstruct the virtual environment, please try to use dependencies_minimum.yml
```
conda env create -f dependencies_minimum.yml
conda activate arms_race_minimum
```

## 2. Preprocessing
Download raw data used for analysis and preprocess them.

The following processes are contained. 

1. Download raw and processed data from Box.
2. Download genome, annotation, and epigenetic state in hESC from public databases.
3. Identification of KRAB-ZFP targets.
4. Identification of escape candidates.
5. Identification of binding sites of KRAB-ZFPs.
6. Obtaining the distance between LTR7_HERVH and TSS.

The preprocessing takes approximately 30 minutes.

```
cd script
python download_rawdata.py
python preprocessing.py
```

if you want the additional data (Phylogenetic tree, distance matrix, MSA, consensus sequence), please run the following script.

Since the additional data is 30GB, please be cautious when downloading it.

```
python download_additional_data.py
```

## 3. Analysis
Figures and Supplementary Figures described in the paper can be reproduced by running the code in the Jupyter notebooks in the notebook/ directory.  

## 4. Citation
**Landscape of evolutionary arms races between transposable elements and KRAB-ZFP family**<br>
Masato Kosuge, Jumpei Ito, Michiaki Hamada. Scientific Reports, (2024)<br>
[https://doi.org/10.1038/s41598-024-73752-7](https://doi.org/10.1038/s41598-024-73752-7)

## 5. Changelogs
24/09/21 Add the information of the additional data and update the URL.<br>
24/10/08 update citation.
