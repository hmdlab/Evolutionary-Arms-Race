# Evolutionary-Arms-Race
The repository contains scripts and processed data for "Landscape of Evolutionary arms race between Transposable elements and KRAB-ZFP family".

## 1. Dependencies
We have tested our script in the following environments.
```
OS: macOS 12.5.1
CPU: Intel Xeon W
conda: 23.1.0
```
You can reconstruct the virtual environment following these command lines:
```
git clone URL
cd Evolutionary-Arms-Race
conda env create -f environment.yml
conda activate arms_race
```

## 2. Preprocessing
Download raw data used for analysis and preprocess them.

The following processes are contained. 

1. Download raw and processed data from Google Drive.
2. Download genome, annotation, and epigenetic state in hESC from public databases.
3. identification of KRAB-ZFP targets.
4. identification of escape candidates.
5. identification of binding sites of KRAB-ZFPs.
6. obtaining the distance between LTR7_HERVH and TSS.

The preprocessing takes approximately 30 minutes.

```
python script/download_rawdata.py
python script/preprocessing.py
```

## 3. Analysis
Figures and Supplementary Figures described in the paper can be reproduced by running the code in the Jupyter notebooks in the notebook/ directory.

## 4. Citation
**Landscape of Evolutionary arms race between Transposable elements and KRAB-ZFP family.**

Masato Kosuge, Jumpei Ito, Michiaki Hamada. 2024
(URL)
