19-Oct-2025

-----
System requirements:

python_v3.9.6
imblearn            0.12.4
joblib              1.3.1
numpy               1.26.3
pandas              2.2.0
sklearn             1.1.3
Tested on: macOS-13.2-arm64-arm-64bit
No non-standard hardware. 

-----
Installation guide: 

Installation not required. 

-----
Demo: 

Demo data is provided which represents a subset of 200 cells where the tumour-reactivity 
of their gamma delta TCRs is known. Expected output is 1) trained model and 2) per barcode prediction scores determined by the trained model. Expected runtime on this demo data is less than 2 minutes. 

meta_discovery_tested.tsv - seurat formatted metadata of these 200 demo cells
PreGame_testing_matrix.feather - expression values for all PreGame features for these 200 demo cells (see Extended Data Figure 6b)
PreGame_testing_labels.tsv - randomized tumor-reactivity labels for these 200 demo cells

-----
Instructions for use: 

Save both scripts and demo data, replace generic '/path/to/'s in each script with the path to where the scripts and demo data was saved. Run code presented in 1_Train_PreGame.ipynb to train and save a prediction model. Then run code presented in 2_Run_PreGame.ipynb to use the model previously generated to calculate prediction scores for these 200 demo cells. The code in both scripts is exactly that which we used to train PreGame, while the demo data contains only a 200 cell subset, meaning exact reproduction of PreGame and related results in the manuscript is not possible. The full training data and/or PreGame models can be provided upon reviewer request allowing for reproduction. 


