# CS466-IDPP

This code is for our final project of CS466 "Introduction to Bioinformatics" in UIUC. We focus on the Incomplete Directed Perfect Phylogeny Problem which handling incomplete phylogenetic data arises whenever some of the data is missing. The question is whether one can complete the missing states in a way that admits a perfect phylogeny. We implement [Pe'er's algorithm](https://doi.org/10.1007/3-540-45123-4_14) and a brute-force method based on the Gusfield algorithm with Python3.7.

## Usage
The input data is a tab-separated matrix stored in .txt files in \test_data folder.

To run the code, execute the following command:
```shell
python3 run.py test_data/test_matrix.txt outfilename --plot
```
It will output a .dot file with the name "outfilename".dot

## Visualize Tree
To draw the graph encoded in outfilename.dot, [Graphviz](https://www.graphviz.org/) needs to be installed.
Then conduct the following command to get the plot:
```shell
dot -Tpng outfilename.dot -o output_tree.png
```
## Experiment
Out experiments shown in our report are implemented in experiments.ipynb to verify correctness and compare the speed between 2 algorithms.
