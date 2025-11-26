This file shows the code to produce figure 4, figure 5, figure6, and figure 7 in the main part of the paper.

The scripts should run in the following order:

1. Algorithm_comparisons1.R; Comparisons for the Louvain clustering and GMM, export V for VB-GMM in Python script.
2. KSneuro.ipynb; VB-GMM clustering
3. Algorithm_comparisons1.R; import the posterior obtained from step 2 to visualise VB-GMM clustering results and make a comparison with the Louvain clustering and GMM.