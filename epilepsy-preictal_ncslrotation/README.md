Epilepsy: Analyzing Network Dynamics of the Brain and Epileptic Seizure Onset Zones
Goal: Determine existence of a preictal state using statistics and computational machine learning.

Main Contributor: Adam Li

Introduction:
Originally worked by Samuel P. Burns and Sabato Santaniello under Sri Sarma at Johns Hopkins Institute for Computational Medicine in collaboration with the Johns Hopkins Hospital. They published in PNAS, "Network dynamics of the brain and influence of the epileptic seizure onset zone", (Burns S, Santaniello S), PNAS, October 7, 2014. 

Their foundings focused in on the seizure zone of various clinical patients, such as: 
  1) complex partial seizure (CPS)
  2) secondarily generalized tonic-clonic seizure (GTC)
  3) partial seizures (PS)

Code for the epilepsy project. Goal is to determine existence of a pre-ictal state. Code to run:

1) Computation of interictal power spectrum using windows of Fourier Transform (FT) - IImspec2.m
2) Computation of seizure power spectrum and R-spectrum (used to determine optimal frequency band to analyze based on variability) - SZmspec2.m
3) Takes generated power spectrum and compute array of eigenvector centralities (EVCs) for time windows. - MakeQ2.m
4) Z-score the avarage connectivity and saves the standard deviation of adjacent matrix for patient - ZscoreADJ2.m
5) Compute clusters of EVCs using K-means and Gap Statistic - ClusterII2.m and GapStat.m

Future Work: 
  1) Develop IPython notebook for all sections of results/methods
  2) Use coherence matrices model to define a fragility metric of the network
  
  
  
