
Influenza stem binding
======================

This example will compute all descriptors for a set of 350 protein-protein
complexes using all four modeling protocols.  The computed descriptors are used
for some statistics regarding cross-correlation among descriptors, principal
component analysis, and correlations between modeling protocols.  Finally, a
subset of descriptors are used to train a simple Random Forest machine learning
model, and some statistics are computed.

Scripts in this example
-----------------------

Computing the descriptors:
 - `config-ppdx.ini` Configuration file for ppdx
 - `compute_descriptors.py` Main script to compute all descriptors
 - `clean.py` Clean unnecessary files from the `models` directory
 - `compute-averages.py` Compute average and median for each descriptor

Making analysis on the descriptors:
 - `convergence_nmodels.py` Study the convergence of the descriptors as a function of how many models are generated
 - `crosscorr.py` Cross-correlation matrix between all descriptors
 - `protocol_correlation.py` Correlation between modeling protocols
 - `timing.py` Timing statistics on making models and computing descriptors

Fitting a Random Forest model:
 - `rf-fit.py` Generate the model
 - `rf-stats.py` Compute some basic statistics

