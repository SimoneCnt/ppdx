
# HIV breadth example

This is an example on how to compute the breadth of the VRC01 and VRC01GL
antibodies using `ppdg`.  The breadth is computed as the fraction of antigens in
the Seaman panel for which the binding affinity is stronger than a given cutoff.
The details of the scoring function are stored in a `pkl` file.

To run this example script:
 1. Install ppdg 
 2. Edit the `ppdg-configure.ini` file with the correct paths
 3. Run the `hiv-breadth.py` file

`ppdg` keeps all intermediate and computed data stored on disk, so a second run is
much faster than the first since all structures and descriptors are already
computed.

The output plot should be similar to `histogram_breadth_ref.png`. The process is
stochastic, so slightly different results are expected.

The computed binding affinities and breadth values should be similar to:

    Binding affinity VRC01       :   -9.843
    Binding affinity VRC01GL     :   -9.062
    Computed breadth for VRC01   :    0.623
    Computed breadth for VRC01GL :    0.028

