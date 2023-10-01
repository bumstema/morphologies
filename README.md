# morphologies
Generate condensed planar distributions of arbitrary shaped 2d particles using Monte Carlo steric interactions.



![resized_](https://github.com/bumstema/morphologies/assets/25807978/99e5a1d0-6f0c-4f11-958a-218825526d18)


Example simulation of two ellipse, with one (blue) having thin hair-like pretursions on its surface. Hexagonal periodic boundary conditions.

Written as part of PhD thesis: [Simulating Self-Assembly of Organic Molecules & Classifying Intermolecular Dispersion](https://macsphere.mcmaster.ca/handle/11375/22038)

Use >make while in "1_compile" to create the executable files for your computer.  They will be put in folder: "2_run".  Note: requires "Boost" c-library to be installed.

Set up the simulation details within the files inside folder: "2_run".  Link the polygons files you wish to simulate using their system paths in "3_data".  Submit batches of multiple simulations to run on machines using sqsub with batch_submit.sh. 

Simulation final configurations will be placed into the folder: "3_data/data".

Examples of analysing the planar distributions can be seen with tools found in "4_analyse".


