### Scalable Temporal Motif Densest Subnetwork Discovery 

This respository contains a C++ implementation of the ALDENTE algorithms (and the baselines used for comparison) for our [KDD2024 paper](https://arxiv.org/abs/2406.10608) (Sarpe, Vandin, and Gionis).

## For Linux

Here is an example of how to use the code:
```
# Install hdf5 and relative dependencies
sudo apt-get install libhdf5-serial-dev
# Build the code
make
cd examples/densemotifs

# Run the suite of algorithms
./densemotifs to execute the suite of algorithms that takes the following paratemeters:
-i:<input_network.txt> in the format "src dst timestamp" for each line of the network see testsgraph.txt as an example
-m:<temporal_motif.txt> in the format "src dst ordering" see 4cycle.txt as an example
-delta:<delta> value of the parameter delta
-type: 1 for 1/k-appx, 2 for Batch, 3 for ProbPeel, 4 for the exact-flow based algorithms, 5 for HybridPeel, other types are for different experiments (see densemotifs.cpp)
-epsilon: value of xi for type 2 or 3, 5
-samp: number of samples if type 3 or 5 is used
-ithyb: parameter J for hybridpeel
-lambda: value of the exponential decay, set to 2 if $\tau_c$ is used
```
## CONTAINER
We now provide also a practical container that can be build trough [Apptainer](https://apptainer.org/) to avoid issues when building the code.

To build the image just make sure that you have installed Apptainer on your system and that it is properly working, then simply run the following commands.

```
# Build the base image of the OS
apptainer build snap-based.sif snap-based.def

# Build the executable that contains all ALDENTE algorithms and baselines
apptainer build ALDENTE.sif ALDENTE.def

# Run app
apptainer run --app densemotifs ALDENTE.sif

# Ask help
apptainer run-help --app densemotifs ALDENTE.sif
```

If you are not practical with Apptainer, you can check [this very simple guide](https://www.dei.unipd.it/~ceccarello/posts/apptainer-devenv/) provided by Prof. Matteo Ceccarello.

## Notes

The code was tested on versions of Ubuntu with OS > 18.04

If you experience problems with the JSON library, you can try to install it system-wide by following the guide:https://github.com/nlohmann/json

No issues should occur when building the image trough Apptainer
