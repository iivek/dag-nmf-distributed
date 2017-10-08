ynopsis

An implementation of the learning algorithm from https://github.com/iivek/dag-constrained-nmf in an older version of Graphlab, by team from Carnegie Mellon University.
Hyperparameter optimization has not been implemented.
For more information, please visit the aforementioned link.

## Usage

Currently, both the input and output data are in form of .mat files and MATLAB API is used for I/O. An example of how to format the data can be bound under /resources. How the data and parameters are fed to the algorithm and how the algorithm is run can be seen in main.cpp.

## Dependencies

No CMake is available for the time being.
The package depends on CMU Graphlab, available here https://github.com/iivek/graphlab-cmu-mirror.
MATLAB API is required to read from and write to .mat files.
Also, Atlas, Blas, Boost.

## To be done
CMake
Command line interface.
Provide other means of input/output (currenlty it goes through .mat files exclusively)

## References

If you have used this code, please cite the following repo:
I.Ivek, "DAG-constrained Variational Bayesian NMF." [Online]. Available: https://github.com/iivek/dag-constrained-nmf

## License

Published under GPL-3.0 License.
