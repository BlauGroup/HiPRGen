<img src="./logo.png">

HiPRGen (**Hi**_gh_ **P**_erformance_ **R**_eaction_ **Gen**_eration_) is a python module for constructing reaction networks via exhaustive reaction enumeration and filtering decision trees with the capacity to be applied to systems with hundreds of billions of possible reactions. HiPRGen is built on top of [MPI4py](https://mpi4py.readthedocs.io/en/stable/) which facilitates multi-node parallelism.

### Installation

HiPRGen depends on `pymatgen`, `openbabel`, `pygraphviz`, `pycairo` and `mpi4py`. In our experience, the Conda version of MPI4py does not work consistently, so we use the [nix package manager](https://nixos.org/) to get HiPRGen running on a wide range of systems.  Instructions for installing nix can be found [here](https://nixos.org/download.html).

The whole process looks like this:
```
# The first step requires sudo to create the directory /nix as root.
# Run the NixOS install script below and follow the prompts.
# Note: On Linux, instructions for uninstalling nix can be found with a quick
# web search. On MacOS, uninstalling can be accompished with this script:
# https://gist.github.com/expelledboy/c00aebb004b178cf78b2c9b344526ff6

sh <(curl -L https://nixos.org/nix/install) --daemon

# If you have an M1 Mac, you also need to force nix to use x86 binaries
# since some of our dependencies don't have native arm binaries.
# Uncomment and run the following two lines if you have an M1 Mac:

# mkdir -p ~/.config/nix
# echo "system = x86_64-darwin" > ~/.config/nix/nix.conf


# Close your existing terminal and open a new one, then run:

git clone https://github.com/BlauGroup/HiPRGen
cd HiPRGen
nix-shell
```

HiPRGen is supported for MacOS and Linux and has been tested on MacOS 11.6 and 12.0.1 as well as Ubuntu 21.10. Installation should take less than five minutes.


### Running on the LRC cluster

On the LRC cluster, an environment where HiPRGen can be run is set up as follows:

```
module load python/3.9.12
conda init

logout, log back in

module load python/3.8.2-gcc
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

pip3 install --user mpi4py
conda create -n HiPRGen_RNMC python=3.8
conda activate HiPRGen_RNMC
conda install -c conda-forge openbabel pygraphviz pymatgen pycairo

git clone https://github.com/BlauGroup/RNMC.git
cd RNMC
module load gsl
CXX=g++ ./build.sh
export PATH=$PATH:$PROJ/RNMC/build


can pick up from reloading the environment:

conda activate HiPRGen_RNMC
module load python/3.8.2-gcc
module load gcc/7.4.0 
module load openmpi/4.0.1-gcc
module load gsl
export PATH=$PATH:$PROJ/RNMC/build
```

### Tests

Once you are in an environment where HiPRGen is installed, the tests can be run with `python test.py 4`. This will run the tests using 4 threads, though you could use as many threads as your machine allows to speed up the execution. Running the tests will populate working directories in `scratch`. Note that `test.py` is heavily commented to explain how to use HiPRGen. With at least 4 threads, the tests should take less than five minutes to run. Along with a variety of other information, the following lines will be printed to standard output to confirm that the tests have passed:

```
mg_test: correct number of species
mg_test: correct number of reactions
li_test: correct number of species
li_test: correct number of reactions
```

Once the tests have finished, you can run `python -i repl.py` and inspect the `network_loader` object, which contains all of the data associated with the test Lithium / Ethylene Carbonate network after running 1000 trajectories. Additionally, HiPRGen has a report generation system for visualizing results. For example, in `scratch/li_test`, run `pdflatex LEDC_pathways.tex` to generate a PDF of the top pathways to Lithium Ethylene Dicarbonate (LEDC) in the test Lithium / Ethylene Carbonate network. Explanation of other types of reports and the commands to generate them are given in `test.py`.


### Design

- Species filtering: This phase loads a JSON generated from our database, constructs molecule entries, filters them by isomorphism, and then runs each molecule through a handcrafted decision tree in `species_questions.py`. The resulting list is then pickled for loading in other phases. The reason we use pickle here instead of JSON is that some of the species questions append non-trivial data structures to the molecule entries which get mangled when serialized to JSON.

- Bucketing: Now we loop through pairs (A,B) where A and B are molecules in the saved pickle and group them by atom counts. These groups are stored in a bucket database.

- Reaction filtering + network generation: This is where MPI is used. The program launches a dispatcher process and many filter processes. The filter processes request buckets from the dispatcher, generate all possible reactions from each bucket, run those reactions through a decision tree from `reaction_questions.py`, and then sends the reactions which pass the decision tree back to the dispatcher as they are generated. The dispatcher writes the reactions sent back from the filter processes into the reaction network database.

- Simulation: Once the reaction network database has been generated, it is provided as an input to [RNMC](https://github.com/BlauGroup/RNMC) which runs simulations and writes them into the reaction network database. This is much more well-suited to Lustre filesystems than an approach involving writing each trajectory to an independent file.

- Analysis: HiPRGen also has important primitives for useful analysis. The ReportGenerator class in `report_generator.py` facilitates the construction of a variety of useful PDFs via functions in `mc_analysis.py`, and the NetworkLoader class in `network_loader.py` allows for straightforward interrogation of the network and trajectories while abstracting away the fact that they are stored in a sqlite db.

The network loader is a great place to start using the codebase and is run as follows:

```
# run from the root directory of HiPRGen after running the tests
from HiPRGen.network_loader import *

network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite',
)
```
