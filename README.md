<img src="./logo.png">

HiPRGen is a python module for constructing reaction networks by running hundreds of billions of reactions through a decision tree. HiPRGen is built intop off [pymatgen](https://pymatgen.org/) and [MPI4py](https://pymatgen.org/) which facilitates multi-node parallelism.

### Installation

The most consistent way to get HiPRGen running (especially on macos, where the conda version of MPI doesn't work consistently) is using the [nix package manager](https://nixos.org/). Instructions for installing nix can be found [here](https://nixos.org/download.html). HiPRGen requires nix version >= 2.4. If you have an older version installed, upgrade instructions can be found [here](https://nixos.org/manual/nix/unstable/installation/upgrading.html).

If you would prefer to use conda, the dependencies are `pymatgen`, `openbabel`, `pygraphviz`, `pycairo` and `mpi4py`. Create a conda environment where these are installed and then run `pip install -e .` from the HiPRGen directory. Again, in our experience, the conda version of MPI does not work well.

The whole process looks like this:
```
# this first step requires sudo to create the directory as root /nix with permissions 755.
sh <(curl -L https://nixos.org/nix/install) --daemon

## if you have an M1 mac, you also need to force nix to use x86 emulation mode
## some of our dependencies don't have native arm versions.
## uncomment the following lines if you have an M1 mac.
# mkdir -p ~/.config/nix
# echo "system = x86_64-darwin" > ~/.config/nix/nix.conf


# close your existing terminal and open a new one
git clone https://github.com/BlauGroup/HiPRGen
cd HiPRGen
nix-shell
```

### Tests

Once you are in an environment where HiPRGen is installed, the tests can be run with `python test.py 4`. This will run the tests using 4 threads (you should use as many as your machine has). The working directories for the tests are in `scratch`. After the tests have finished, run `python -i repl.py` and inspect the `network_loader` object. This is the test Lithium / Ethylene Carbonate network.

### Design

- species filtering: This phase loads a json generated from the database, generates a bunch of mol_entries, filters them by isomorphism and then runs each molecule through a hand crafted decision tree from `species_questions.py`. The resulting list is then pickled for loading in other phases. The reason we use pickle here instead of json is that some of the species questions append non trivial data structures to the mol entries which get mangled when serialized to json.

- bucketing: Now we loop through pairs (A,B) where A and B are mols in the saved pickle and group them by atom counts. These groups are stored in a bucket database.

- reaction filtering + network generation: This is where MPI comes in. The program launches a dispatcher process and tons of filter processes. The filter processes request buckets from the dispatcher, generate all possible reactions from that bucket, run them through a decision tree from `reaction_questions.py` and then send the ones which pass back as they are generated. The dispatcher writes the reactions sent back from the filter processes into the reaction network database.

- simulation: once the reaction network database has been generated, it is fed into [RNMC](https://github.com/BlauGroup/RNMC) which runs simulations and writes them into the reaction network database. This is much more Lustre friendly than the previous approach which was writing each trajectory to an independent file. It also greatly cleaned up the RNMC code so now there are no superfluous code paths

- analysis: HiPRGen also has some primitives for analysis which are useful. There is a much improved report generator in `report_generator.py` and a network loader in `network_loader.py` which completely abstracts away the fact that the network and trajectories are stored in a sqlite db.

The network loader is a great place to start using the codebase and is run as follows:

```
from HiPRGen.network_loader import *

    network_loader = NetworkLoader(
        './rn.sqlite',
        './mol_entries.pickle'
        )
```
