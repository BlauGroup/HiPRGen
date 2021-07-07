# HiPRGen

HiPRGen is a python module for constructing reaction networks by running hundreds of billions of reactions through a decision tree. HiPRGen is built intop off [pymatgen](https://pymatgen.org/) and [MPI4py](https://pymatgen.org/) which facilitates multi-node parallelism.

- species filtering: This phase loads a json generated from the database, generates a bunch of mol_entries, filters them by isomorphism and then runs each molecule through a hand crafted decision tree currently called `standard_mol_decision_tree` in `species_questions.py`. The resulting list is then pickled for loading in other phases. The reason we use pickle here instead of json is that some of the species questions append non trivial data structures to the mol entries which get mangled when serialized to json, but i consider this a bug and would definately prefer to fix this and use json

- bucketing: Now we loop through pairs (A,B) where A and B are mols in the saved pickle and group them by atom counts. These groups are stored in a bucket database.

- reaction filtering + network generation: This is where MPI comes in. The program launches a dispatcher process and tons of filter processes. The filter processes request buckets from the dispatcher, generate all possible reactions from that bucket, run them through the decision tree from `reaction_questions.py` and then send the ones which pass back as they are generated. The dispatcher writes the reactions sent back from the filter processes into the reaction network database.

- simulation: once the reaction network database has been generated, it is fed into RNMC which runs simulations and writes them into the reaction network database. This is much more Lustre friendly than the previous approach which was writing each trajectory to an independent file. It also greatly cleaned up the RNMC code so now there are no superfluous code paths

- analysis: HiPRGen also has some primitives for analysis which are useful. There is a much improved report generator in `report_generator.py` and a network loader in `network_loader.py` which completely abstracts away the fact that the network and trajectories are stored in a sqlite db
