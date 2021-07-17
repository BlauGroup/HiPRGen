FROM sleak75/conda-mpi4py-haswell:latest
SHELL ["/bin/bash", "-c"]
WORKDIR /app

RUN conda install -c conda-forge pymatgen openbabel pygraphviz

# do this to reduce image size:
RUN conda clean -a

RUN /sbin/ldconfig