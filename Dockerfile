FROM conda/miniconda3:latest

COPY * /opt

RUN cd /opt
RUN /opt/conda/bin/conda env create -f environment.yml && conda clean -a
RUN /opt/conda/bin/conda env export --name plasmidID > plasmidID.yml
ENV PATH /opt/conda/envs/plasmidID/bin:$PATH
ENV PATH /opt/bin:/opt:$PATH
