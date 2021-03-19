FROM continuumio/miniconda3:latest

RUN mkdir /opt/plasmidID
COPY * /opt/plasmidID/

RUN cd /opt/plasmidID
RUN /opt/conda/bin/conda env create -f environment.yml && /opt/conda/bin/conda clean -a
RUN /opt/conda/bin/conda env export --name plasmidID > plasmidID.yml
ENV PATH /opt/conda/envs/plasmidID/bin:$PATH
ENV PATH /opt/plasmidID/bin:/opt/plasmidID:$PATH
