FROM conda/miniconda3:latest

COPY environment.yml /

RUN /opt/conda/bin/conda env create -f environment.yml && conda clean -a
RUN /opt/conda/bin/conda env export --name plasmidID > plasmidID.yml
ENV PATH /opt/conda/envs/plasmidID/bin:$PATH
