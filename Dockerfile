FROM continuumio/miniconda3:latest

RUN mkdir /opt/plasmidID/
ADD bin /opt/plasmidID/bin
ADD config_files /opt/plasmidID/config_files
ADD databases /opt/plasmidID/databases
ADD documents /opt/plasmidID/documents
ADD img /opt/plasmidID/img
ADD test /opt/plasmidID/test
ADD plasmidID /opt/plasmidID/
ADD environment.yml /opt/plasmidID/
ADD CHANGELOG.md /opt/plasmidID/
ADD LICENSE /opt/plasmidID/

RUN cd /opt/plasmidID
RUN /opt/conda/bin/conda env create -f /opt/plasmidID/environment.yml && /opt/conda/bin/conda clean -a
RUN /opt/conda/bin/conda env export --name plasmidID > plasmidID.yml
ENV PATH /opt/conda/envs/plasmidID/bin:$PATH
ENV PATH /opt/plasmidID/bin:/opt/plasmidID:$PATH
