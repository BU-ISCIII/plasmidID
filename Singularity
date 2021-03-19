Bootstrap: docker
From: conda/miniconda3:latest

%files
    environment.yml

%post
    /opt/conda/bin/conda env create -f environment.yml && conda clean -a
    /opt/conda/bin/conda env export --name plasmidID > plasmidID.yml

%runscript
    exec /opt/conda/envs/$(head -n 1 environment.yml | cut -f 2 -d ' ')/bin/"$@"
