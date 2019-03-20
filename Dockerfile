FROM centos:latest

COPY ./scif_app_recipes  /opt/scif_app_recipes

RUN echo "Install basic development tools" && \
    yum -y groupinstall "Development Tools" && \
    yum -y update && yum -y install wget curl && \
    echo "Install python2.7 setuptools and pip" && \
    yum -y install python-setuptools && \
    easy_install pip && \
    echo "Installing SCI-F" && \
    pip install scif ipython && \
    echo "Installing plasmidID app" && \
    scif install /opt/scif_app_recipes/plasmidid_v1.4.2_centos7.scif

# Include ENV variables
ENV LC_ALL=en_US.UTF-8
ENV PATH=$PATH:/scif/apps/hmmer3/bin
ENV PATH=$PATH:/scif/apps/prodigal/bin
ENV PATH=$PATH:/scif/apps/circos/bin
ENV PATH=$PATH:/scif/apps/minced/bin
ENV PATH=$PATH:/scif/apps/aragorn/bin
ENV PATH=$PATH:/scif/apps/trimmomatic/bin
ENV PATH=$PATH:/scif/apps/prokka/bin
ENV PATH=$PATH:/scif/apps/barrnap/bin
ENV PATH=$PATH:/scif/apps/tbl2asn/bin
ENV PATH=$PATH:/scif/apps/cdhit/bin
ENV PATH=$PATH:/scif/apps/bedtools/bin
ENV PATH=$PATH:/scif/apps/spades/bin
ENV PATH=$PATH:/scif/apps/ncbiblast/bin
ENV PATH=$PATH:/scif/apps/bowtie2/bin
ENV PATH=$PATH:/scif/apps/plasmidid/bin
ENV PATH=$PATH:/scif/apps/samtools/bin
ENV PATH=$PATH:/scif/apps/gcc/bin
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/hmmer3/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/prodigal/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/circos/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/minced/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/aragorn/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/trimmomatic/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/prokka/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/barrnap/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/tbl2asn/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/cdhit/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/bedtools/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/spades/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/ncbiblast/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/bowtie2/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/plasmidid/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/samtools/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scif/apps/gcc/lib

# Include them also in /etc/bashrc
RUN echo "export LC_ALL=en_US.UTF-8" >> /etc/bashrc
RUN find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\$PATH:$in" >> /etc/bashrc;done
RUN if [ -z "${LD_LIBRARY_PATH-}" ]; then echo "export LD_LIBRARY_PATH=/usr/local/lib" >> /etc/bashrc;fi
RUN find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$in" >> /etc/bashrc;done
