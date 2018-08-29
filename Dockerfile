FROM centos:latest

COPY ./scif_app_recipes/*  /opt/

RUN echo "Install basic development tools" && \
    yum -y groupinstall "Development Tools" && \
    yum -y update && yum -y install wget curl && \
    echo "Install python2.7 setuptools and pip" && \
    yum -y install python-setuptools && \
    easy_install pip && \
    echo "Installing SCI-F" && \
    pip install scif ipython && \
    echo "Installing plasmidID app" && \
    #scif install /opt/samtools_v1.9_centos7.scif
    scif install /opt/plasmidid_vdev_centos7.scif
#ENV find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "PATH=\${PATH}:$in";done | tr '\n' ' '


ENTRYPOINT ["/opt/docker-entrypoint.sh"]
CMD ["plasmidID.sh"]
