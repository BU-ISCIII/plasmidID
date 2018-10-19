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
    scif install /opt/plasmidid_v1.4.1_centos7.scif
## Test

ENTRYPOINT ["/opt/docker-entrypoint.sh"]
CMD ["plasmidID.sh"]
