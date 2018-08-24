FROM centos

# Dependencies
RUN	yum -y groupinstall "Development Tools"
RUN yum -y update && yum -y install wget curl
RUN yum -y install python-setuptools
RUN easy_install pip

# Install scif from pypi
RUN pip install scif

# Install the filesystem from the recipe
ADD apps_recipes/*.scif /opt/
RUN scif install /opt/bwa_v0.7.17_centos7.scif
ENV PATH=${PATH}:/scif/apps/bwa/bin
RUN scif install /opt/samtools_v1.9_centos7.scif
ENV PATH=${PATH}:/scif/apps/samtools/bin
RUN scif install /opt/bcftools_v1.9_centos7.scif
ENV PATH=${PATH}:/scif/apps/bcftools/bin

# SciF Entrypoint
# Disabled because of compatibility with nextflow.
#ENTRYPOINT ["scif"]
CMD scif
