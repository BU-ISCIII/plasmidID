Bootstrap: docker
From: centos:latest

%files
    ./apps_recipes/bwa_v0.7.17_centos7.scif /opt
	./apps_recipes/samtools_v1.9_centos7.scif /opt
	./apps_recipes/bcftools_v1.9_centos7.scif /opt
%post
	echo "Install basic development tools"
	yum -y groupinstall "Development Tools"
	yum -y update && yum -y install wget curl

	echo "Install python2.7 setuptools and pip"
	yum -y install python-setuptools
	easy_install pip

	echo "Installing SCI-F"
    pip install scif

	echo "Installing bwa app"
    scif install /opt/bwa_v0.7.17_centos7.scif
    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It will be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    echo 'export PATH=${PATH}:/scif/apps/bwa/bin' >> $SINGULARITY_ENVIRONMENT

	echo "Installing samtools app"
    scif install /opt/samtools_v1.9_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/samtools/bin' >> $SINGULARITY_ENVIRONMENT

	echo "Installing bcftools app"
    scif install /opt/bcftools_v1.9_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/bcftools/bin' >> $SINGULARITY_ENVIRONMENT

%runscript
    exec scif "$@"
