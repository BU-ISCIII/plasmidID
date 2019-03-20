Bootstrap: docker
From: centos:latest

%files
	./scif_app_recipes/ /opt/
%post
	echo "Install basic development tools"
	yum -y groupinstall "Development Tools"
	yum -y update && yum -y install wget curl

	echo "Install python2.7 setuptools and pip"
	yum -y install python-setuptools
	easy_install pip

	echo "Installing SCI-F"
    pip install scif

    echo "Installing plasmidID app"
    scif install /opt/scif_app_recipes/plasmidid_v1.4.2_centos7.scif

    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It would be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    echo "PlasmidID Done"
    echo "export LC_ALL=en_US.UTF-8" >> $SINGULARITY_ENVIRONMENT
    find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\${PATH}:$in" >> $SINGULARITY_ENVIRONMENT;done
    find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:$in" >> $SINGULARITY_ENVIRONMENT;done

%runscript
    exec scif "$@"
