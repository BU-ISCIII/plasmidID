Bootstrap: docker
From: centos:latest

%files
    ./scif_app_recipes/trimmomatic_v0.38_centos7.scif /opt
	./scif_app_recipes/spades_v3.12.0_centos7.scif /opt
	./scif_app_recipes/ncbiblast_v2.7.1_centos7.scif /opt
	./scif_app_recipes/bedtools_v2.27_centos7.scif /opt
	./scif_app_recipes/bowtie2_v2.3.4.2_centos7.scif /opt
	./scif_app_recipes/samtools_v1.9_centos7.scif /opt
	./scif_app_recipes/prokka_v1.13_centos7.scif /opt
	./scif_app_recipes/cdhit_v4.6.6_centos7.scif /opt
	./scif_app_recipes/circos_v0.69.6_centos7.scif /opt
	./scif_app_recipes/plasmidid_v1.3.0_centos7.scif /opt

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
    scif install /opt/plasmidid_v1.3.0_centos7.scif

    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It would be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\${PATH}:$in >> $SINGULARITY_ENVIRONMENT";done


%runscript
    exec scif "$@"
