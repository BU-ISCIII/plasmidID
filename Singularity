Bootstrap: docker
From: centos:latest

%files
    ./scif_app_recipes/trimmomatic_v0.33_centos7.scif /opt
	./scif_app_recipes/spades_v3.12.0_centos7.scif /opt
	./scif_app_recipes/ncbiblast_v2.7.1_centos7.scif /opt
	./scif_app_recipes/bedtools_v2.27_centos7.scif /opt
	./scif_app_recipes/bowtie2_v2.3.4.2_centos7.scif /opt
	./scif_app_recipes/samtools_v1.9_centos7.scif /opt
	./scif_app_recipes/prokka_v1.13_centos7.scif /opt
	./scif_app_recipes/cdhit_v4.6.6_centos7.scif /opt
	./scif_app_recipes/circos_v0.69.6_centos7.scif /opt

%post
	echo "Install basic development tools"
	yum -y groupinstall "Development Tools"
	yum -y update && yum -y install wget curl

	echo "Install python2.7 setuptools and pip"
	yum -y install python-setuptools
	easy_install pip

	echo "Installing SCI-F"
    pip install scif

	echo "Installing trimmomatic app"
    scif install /opt/trimmomatic_v0.33_centos7.scif
    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It will be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    echo 'export PATH=${PATH}:/scif/apps/trimmomatic/bin' >> $SINGULARITY_ENVIRONMENT

	echo "Installing samtools app"
    scif install /opt/samtools_v1.9_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/samtools/bin' >> $SINGULARITY_ENVIRONMENT

	echo "Installing spades app"
    scif install /opt/spades_v3.12.0_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/spades/bin' >> $SINGULARITY_ENVIRONMENT

    echo "Installing NCBI-BLAST + app"
    scif install /opt/ncbiblast_v2.7.1_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/ncbiblast/bin' >> $SINGULARITY_ENVIRONMENT

    echo "Installing bedtools app"
    scif install /opt/bedtools_v2.27_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/bedtools/bin' >> $SINGULARITY_ENVIRONMENT

    echo "Installing bowtie2 app"
    scif install /opt/bowtie2_v2.3.4.2_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/bowtie2/bin' >> $SINGULARITY_ENVIRONMENT

    echo "Installing prokka app"
    scif install /opt/prokka_v1.13_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/prokka/bin' >> $SINGULARITY_ENVIRONMENT

    echo "Installing cdhit app"
    scif install /opt/cdhit_v4.6.6_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/cdhit/bin' >> $SINGULARITY_ENVIRONMENT

    echo "Installing circos app"
    scif install /opt/circos_v0.69.6_centos7.scif
    echo 'export PATH=${PATH}:/scif/apps/circos/bin' >> $SINGULARITY_ENVIRONMENT

%runscript
    exec scif "$@"
