# Trimmomatic
- wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
- unzip Trimmomatic-0.38.zip
- copy to /opt/Trimmomatic or use trimmomatic-dir PATH/TO/Trimmomatic-0.38

# SPAdes

- wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz
- tar -xzf SPAdes-3.12.0-Linux.tar.gz
- Add to PATH SPAdes-3.12.0-Linux/bin/

# Blast+

- sudo apt-get install ncbi-blast+

# Bowtie2

- sudo apt install bowtie2

# Cd-hit-est

- sudo apt-get install cd-hit

# Bedtools

- sudo apt install bedtools

# Prokka

- sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
- sudo cpan Bio::Perl
- git clone https://github.com/tseemann/prokka.git $HOME/prokka
- $HOME/prokka/bin/prokka --setupdb
- Add $HOME/prokka/bin/ to PATH

# Circos


- wget http://www.circos.ca/distribution/circos-0.69-6.tgz
- tar xvfz circos-0.69-6.tgz
- sudo apt-get -y install libgd2-xpm-dev
- Add circos-0.69-6.tgz/bin to PATH
- sudo sed -i 's/max_points_per_track = 25000/max_points_per_track = 20000000/g' /opt/circos-0.69-6/etc/housekeeping.conf






##g++
- sudo apt-get install build-essential
##libz.h
- sudo apt-get install libz-dev
##circos dependencies
- sudo apt install circos
