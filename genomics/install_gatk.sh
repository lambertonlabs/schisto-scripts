#installing gatk on a CLIMB instance
#if you are running a GVL, htslib is installed already, so install if not

###########################################################################
#git-lfs
###########################################################################
#install add-apt repository
sudo apt-get install software-properties-common

#update git
sudo add-apt-repository ppa:git-core/ppa

#call apt-get update to get git-lfs
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash

#install git lfs
sudo apt-get install git-lfs

#final install and initialise
git lfs install
###########################################################################
#GATKtools - SNP calling
###########################################################################
#navigate to https://github.com/broadinstitute/gatk/#building for more info
#clone the up-to-date GATK repository
git clone https://github.com/broadinstitute/gatk.git

#go to GATK
cd gatk

#build .jar
./gradlew bundle

###########################################################################
#Picardtools - BAM file fannying about
###########################################################################
#go to https://github.com/broadinstitute/picard
#clone picardtools
git clone https://github.com/broadinstitute/picard.git
cd picard/

#build .jar
./gradlew shadowJar

###########################################################################
#Freebayes - SNP calling
###########################################################################
#go to https://github.com/ekg/freebayes
#clone freebayes

git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes

#make in local
make

#make in usr/bin/lib
sudo make install

###########################################################################
#vcflib - snp filtering and summary statistics
###########################################################################
#clone and build
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make

###########################################################################
#vcftools - snp filtering and summary statistics
###########################################################################
#go to https://github.com/ekg/freebayes
#clone and build
git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure
make
make install
