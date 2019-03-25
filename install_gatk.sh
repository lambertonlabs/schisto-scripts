#installing gatk on a CLIMB instance
#if you are running a GVL, htslib is installed already, so install if not

#install git-lfs

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

#navigate to https://github.com/broadinstitute/gatk/#building for more info
#clone the up-to-date GATK repository
git clone https://github.com/broadinstitute/gatk.git

#go to GATK
cd gatk

#build .jar
./gradlew bundle
