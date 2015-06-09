#This docker file will create a docker container for the CommonKmers method

FROM cami/profiling

MAINTAINER David Koslicki version: 0.1

#Install Jellyfish
RUN apt-get -y update && \ 
 apt-get install -y wget && \ 
 mkdir jellyfish && \ 
 cd jellyfish && \ 
 wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.0/jellyfish-2.2.0.tar.gz && tar -xzf jellyfish-2.2.0.tar.gz && \
 apt-get install -y g++ && \ 
 apt-get install -y build-essential  && \ 
 cd /jellyfish/jellyfish-2.2.0/ && \ 
 ./configure && \ 
 make

#Install Julia, force zero exit code for apt-get update due to outdated package called "wheezy" (who know what that is!)
RUN apt-get install -y software-properties-common python-software-properties  && \ 
 add-apt-repository -y ppa:staticfloat/juliareleases && \ 
 add-apt-repository -y ppa:staticfloat/julia-deps
RUN apt-get -y update || echo "ok" 
RUN apt-get install -y julia && \ 
 apt-get install -y hdf5-tools && \
 julia -e 'Pkg.add("HDF5"); Pkg.add("ArgParse");'

#Get the CommonKmers code, Install query per sequence
RUN apt-get install -y git && \
 git clone https://github.com/dkoslicki/CommonKmers.git && \
 cd /CommonKmers/src/QueryPerSeq/ && \
 g++ -I /jellyfish/jellyfish-2.2.0/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.0/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.0/.libs query_per_sequence.cc sequence_mers.hpp -o query_per_sequence && \
 cp /CommonKmers/src/QueryPerSeq/query_per_sequence /usr/local/bin

#Install gzip
RUN apt-get install -y gzip

#Install fastx toolkit
RUN mkdir Fastx && \
 cd Fastx && \
 wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
 tar -jxf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2

#Create tasks
RUN cp /CommonKmers/Docker/default /dckr/etc/tasks.d/default && \
 cp /CommonKmers/Docker/sensitive /dckr/etc/tasks.d/sensitive && \
 cp /CommonKmers/Docker/specific /dckr/etc/tasks.d/specific && \
 echo -e "\n \n NOTE: Be sure to download the training data at http://www.math.oregonstate.edu/~koslickd/CommonKmersData.tar.gz and uncompress it. The resulting directory should be passed to docker at the location CONT_DATABASES_DIR via: -v /path/to/uncompressed/Data:/dckr/mnt/camiref/CommonKmersData:ro"
