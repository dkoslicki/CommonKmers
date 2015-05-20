#This docker file will create a docker container for the CommonKmers method

FROM cami/profiling

MAINTAINER David Koslicki version: 0.1

#Install Jellyfish
RUN apt-get -y update
RUN apt-get install -y wget
RUN mkdir jellyfish
RUN cd jellyfish && wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.0/jellyfish-2.2.0.tar.gz && tar -xzf jellyfish-2.2.0.tar.gz
RUN apt-get install -y g++
RUN apt-get install -y build-essential
RUN cd /jellyfish/jellyfish-2.2.0/ && ./configure && make

#Install Julia
RUN apt-get install -y software-properties-common python-software-properties
RUN add-apt-repository -y ppa:staticfloat/juliareleases
RUN add-apt-repository -y ppa:staticfloat/julia-deps
#force zero exit code due to outdated package called wheezy
RUN apt-get -y update || echo "ok" 
RUN apt-get install -y julia
RUN apt-get install -y hdf5-tools

#Get the CommonKmers code
RUN apt-get install -y git
RUN git clone https://github.com/dkoslicki/CommonKmers.git

#Install Julia packages
RUN julia -e 'Pkg.add("HDF5"); Pkg.add("ArgParse");'

#Install query per sequence
RUN cd /CommonKmers/src/QueryPerSeq/ && g++ -I /jellyfish/jellyfish-2.2.0/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.0/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.0/.libs query_per_sequence.cc sequence_mers.hpp -o query_per_sequence
RUN cp /CommonKmers/src/QueryPerSeq/query_per_sequence /usr/local/bin

#Install gzip
RUN apt-get install -y gzip

#Create tasks
RUN cp /CommonKmers/Docker/default /dckr/etc/tasks.d/default
RUN cp /CommonKmers/Docker/sensitive /dckr/etc/tasks.d/sensitive
RUN cp /CommonKmers/Docker/specific /dckr/etc/tasks.d/specific

#Install CommonKmers database
RUN echo "Be sure to download the training data at http://www.math.oregonstate.edu/~koslickd/CommonKmersData.tar.gz and uncompress it. The resulting directory should be passed to docker at the location CONT_DATABASES_DIR via: -v /path/to/uncompressed/Data:/dckr/mnt/camiref/CommonKmersData:ro"