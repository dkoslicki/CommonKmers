# CommonKmers #
This is a work in progress, and not yet ready for use.

This is the source code for CommonKmers.

## What is CommonKmers? ##
CommonKmers is a k-mer based bacterial community reconstruction technique that utilizes sparsity promoting ideas from the field of compressed sensing to reconstruct the composition of a bacterial community. This method allows for strain-level abundance estimation, and can quantify the evolutionary distance between organisms in the sample and in the training database (thereby allowing for successful classification even with incomplete training data).


## How Do I Install CommonKmers? ##
###Build from source###
You will need the Kmer counting tool Jellyfish to be installed. Please see [the Jellyfish installation page](http://www.genome.umd.edu/jellyfish.html) for installation directions.

You will need to download [this data repository](http://www.math.oregonstate.edu/~koslickd/CommonKmersData.tar.gz). Place this file in the CommonKmers directory, and then extract using ``tar -xf CommonKmersData.tar.gz``. This folder contains all the default training data.

Please refer to [the Julia installation page](http://julialang.org/downloads/) to install Julia.
You will need to add the HDF5 and ArgParse packages. These can be added using `Pkg.add("HDF5")` and `Pkg.add("ArgParse")`.

You will also need to compile the ``query_per_sequence`` code using a command such as:
```bash
g++ -I /jellyfish/jellyfish-2.2.0/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.0/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.0/.libs query_per_sequence.cc sequence_mers.hpp -o query_per_sequence
```

###Using Docker###
A Dockerfile is included in this repository. See the [Docker homepage](https://www.docker.com/) for more information.

You can build the docker image by cloning the repository, starting Docker, and then in the ``CommonKmers/Docker`` folder, using the command:
```bash 
docker build -t username/imagename .
```


## Running the program ##
###From the command line###
To classify a sample using the Julia version, use the ``ClassifyFull.jl`` command located in ``CommonKmers/src/Julia``. An example of running the program in the sensitive mode using 48 threads and a minimum quality score (for kmers to be counted) of C (phred33 ascii code 35) is given by

```julia
julia -p 48 ClassifyFull.jl -d /path/to/CommonKmersData/ -o /path/to/output/dir/ -i /path/to/input/file.fastq -Q C -k sensitive -j /path/to/./jellyfish
```

Only FASTQ files are acceptable input.

###Using Docker###
To run the tool from docker, mount the appropriate folders and run using the following command:
```bash
docker run --rm -e "QUALITY=C" -e "DCKR_THREADS=48" -v /path/to/CommonKmersData:/dckr/mnt/camiref/CommonKmersData:ro -v /path/to/Output:/dckr/mnt/output:rw -v /path/to/Input:/dckr/mnt/input:ro -t username/imagename [type]
```
In the input folder must be a collection of gzipped FASTQ (not FASTA) files, as well as a file (called ``sample.fq.gz.list`` (given by the docker image environmental variable ``$CONT_FASTQ_FILE_LISTING``) listing the files on which to run the tool.
Here ``[type]`` is one of ``default, sensitive, specific``.
The ``--rm`` flag deletes temporary files after exit (otherwise they might persist in ``/var/lib/docker/volumes`` or the like).
If the environmental variable ``QUALITY`` is not passed to docker (via ``-e QUALITY=<ascii character>``, a default value of "C" will be used. 



## Output format ##
The output format complies with the [CAMI format](https://github.com/CAMI-challenge/contest_information/blob/master/file_formats/CAMI_TP_specification.mkd).
The docker complies with the [Bioboxes profiling format 0.9](https://github.com/bioboxes/rfc/tree/master/data-format).

## Recommendations ##
I recommend using a quality score roughly equal to the average first quartile quality score in the file. This can be found with the following commands:
```bash
#Convert non ACTGN characters to N
awk '{if(NR%4==2){gsub(/[^ACGT]/,"N");print $0}else{print $0}}' input.fq > input_ACTGN.fq 
#Compute average first quartile quality score
/Fastx/bin/./fastx_quality_stats -i input_ACTGN.fq -Q33 | cut -f7 | sed -n '1!p' | awk '{a+=$1} END{print a/NR}' | awk '{printf "%.0f",$1}'
```
The FastX toolbox can be downloaded [here](http://hannonlab.cshl.edu/fastx_toolkit/).
## Contact ##
For issues with this software, contact david.koslicki@math.oregonstate.edu


## License ##
This project is released under the GPL-3 License. Please view the [LICENSE](LICENSE)
file for more details.


## Contributors ##
+ David Koslicki
+ Daniel Alonso Alemany
+ Daniel Falush
+ Nam Nguyen
