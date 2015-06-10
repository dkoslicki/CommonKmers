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
```g++ -I /jellyfish/jellyfish-2.2.0/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.0/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.0/.libs query_per_sequence.cc sequence_mers.hpp -o query_per_sequence```

###Using Docker###
A Dockerfile is included in this repository. See the [Docker homepage](https://www.docker.com/) for more information.

You can build the docker image by cloning the repository, starting Docker, and then in the ``CommonKmers/Docker`` folder, using the command:
```docker build -t username/imagename .```


## Running the program ##
###From the command line###
To classify a sample using the Julia version, use the ``ClassifyFull.jl`` command located in ``CommonKmers/src/Julia``. An example of running the program in the sensitive mode using 48 threads is given by

```julia -p 48 ClassifyFull.jl -d /path/to/CommonKmersData/ -o /path/to/output/dir/ -i /path/to/input/file.fa -k sensitive -j /path/to/./jellyfish```

Both FASTA and FASTQ files are acceptable input.

###Using Docker###
To run the tool from docker, mount the appropriate folders and run using the following command:
```docker run -v /path/to/CommonKmersDataTest:/dckr/mnt/camiref/CommonKmersData:ro -v /path/to/Output:/dckr/mnt/output:rw -v /path/to/Input:/dckr/mnt/input:ro -t username/imagename [type]```
where ``[type]`` is one of ``default, sensitive, specific``.
In the input folder must be a collection of gzipped FASTQ (not FASTA) files, as well as a file (called ``sample.fq.gz.list`` (given by the docker image environmental variable ``$CONT_FASTQ_FILE_LISTING``) listing the files on which to run the tool.

## Output format ##
The output format complies with the [CAMI format](https://github.com/CAMI-challenge/contest_information/blob/master/file_formats/CAMI_TP_specification.mkd).
The docker complies with the [Bioboxes profiling format 0.9](https://github.com/bioboxes/rfc/tree/master/data-format).

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
