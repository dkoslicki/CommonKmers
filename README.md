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
####From the command line####
To classify a sample using the Julia version, use the ``ClassifyFull.jl`` command located in ``CommonKmers/src/Julia``. An example of running the program in the sensitive mode using 48 threads and a minimum quality score (for kmers to be counted) of C (phred33 ascii code 35) is given by

```julia
julia -p 48 Classify.jl -d /path/to/CommonKmersData/ -o /path/to/output/file.profile -i /path/to/input/file.fastq -Q C -k sensitive -j /path/to/./jellyfish -q /path/to/./query_per_sequence
```

Only FASTQ files are acceptable input.

####Using Docker####
To run the tool from docker, mount the appropriate folders and run using the following command:
```bash
docker run --rm -e "QUALITY=C" -e "DCKR_THREADS=48" -v /path/to/CommonKmersData:/dckr/mnt/camiref/CommonKmersData:ro -v /path/to/Output:/dckr/mnt/output:rw -v /path/to/Input:/dckr/mnt/input:ro -t username/imagename [type]
```
In the input folder must be a collection of gzipped FASTQ (not FASTA) files, as well as a file (called ``sample.fq.gz.list`` (given by the docker image environmental variable ``$CONT_FASTQ_FILE_LISTING``) listing the files on which to run the tool.
Here ``[type]`` is one of ``default, sensitive, specific``.
The ``--rm`` flag deletes temporary files after exit (otherwise they might persist in ``/var/lib/docker/volumes`` or the like).
If the environmental variable ``QUALITY`` is not passed to docker (via ``-e QUALITY=<ascii character>``), a default value of "C" will be used. 



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


## Custom Training Databases ##
If you wish to use a custom training database, the following steps must be performed:

1. Create a directory to contain the training data (called ``CommonKmerTrainingData`` below).
2. Create an acceptable taxonomy for the training genomes, and place it in the ``CommonKmerTrainingData`` folder.
3. Create a file consisting of the full paths of the training genomes, and save this to a file (for example, ``FileNames.txt``).
4. Compile the code contained in ``CommonKmers/src/CountInFile/``.
5. Run the script ``Train.jl``.

####Creating custom taxonomy####
For each genome in ``FileNames.txt`` (and in the same order), a taxonomy file must be created. This file MUST be a newline delimitated file with each line having the following format:
```bash
<organismName>\t<TaxID>\t<TaxPath>
```

``<organismName>`` must be a unique identifier for each genome.

``<TaxID>`` must be a unique TaxID for each genome

``<TaxPath>`` must be a pipe delimitated list that gives the taxonomy of the given organism. The format is: 

```
k__<KingdomTaxID>_<KingdomName>|p__<PhylumTaxID>_<PhylumName>|c__<ClassTaxID>_<ClassName>|o__<OrderTaxID>_<OrderName>|f__<FamilyTaxID>_<FamilyName>|g__<GenusTaxID>_<GenusName>|s__<SpeciesTaxID>_<SpeciesName>|t__<StrainTaxID>_<StrainName>
```

The taxonomy is only required at the kingdom level, with lower levels being optional.

An example line is as follows:

```
1184607_Austwickia_chelonae_NBRC_105200	1184607	k__2_Bacteria|p__201174_Actinobacteria|c__1760_Actinobacteria|o__2037_Actinomycetales|f__85018_Dermatophilaceae|g__1184606_Austwickia|s__100225_Austwickia_chelonae|t__1184607_Austwickia_chelonae_NBRC_105200
```

For your convenience, the script ``CommonKmers/src/Taxonomy/generate_taxonomy_taxid.py`` generates such a taxonomy using the NCBI taxonomy. This file must be placed in the ``CommonKmerTrainingData`` folder.

####Compile the ``count_in_file`` code####
The ``/CommonKmers/src/CountInFile/count_in_file.cc`` code can be compiled using a command like:

```bash
g++ -I /jellyfish/jellyfish-2.2.0/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.0/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.0/.libs count_in_file.cc -o count_in_file
```

####Run the script ``Train.jl``####
The script ``Train.jl`` can be called using a command such as:
```bash
julia -p 48 Train.jl -i FullFileNames.txt -o /path/to/output/CommonKmerTrainingData/ -b /path/to/./bcalm -r /path/to/fast/IO/device/ -j /path/to/jellyfish -c /path/to./count_in_file -s 500 -t 20
```

The option ``-s`` specifies how many training genomes at a time are held in memory. Increasing/decreasing this increases/decreases the amount of RAM used.

The option ``-t`` specifies how many jellyfish instances are created. Too many will cause disk thrashing (default is 20).

The option ``-r`` specifics the location of a temporary folder on a fast IO device. Unfortunately Bcalm uses a considerable amount of file IO, and so a fast storage device is required. Note that you can create a RAM disk using a command like:

```bash
mkdir /tmp/ramdisk; chmod 777 /tmp/ramdisk
mount -t tmpfs -o size=100G tmpfs /tmp/ramdisk/
```
####Run the ``Classify.jl`` script####
You can now run the ``Classify.jl`` script as before, but this time utilizing the directory ``CommonKmerTrainingData`` for the option ``-d``.
## Contact ##
For issues with this software, contact david.koslicki@math.oregonstate.edu

## License ##
This project is released under the GPL-3 License. Please view the [LICENSE](LICENSE)
file for more details.


## Contributors ##
+ David Koslicki (all main source code, unless otherwise noted)
+ Daniel Alonso Alemany (early version of the software)
+ Daniel Falush
+ Nam Nguyen
