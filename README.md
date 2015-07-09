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
julia -p 48 ClassifyFull.jl -d /path/to/CommonKmersData/ -o /path/to/output/file.profile -i /path/to/input/file.fastq -Q C -k sensitive -j /path/to/./jellyfish -q /path/to/./query_per_sequence
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

1. Create an acceptable taxonomy for the training databases.
2. Create 30mer and 50mer jellyfish files for each training genome.
3. Create common kmer matrices for the 30mers and 50mers.
4. Create bcalms for each 30mer jellyfish file.
5. Run CommonKmers using the custom training data.

All of these files must be placed in a directory (for example, called ``CommonKmerTrainingData`` below).

Before getting started, create a file consisting of the base names of each of the training genomes, and save this to a file (for example, ``fileNames.txt``).

####1. Creating custom taxonomy####
For each genome in ``fileNames.txt`` (and in the same order), a taxonomy file must be created. This file MUST be a newline delimitated file with each line having the following format:
```bash
<organismName>\t<TaxID>\t<TaxPath>
```

``<organismName>`` must be a unique identifier for each genome.

``<TaxID>`` must be a unique TaxID for each genome

``<TaxPath>`` must be a pipe delimitated list that gives the taxonomy of the given organism. The format is: 

```
k__<KingdomTaxID>_<KingdomName>|p__<PhylumTaxID>_<PhylumName>|c__<ClassTaxID>_<ClassName>|o__<OrderTaxID>_<OrderName>|f__<FamilyTaxID>_<FamilyName>|g__<GenusTaxID>_<GenusName>|s__<SpeciesTaxID>_<SpeciesName>|t__<StrainTaxID>_<StrainName>
```

An example line is as follows:

```
1184607_Austwickia_chelonae_NBRC_105200	1184607	k__2_Bacteria|p__201174_Actinobacteria|c__1760_Actinobacteria|o__2037_Actinomycetales|f__85018_Dermatophilaceae|g__1184606_Austwickia|s__100225_Austwickia_chelonae|t__1184607_Austwickia_chelonae_NBRC_105200
```

For your convenience, the script ``CommonKmers/src/Taxonomy/generate_taxonomy_taxid.py`` generates such a taxonomy using the NCBI taxonomy. This file must be placed in the ``CommonKmerTrainingData`` folder.

####2. Create 30mer and 50mer jellyfish files####
For each genome in ``fileNames.txt``, 30mer and 50mer jellyfish files must be created. And example command to do this is:

```bash
cat fileNames.txt | xargs -I{} -P <num_threads> /path/to/jellyfish count {} -m <kmer_size> -t 1 -s 100M -C -o /counts/{}-30mers.jf
```

The resulting jellyfish files MUST begin with the corresponding ``fileNames.txt`` name, and end in ``-xmers.jf`` with x=30 or x=50. For example, a file might be ``G000022605.fna-30mers.bcalm.fa``.

####3. Create common kmer matrices####
First, compile the ``/CommonKmers/src/CountInFile/count_in_file.cc`` code using a command like:

```bash
g++ -I /jellyfish/jellyfish-2.2.0/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.0/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.0/.libs count_in_file.cc -o count_in_file
```

Next, for the common kmer matrix using:

```bash
	julia -p <NumThreads> CommonKmers/src/Julia/FormTrainingMatrix.jl -i <JellyfishFiles> -o <OutFile> -c <count_in_file_binary> -s <chunk_size>
```

where:

``<NumThreads>`` is the number of threads to run.

``<JellyfishFiles>`` is a list of the full paths to the 30mers or 50mers jellyfish files (in the same order as ``fileNames.txt``).

``<Outfile>`` is the output 30mer common kmer matrix or 50mer common kmer matrix. You will need to place it in the ``CommonKmerTrainingData`` directory.

``<count_in_file_binary>`` is the location of the ``count_in_file_binary`` binary created previously.

``<chunk_size>`` is optional, but depends on the amount of RAM available. I have found that a value of approximately 700 is acceptable with 256GB of RAM. Specifying this is optional.

Note that you will need to do this for both the 30mers and the 50mers (so two ``<JellyfishFiles>`` will be needed, and two ``<Outfiles>`` will be created.
The resulting files must be placed in the ``CommonKmerTrainingData`` folder.

####4. Create bcalms for each 30mer jellyfish file####
For each 30mer jellyfish file, you will need to create a Bcalm file. The Bcalm source code can [be found here](https://github.com/Malfoy/bcalm).

The bcalm files can then be formed using:

```bash
julia -p <NumThreads> CommonKmers/src/Julia/FormBcalms.jl -i <FileNames> -c <LocationOf30merJellyfishFiles> -o <OutputFolder> -b <BcalmBinaryLocation> -j <JellyfishBinaryLocation> -r <RamdiskOrSSDLocation>
```

where:

``<NumThreads>`` is the number of threads to run.

``<FileNames>`` is the file ``fileNames.txt`` referred to above (base names of each of the training genomes).

``<LocationOf30merJellyfishFiles>`` is the location of the 30mer jellyfish files. Recall that these MUST be named like: ``<LocationOf30merJellyfishFiles>/<FileName>-30mers.jf``.

``<OutputFolder>`` MUST be the directory ``CommonKmersData/Bcalms/``.

``<BcalmBinaryLocation>`` is the location of the previously compiled code (eg. ``path/to/./count_in_file``).

``<JellyfishBinaryLocation>`` is the location of the jellyfish binary (eg. ``path/to/bin/jellyfish``).

``<RamdiskOrSSDLocation>`` is the location of a RAM disk or SSD (or other fast storage device. Unfortunately Bcalm uses a considerable amount of file IO, and so a fast storage device is required. Note that you can create a RAM disk using a command like:

```bash
mkdir /tmp/ramdisk; chmod 777 /tmp/ramdisk
mount -t tmpfs -o size=100G tmpfs /tmp/ramdisk/
```

####5. Run CommonKmers using the custom training data####
After the above steps are completed, you can utilize the custom training data by calling the ``CommonKmers/src/Julia/ClassifyFull.jl`` script. An example follows:

```bash
julia -p 48 ClassifyFull.jl -d <CustomDataPath> -o /path/to/output/file.profile -i /path/to/input/file.fastq -Q C -k sensitive -j /path/to/./jellyfish -q /path/to/./query_per_sequence --common_kmer_30_filename CommonKmers-30mers.h5 --common_kmer_50_filename CommonKmers-50mers.h5 --FileNames fileNames.txt --Taxonomy Taxonomy.txt
```

Where the files ``CommonKmers-30mers.h5``, ``CommonKmers-50mers.h5``, ``fileNames.txt``, and ``Taxonomy.txt`` (along with the folder ``Bcalms``) exist in the folder ``<CustomDataPath>`` (referred to above as the ``CommonKmerTrainingData``).

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
