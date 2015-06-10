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


## Running the program ##
To classify a sample using the Julia version, use the ``Classify.jl`` command located in ``CommonKmers/src/Julia``. An example of running the program in the sensitive mode using 48 threads is given by

``julia -p 48 Classify.jl -d ../Data/ -o /path/to/output/dir/ -i /path/to/input/file.fa -k sensitive -j /path/to/./jellyfish``

Both FASTA and FASTQ files are acceptable input.


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
