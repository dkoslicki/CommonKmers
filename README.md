# CommonKmers #
This is the source code for CommonKmers.

## What is CommonKmers? ##
CommonKmers is a k-mer based bacterial community reconstruction technique that utilizes sparsity promoting ideas from the field of compressed sensing to reconstruct the composition of a bacterial community. This method allows for strain-level abundance estimation, and can quantify the evolutionary distance between organisms in the sample and in the training database (thereby allowing for successful classification even with incomplete training data).


## How Do I Install CommonKmers? ##
There are two versions available, one written in Python and one in Julia. The Julia version is nearly two orders of magnitude faster than the python version.

Both versions require the Kmer counting tool Jellyfish to be installed. Please see [the Jellyfish installation page](http://www.genome.umd.edu/jellyfish.html) for installation directions.

#### Julia version ####
Please refer to [the Julia installation page](http://julialang.org/downloads/) to install Julia.
You will need to add the HDF5 and ArgParse packages. These can be added using `Pkg.add("HDF5")` and `Pkg.add("ArgParse")`.

##### Running the Julia version #####
To classify a sample using the Julia version, use the commands

##### Retraining the Julia version #####
To train on a custom database, use the functions ` FormCountsAndKCounts.jl` and `FormCommonKmerMatrix.jl`. For example, say the directory `~/TrainingSequences` contains your training sequences, the file `~/TrainingSequencesFileNames.txt` contains the file names of these sequences, and the Jellyfish binary is located in `~/jellyfish-2.1.1/bin/`. Then you can retrain the method by using:
```
mkdir Counts20
julia -p 10 FormCountsAndKCounts.jl -o ~/Counts20 -f ~/TrainingSequencesFileNames.txt -s ~/TrainingSequences/ -k 20 -j ~/jellyfish-2.1.1/bin/./jellyfish
julia -p 45 FormCommonKmerMatrix.jl -c ~/Counts20 -k 20 -s 100 -f ~/TrainingSequencesFileNames.txt -o ~/CommonKmerMatrix-20mers.h5
```



#### Python version ####
##### Running the Python version #####

##### Retraining the Python version #####
Usual workflow:
First split the Fasta files containing the training sequences (in this example: `./sequences/multiple-fasta-*.fa`) into Fasta files with exactly one sequence:
```
$ common-kmers.py split --prefix ./sequences/single-fasta_ ./sequences/multiple-fasta*
```

Then count the k-mers:
```
$ common-kmers.py count --kmer-length 7 --prefix ./counts/counts7_ ./sequences/single-fasta*
```

Finally, build the matrix with the common k-mers:
```
$ common-kmers.py common --output ./matrix7 ./counts/counts7_*
```

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
