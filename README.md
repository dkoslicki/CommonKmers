# CommonKmers #
This is the source code for CommonKmers.

## What is CommonKmers? ##
CommonKmers is a k-mer based bacterial community reconstruction technique that utilizes sparsity promoting ideas from the field of compressed sensing to reconstruct the composition of a bacterial community. This method allows for strain-level abundance estimation, and can quantify the evolutionary distance between organisms in the sample and in the training database (thereby allowing for successful classification even with incomplete training data).


## How Do I Install CommonKmers? ##
There are two versions available, one written in Python and one in Julia. The Julia version is nearly two orders of magnitude faster than the python version.

#### Python version ####
##### Running the Julia version #####
To classify a sample using the Python version, use the command

#### Julia version ####
Please refer to [the Julia installation page](http://julialang.org/downloads/) to install Julia.
You will need to add the HDF5 and ArgParse packages. These can be added using `Pkg.add("HDF5")` and `Pkg.add("ArgParse")`.

##### Running the Julia version #####
To classify a sample using the Julia version, use the commands

##### Retraining the Julia version #####
To train on a custom database, use the functions ` FormCountsAndKCounts.jl` and `FormCommonKmerMatrix.jl`. For example, say the directory `~/TrainingSequences` contains your training sequences, and the file `TrainingSequencesFileNames.txt` contains the file names of these sequences. Then you can train the method by using:
```
mkdir Counts20
julia -p 10 FormCountsAndKCounts.jl -o ~/ -f Test/testFileNames.txt -s ../microbes/ -k 21 -j /home/pi/koslickd/jellyfish-2.1.1/bin/./jellyfish

```





Please read the directions on the [installation page](doc/install.markdown).


## How Do I use Quikr? ##
We have several ways to use quikr. Quikr is first and formost a command
line utility, but we also provide python and matlab scripts.

+ [Command Line Utilities](doc/cli.markdown)
+ [Matlab documentation](doc/matlab.markdown)
+ [Python documentation](doc/python.markdown)


## Contact ##
For issues with the quikr software, contact gailro@gmail.com


## License ##
The Quikr project is released under the GPL-3 License. Please view the [LICENSE](LICENSE)
file for more details.


## Contributors ##
+ David Koslicki
+Daniel Alonso Alemany
+Daniel Falush
+Nam Nguyen
