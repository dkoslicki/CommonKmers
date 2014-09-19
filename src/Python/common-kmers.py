#! /usr/bin/env python

import argparse, h5py, numpy, os, os.path, subprocess, sys
from Bio import SeqIO
from multiprocessing import Pool

def main():
    parser = argparse.ArgumentParser(description = 'Count occurrences of \
common k-mers between pairs of sequences.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog = 'Usual workflow:\n\
===============\n\
First split the Fasta files containing the sequences (in this example: \n\
"./sequences/multiple-fasta-*.fa") into Fasta files with exactly one sequence:\n\n\
$ common-kmers.py split --prefix ./sequences/single-fasta_ ./sequences/multiple-fasta*\n\n\
\
Then count the k-mers:\n\n\
$ common-kmers.py count --kmer-length 7 --prefix ./counts/counts7_ ./sequences/single-fasta*\n\n\
\
Finally, build the matrix with the common k-mers:\n\n\
$ common-kmers.py common --output ./matrix7 ./counts/counts7_*')

    subparsers = parser.add_subparsers(title = 'subcommands',
                                       help = 'additional help')

    # Split sequences
    parser_split = subparsers.add_parser('split', description = 'Split the input \
                   fasta files such that each output file only contains a sequence.')
    parser_split.add_argument('files', metavar = 'file', nargs='+',
                              help = 'Input fasta file with genome sequences.')
    parser_split.add_argument('-p', '--prefix', default = './seq_',
                              help = 'Prefix of the generated files (default "./seq_"). \
                              The generated files will be PREFIX_0.fa, PREFIX_1.fa...')
    parser_split.set_defaults(func = split_fasta)


    # Count k-mers
    parser_count = subparsers.add_parser('count', description = 'Count k-mers of a given length.')
    parser_count.add_argument('files', nargs='+', metavar = 'file',
                              help = 'Input fasta file with exactly one sequence.')
    parser_count.add_argument('-k', '--kmer-length', type = positive_int, required = True,
                              help = 'Length of the k-mers.')
    parser_count.add_argument('-p', '--prefix', default = './seq_',
                              help = 'Prefix for the generated files (default: "./seq_"). \
                              The generated files will be PREFIX_0.kcount, PREFIX_1.kcount...')
    parser_count.add_argument('-s', '--initial-hash-size', default = '100M',
                              help = 'Initial hash size for Jellyfish (default: 100M).')
    parser_count.add_argument('-t', '--num-threads', type = positive_int, default = 1,
                              help='Number of Jellyfish instances to run in parallel (default: 1).')
    parser_count.add_argument('-j', '--jellyfish-threads', type = positive_int, default = 1,
                              help='Number of threads used by each Jellyfish (default 1).')
    parser_count.add_argument('-b', '--jellyfish-binary', default = 'jellyfish',
                              help='Path to the jellyfish binary (default= "jellyfish").')
    parser_count.set_defaults(func = count_kmers)


    # Common k-mers
    parser_common = subparsers.add_parser('common', description = 'Count common occurrences \
                    of k-mers for pairs of sequences.')
    parser_common.add_argument('files', nargs='+', metavar = 'file',
                               help = 'File obtained with "common-kmers.py count".')
    parser_common.add_argument('-o', '--output', help = 'Output file.')
    parser_common.add_argument('-t', '--num-threads', type = positive_int, default = 1,
                              help = 'Number of threads used when building the matrix (default: 1).')
    parser_common.set_defaults(func = common_kmers)

    args = parser.parse_args()
    args.func(args)


def split_fasta(args):
    num = 0
    for f in args.files:
        with open(f, 'rU') as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                name = args.prefix + '_' + str(num) + '.fa'
                SeqIO.write(seq, name, 'fasta')
                num += 1

def positive_int(string):
    try:
        value = int(string)
        if (value < 1):
            raise Exception
    except:
        msg = "%r is not a positive integer." % string
        raise argparse.ArgumentTypeError(msg)
    return value


def count_kmers(args):
    kmer_args = { 'kmer_length' : args.kmer_length,
                  'hash_size' : args.initial_hash_size,
                  'jellyfish_threads' : args.jellyfish_threads,
                  'jellyfish_binary' : args.jellyfish_binary,
                  'prefix' : args.prefix }

    n = len(args.files)
    #with_args = [(kmer_args, i, args.files[i]) for i in range(n)]
    with_args = [(kmer_args, f) for f in args.files]
    pool = Pool(processes = args.num_threads)
    results = pool.map(count_kmers_wrapper, with_args)
    if False in results:
        sys.stderr.write('There was a problem counting kmers.\n')
        exit(1)


def count_kmers_wrapper(arg):
    "In Python 2 there seems to be no way to forward the Exception\
    from a Process in a Pool to its parent. Print to stderr and\
    return a Boolean instead."
    try:
        count_kmers_worker(arg)
    except Exception as e:
        sys.stderr.write(str(e) + '\n')
        return False
    return True

def count_kmers_worker(arg):
    args, fasta_file = arg
    i = get_number(fasta_file)
    iden = None
    with open(fasta_file, 'rU') as handle:
        for seq in SeqIO.parse(handle, 'fasta'):
            if (iden is None):
                iden = seq.id
            else:
                raise Exception('There are multiple sequences in file ' + fasta_file)
    if iden is None:
        raise Exception('There is no sequence in file ' + fasta_file)

    jf_file = args['prefix'] + str(i) + '.jf'
    dump_file = args['prefix'] + str(i) + '.djf'
    integer_file = args['prefix'] + str(i) + '.kcount'

    count = args['jellyfish_binary'] + ' count -m ' + str(args['kmer_length']) +\
            ' -s ' + args['hash_size'] + ' -t ' + str(args['jellyfish_threads']) +\
            ' -o ' + jf_file + ' ' + fasta_file

    subprocess.check_call(count, shell = True)

    dump = args['jellyfish_binary'] + ' dump -c -o ' + dump_file + ' ' + jf_file
    try:
        subprocess.check_call(dump, shell = True)
    finally:
        None
        #os.remove(jf_file)

    try:
        kmer2int(dump_file, integer_file, identifier = iden, kmer_length = args['kmer_length'])
    finally:
        None
        #os.remove(dump_file)

def get_number(filename):
    "Return the number at the end of the filename (before the extension)."
    base = os.path.splitext(os.path.basename(filename))[0]
    for i in range(len(base) - 1, -1, -1):
        if (not base[i].isdigit()):
            return base[i+1:]

def kmer2int(dump_file, integer_file, identifier, kmer_length):
    """Encode the kmers as Ints and store them in a h5py file, with the
    structure:

    Attributes
    ----------
    '/sequence_id'
    '/kmer_length'
    '/total_count': number of k-mers in the sequence

    Data
    ----
    '/kmers': (n,2) array where first column contains the
    kmers (encoded as Ints) and second column is their count.

    """
    kmers = []
    total = 0
    with open(dump_file, "r") as f:
        for line in f:
            (kmer, _, count) = line.partition(" ")
            i = encode(kmer)
            icount = int(count)
            kmers.append((i, icount))
            total += icount
        kmers.sort()
        array = numpy.array(kmers)

        with h5py.File(integer_file, 'w') as outfile:
            outfile.attrs['sequence_id'] = identifier
            outfile.attrs['kmer_length'] = kmer_length
            outfile.attrs['total_count'] = total
            dset = outfile.create_dataset('kmers', data = array)


def encode(kmer):
    """Encode a kmer to an Int, using 2 bits per character."""
    i = 0
    bits = 0
    for c in kmer:
        if c == 'A':
            bits = 0
        elif c == 'C':
            bits = 1
        elif c == 'T':
            bits = 2
        elif c == 'G':
            bits = 3
        else:
            raise Exception('Found non ACTG character in ' + kmer)
        i = (i << 2) + bits
    return i

def decode(num, k):
    """Given an Int num and the length of the k-mer it encodes, decode num into a string."""
    string = ''
    for i in range(k):
        bits = num % 4
        c = ''
        if bits == 0:
            c = 'A'
        elif bits == 1:
            c = 'C'
        elif bits == 2:
            c = 'T'
        else:
            c = 'G'
        string = c + string
        num = num >> 2
    return string

def common_kmers(args):
    """ Form a matrix where the (i,j)th entry is the sum, over all the
    kmers common to sequence i and j, of the occurrence of said kmer
    in sequence j. Store it in a h5py file with the structure:

    Attributes
    ----------
    'kmer_length'

    Data
    ----
    '/sequence_ids': sequence identifiers
    '/common_kmers': matrix

    """
    n = len(args.files)
    ids = []
    klength = None
    for f in args.files:
        with h5py.File(f, 'r') as hFile:
            k = hFile.attrs['kmer_length']
            if klength is None:
                klength = k
            elif klength != k:
                raise Exception('Different k-mer lengths: %s (%d), %s (%d).' %\
                                (args.files[0], klength, f, k))
            ids.append(hFile.attrs['sequence_id'])
    array_ids = numpy.array(ids)

    upper_diagonal_files = [args.files[i : ] for i in range(n)]
    pool = Pool(processes = args.num_threads)
    rows_and_columns = pool.map(common_kmers_worker, upper_diagonal_files)
    matrix = numpy.empty([n, n], dtype = numpy.int)
    for i in range(n):
        matrix[i, i:] = rows_and_columns[i][0] # row
        matrix[i:, i] = rows_and_columns[i][1] # column

    with h5py.File(args.output, 'w') as outfile:
        outfile.attrs['kmer_length'] = klength
        outfile.create_dataset('sequence_ids', data = array_ids)
        outfile.create_dataset('common_kmers', data = matrix)


def common_kmers_worker(files):
    """For the matrix n * n, given the files [i..n] compute row i and column i
       (only the slices [(i, i)..(i,n)] and [(i,i)..(n,i)]
    """

    with h5py.File(files[0], 'r') as iFile:
        iKmers = iFile['kmers'].value
        irow = [iFile.attrs['total_count']]
        icolumn = [iFile.attrs['total_count']]
        for jName in files[1 : ]:
            with h5py.File(jName, 'r') as jFile:
                jKmers = jFile['kmers'].value
                iTotal, jTotal = countPair(iKmers, jKmers)
                irow.append(jTotal)
                icolumn.append(iTotal)
        return (irow, icolumn)


def countPair(iKmers, jKmers):
    iMax = len(iKmers)
    jMax = len(jKmers)
    i = 0
    j = 0
    iTotal = 0
    jTotal = 0

    while (i < iMax) and (j < jMax):
        if iKmers[i,0] == jKmers[j,0]:
            iTotal += iKmers[i,1]
            jTotal += jKmers[j,1]
            i += 1
            j += 1
        elif iKmers[i,0] < jKmers[j,0]:
            i += 1
        else:
            j += 1

    return (iTotal, jTotal)

if __name__ == '__main__':
    main()