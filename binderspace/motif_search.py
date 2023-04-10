import getopt
import pandas as pd
import re
import multiprocessing as mp
import time
import sys
import os
import operator

arguments = sys.argv[1:]
posfile = -1  # the input of positive sequence file (required)
negfile = -1  # the input of negative sequence file (optional)
motif_file = 'motifs'  # the output file name (optional)
minfreq = -1  # the minimum frequency of motifs (optional)
minlength = -1  # the minimum length of the motif (optional)
maxlength = -1  # the maximum length of the motif (optional)
maxgaps = 0  # the maximum number of gaps (optional)
maxgaplength = 1  # the maximum length of gaps (optional)
processor = 20  # processors used for the multiprocessing running (optional)
category = 'amino_acids'  # search motifs for dna or amino_acids (required)


def process_options(args):
    global posfile, negfile, motif_file, minfreq, minlength, maxlength, maxgaps, maxgaplength, processor, category
    opts, args = getopt.getopt(args, "i:n:o:f:m:l:g:a:p:c:")
    for opt, arg in opts:
        if opt in ['-i']:
            if not os.path.exists(arg):
                sys.exit(f"File {arg} (option -i) not found.\n")
            posfile = arg
        elif opt in ['-n']:
            if not os.path.exists(arg):
                sys.exit(f"File {arg} (option -p) not found.\n")
            negfile = arg
        elif opt in ['-o']:
            motif_file = arg
        elif opt in ['-f']:
            minfreq = float(arg)
        elif opt in ['-m']:
            minlength = int(arg)
        elif opt in ['-l']:
            maxlength = int(arg)
        elif opt in ['-g']:
            maxgaps = int(arg)
        elif opt in ['-a']:
            maxgaplength = int(arg)
        elif opt in ['-p']:
            processor = int(arg)
        elif opt in ['-c']:
            category = arg
        else:
            sys.exit(f"Unknown option{opt}. Run the program without any options to see its usage.\n")


def check_arguments():
    if posfile == -1:
        sys.exit('No csv file with (positive) sequences given. Please set option -i.\n')


def scan(posset, cand_dict):
    """input a dictionary with candidates and indexes, out put a new dictionary with candidates and indexed {candidates:
    list of the indexes}"""
    cand_new_dict = {}
    for key, value in cand_dict.items():
        cand = key
        index_pool = value
        new_cands = add_elements(cand, index_pool)
        for candstring in new_cands:
            candstring_new = candstring.replace('*', gap)
            new_index = []
            count = 0
            for i in index_pool:
                if re.search(candstring_new, posset[i]):
                    count += 1
                    new_index.append(i)
            if count >= min_freq:
                cand_new_dict[candstring] = new_index
    return cand_new_dict


def add_elements(candidate, index):
    """add one base/gap to the old candidate and generate a list of new candidates"""
    new_cands = []
    n_gaps = number_gaps(candidate)
    if len(candidate) < max_length and len(index) >= min_freq:
        new_cands.extend(candidate + base for base in bases)
        if max_gaps > 0 and n_gaps < max_gaps:
            new_cands.append(candidate + '*')
    return new_cands


def number_gaps(candidate):
    """check the number of gaps in the candidates"""
    n_gaps = 0
    for base in candidate:
        if base == '*':
            n_gaps += 1
    return n_gaps


def prune(cand_dict):
    """ filter out motifs look like X****, X**, XX**"""
    dict_new = {}
    for k, v in cand_dict.items():
        candstring = k
        temp = list(candstring)
        while temp[-1] == '*':
            temp.pop()
        if len(temp) >= min_length:
            dict_new[k] = v
    dict_new2 = {}
    for k in sorted(dict_new, key=lambda k: len(dict_new[k]), reverse=True):
        dict_new2[k] = dict_new[k]

    return dict_new2


def prepare_cand_dict0(posset, bases):
    """prepare the input for multiple processing, for the dna sequences, scan twice to get more candidates,
    and fully use the processors"""
    index0 = list(range(0, len(posset)))
    cand_dict = {}
    for base in bases:
        new_index = []
        count = 0
        for i in index0:
            if re.search(base, posset[i]):
                count += 1
                new_index.append(i)
        if count >= min_freq:
            cand_dict[base] = new_index

    cand_dict_list = []
    if category == 'dna':
        cand_dict_2 = scan(posset, cand_dict)
        for k, v in cand_dict_2.items():
            cand_dict_list.append({k: v})
    else:
        for k, v in cand_dict.items():
            cand_dict_list.append({k: v})

    return cand_dict_list


#
def add_neg(cand_dict):
    """check negative occurrence, and add the occurrence index after the positive index"""
    new_can_dict = {}
    for k, v in cand_dict.items():
        cand_str = k.replace('*', gap)
        index_list = []
        for i in range(0, len(negset)):
            if re.search(cand_str, negset[i]):
                index_list.append(i)
        new_can_dict[k] = [v, index_list]
    return new_can_dict


def rescan(re_dict, cand_dict, posset):
    """recall the scan function until get the defined length motifs and append all the {motifs:index} in a dictionary"""
    a = [{}] * max_length
    a[0] = cand_dict
    for i in range(1, max_length):
        a[i] = scan(posset, a[i - 1])
        if negfile == -1:
            re_dict.update(prune(a[i]))
        else:
            re_dict.update(add_neg(prune(a[i])))
    return re_dict


def generate_df(re_dict):
    """convert the dictionary result to a dataframe"""
    print("motif search finished, start generating dataframe")
    data = []
    if negfile == -1:
        for i in re_dict.keys():
            pos_perc = (len(re_dict[i]) / len(posset)) * 100
            data.append([i, len(re_dict[i]), pos_perc, re_dict[i]])
            df = pd.DataFrame(data, columns=['motifs', 'positive_occurrence', 'pos_percentage (%)', 'index'])
            df2 = df.sort_values('pos_percentage', ascending=False)
            df2.drop(['index'], axis=1, inplace=True)

    else:
        for i in re_dict.keys():
            pos_perc = (len(re_dict[i][0]) / len(posset)) * 100
            neg_perc = (len(re_dict[i][1]) / len(negset)) * 100
            data_chi = np.array([[len(re_dict[i][0]),len(re_dict[i][1])],
                                 [len(posset)-len(re_dict[i][0]),len(negset)-len(re_dict[i][1])]])
            chi2, p_value, degrees_of_freedom, expected_values = chi2_contingency(data_chi)
            
            if pos_perc > neg_perc:
                data.append(
                    [i, len(re_dict[i][0]), len(re_dict[i][1]), pos_perc, neg_perc, pos_perc - neg_perc, re_dict[i][0], re_dict[i][1],round(chi2, 2),p_value])
        df = pd.DataFrame(data, columns=['motifs', 'positive_occurrence', 'negative_occurrence ', 'pos_percentage (%)',
                                         'neg_percentage (%)', 'difference (%)', 'index_pos', 'index_neg','chi2','p-value'])
        df2 = df.sort_values('difference (%)', ascending=False)
        df2.drop(['index_pos','index_neg'], axis=1, inplace=True)

    return df2


def print_csv(dataframe):
    print("printing the motifs into different csv files regarding to their length")
    for i in range(min_length, max_length + 1):
        df_sub = dataframe[dataframe.motifs.str.len() == i]
        if not df_sub.empty:
            df_sub.to_csv('{}_length {}.csv'.format(motif_file, i), index=False)


def process_arguments():
    """read input file, process the default parameters"""
    global posset, negset, min_freq, max_length, min_length, max_gaps, max_gaplength, n_cores, bases, gap
    process_options(arguments)
    check_arguments()

    df_pos = pd.read_csv(posfile)
    posset = list(df_pos['Sequence'])
    print('{} positive sequences read'.format(len(posset)))

    if negfile == -1:
        print('no negative sequences read')
    else:
        df_neg = pd.read_csv(negfile)
        negset = list(df_neg['Sequence'])
        print('{} negative sequences read'.format(len(negset)))
    n_cores = processor
    print("{} cores used for running this task".format(n_cores))
    if len(posset) > 1000 and minfreq == -1:
        min_freq = int(0.001 * len(posset))

    elif len(posset) < 1000 and minfreq == -1:
        print('Error: need to give frequency value (-f option) if the sequence number less than 1000')
        exit()
    else:
        min_freq = int(minfreq * len(posset))
    print("minimal frequency {}".format(min_freq))
    if maxlength == -1:
        max_length = len(posset[0])
    else:
        max_length = maxlength
    print("maximum length {}".format(max_length))
    if minlength == -1:
        min_length = 3
    else:
        min_length = minlength
    print("minimal length {}".format(min_length))

    max_gaps = maxgaps
    print("maximum number of gaps {}".format(max_gaps))

    max_gaplength = maxgaplength
    print("maximum length of the gaps {}".format(max_gaplength))

    if category == 'dna':
        bases = ['A', 'T', 'G', 'C']
        gap = '[ATGC]'
    elif category == 'amino_acids':
        bases = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',
                 'S', 'T', 'V', 'W', 'Y']
        gap = '[ACDEFGHIKLMNPQRSTVWY]'
    print("bases: {}".format(bases))


def main():
    start_scan = time.time()
    if len(arguments) <= 0:
        print("To run the program, please define at least the required options:\n\n")
        print("  -i posfile          # the file with positive sequences (REQUIRED)\n")
        print("  -n negfile          # the file with negative sequences \n")
        print("  -o motif_file		# the name of the output file (default: 'motifs')\n")
        print("  -f minfreq         # the minimal frequency (absolute number) for the positive sequences "
              "(default: 0.01*(positive sequences number))\n")
        print("  -m minlength          # the minimal length of the motifs (default: 3) \n")
        print("  -l maxlength        # the maximal length of the motifs (default: the length of the "
              "positive sequences)\n")
        print("  -c category        # choose from 'dna' or 'amino_acids',(default: 'amino_acids')\n")
        print("  -g maxgaps	     # maximal number of gaps (default: 0)\n")
        print("  -a maxgaplength    # maximal gap length (default: 1)\n")
        print("  -p processor      # number of processors used for running the program (default: 20)\n\n")
        print("See manual for more information on the options.\n\n")
        sys.exit()
    else:
        process_arguments()
        print("start searching ...")
        manager = mp.Manager()
        return_dict = manager.dict()

        pool = mp.Pool(n_cores)

        cand_dict_list = prepare_cand_dict0(posset, bases)
        for cand_dict in cand_dict_list:
            pool.apply_async(rescan, args=(return_dict, cand_dict, posset))

        pool.close()
        pool.join()
        print("time for motif search: ", time.time() - start_scan)

        df = generate_df(return_dict)
        print_csv(df)


if __name__ == '__main__':
    start = time.time()
    main()
    print("total time: ", time.time() - start)

