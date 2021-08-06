import pandas as pd
from operator import itemgetter


class Sequence:
    def __init__(self, df, sequence_column='Sequence', class_column='class'):
        self.df = df  # initializing the DataFrame
        self.sequence_column = sequence_column  # the column which contains the Sequences
        self.class_column = class_column

    """"" 
    removing residues which are repeating in the specific position along all sequences
    example: ACGTTACT
    """""
    def DropRepeating(self):
        position = []
        for i in range(len(self.df[self.sequence_column][0])):
            if len((self.df[self.sequence_column].astype(str).str[i]).unique()) != 1:
                position.append(i)
        self.df[self.sequence_column] = self.df[self.sequence_column].apply(lambda x: ''.join(itemgetter(*position)(x)))
        return self.df
    #
    # # method to find the amino acids that appear a long the dataset
    # def AminoOccurrence(self):
    #     amino_dict = {}
    #     amino_list = list('ACDEFGHIKLMNPQRSTVWY')
    #     for amino in amino_list:
    #         amino_dict[amino] = df.Sequence.str.count(amino).sum()
    #     return amino_dict
