from jakomics.file import FILE
import pandas as pd


class TABLE(FILE):

    '''
    Initiated either as a file or directly from a pandas dataframe.

    If passed as a file, a dataframe gets made, plus access to all of the FILE
    class options are available.

    If passed as a dataframe, only df specific operations are available, and not
    any FILE class methods.
    '''

    def __str__(self):
        return "<JAKomics TABLE class>"

    def __init__(self, df):
        if isinstance(df, pd.DataFrame):
            self.df = df
        else:
            self.file = FILE(df)
            self.df = pd.read_csv(self.file.file_path, sep="\t")

    def column_value_counts(self, column):
        '''
        takes a column, finds and tabulates the unique values. Returns a series
        '''
        df = self.df[column].value_counts()
        return df


if __name__ == "__main__":
    t = TABLE("/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/ARW/data/nrMAGs/faa/mPt_1.faa.dbcan8.txt")
    d = t.column_value_counts('SUBSTRATE')
    print(d)
    print(t.file)
