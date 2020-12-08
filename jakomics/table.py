from jakomics.file import FILE
import pandas as pd
from natsort import natsorted


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
            self.df = pd.read_csv(self.file.file_path, sep="\t", comment="#")

    def column_value_counts(self, column, remove_duplicates=False):
        '''
        takes a column, finds and tabulates the unique values. Returns a series
        '''
        df = self.df

        if remove_duplicates:
            # https://stackoverflow.com/questions/32093829/remove-duplicates-from-dataframe-based-on-two-columns-a-b-keeping-row-with-max
            df = df.sort_values('HMM_COVERAGE').drop_duplicates(
                subset=['LOCUS', 'HMM'], keep='last')

        df = df[column].value_counts()
        return df


def merge_value_counts(file_list, column, remove_duplicates=False):
    '''
    Take a list of files and counts/merges a column, returning a df
    '''
    shared_df = pd.DataFrame()

    for f in natsorted(file_list):
        df = TABLE(f)
        s = df.column_value_counts(column, remove_duplicates)

        if hasattr(df.file, 'short_name'):
            s = s.rename(df.file.short_name)
        else:
            s = s.rename(os.path.basename(f))

        shared_df = pd.merge(shared_df, s, how='outer', left_index=True, right_index=True)
        shared_df = shared_df.fillna(0)

    return shared_df


if __name__ == "__main__":
    file_list = [
        '/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/ARW/data/nrMAGs/faa/mPt_15.dbcan8.txt']
    # ,                 '/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/ARW/data/nrMAGs/faa/mPt_16.dbcan8.txt', '/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/ARW/data/nrMAGs/faa/mPt_17.dbcan8.txt']

    a = merge_value_counts(file_list, 'HMM', remove_duplicates=True)
    print(a)
