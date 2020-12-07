from jakomics.file import FILE


class IMAGE(FILE):

    def __str__(self):
        return "<JAKomics IMAGE class>"

    def edge_crop(self, dims):
        # dims are l,r,t,b

        pass


if __name__ == "__main__":
    t = IMAGE("/Users/kimbrel1/Dropbox/LLNL/Projects/BlueCarbon/analysis/images/nanosims/glycolate_1through18@_1_T_14N 12C.tif")
    t.view()
