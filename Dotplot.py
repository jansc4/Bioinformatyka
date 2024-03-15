import numpy as np


class Dotplot:
    def __init__(self, array):
        self.array = array
        self.size = len(array)

    def __getitem__(self, index):
        n, m = index
        return self.array[n][m]

    def __setitem__(self, index, val):
        n, m = index
        self.array[n][m] = val

    def __str__(self):
        return str(self.array)

    def saveTxt(self, filename):
        np.savetxt(filename, self.array, fmt="%d")

# %%
