import numpy as np
import matplotlib.pyplot as plt
from app.Classes.Sequence import Sequence

"""
Dotplot class contain methods for plotting sequences.
    @:param Sequence: Sequence 1 to be compared
    @:param Sequence: Sequence 2 to be compared
"""


class Dotplot:
    def __init__(self, sequence1, sequence2):
        self.seq1 = sequence1.getSequence()
        self.seq2 = sequence2.getSequence()
        self.seq1Id = sequence1.seq1.get_id()
        self.seq2Id = sequence2.seq2.get_id()
        self.seq1Name = sequence1.get_name()
        self.seq2Name = sequence2.get_name()
        self.array = self._array()
        self.size = self.array.shape

    def _array(self):
        len_seq1 = len(self.seq1)
        len_seq2 = len(self.seq2)
        dp = np.zeros((len_seq1, len_seq2), dtype=int)

        for i in range(len_seq1):
            for j in range(len_seq2):
                if self.seq1[i] == self.seq2[j]:
                    dp[i, j] = 1
        return dp

    def __getitem__(self, index):
        return self.array[index]

    def __setitem__(self, index, val):
        self.array[index] = val

    def __str__(self):
        return str(self.array)

    def graphic(self, filename, title="Dotplot Matrix", window=None, threshold=None):
        plt.imshow(self.array, cmap='binary', interpolation='nearest')
        plt.title(title)
        plt.xlabel(f"{self.seq2Name} {self.seq2Id}")
        plt.ylabel(f"{self.seq1Name} {self.seq1Id}")
        plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')

        # Dodanie skali
        plt.colorbar(label='Intensity')

        # Dodanie legendy z parametrami filtracji
        if window is not None and threshold is not None:
            plt.text(0.9, -0.1, f"Window: {window}, Threshold: {threshold}",
                     horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

        plt.savefig(filename)
        plt.show()

    def filter(self, window, threshold):
        wynik = np.zeros_like(self.array)
        high, width = self.array.shape

        for i in range(high - window + 1):
            for j in range(width - window + 1):
                sum_diag = np.sum(np.diagonal(self.array[i:i+window, j:j+window]))
                if sum_diag >= threshold:
                    wynik[i, j] = 1

        self.array = wynik

    def save_txt(self, filename):
        np.savetxt(filename, self.array, fmt="%d")
