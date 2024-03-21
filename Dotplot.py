import numpy as np
import matplotlib.pyplot as plt

class Dotplot:
    def __init__(self, array, seq1Id="1", seq2Id="2", seq1Name="Seq 1", seq2Name="Seq 2"):
        self.array = array
        self.size = array.shape
        self.seq1Id = seq1Id
        self.seq2Id = seq2Id
        self.seq1Name = seq1Name
        self.seq2Name = seq2Name

    @classmethod
    def from_sequences(cls, seq1, seq2):
        len_seq1 = len(seq1.sequence)
        len_seq2 = len(seq2.sequence)

        dp = np.zeros((len_seq1, len_seq2), dtype=int)

        for i in range(len_seq1):
            for j in range(len_seq2):
                if seq1.sequence[i] == seq2.sequence[j]:
                    dp[i, j] = 1

        return cls(dp, seq1.id, seq2.id, seq1.name, seq2.name)

    def __getitem__(self, index):
        return self.array[index]

    def __setitem__(self, index, val):
        self.array[index] = val

    def __str__(self):
        return str(self.array)

    def save_txt(self, filename):
        np.savetxt(filename, self.array, fmt="%d")

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

        return Dotplot(wynik, self.seq1Id, self.seq2Id, self.seq1Name, self.seq2Name)

