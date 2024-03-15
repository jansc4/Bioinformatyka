import numpy as np
import matplotlib.pyplot as plt

class Dotplot:
    def __init__(self, array):
        self.array = array
        self.size = len(array)
        self.seq1Id = "1"
        self.seq2Id = "2"
        self.seq1Name = "Seq 1"
        self.seq2Name = "Seq 2"

    @classmethod
    def from_sequences(cls, seq1, seq2):
        len_seq1 = len(seq1.sequence)
        len_seq2 = len(seq2.sequence)

        dp = np.zeros((len_seq1, len_seq2), dtype=int)

        for i in range(len_seq1):
            for j in range(len_seq2):
                if seq1.sequence[i] == seq2.sequence[j]:
                    dp[i, j] = 1

        # Utwórz nową instancję klasy Dotplot i przekaż dp jako argument konstruktora
        instance = cls(dp)

        # Ustaw atrybuty instancji na podstawie przekazanych sekwencji
        instance.seq1Id = seq1.id
        instance.seq2Id = seq2.id
        instance.seq1Name = seq1.name
        instance.seq2Name = seq2.name

        return instance

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

    def graphic(self, title="Macierz kropkowa"):
        plt.imshow(self.array, cmap='binary', interpolation='nearest')
        plt.title(title)
        plt.xlabel(self.seq2Name + " " + self.seq2Id)
        plt.ylabel(self.seq1Name + " " + self.seq1Id)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
        plt.show()