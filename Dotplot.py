import numpy as np
import matplotlib.pyplot as plt

class Dotplot:
    def __init__(self, array, seq1, seq2, seq1Id="1", seq2Id="2", seq1Name="Seq 1", seq2Name="Seq 2"):
        self.array = array
        self.size = array.shape
        self.seq1Id = seq1Id
        self.seq2Id = seq2Id
        self.seq1Name = seq1Name
        self.seq2Name = seq2Name
        self.seq1 = seq1
        self.seq2 = seq2

    @classmethod
    def from_sequences(cls, seq1, seq2):
        len_seq1 = len(seq1.sequence)
        len_seq2 = len(seq2.sequence)

        dp = np.zeros((len_seq1, len_seq2), dtype=int)

        for i in range(len_seq1):
            for j in range(len_seq2):
                if seq1.sequence[i] == seq2.sequence[j]:
                    dp[i, j] = 1

        return cls(dp, seq1, seq2, seq1.id, seq2.id, seq1.name, seq2.name)

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

    def needleman_wunsch_support(self, match, mismatch, gap):
        high, width = self.array.shape
        wynik = np.zeros((high + 1, width + 1))
        supprot_matrix = np.zeros((high, width))


        for i in range(1, width + 1):
            wynik[i][0] = wynik[i-1][0]  + gap

        for i in range(1, high + 1):
            wynik[0][i] = wynik[0][i-1] + gap

        #print(wynik)

        for i in range(1, high + 1):
            for j in range(1, width + 1):
                match_score = 0
                match_score = match if self.array[i - 1, j - 1] else mismatch  # Indeksowanie od zera
                U = wynik[i - 1, j] + gap
                L = wynik[i, j - 1] + gap
                D = wynik[i - 1, j - 1] + match_score
                wynik[i, j] = max(U, L, D)
                if wynik[i, j] == U:
                    supprot_matrix[i - 1, j - 1] = 1
                elif wynik[i, j] == L:
                    supprot_matrix[i - 1, j - 1] = 2
                elif wynik[i, j] == D:
                    supprot_matrix[i - 1, j - 1] = 0

         #zapis do pliku dla widocznosci
        np.savetxt("needleman", wynik, fmt="%d")
        print(wynik)

        '''
        D = 0       do algorytmu i-1,j-1 przepisz literkę
        U = 1       do algorytmu i-1, j wstaw gap w gornej sekwencji
        L = 2       do algorytmu i, j-1 wstaw gap w lewej sekwencji
        
        '''
        print("####################")
        print(supprot_matrix)
        np.savetxt("support_matrix", supprot_matrix, fmt="%d")

        return supprot_matrix, wynik

    def needleman_wunsch(self, match, mismatch, gap):
        seq1 = ""
        seq2 = ""

        supprot_matrix, wynik = self.needleman_wunsch_support(match, mismatch, gap)
        high, width = supprot_matrix.shape
        i = high -1
        j = width -1
        while i and j >= 0:

            if supprot_matrix[i, j] == 0:
                seq1 += self.seq1.sequence[i]
                seq2 += self.seq2.sequence[j]
                i = i - 1
                j = j - 1
            elif supprot_matrix[i, j] == 1:
                seq1 += self.seq1.sequence[i]
                seq2 += "_"
                i = i - 1
            elif supprot_matrix[i, j] == 2:
                seq1 += "_"
                seq2 += self.seq1.sequence[j]
                j = j - 1

        reversed_seq1 = ''.join(reversed(seq1))
        reversed_seq2 = ''.join(reversed(seq2))
        return reversed_seq1, reversed_seq2


