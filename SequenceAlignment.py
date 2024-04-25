import numpy as np

class SequenceAlignment:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.array = np.zeros((len(seq1), len(seq2)))

    def needleman_wunsch_support(self, match, mismatch, gap):
        high, width = self.array.shape
        wynik = np.zeros((high + 1, width + 1))
        support_matrix = np.zeros((high, width))

        for i in range(1, high + 1):
            wynik[i][0] = wynik[i-1][0] + gap

        for i in range(1, width + 1):
            wynik[0][i] = wynik[0][i-1] + gap

        for i in range(1, high + 1):
            for j in range(1, width + 1):
                match_score = match if self.seq1[i - 1] == self.seq2[j - 1] else mismatch
                U = wynik[i - 1, j] + gap
                L = wynik[i, j - 1] + gap
                D = wynik[i - 1, j - 1] + match_score
                wynik[i, j] = max(U, L, D)
                if wynik[i, j] == U:
                    support_matrix[i - 1, j - 1] = 1
                elif wynik[i, j] == L:
                    support_matrix[i - 1, j - 1] = 2
                else:
                    support_matrix[i - 1, j - 1] = 0

        np.savetxt("needleman", wynik, fmt="%d")
        print(wynik)
        print("####################")
        print(support_matrix)
        np.savetxt("support_matrix", support_matrix, fmt="%d")

        return support_matrix, wynik

    def needleman_wunsch(self, match, mismatch, gap):
        seq1 = ""
        seq2 = ""

        support_matrix, wynik = self.needleman_wunsch_support(match, mismatch, gap)
        high, width = support_matrix.shape
        i, j = high - 1, width - 1

        while i >= 0 and j >= 0:
            if support_matrix[i, j] == 0:
                seq1 += self.seq1[i]
                seq2 += self.seq2[j]
                i -= 1
                j -= 1
            elif support_matrix[i, j] == 1:
                seq1 += self.seq1[i]
                seq2 += "_"
                i -= 1
            elif support_matrix[i, j] == 2:
                seq1 += "_"
                seq2 += self.seq2[j]
                j -= 1

        return seq1[::-1], seq2[::-1]  # Reversed sequences are returned