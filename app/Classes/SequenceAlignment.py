import numpy as np
from app.Classes.Sequence import Sequence
import matplotlib.pyplot as plt
"""
SequenceAlignment class contain methods to align and compare sequences.
    @:param Sequence: Sequence 1 to be compared
    @:param Sequence: Sequence 2 to be compared
"""


class SequenceAlignment:
    def __init__(self, sequence1, sequence2):
        self.seq1 = sequence1.get_sequence()
        self.seq2 = sequence2.get_sequence()
        self.seq1Id = sequence1.get_id()
        self.seq2Id = sequence2.get_id()
        self.seq1Name = sequence1.get_name()
        self.seq2Name = sequence2.get_name()
        self.array = self._array()

    def _array(self):
        len_seq1 = len(self.seq1)
        len_seq2 = len(self.seq2)
        return np.zeros((len_seq1, len_seq2), dtype=int)

    def _needleman_wunsch_support(self, match, mismatch, gap):
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
                upper = wynik[i - 1, j] + gap
                left = wynik[i, j - 1] + gap
                diagonal = wynik[i - 1, j - 1] + match_score
                wynik[i, j] = max(upper, left, diagonal)
                if wynik[i, j] == upper:
                    support_matrix[i - 1, j - 1] = 1
                elif wynik[i, j] == left:
                    support_matrix[i - 1, j - 1] = 2
                else:
                    support_matrix[i - 1, j - 1] = 0

        np.savetxt(r"D:\bioinformatyka\Lista 1\app\Output\needleman", wynik, fmt="%d")
        #print(wynik)
        #print("####################")
        #print(support_matrix)
        np.savetxt(r"D:\bioinformatyka\Lista 1\app\Output\needleman_support_matrix", support_matrix, fmt="%d")

        return support_matrix, wynik

    def needleman_wunsch(self, match, mismatch, gap, add_data = True):
        seq1 = ""
        seq2 = ""

        support_matrix, wynik = self._needleman_wunsch_support(match, mismatch, gap)
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

                # Odwróć sekwencje, aby zachować kolejność
        aligned_seq1 = seq1[::-1]
        aligned_seq2 = seq2[::-1]

        # Oblicz punktację dopasowania, długość dopasowania, procent identycznych pozycji i przerw
        alignment_score = 0
        identical_positions = 0
        gaps = 0
        for i in range(len(aligned_seq1)):
            if aligned_seq1[i] == aligned_seq2[i]:
                alignment_score += match
                identical_positions += 1
            elif aligned_seq1[i] == '_' or aligned_seq2[i] == '_':
                alignment_score += gap
                gaps += 1
            else:
                alignment_score += mismatch

        alignment_length = len(aligned_seq1)
        identical_positions_percent = (identical_positions / alignment_length) * 100
        gaps_percent = (gaps / alignment_length) * 100
        if add_data:
            # Wyświetl parametry programu i wyniki
            print("Parametry programu dla needleman_wunsch:")
            print("Match:", match)
            print("Mismatch:", mismatch)
            print("Gap:", gap)
            print("Wynik dopasowania:", alignment_score)
            print("Długość dopasowania:", alignment_length)
            print("Procent identycznych pozycji:", identical_positions_percent)
            print("Procent przerw:", gaps_percent)

            with open(r"D:\bioinformatyka\Lista 1\app\Output\parameters_needleman.txt", 'w') as f:
                f.write("Sekwencja 1: {}\n".format(aligned_seq1))
                f.write("Sekwencja 2: {}\n".format(aligned_seq2))
                f.write("Parametry programu:\n")
                f.write("Match: {}\n".format(match))
                f.write("Mismatch: {}\n".format(mismatch))
                f.write("Gap: {}\n".format(gap))
                f.write("Wynik dopasowania: {}\n".format(alignment_score))
                f.write("Długość dopasowania: {}\n".format(alignment_length))
                f.write("Procent identycznych pozycji: {:.2f}%\n".format(identical_positions_percent))
                f.write("Procent przerw: {:.2f}%\n".format(gaps_percent))

        # Zwróć dopasowane sekwencje oraz macierz wsparcia
        return (Sequence(aligned_seq1, self.seq1Name, self.seq1Id), Sequence(aligned_seq2, self.seq2Name, self.seq2Id),
                support_matrix, wynik)

    def _smith_waterman_support(self, match, mismatch, gap):
        high, width = self.array.shape
        wynik = np.zeros((high + 1, width + 1))
        support_matrix = np.zeros((high, width))

        for i in range(1, high + 1):
            wynik[i][0] = 0  # W algorytmie Smitha-Watermana zaczynamy od zera w każdym wierszu i kolumnie.

        for i in range(1, width + 1):
            wynik[0][i] = 0

        max_score = 0  # Zmienna do śledzenia maksymalnej wartości wyniku w macierzy.

        for i in range(1, high + 1):
            for j in range(1, width + 1):
                match_score = match if self.seq1[i - 1] == self.seq2[j - 1] else mismatch
                upper = max(wynik[i - 1, j] + gap, 0)  # Ustawiamy na 0, jeśli wynik jest mniejszy od zera.
                left = max(wynik[i, j - 1] + gap, 0)
                diagonal = max(wynik[i - 1, j - 1] + match_score, 0)
                wynik[i, j] = max(upper, left, diagonal)
                if wynik[i, j] > max_score:  # Aktualizujemy maksymalną wartość wyniku.
                    max_score = wynik[i, j]
                    max_i, max_j = i, j

                # Wprowadzamy zmiany w obsłudze macierzy wsparcia, jeśli wynik jest ujemny.
                if wynik[i, j] <= 0:
                    support_matrix[i - 1, j - 1] = -1
                elif wynik[i, j] == upper:
                    support_matrix[i - 1, j - 1] = 1
                elif wynik[i, j] == left:
                    support_matrix[i - 1, j - 1] = 2
                else:
                    support_matrix[i - 1, j - 1] = 0

        np.savetxt(r"D:\bioinformatyka\Lista 1\app\Output\smith-wateman", wynik, fmt="%d")
        #print(wynik)
        print("####################")
        #print(support_matrix)
        np.savetxt(r"D:\bioinformatyka\Lista 1\app\Output\smith_support_matrix", support_matrix, fmt="%d")

        return support_matrix, wynik, max_i, max_j  # Zwracamy także pozycję maksymalnego wyniku.

    def smith_waterman(self, match, mismatch, gap):
        seq1 = ""
        seq2 = ""

        support_matrix, wynik, max_i, max_j = self._smith_waterman_support(match, mismatch, gap)

        # Początek ścieżki - największy wynik w macierzy wynikowej.
        i, j = max_i, max_j

        while i > 0 and j > 0 and wynik[
            i, j] > 0:  # Przerywamy, gdy dojdziemy do krawędzi macierzy lub wynik stanie się ujemny.
            if support_matrix[i - 1, j - 1] == 0:  # Diagonalna, kontynuuj dopasowanie.
                seq1 += self.seq1[i - 1]
                seq2 += self.seq2[j - 1]
                i -= 1
                j -= 1
            elif support_matrix[i - 1, j - 1] == 1:  # Górna, wstawienie przerwy w pierwszej sekwencji.
                seq1 += self.seq1[i - 1]
                seq2 += "_"
                i -= 1
            elif support_matrix[i - 1, j - 1] == 2:  # Lewa, wstawienie przerwy w drugiej sekwencji.
                seq1 += "_"
                seq2 += self.seq2[j - 1]
                j -= 1

                # Odwróć sekwencje, aby zachować kolejność
        aligned_seq1 = seq1[::-1]
        aligned_seq2 = seq2[::-1]

        # Oblicz punktację dopasowania, długość dopasowania, procent identycznych pozycji i przerw
        alignment_score = 0
        identical_positions = 0
        gaps = 0
        for i in range(len(aligned_seq1)):
            if aligned_seq1[i] == aligned_seq2[i]:
                alignment_score += match
                identical_positions += 1
            elif aligned_seq1[i] == '_' or aligned_seq2[i] == '_':
                alignment_score += gap
                gaps += 1
            else:
                alignment_score += mismatch

        alignment_length = len(aligned_seq1)
        identical_positions_percent = (identical_positions / alignment_length) * 100
        gaps_percent = (gaps / alignment_length) * 100

        # Wyświetl parametry programu i wyniki
        print("Parametry programu dla smith_waterman:")
        print("Match:", match)
        print("Mismatch:", mismatch)
        print("Gap:", gap)
        print("Wynik dopasowania:", alignment_score)
        print("Długość dopasowania:", alignment_length)
        print("Procent identycznych pozycji:", identical_positions_percent)
        print("Procent przerw:", gaps_percent)

        with open(r"D:\bioinformatyka\Lista 1\app\Output\parameters_smith.txt", 'w') as f:
            f.write("Sekwencja 1: {}\n".format(aligned_seq1))
            f.write("Sekwencja 2: {}\n".format(aligned_seq2))
            f.write("Parametry programu:\n")
            f.write("Match: {}\n".format(match))
            f.write("Mismatch: {}\n".format(mismatch))
            f.write("Gap: {}\n".format(gap))
            f.write("Wynik dopasowania: {}\n".format(alignment_score))
            f.write("Długość dopasowania: {}\n".format(alignment_length))
            f.write("Procent identycznych pozycji: {:.2f}%\n".format(identical_positions_percent))
            f.write("Procent przerw: {:.2f}%\n".format(gaps_percent))

        # Zwróć dopasowane sekwencje oraz macierz wsparcia
        return aligned_seq1, aligned_seq2, support_matrix
    def sequence_path(self, support_matrix, filename):
        plt.imshow(support_matrix, cmap='viridis', interpolation='nearest')
        plt.colorbar()

        # Dodanie ścieżki
        for i in range(support_matrix.shape[0]):
            for j in range(support_matrix.shape[1]):
                plt.text(j, i, int(support_matrix[i, j]), ha='center', va='center', color='white')

        plt.xlabel(self.seq1Id)
        plt.ylabel(self.seq2Id)
        plt.title('Macierz z ścieżką')
        plt.savefig("D:\\bioinformatyka\\Lista 1\\app\\Output\\"+filename+".png")
        plt.show()
