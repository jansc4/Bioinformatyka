# -*- coding: utf-8 -*-
from app.Classes.Alignment import Alignment
from app.Classes.Sequence import Sequence
from app.Classes.SequenceAlignment import SequenceAlignment
import numpy as np



class MultipleSequenceAlignment:
    def __init__(self, match, mismatch, gap, *sequences):
        self.sequences = sequences
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.high = self.width = len(sequences)
        self.score_matrix = np.zeros((self.high, self.width + 1))
        self.alignments = []

    def _sum_alignment_scores(self):
        for i in range(self.high):
            sum = 0
            for j in range(self.width):
                sum += self.score_matrix[i, j]
            self.score_matrix[i, -1] = sum

    def _msa_sequences(self):

        for i, seq1 in enumerate(self.sequences):

            for j, seq2 in enumerate(self.sequences):
                if i != j:  # Unikamy porównywania sekwencji z samymi sobą

                    if i < j:  # Unikamy powtarzających się porównań, ponieważ score_matrix jest symetryczna
                        aligned_seq1, aligned_seq2, support_matrix, wynik = SequenceAlignment(seq1,
                                                                                              seq2).needleman_wunsch(
                            self.match, self.mismatch, self.gap, False)
                        # Obliczamy punktację końcową dopasowania
                        alignment_score = wynik[-1, -1]  # Ostatni element macierzy to końcowy wynik dopasowania
                        # Zapisujemy wynik punktacji w odpowiedniej komórce macierzy wynikowej
                        self.score_matrix[i, j] = alignment_score
                        # Ponieważ score_matrix jest symetryczna, zapisujemy także wynik w odpowiednim miejscu symetrycznym
                        self.score_matrix[j, i] = alignment_score
                        temp = []
                        temp.append(aligned_seq1)
                        temp.append(aligned_seq2)
                        self.alignments.append(temp)

        self._sum_alignment_scores()
        #print(self.score_matrix)
        return self.alignments, self.score_matrix

    def _msa_center_seq(self):
        alignments, score_matrix = self._msa_sequences()
        score = []
        for i in range(self.score_matrix.shape[0]):
            score.append(self.score_matrix[i, -1])

        max_score_index = np.argmax(score)  # Znajdź indeks sekwencji z najwyższym wynikiem
        max_score_sequence_id = self.sequences[max_score_index].get_id()  # ID sekwencji z najwyższym wynikiem
        return max_score_sequence_id

    def msa(self, alignment_id=0):
        max_score_sequence_id = self._msa_center_seq()
        msa_alignment = Alignment(max_score_sequence_id, alignment_id)

        for al in self.alignments:
            if al[0].get_id() == max_score_sequence_id or al[1].get_id() == max_score_sequence_id:
                for i in al:
                    if not msa_alignment.alignment_check(i):
                        msa_alignment.add_alignment(i)

        print(msa_alignment.alignment())
        print(msa_alignment)

        # for al in self.alignments:
        #     if al[0].get_id() or al[1].get_id() == max_score_sequence_id:
