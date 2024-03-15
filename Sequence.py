import numpy as np


class Sequence:
    def __init__(self, seq, name="default", id=000000):
        self.sequence = seq
        self.length = len(seq)
        self.name = name
        self.id = id

    def __str__(self):
        return self.sequence

    def get_length(self):
        return self.length

    def dotplot(self, seq2):
        dp = np.zeros((self.length, seq2.get_length()))
        for i in range(self.length):
            for j in range(seq2.get_length()):
                if self.sequence[i] == seq2.sequence[j]:
                    dp[i, j] = 1
                else:
                    dp[i, j] = 0
        return dp


# %%
import requests