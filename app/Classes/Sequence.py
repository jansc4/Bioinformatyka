"""
Sequence Class contains DNA sequence and associated metadata
"""
class Sequence:
    def __init__(self, seq, name="default", id=000000):
        self.sequence = seq
        self.length = len(seq)
        self.name = name
        self.id = id

    def __str__(self):
        return "Name="+ self.name+" Id="+ str(self.id) + " Sequence="+ self.sequence

    def get_length(self):
        return self.length
    def get_name(self):
        return self.name
    def get_id(self):
        return self.id
    def get_sequence(self):
        return self.sequence