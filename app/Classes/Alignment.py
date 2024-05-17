class Alignment:
    def __init__(self, id = 00, *alignments):
        self.id = id
        self.alignments = alignments

    def __str__(self):
        for alignment in self.alignments:
            print(alignment)

    def get_id(self):
        return self.id

    def get_alignments(self):
        return self.alignments

    def set_id(self, id):
        self.id = id

    def set_alignments(self, alignments):
        self.alignments = alignments