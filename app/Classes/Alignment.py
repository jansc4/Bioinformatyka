class Alignment:
    def __init__(self, c_seq_id, id=0, *alignments):
        self.id = id
        self.central_seq_id = c_seq_id
        self.alignments = list(alignments)
        self.output_alignments = []

    @staticmethod
    def insert_char_at(alignment, index, char):
        if index < 0 or index > len(alignment):
            raise ValueError("Index is out of bounds")
        return alignment[:index +1] + char + alignment[index+1:]

    def __str__(self):
        alignment_str = f"ID: {self.id}\nCentral sequence ID: {self.central_seq_id}\nAlignments:\n"
        for alignment in self.alignments:
            alignment_str += str(alignment) + "\n"
        return alignment_str

    def get_id(self):
        return self.id

    def get_alignments(self):
        return self.alignments

    def set_id(self, id):
        self.id = id

    def set_alignments(self, alignments):
        self.alignments = list(alignments)

    def add_alignment(self, alignment):
        self.alignments.append(alignment)

    def del_alignment(self, alignment_id):
        for al in self.alignments:
            if al.get_id() == alignment_id:
                self.alignments.remove(al)
                return True
        return False

    def alignment_check(self, alignment):
        for al in self.alignments:
            if (al.get_id() == alignment.get_id() and
                    al.get_name() == alignment.get_name() and
                    al.get_sequence() == alignment.get_sequence()):
                return True
        return False

    @staticmethod
    def find_all_indices(string, char):
        indices = []
        index = string.find(char)
        while index != -1:
            indices.append(index)
            index = string.find(char, index + 1)
        return indices
    def alignment(self):
        cent = None
        for al in self.alignments:
            if al.get_id() == self.central_seq_id and cent is None:
                cent = al.get_sequence()
            elif al.get_id() == self.central_seq_id:
                temp = al.get_sequence()
                indices = self.find_all_indices(temp, '_')
                for ind in indices:
                    cent = self.insert_char_at(cent, ind, '_')
                    self.output_alignments.append(cent)
            if al.get_id() != self.central_seq_id:
                temp = al.get_sequence()
                indices = self.find_all_indices(cent, '_')
                for ind in indices:
                    temp = self.insert_char_at(cent, ind, '_')
                    self.output_alignments.append(temp)
        return self.output_alignments
