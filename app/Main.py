from app.UI.View import start_menu, display_sequence_input_prompt, display_sequence_retrieval_prompt, display_error
from app.Service.MainService import seq_constructor, get_uniprot_sequence, read_sequence_from_file
from app.Classes.MultipleSequenceAlignment import MultipleSequenceAlignment
from app.Classes.Sequence import Sequence

def start():
    while True:
        choice = start_menu()
        if choice == "1":
            name = display_sequence_input_prompt()
            seq = read_sequence_from_file(name)
            if seq:
                return seq
        elif choice == "2":
            uniprot_id = display_sequence_retrieval_prompt()
            seq = get_uniprot_sequence(uniprot_id)
            if seq:
                return seq
        else:
            display_error("Błędna opcja")


def main():
    from app.Classes.Sequence import Sequence
    from app.Classes.SequenceAlignment import SequenceAlignment

    # Tworzymy obiekty sekwencji
    seq1 = Sequence("CATWALK", "seq1_id")
    seq2 = Sequence("CQWARD", "seq2_id")
    seq3 = Sequence("CADARDC", "seq3_id")

    # Tworzymy obiekt klasy SequenceAlignment
    #alignment = SequenceAlignment(seq1, seq2)

    # Parametry do algorytmu
    match_score = 2
    mismatch_score = -1
    gap_penalty = -2

    # Przeprowadzamy lokalne dopasowanie za pomocą algorytmu Smitha-Watermana
    #aligned_seq1, aligned_seq2, support_matrix = alignment.smith_waterman(match_score, mismatch_score, gap_penalty)

    # Wyświetlamy rezultaty
    #print("Aligned Sequence 1:", aligned_seq1)
    #print("Aligned Sequence 2:", aligned_seq2)

    # Wyświetlamy graficzną reprezentację macierzy z ścieżką
    #alignment.sequence_path(support_matrix, "smith_path")
    #print("#######################")
    # Przeprowadzamy globalne dopasowanie za pomocą algorytmu Needleman-Wunsch
    #aligned_seq1_1, aligned_seq2_1, support_matrix_1 = alignment.needleman_wunsch(match_score, mismatch_score, gap_penalty)

    # Wyświetlamy rezultaty
    #print("Aligned Sequence 1:", aligned_seq1_1)
    #print("Aligned Sequence 2:", aligned_seq2_1)

    # Wyświetlamy graficzną reprezentację macierzy z ścieżką
    #alignment.sequence_path(support_matrix_1, "neddleman_path")
    alignments, score_matrix = MultipleSequenceAlignment(match_score, mismatch_score, gap_penalty,seq1, seq2, seq3)._msa_center_seq()

    for a in alignments:
        print("[")
        for i in a:
            print(i)
        print("]")
if __name__ == "__main__":
    main()
