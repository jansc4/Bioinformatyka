def start_menu():
    print("Wybierz metodę wprowadzania sekwencji: plik FASTA(1) lub ID Uniprot(2)")
    return input()


def display_error(message):
    print(message)


def display_sequence_input_prompt():
    print("Podaj nazwę pliku")
    return input()


def display_sequence_retrieval_prompt():
    print("Podaj ID Uniprot")
    return input()
