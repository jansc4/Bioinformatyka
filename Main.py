import requests
from Sequence import Sequence
from Dotplot import Dotplot
import re


def seqConstructor(fasta_data):
    lines = fasta_data.split('\n')

    # Sprawdź, czy sekwencja zawiera nagłówek
    if not lines[0].startswith('>'):
        print("Błąd: Brak nagłówka w sekwencji FASTA.")
        return None

    seq = ''.join(lines[1:])

    # Sprawdź, czy sekwencja zawiera nielegalne znaki
    basic_amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                         'Y'}

    if not set(seq).issubset(basic_amino_acids):
        print("Błąd: Sekwencja zawiera nielegalne znaki.")
        return None

    match = re.match(r">sp\|(\w+)\|(\w+)OS", lines[0])
    id = None
    name = None
    if match:
        id = match.group(1)
        name = match.group(2)

    protein = Sequence(seq, id, name)
    return protein


def get_uniprot_sequence(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Sprawdź, czy odpowiedź nie zawiera błędu HTTP
        fasta_data = response.text
        seq = seqConstructor(fasta_data)
        return seq
    except requests.exceptions.HTTPError as e:
        print(f"Wystąpił błąd HTTP: {e}")
        return None
    except Exception as e:
        print(f"Wystąpił inny błąd: {e}")
        return None


def start():
    while True:
        print("Wybierz metodę wprowadzania sekwencji: plik FASTA(1) lub ID Uniprot(2)")
        choice = input()
        if choice == "1":
            print("Podaj nazwę pliku")
            name = input()
            seq = ""  # Inicjalizacja zmiennej seq
            try:
                with open(name, "r") as f:  # Użyj 'with' do otwarcia pliku, aby uniknąć zapominania o zamknięciu pliku
                    read = f.read()
                    seq = seqConstructor(read)
            except IOError as e:
                print("Wystąpił wyjątek:", e)
            return seq
        elif choice == "2":
            print("Podaj ID Uniprot")
            uniprot_id = input()
            seq = get_uniprot_sequence(uniprot_id)
            if seq:
                return seq
            else:
                print("Nie udało się pobrać sekwencji dla podanego ID Uniprot.")
        else:
            print("Błędna opcja")

def filter(dotplot, threshold, window):


    return dotplot

print("Wybierz pierwszą sekwencję:")
#seq1 = start()
print("Wybierz drugą sekwencję:")
#seq2 = start()
seq1 = Sequence("MVHLTDAEKSAVSCLWAKVNPDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPKVKAHGKKVITAFNEGLKNLDNLKGTFASLSELHCDKLHVDPENFRLLGNAIVIVLGHHLGKDFTPAAQAAFQKVVAGVATALAHKYH", "human", "P02089")
seq2 = Sequence("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", "mouse", "P68871")
dp = Dotplot.from_sequences(seq1, seq2)

#print(seq1.name)
dp.graphic("test")
#dp.saveTxt("dotplot.txt")