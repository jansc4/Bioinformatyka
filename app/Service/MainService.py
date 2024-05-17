import requests
import re
from app.Classes.Sequence import Sequence
from app.UI.View import display_error


def seq_constructor(fasta_data):
    """
    Function construct Sequence object from fasta data
    :type fasta_data: string
    """
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

    return Sequence(seq, id, name)


def get_uniprot_sequence(uniprot_id):
    """
    Function get uniprot sequence from uniprot id
    :param uniprot_id:
    :return: Sequence object
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Sprawdź, czy odpowiedź nie zawiera błędu HTTP
        fasta_data = response.text
        seq = seq_constructor(fasta_data)
        return seq
    except requests.exceptions.HTTPError as e:
        print(f"Wystąpił błąd HTTP: {e}")
        return None
    except Exception as e:
        print(f"Wystąpił inny błąd: {e}")
        return None


def read_sequence_from_file(name):
    try:
        with open(name, "r") as f:
            read = f.read()
            seq = seq_constructor(read)
            return seq
    except IOError as e:
        display_error(f"Wystąpił wyjątek: {e}")
        return None
