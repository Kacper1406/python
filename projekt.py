import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Stałe konfiguracyjne programu ---
FASTA_FILENAME = "sekwencje.txt"
MIN_SEQUENCE_LENGTH = 50
MAX_SEQUENCE_LENGTH = 199

# Liczba jawnie dodanych duplikatów dla celów testowych
NUM_DUPLICATES_TO_ADD = 6
# Liczba jawnie dodanych niepoprawnych sekwencji
NUM_INVALID_TO_ADD = 6

# Prawdopodobieństwo (od 0.0 do 1.0), że "normalna" sekwencja
# również będzie zawierała przypadkowe nieprawidłowe znaki
PROB_INCLUDE_RANDOM_INVALID_CHAR = 0.25

# Zakres losowej liczby sekwencji generowanych w pliku FASTA
MIN_RANDOM_SEQUENCES = 31
MAX_RANDOM_SEQUENCES = 50

# Wymagane minimum sekwencji w pliku dla komunikatu ostrzegawczego
REQUIRED_MIN_SEQUENCES_FOR_MESSAGE = 30

# Liczba sekwencji do wyświetlenia na początku
DISPLAY_INITIAL_SEQUENCES_COUNT = 5
# Liczba sekwencji do wyświetlenia szczegółowych informacji
DISPLAY_USER_INFO_SEQUENCES_COUNT = 3


class DNASequence:
    """
    Reprezentuje pojedynczą sekwencję DNA z jej nazwą i danymi nukleotydowymi.

    Atrybuty:
        name (str): Nazwa sekwencji (nagłówek FASTA).
        sequence (str): Ciąg nukleotydów (ATCG).
    """

    def __init__(self, name: str, sequence: str):
        """
        Inicjalizuje obiekt DNASequence.

        Args:
            name (str): Nazwa sekwencji.
            sequence (str): Ciąg znaków reprezentujących sekwencję DNA.
        """
        self.name = name.strip()
        self.sequence = sequence.strip().upper()

    def __str__(self) -> str:
        """
        Zwraca czytelną dla człowieka reprezentację obiektu.

        Returns:
            str: Sformatowany ciąg zawierający nazwę i sekwencję.
        """
        return f"Nazwa: {self.name}\nSekwencja: {self.sequence}"

    def get_length(self) -> int:
        """
        Oblicza długość sekwencji DNA.

        Returns:
            int: Długość sekwencji w nukleotydach.
        """
        return len(self.sequence)

    def calculate_gc_content(self) -> float:
        """
        Oblicza procentową zawartość par G-C w sekwencji.

        Returns:
            float: Zawartość GC w procentach. Zwraca 0.0 dla pustej sekwencji.
        """
        if not self.sequence:
            return 0.0
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / len(self.sequence)) * 100

    def is_valid(self) -> bool:
        """
        Sprawdza, czy sekwencja zawiera tylko prawidłowe nukleotydy (A, T, C, G).

        Returns:
            bool: True, jeśli sekwencja jest poprawna; False w przeciwnym razie.
        """
        valid_nucleotides = set("ATCG")
        return all(char in valid_nucleotides for char in self.sequence)

    def get_sequence_type(self) -> str:
        """
        Określa typ sekwencji na podstawie zawartości GC.

        Returns:
            str: "GC-rich" (>60% GC), "AT-rich" (<40% GC) lub "Standard".
        """
        gc_content = self.calculate_gc_content()
        if gc_content > 60:
            return "GC-rich"
        elif gc_content < 40:
            return "AT-rich"
        else:
            return "Standard"


def generate_random_dna_sequence(
    length: int, force_invalid_char: bool = False
) -> str:
    """
    Generuje losową sekwencję DNA.

    Args:
        length (int): Pożądana długość sekwencji.
        force_invalid_char (bool): Jeśli True, gwarantuje, że w sekwencji
                                   znajdzie się co najmniej jeden nieprawidłowy
                                   nukleotyd (N, X, Z). Domyślnie False.

    Returns:
        str: Wygenerowana sekwencja DNA.
    """
    nucleotides = ['A', 'T', 'C', 'G']

    if force_invalid_char:
        invalid_chars = ['N', 'X', 'Z']
        if length == 0:
            return ""
        insert_pos = random.randint(0, length - 1)

        seq_list = [random.choice(nucleotides) for _ in range(length)]

        seq_list[insert_pos] = random.choice(invalid_chars)
        return ''.join(seq_list)
    else:
        if random.random() < PROB_INCLUDE_RANDOM_INVALID_CHAR:
            nucleotides_with_errors = nucleotides + ['N', 'X', 'Z']
            return ''.join(random.choice(nucleotides_with_errors)
                           for _ in range(length))
        else:
            return ''.join(random.choice(nucleotides) for _ in range(length))


def create_fasta_file(
    filename: str = FASTA_FILENAME,
    num_total_sequences_to_generate: int = MIN_RANDOM_SEQUENCES
):
    """
    Tworzy lub nadpisuje plik FASTA z losowymi sekwencjami DNA,
    w tym z duplikatami i wadliwymi sekwencjami.

    Args:
        filename (str): Nazwa pliku FASTA do utworzenia.
        num_total_sequences_to_generate (int): Całkowita liczba sekwencji
                                               do wygenerowania w pliku.
    """
    print(f"Tworzę plik '{filename}' z przykładowymi sekwencjami "
          f"(łącznie {num_total_sequences_to_generate} wpisów)...")
    with open(filename, 'w') as f:
        actual_duplicates_to_add = min(NUM_DUPLICATES_TO_ADD,
                                       num_total_sequences_to_generate // 4)
        actual_invalid_to_add = min(NUM_INVALID_TO_ADD,
                                    num_total_sequences_to_generate // 4)

        duplicate_base_seq_data = None
        duplicate_base_seq_name_index = -1

        num_base_sequences = (num_total_sequences_to_generate
                              - actual_duplicates_to_add
                              - actual_invalid_to_add)

        if actual_duplicates_to_add > 0 and num_base_sequences < 1:
            num_base_sequences = 1
        if num_base_sequences < 0:
            num_base_sequences = 0

        if actual_duplicates_to_add > 0 and num_base_sequences > 0:
            duplicate_base_seq_name_index = random.randint(1, num_base_sequences)
            length_for_duplicate_base = random.randint(MIN_SEQUENCE_LENGTH,
                                                       MAX_SEQUENCE_LENGTH)
            duplicate_base_seq_data = generate_random_dna_sequence(
                length_for_duplicate_base, force_invalid_char=False
            )

        base_sequence_counter = 0
        for i in range(num_base_sequences):
            base_sequence_counter += 1
            length = random.randint(MIN_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH)

            sequence_name = f"Sekwencja_{base_sequence_counter}"

            if duplicate_base_seq_name_index == base_sequence_counter:
                sequence_data = duplicate_base_seq_data
            else:
                sequence_data = generate_random_dna_sequence(
                    length, force_invalid_char=False
                )

            f.write(f">{sequence_name}\n{sequence_data}\n")

        if actual_duplicates_to_add > 0:
            if duplicate_base_seq_data is None:
                length_for_duplicate_base = random.randint(MIN_SEQUENCE_LENGTH,
                                                           MAX_SEQUENCE_LENGTH)
                duplicate_base_seq_data = generate_random_dna_sequence(
                    length_for_duplicate_base, force_invalid_char=False
                )

                if num_base_sequences == 0:
                    f.write(f">Sekwencja_Duplikat_Baza\n"
                            f"{duplicate_base_seq_data}\n")

            for i in range(1, actual_duplicates_to_add + 1):
                f.write(f">Duplikat_A_{i}\n{duplicate_base_seq_data}\n")

        for i in range(actual_invalid_to_add):
            length = random.randint(MIN_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH)
            f.write(f">Niepoprawna_{i+1}\n"
                    f"{generate_random_dna_sequence(length, True)}\n")

    print(f"Plik '{filename}' został utworzony pomyślnie.")


def read_fasta_file(filename: str) -> dict[str, DNASequence] | None:
    """
    Odczytuje sekwencje DNA z pliku FASTA.

    Args:
        filename (str): Nazwa pliku FASTA do odczytania.

    Returns:
        dict[str, DNASequence] | None: Słownik zawierający nazwy sekwencji
        jako klucze i obiekty DNASequence jako wartości. Zwraca None w
        przypadku krytycznego błędu (np. plik nie istnieje) lub pusty słownik,
        jeśli plik jest pusty lub nie zawiera poprawnych wpisów.
    """
    sequences = {}
    current_name = None
    current_sequence_lines = []
    line_number = 0

    try:
        if not os.path.exists(filename):
            print(f"BŁĄD: Plik '{filename}' nie został znaleziony przed "
                  "odczytem.")
            return None

        if os.path.getsize(filename) == 0:
            print(f"OSTRZEŻENIE: Plik '{filename}' jest pusty.")
            return {}

        with open(filename, 'r') as f:
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    if current_name is not None:
                        if not current_sequence_lines:
                            print(f"OSTRZEŻENIE: Brak sekwencji dla nagłówka "
                                  f"'{current_name}' przed linią "
                                  f"{line_number}. Pomijam ten wpis.")
                        else:
                            sequences[current_name] = DNASequence(
                                current_name, "".join(current_sequence_lines))
                    current_name = line[1:]
                    if not current_name:
                        print(f"OSTRZEŻENIE: Pusty nagłówek w linii "
                              f"{line_number}. Ten wpis może zostać pominięty "
                              "lub źle zinterpretowany.")
                    current_sequence_lines = []
                else:
                    if current_name is None:
                        print(f"OSTRZEŻENIE: Znaleziono linię sekwencji bez "
                              "poprzedzającego nagłówka w linii "
                              f"{line_number}: '{line}'. Pomijam.")
                        continue
                    current_sequence_lines.append(line)

            if current_name is not None:
                if not current_sequence_lines:
                    print(f"OSTRZEŻENIE: Brak sekwencji dla ostatniego "
                          f"nagłówka '{current_name}'. Pomijam ten wpis.")
                else:
                    sequences[current_name] = DNASequence(
                        current_name, "".join(current_sequence_lines))
            elif not sequences:
                print(f"OSTRZEŻENIE: Plik '{filename}' nie zawiera żadnych "
                      "poprawnych wpisów FASTA.")

    except FileNotFoundError:
        print(f"BŁĄD: Plik '{filename}' nie został znaleziony. "
              "Upewnij się, że jest w tym samym katalogu co skrypt.")
        return None
    except Exception as e:
        print(f"Wystąpił błąd podczas odczytu pliku '{filename}': {e}")
        return None
    return sequences


def add_user_sequence(
    sequences_dict: dict[str, DNASequence], filename: str = FASTA_FILENAME
) -> dict[str, DNASequence]:
    """
    Prosi użytkownika o podanie własnej sekwencji DNA i dodaje ją
    do pamięci programu oraz do pliku FASTA.

    Args:
        sequences_dict (dict): Słownik istniejących sekwencji DNA.
        filename (str): Nazwa pliku FASTA, do którego ma być dopisana sekwencja.

    Returns:
        dict[str, DNASequence]: Zaktualizowany słownik sekwencji DNA.
    """
    print("\n--- Dodaj swoją sekwencję ---")
    name = input("Podaj nazwę dla Twojej sekwencji (np. Moja_Sekwencja): ").strip()
    if not name:
        print("Nazwa sekwencji nie może być pusta. Anulowano dodawanie "
              "sekwencji.")
        return sequences_dict

    sequence_data = input("Wklej sekwencję DNA (tylko A, T, C, G): ").strip().upper()

    if not sequence_data:
        print("Sekwencja nie może być pusta. Anulowano dodawanie sekwencji.")
        return sequences_dict

    temp_seq_obj = DNASequence(name, sequence_data)
    if not temp_seq_obj.is_valid():
        print("UWAGA: Twoja sekwencja zawiera nieprawidłowe znaki. Zostanie "
              "oznaczona jako niepoprawna.")

    if name in sequences_dict:
        print(f"UWAGA: Sekwencja o nazwie '{name}' już istnieje w pamięci. "
              "Zostanie nadpisana.")

    sequences_dict[name] = temp_seq_obj
    print(f"Sekwencja '{name}' została dodana do pamięci programu.")

    try:
        with open(filename, 'a') as f:
            f.write(f">{name}\n{sequence_data}\n")
        print(f"Sekwencja '{name}' została trwale dopisana do pliku "
              f"'{filename}'.")
    except IOError as e:
        print(f"BŁĄD: Nie udało się zapisać sekwencji do pliku '{filename}': "
              f"{e}")
    except Exception as e:
        print(f"Wystąpił nieoczekiwany błąd podczas zapisu do pliku: {e}")

    return sequences_dict


def display_sequence_info(seq_object: DNASequence):
    """
    Wyświetla szczegółowe informacje o pojedynczej sekwencji DNA.

    Args:
        seq_object (DNASequence): Obiekt sekwencji DNA do wyświetlenia.
    """
    print(f"\n--- Informacje o sekwencji: {seq_object.name} ---")
    print(f"  Długość: {seq_object.get_length()} nukleotydów")
    print(f"  Zawartość GC: {seq_object.calculate_gc_content():.2f}%")
    print(f"  Poprawna: {'Tak' if seq_object.is_valid() else 'Nie (zawiera niepoprawne nukleotydy)'}")
    print(f"  Typ: {seq_object.get_sequence_type()}")
    print(f"  Sekwencja (pierwsze 30 znaków): {seq_object.sequence[:30]}...")


def process_sequences(
    sequences_dict: dict[str, DNASequence]
) -> pd.DataFrame:
    """
    Przetwarza sekwencje DNA, usuwając wadliwe i zduplikowane wpisy,
    a następnie konwertuje pozostałe dane do DataFrame.

    Args:
        sequences_dict (dict): Słownik sekwencji DNA do przetworzenia.

    Returns:
        pd.DataFrame: DataFrame zawierający oczyszczone i przetworzone
                      sekwencje. Pusty DataFrame, jeśli nie ma danych
                      do przetworzenia.
    """
    if not sequences_dict:
        print("Brak sekwencji do przetworzenia.")
        return pd.DataFrame()

    print("\n--- Przetwarzanie sekwencji ---")
    processed_data = []
    seen_sequences_content = set()
    initial_count = len(sequences_dict)
    removed_duplicates = 0
    removed_invalid = 0

    for name, seq_obj in sequences_dict.items():
        if not seq_obj.is_valid():
            print(f"  Usuwam wadliwą sekwencję: {name} "
                  "(zawiera niepoprawne nukleotydy)")
            removed_invalid += 1
            continue

        if seq_obj.sequence in seen_sequences_content:
            print(f"  Usuwam duplikat: {name} (treść sekwencji już istnieje)")
            removed_duplicates += 1
            continue

        seen_sequences_content.add(seq_obj.sequence)
        processed_data.append({
            "Nazwa": seq_obj.name,
            "Sekwencja": seq_obj.sequence,
            "Dlugosc": seq_obj.get_length(),
            "GC_Zawartosc": seq_obj.calculate_gc_content(),
            "Typ_Sekwencji": seq_obj.get_sequence_type()
        })

    print(f"  Usunięto duplikatów: {removed_duplicates}")
    print(f"  Usunięto wadliwych sekwencji: {removed_invalid}")
    print(f"  Pozostało sekwencji po oczyszczeniu: {len(processed_data)} "
          f"z {initial_count} początkowych.")

    data_frame = pd.DataFrame(processed_data)
    return data_frame


def visualize_data(dataframe: pd.DataFrame):
    """
    Generuje i wyświetla różne wizualizacje danych z DataFrame.

    Args:
        dataframe (pd.DataFrame): DataFrame zawierający przetworzone dane
                                  sekwencji DNA.
    """
    if dataframe.empty:
        print("Brak danych do wizualizacji.")
        return

    print("\n--- Tworzenie wizualizacji ---")

    # Wykres 1: Rozkład długości sekwencji
    plt.figure(figsize=(10, 6))
    plt.hist(dataframe['Dlugosc'], bins=10, edgecolor='black')
    plt.title('Rozkład długości sekwencji DNA')
    plt.xlabel('Długość sekwencji (nukleotydy)')
    plt.ylabel('Liczba sekwencji')
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()

    # Wykres 2: Rozkład zawartości GC
    plt.figure(figsize=(10, 6))
    plt.hist(dataframe['GC_Zawartosc'], bins=10,
             edgecolor='black', color='lightgreen')
    plt.title('Rozkład zawartości GC w sekwencjach DNA')
    plt.xlabel('Zawartość GC (%)')
    plt.ylabel('Liczba sekwencji')
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()

    # Wykres 3: Długość sekwencji vs. Zawartość GC (scatter plot)
    plt.figure(figsize=(10, 6))
    plt.scatter(dataframe['Dlugosc'], dataframe['GC_Zawartosc'],
                alpha=0.7, color='skyblue')
    plt.title('Długość sekwencji vs. Zawartość GC')
    plt.xlabel('Długość sekwencji (nukleotydy)')
    plt.ylabel('Zawartość GC (%)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Wykres 4: Liczba sekwencji według typu (bar plot)
    plt.figure(figsize=(8, 5))
    dataframe['Typ_Sekwencji'].value_counts().plot(
        kind='bar', color=['purple', 'orange', 'cyan'])
    plt.title('Liczba sekwencji według typu')
    plt.xlabel('Typ sekwencji')
    plt.ylabel('Liczba sekwencji')
    plt.xticks(rotation=45)
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()


def main():
    """
    Główna funkcja programu, orkiestrująca całym procesem:
    generowaniem pliku, odczytem, przetwarzaniem, wizualizacją
    i interakcją z użytkownikiem.
    """
    target_num_sequences_for_file = random.randint(MIN_RANDOM_SEQUENCES,
                                                   MAX_RANDOM_SEQUENCES)

    if not os.path.exists(FASTA_FILENAME):
        create_fasta_file(FASTA_FILENAME,
                          num_total_sequences_to_generate=
                          target_num_sequences_for_file)
        print(f"Plik '{FASTA_FILENAME}' nie istniał, został utworzony "
              "z losową liczbą sekwencji "
              f"({target_num_sequences_for_file}).")

    while True:
        all_sequences_dict = read_fasta_file(FASTA_FILENAME)

        if all_sequences_dict is None:
            print("Nie można kontynuować z powodu krytycznego błędu odczytu "
                  "pliku. Zakończono program.")
            return

        num_entries_in_file = len(all_sequences_dict)

        if (num_entries_in_file < REQUIRED_MIN_SEQUENCES_FOR_MESSAGE
                or num_entries_in_file > MAX_RANDOM_SEQUENCES):
            print(f"\nOSTRZEŻENIE: Plik '{FASTA_FILENAME}' zawiera "
                  f"{num_entries_in_file} sekwencji. Wymagane minimum to "
                  f"{REQUIRED_MIN_SEQUENCES_FOR_MESSAGE} sekwencji.")
        else:
            print(f"\nPlik '{FASTA_FILENAME}' zawiera {num_entries_in_file} "
                  "sekwencji. Spełnia wymagane minimum "
                  f"{REQUIRED_MIN_SEQUENCES_FOR_MESSAGE} sekwencji.")

        user_choice = input("Czy chcesz nadpisać obecny plik i wygenerować "
                            "nowe sekwencje z losową liczbą? (tak/nie): ").lower()
        if user_choice == 'tak':
            target_num_sequences_for_file = random.randint(
                MIN_RANDOM_SEQUENCES, MAX_RANDOM_SEQUENCES)
            create_fasta_file(FASTA_FILENAME,
                              num_total_sequences_to_generate=
                              target_num_sequences_for_file)
            continue
        else:
            print(f"Kontynuuję pracę z obecnym stanem pliku "
                  f"({num_entries_in_file} sekwencji).")
            break

    print("\n--- Odczytano sekwencje z pliku ---")
    if all_sequences_dict:
        for name, seq_obj in \
                list(all_sequences_dict.items())[:DISPLAY_INITIAL_SEQUENCES_COUNT]:
            print(f"  {name}: {seq_obj.get_length()} nt, "
                  f"GC: {seq_obj.calculate_gc_content():.2f}%")
        if len(all_sequences_dict) > DISPLAY_INITIAL_SEQUENCES_COUNT:
            print(f"  ...i jeszcze "
                  f"{len(all_sequences_dict) - DISPLAY_INITIAL_SEQUENCES_COUNT} "
                  "innych sekwencji.")
    else:
        print("Brak sekwencji do wyświetlenia z pliku.")

    add_more = input("\nCzy chcesz dodać własną sekwencję? (tak/nie): ").lower()
    if add_more == 'tak':
        all_sequences_dict = add_user_sequence(all_sequences_dict,
                                               FASTA_FILENAME)

    print("\n--- Wyświetlam informacje o kilku sekwencjach (przed czyszczeniem) ---")
    if all_sequences_dict:
        count = 0
        for name, seq_obj in all_sequences_dict.items():
            display_sequence_info(seq_obj)
            count += 1
            if count >= DISPLAY_USER_INFO_SEQUENCES_COUNT:
                break
    else:
        print("Brak sekwencji do wyświetlenia informacji.")

    cleaned_df = process_sequences(all_sequences_dict)

    # Wyświetlenie końcowego DataFrame
    if not cleaned_df.empty:
        print("\n--- Ostateczny DataFrame z czystymi sekwencjami ---")
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        print(cleaned_df)
        pd.reset_option('display.max_rows')
        pd.reset_option('display.max_columns')
        pd.reset_option('display.width')
    else:
        print("DataFrame jest pusty po czyszczeniu.")

    visualize_data(cleaned_df)

    print("\nProgram zakończył działanie. Miłego dnia!")


if __name__ == "__main__":
    main()
