import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

NAZWA_PLIKU_FASTA = "sekwencje.txt"
MIN_DLUGOSC_SEKWENCJI = 50
MAX_DLUGOSC_SEKWENCJI = 199

LICZBA_DUPLIKATOW_DO_DODANIA = 6
LICZBA_NIEPOPRAWNYCH_DO_DODANIA = 6

PRAWDOPODOBIENSTWO_LOSOWYCH_NIEPRAWIDLOWYCH_ZNAKOW = 0.25

MIN_LOSOWYCH_SEKWENCJI = 31
MAX_LOSOWYCH_SEKWENCJI = 50

WYMAGANE_MIN_SEKWENCJI_DLA_OSTRZEZENIA = 30

LICZBA_WYSWIETLANYCH_POCZATKOWYCH_SEKWENCJI = 5
LICZBA_WYSWIETLANYCH_INFO_SEKWENCJI_OD_UZYTKOWNIKA = 3


class SekwencjaDNA:
    """
    Reprezentuje pojedynczą sekwencję DNA z jej nazwą i danymi nukleotydowymi.

    Atrybuty:
        nazwa (str): Nazwa sekwencji (nagłówek FASTA).
        sekwencja (str): Ciąg nukleotydów (ATCG).
    """

    def __init__(self, nazwa: str, sekwencja: str):
        """
        Inicjalizuje obiekt SekwencjaDNA.

        Args:
            nazwa (str): Nazwa sekwencji.
            sekwencja (str): Ciąg znaków reprezentujących sekwencję DNA.
        """
        self.nazwa = nazwa.strip()
        self.sekwencja = sekwencja.strip().upper()

    def __str__(self) -> str:
        """
        Zwraca czytelną dla człowieka reprezentację obiektu.

        Returns:
            str: Sformatowany ciąg zawierający nazwę i sekwencję.
        """
        return f"Nazwa: {self.nazwa}\nSekwencja: {self.sekwencja}"

    def pobierz_dlugosc(self) -> int:
        """
        Oblicza długość sekwencji DNA.

        Returns:
            int: Długość sekwencji w nukleotydach.
        """
        return len(self.sekwencja)

    def oblicz_zawartosc_gc(self) -> float:
        """
        Oblicza procentową zawartość par G-C w sekwencji.

        Returns:
            float: Zawartość GC w procentach. Zwraca 0.0 dla pustej sekwencji.
        """
        if not self.sekwencja:
            return 0.0
        ilosc_gc = self.sekwencja.count('G') + self.sekwencja.count('C')
        return (ilosc_gc / len(self.sekwencja)) * 100

    def jest_poprawna(self) -> bool:
        """
        Sprawdza, czy sekwencja zawiera tylko prawidłowe nukleotydy (A, T, C, G).

        Returns:
            bool: True, jeśli sekwencja jest poprawna; False w przeciwnym razie.
        """
        prawidlowe_nukleotydy = set("ATCG")
        return all(char in prawidlowe_nukleotydy for char in self.sekwencja)

    def pobierz_typ_sekwencji(self) -> str:
        """
        Określa typ sekwencji na podstawie zawartości GC.

        Returns:
            str: "Bogata_w_GC" (>60% GC), "Bogata_w_AT" (<40% GC) lub "Standardowa".
        """
        zawartosc_gc = self.oblicz_zawartosc_gc()
        if zawartosc_gc > 60:
            return "Bogata_w_GC"
        elif zawartosc_gc < 40:
            return "Bogata_w_AT"
        else:
            return "Standardowa"


def generuj_losowa_sekwencje_dna(
    dlugosc: int, wymus_niepoprawny_znak: bool = False
) -> str:
    """
    Generuje losową sekwencję DNA.

    Args:
        dlugosc (int): Pożądana długość sekwencji.
        wymus_niepoprawny_znak (bool): Jeśli True, gwarantuje, że w sekwencji
                                       znajdzie się co najmniej jeden nieprawidłowy
                                       nukleotyd (N, X, Z). Domyślnie False.

    Returns:
        str: Wygenerowana sekwencja DNA.
    """
    nukleotydy = ['A', 'T', 'C', 'G']

    if wymus_niepoprawny_znak:
        niepoprawne_znaki = ['N', 'X', 'Z']
        if dlugosc == 0:
            return ""
        pozycja_wstawienia = random.randint(0, dlugosc - 1)

        lista_sekwencji = [random.choice(nukleotydy) for _ in range(dlugosc)]

        lista_sekwencji[pozycja_wstawienia] = random.choice(niepoprawne_znaki)
        return ''.join(lista_sekwencji)
    else:
        if random.random() < PRAWDOPODOBIENSTWO_LOSOWYCH_NIEPRAWIDLOWYCH_ZNAKOW:
            nukleotydy_z_bledami = nukleotydy + ['N', 'X', 'Z']
            return ''.join(random.choice(nukleotydy_z_bledami)
                           for _ in range(dlugosc))
        else:
            return ''.join(random.choice(nukleotydy) for _ in range(dlugosc))


def utworz_plik_fasta(
    nazwa_pliku: str = NAZWA_PLIKU_FASTA,
    liczba_wszystkich_sekwencji_do_wygenerowania: int = MIN_LOSOWYCH_SEKWENCJI
):
    """
    Tworzy lub nadpisuje plik FASTA z losowymi sekwencjami DNA,
    w tym z duplikatami i wadliwymi sekwencjami.

    Args:
        nazwa_pliku (str): Nazwa pliku FASTA do utworzenia.
        liczba_wszystkich_sekwencji_do_wygenerowania (int): Całkowita liczba sekwencji
                                                           do wygenerowania w pliku.
    """
    print(f"Tworzę plik '{nazwa_pliku}' z przykładowymi sekwencjami "
          f"(łącznie {liczba_wszystkich_sekwencji_do_wygenerowania} wpisów)...")
    with open(nazwa_pliku, 'w') as f:
        faktyczna_liczba_duplikatow_do_dodania = min(LICZBA_DUPLIKATOW_DO_DODANIA,
                                                     liczba_wszystkich_sekwencji_do_wygenerowania // 4)
        faktyczna_liczba_niepoprawnych_do_dodania = min(LICZBA_NIEPOPRAWNYCH_DO_DODANIA,
                                                        liczba_wszystkich_sekwencji_do_wygenerowania // 4)

        dane_bazowej_sekwencji_duplikatu = None
        indeks_nazwy_bazowej_sekwencji_duplikatu = -1

        liczba_sekwencji_bazowych = (liczba_wszystkich_sekwencji_do_wygenerowania
                                     - faktyczna_liczba_duplikatow_do_dodania
                                     - faktyczna_liczba_niepoprawnych_do_dodania)

        if faktyczna_liczba_duplikatow_do_dodania > 0 and liczba_sekwencji_bazowych < 1:
            liczba_sekwencji_bazowych = 1
        if liczba_sekwencji_bazowych < 0:
            liczba_sekwencji_bazowych = 0

        if faktyczna_liczba_duplikatow_do_dodania > 0 and liczba_sekwencji_bazowych > 0:
            indeks_nazwy_bazowej_sekwencji_duplikatu = random.randint(1, liczba_sekwencji_bazowych)
            dlugosc_dla_bazowej_sekwencji_duplikatu = random.randint(MIN_DLUGOSC_SEKWENCJI,
                                                                      MAX_DLUGOSC_SEKWENCJI)
            dane_bazowej_sekwencji_duplikatu = generuj_losowa_sekwencje_dna(
                dlugosc_dla_bazowej_sekwencji_duplikatu, wymus_niepoprawny_znak=False
            )

        licznik_sekwencji_bazowych = 0
        for i in range(liczba_sekwencji_bazowych):
            licznik_sekwencji_bazowych += 1
            dlugosc = random.randint(MIN_DLUGOSC_SEKWENCJI, MAX_DLUGOSC_SEKWENCJI)

            nazwa_sekwencji = f"Sekwencja_{licznik_sekwencji_bazowych}"

            if indeks_nazwy_bazowej_sekwencji_duplikatu == licznik_sekwencji_bazowych:
                dane_sekwencji = dane_bazowej_sekwencji_duplikatu
            else:
                dane_sekwencji = generuj_losowa_sekwencje_dna(
                    dlugosc, wymus_niepoprawny_znak=False
                )

            f.write(f">{nazwa_sekwencji}\n{dane_sekwencji}\n")

        if faktyczna_liczba_duplikatow_do_dodania > 0:
            if dane_bazowej_sekwencji_duplikatu is None:
                dlugosc_dla_bazowej_sekwencji_duplikatu = random.randint(MIN_DLUGOSC_SEKWENCJI,
                                                                          MAX_DLUGOSC_SEKWENCJI)
                dane_bazowej_sekwencji_duplikatu = generuj_losowa_sekwencje_dna(
                    dlugosc_dla_bazowej_sekwencji_duplikatu, wymus_niepoprawny_znak=False
                )
                if liczba_sekwencji_bazowych == 0:
                    f.write(f">Sekwencja_Duplikat_Baza\n"
                            f"{dane_bazowej_sekwencji_duplikatu}\n")

            for i in range(1, faktyczna_liczba_duplikatow_do_dodania + 1):
                f.write(f">Duplikat_A_{i}\n{dane_bazowej_sekwencji_duplikatu}\n")

        for i in range(faktyczna_liczba_niepoprawnych_do_dodania):
            dlugosc = random.randint(MIN_DLUGOSC_SEKWENCJI, MAX_DLUGOSC_SEKWENCJI)
            f.write(f">Niepoprawna_{i+1}\n"
                    f"{generuj_losowa_sekwencje_dna(dlugosc, True)}\n")

    print(f"Plik '{nazwa_pliku}' został utworzony pomyślnie.")


def wczytaj_plik_fasta(nazwa_pliku: str) -> dict[str, SekwencjaDNA] | None:
    """
    Odczytuje sekwencje DNA z pliku FASTA.

    Args:
        nazwa_pliku (str): Nazwa pliku FASTA do odczytania.

    Returns:
        dict[str, SekwencjaDNA] | None: Słownik zawierający nazwy sekwencji
        jako klucze i obiekty SekwencjaDNA jako wartości. Zwraca None w
        przypadku krytycznego błędu (np. plik nie istnieje) lub pusty słownik,
        jeśli plik jest pusty lub nie zawiera poprawnych wpisów.
    """
    sekwencje = {}
    aktualna_nazwa = None
    linie_aktualnej_sekwencji = []
    numer_linii = 0

    try:
        if not os.path.exists(nazwa_pliku):
            print(f"BŁĄD: Plik '{nazwa_pliku}' nie został znaleziony przed "
                  "odczytem.")
            return None

        if os.path.getsize(nazwa_pliku) == 0:
            print(f"OSTRZEŻENIE: Plik '{nazwa_pliku}' jest pusty.")
            return {}

        with open(nazwa_pliku, 'r') as f:
            for linia in f:
                numer_linii += 1
                linia = linia.strip()
                if not linia:
                    continue

                if linia.startswith('>'):
                    if aktualna_nazwa is not None:
                        if not linie_aktualnej_sekwencji:
                            print(f"OSTRZEŻENIE: Brak sekwencji dla nagłówka "
                                  f"'{aktualna_nazwa}' przed linią "
                                  f"{numer_linii}. Pomijam ten wpis.")
                        else:
                            sekwencje[aktualna_nazwa] = SekwencjaDNA(
                                aktualna_nazwa, "".join(linie_aktualnej_sekwencji))
                    aktualna_nazwa = linia[1:]
                    if not aktualna_nazwa:
                        print(f"OSTRZEŻENIE: Pusty nagłówek w linii "
                              f"{numer_linii}. Ten wpis może zostać pominięty "
                              "lub źle zinterpretowany.")
                    linie_aktualnej_sekwencji = []
                else:
                    if aktualna_nazwa is None:
                        print(f"OSTRZEŻENIE: Znaleziono linię sekwencji bez "
                              "poprzedzającego nagłówka w linii "
                              f"{numer_linii}: '{linia}'. Pomijam.")
                        continue
                    linie_aktualnej_sekwencji.append(linia)

            if aktualna_nazwa is not None:
                if not linie_aktualnej_sekwencji:
                    print(f"OSTRZEŻENIE: Brak sekwencji dla ostatniego "
                          f"nagłówka '{aktualna_nazwa}'. Pomijam ten wpis.")
                else:
                    sekwencje[aktualna_nazwa] = SekwencjaDNA(
                        aktualna_nazwa, "".join(linie_aktualnej_sekwencji))
            elif not sekwencje:
                print(f"OSTRZEŻENIE: Plik '{nazwa_pliku}' nie zawiera żadnych "
                      "poprawnych wpisów FASTA.")

    except FileNotFoundError:
        print(f"BŁĄD: Plik '{nazwa_pliku}' nie został znaleziony. "
              "Upewnij się, że jest w tym samym katalogu co skrypt.")
        return None
    except Exception as e:
        print(f"Wystąpił błąd podczas odczytu pliku '{nazwa_pliku}': {e}")
        return None
    return sekwencje


def dodaj_sekwencje_uzytkownika(
    slownik_sekwencji: dict[str, SekwencjaDNA], nazwa_pliku: str = NAZWA_PLIKU_FASTA
) -> dict[str, SekwencjaDNA]:
    """
    Prosi użytkownika o podanie własnej sekwencji DNA i dodaje ją
    do pamięci programu oraz do pliku FASTA.

    Args:
        slownik_sekwencji (dict): Słownik istniejących sekwencji DNA.
        nazwa_pliku (str): Nazwa pliku FASTA, do którego ma być dopisana sekwencja.

    Returns:
        dict[str, SekwencjaDNA]: Zaktualizowany słownik sekwencji DNA.
    """
    print("\n--- Dodaj swoją sekwencję ---")
    nazwa = input("Podaj nazwę dla Twojej sekwencji (np. Moja_Sekwencja): ").strip()
    if not nazwa:
        print("Nazwa sekwencji nie może być pusta. Anulowano dodawanie "
              "sekwencji.")
        return slownik_sekwencji

    dane_sekwencji = input("Wklej sekwencję DNA (tylko A, T, C, G): ").strip().upper()

    if not dane_sekwencji:
        print("Sekwencja nie może być pusta. Anulowano dodawanie sekwencji.")
        return slownik_sekwencji

    tymczasowy_obiekt_sekwencji = SekwencjaDNA(nazwa, dane_sekwencji)
    if not tymczasowy_obiekt_sekwencji.jest_poprawna():
        print("UWAGA: Twoja sekwencja zawiera nieprawidłowe znaki. Zostanie "
              "oznaczona jako niepoprawna.")

    if nazwa in slownik_sekwencji:
        print(f"UWAGA: Sekwencja o nazwie '{nazwa}' już istnieje w pamięci. "
              "Zostanie nadpisana.")

    slownik_sekwencji[nazwa] = tymczasowy_obiekt_sekwencji
    print(f"Sekwencja '{nazwa}' została dodana do pamięci programu.")

    try:
        with open(nazwa_pliku, 'a') as f:
            f.write(f">{nazwa}\n{dane_sekwencji}\n")
        print(f"Sekwencja '{nazwa}' została trwale dopisana do pliku "
              f"'{nazwa_pliku}'.")
    except IOError as e:
        print(f"BŁĄD: Nie udało się zapisać sekwencji do pliku '{nazwa_pliku}': "
              f"{e}")
    except Exception as e:
        print(f"Wystąpił nieoczekiwany błąd podczas zapisu do pliku: {e}")

    return slownik_sekwencji


def wyswietl_info_o_sekwencji(obiekt_sekwencji: SekwencjaDNA):
    """
    Wyświetla szczegółowe informacje o pojedynczej sekwencji DNA.

    Args:
        obiekt_sekwencji (SekwencjaDNA): Obiekt sekwencji DNA do wyświetlenia.
    """
    print(f"\n--- Informacje o sekwencji: {obiekt_sekwencji.nazwa} ---")
    print(f"  Długość: {obiekt_sekwencji.pobierz_dlugosc()} nukleotydów")
    print(f"  Zawartość GC: {obiekt_sekwencji.oblicz_zawartosc_gc():.2f}%")
    print(f"  Poprawna: {'Tak' if obiekt_sekwencji.jest_poprawna() else 'Nie (zawiera niepoprawne nukleotydy)'}")
    print(f"  Typ: {obiekt_sekwencji.pobierz_typ_sekwencji()}")
    print(f"  Sekwencja (pierwsze 30 znaków): {obiekt_sekwencji.sekwencja[:30]}...")


def przetworz_sekwencje(
    slownik_sekwencji: dict[str, SekwencjaDNA]
) -> pd.DataFrame:
    """
    Przetwarza sekwencje DNA, usuwając wadliwe i zduplikowane wpisy,
    a następnie konwertuje pozostałe dane do DataFrame.

    Args:
        slownik_sekwencji (dict): Słownik sekwencji DNA do przetworzenia.

    Returns:
        pd.DataFrame: DataFrame zawierający oczyszczone i przetworzone
                      sekwencje. Pusty DataFrame, jeśli nie ma danych
                      do przetworzenia.
    """
    if not slownik_sekwencji:
        print("Brak sekwencji do przetworzenia.")
        return pd.DataFrame()

    print("\n--- Przetwarzanie sekwencji ---")
    przetworzone_dane = []
    widziane_tresci_sekwencji = set()
    poczatkowa_liczba = len(slownik_sekwencji)
    usunietych_duplikatow = 0
    usunietych_niepoprawnych = 0

    for nazwa, obiekt_sekwencji in slownik_sekwencji.items():
        if not obiekt_sekwencji.jest_poprawna():
            print(f"  Usuwam wadliwą sekwencję: {nazwa} "
                  "(zawiera niepoprawne nukleotydy)")
            usunietych_niepoprawnych += 1
            continue

        if obiekt_sekwencji.sekwencja in widziane_tresci_sekwencji:
            print(f"  Usuwam duplikat: {nazwa} (treść sekwencji już istnieje)")
            usunietych_duplikatow += 1
            continue

        widziane_tresci_sekwencji.add(obiekt_sekwencji.sekwencja)
        przetworzone_dane.append({
            "Nazwa": obiekt_sekwencji.nazwa,
            "Sekwencja": obiekt_sekwencji.sekwencja,
            "Dlugosc": obiekt_sekwencji.pobierz_dlugosc(),
            "Zawartosc_GC": obiekt_sekwencji.oblicz_zawartosc_gc(),
            "Typ_Sekwencji": obiekt_sekwencji.pobierz_typ_sekwencji()
        })

    print(f"  Usunięto duplikatów: {usunietych_duplikatow}")
    print(f"  Usunięto wadliwych sekwencji: {usunietych_niepoprawnych}")
    print(f"  Pozostało sekwencji po oczyszczeniu: {len(przetworzone_dane)} "
          f"z {poczatkowa_liczba} początkowych.")

    tabela_danych = pd.DataFrame(przetworzone_dane)
    return tabela_danych


def wizualizuj_dane(tabela_danych: pd.DataFrame):
    """
    Generuje i wyświetla różne wizualizacje danych z DataFrame.

    Args:
        tabela_danych (pd.DataFrame): DataFrame zawierający przetworzone dane
                                      sekwencji DNA.
    """
    if tabela_danych.empty:
        print("Brak danych do wizualizacji.")
        return

    print("\n--- Tworzenie wizualizacji ---")

    plt.figure(figsize=(10, 6))
    plt.hist(tabela_danych['Dlugosc'], bins=10, edgecolor='black')
    plt.title('Rozkład długości sekwencji DNA')
    plt.xlabel('Długość sekwencji (nukleotydy)')
    plt.ylabel('Liczba sekwencji')
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.hist(tabela_danych['Zawartosc_GC'], bins=10,
             edgecolor='black', color='lightgreen')
    plt.title('Rozkład zawartości GC w sekwencjach DNA')
    plt.xlabel('Zawartość GC (%)')
    plt.ylabel('Liczba sekwencji')
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.scatter(tabela_danych['Dlugosc'], tabela_danych['Zawartosc_GC'],
                alpha=0.7, color='skyblue')
    plt.title('Długość sekwencji vs. Zawartość GC')
    plt.xlabel('Długość sekwencji (nukleotydy)')
    plt.ylabel('Zawartość GC (%)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    tabela_danych['Typ_Sekwencji'].value_counts().plot(
        kind='bar', color=['purple', 'orange', 'cyan'])
    plt.title('Liczba sekwencji według typu')
    plt.xlabel('Typ sekwencji')
    plt.ylabel('Liczba sekwencji')
    plt.xticks(rotation=45)
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()


def glowna_funkcja():
    """
    Główna funkcja programu, orkiestrująca całym procesem:
    generowaniem pliku, odczytem, przetwarzaniem, wizualizacją
    i interakcją z użytkownikiem.
    """
    docelowa_liczba_sekwencji_do_pliku = random.randint(MIN_LOSOWYCH_SEKWENCJI,
                                                        MAX_LOSOWYCH_SEKWENCJI)

    if not os.path.exists(NAZWA_PLIKU_FASTA):
        utworz_plik_fasta(NAZWA_PLIKU_FASTA,
                          liczba_wszystkich_sekwencji_do_wygenerowania=
                          docelowa_liczba_sekwencji_do_pliku)
        print(f"Plik '{NAZWA_PLIKU_FASTA}' nie istniał, został utworzony "
              "z losową liczbą sekwencji "
              f"({docelowa_liczba_sekwencji_do_pliku}).")

    while True:
        wszystkie_sekwencje_slownik = wczytaj_plik_fasta(NAZWA_PLIKU_FASTA)

        if wszystkie_sekwencje_slownik is None:
            print("Nie można kontynuować z powodu krytycznego błędu odczytu "
                  "pliku. Zakończono program.")
            return

        liczba_wpisow_w_pliku = len(wszystkie_sekwencje_slownik)

        if (liczba_wpisow_w_pliku < WYMAGANE_MIN_SEKWENCJI_DLA_OSTRZEZENIA
                or liczba_wpisow_w_pliku > MAX_LOSOWYCH_SEKWENCJI):
            print(f"\nOSTRZEŻENIE: Plik '{NAZWA_PLIKU_FASTA}' zawiera "
                  f"{liczba_wpisow_w_pliku} sekwencji. Wymagane minimum to "
                  f"{WYMAGANE_MIN_SEKWENCJI_DLA_OSTRZEZENIA} sekwencji.")
        else:
            print(f"\nPlik '{NAZWA_PLIKU_FASTA}' zawiera {liczba_wpisow_w_pliku} "
                  "sekwencji. Spełnia wymagane minimum "
                  f"{WYMAGANE_MIN_SEKWENCJI_DLA_OSTRZEZENIA} sekwencji.")

        wybor_uzytkownika = input("Czy chcesz nadpisać obecny plik i wygenerować "
                                  "nowe sekwencje z losową liczbą? (tak/nie): ").lower()
        if wybor_uzytkownika == 'tak':
            docelowa_liczba_sekwencji_do_pliku = random.randint(
                MIN_LOSOWYCH_SEKWENCJI, MAX_LOSOWYCH_SEKWENCJI)
            utworz_plik_fasta(NAZWA_PLIKU_FASTA,
                              liczba_wszystkich_sekwencji_do_wygenerowania=
                              docelowa_liczba_sekwencji_do_pliku)
            continue
        else:
            print(f"Kontynuuję pracę z obecnym stanem pliku "
                  f"({liczba_wpisow_w_pliku} sekwencji).")
            break

    print("\n--- Odczytano sekwencje z pliku ---")
    if wszystkie_sekwencje_slownik:
        for nazwa, obiekt_sekwencji in \
                list(wszystkie_sekwencje_slownik.items())[:LICZBA_WYSWIETLANYCH_POCZATKOWYCH_SEKWENCJI]:
            print(f"  {nazwa}: {obiekt_sekwencji.pobierz_dlugosc()} nt, "
                  f"GC: {obiekt_sekwencji.oblicz_zawartosc_gc():.2f}%")
        if len(wszystkie_sekwencje_slownik) > LICZBA_WYSWIETLANYCH_POCZATKOWYCH_SEKWENCJI:
            print(f"  ...i jeszcze "
                  f"{len(wszystkie_sekwencje_slownik) - LICZBA_WYSWIETLANYCH_POCZATKOWYCH_SEKWENCJI} "
                  "innych sekwencji.")
    else:
        print("Brak sekwencji do wyświetlenia z pliku.")

    dodaj_wiecej = input("\nCzy chcesz dodać własną sekwencję? (tak/nie): ").lower()
    if dodaj_wiecej == 'tak':
        wszystkie_sekwencje_slownik = dodaj_sekwencje_uzytkownika(wszystkie_sekwencje_slownik,
                                                                  NAZWA_PLIKU_FASTA)

    print("\n--- Wyświetlam informacje o kilku sekwencjach (przed czyszczeniem) ---")
    if wszystkie_sekwencje_slownik:
        licznik = 0
        for nazwa, obiekt_sekwencji in wszystkie_sekwencje_slownik.items():
            wyswietl_info_o_sekwencji(obiekt_sekwencji)
            licznik += 1
            if licznik >= LICZBA_WYSWIETLANYCH_INFO_SEKWENCJI_OD_UZYTKOWNIKA:
                break
    else:
        print("Brak sekwencji do wyświetlenia informacji.")

    wyczyszczona_tabela_danych = przetworz_sekwencje(wszystkie_sekwencje_slownik)

    if not wyczyszczona_tabela_danych.empty:
        print("\n--- Ostateczna Tabela Danych z czystymi sekwencjami ---")
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        print(wyczyszczona_tabela_danych)
        pd.reset_option('display.max_rows')
        pd.reset_option('display.max_columns')
        pd.reset_option('display.width')
    else:
        print("Tabela Danych jest pusta po czyszczeniu.")

    wizualizuj_dane(wyczyszczona_tabela_danych)

    print("\nProgram zakończył działanie. Miłego dnia!")


if __name__ == "__main__":
    glowna_funkcja()
