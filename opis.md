#Co ten program robi?#

1.  **Tworzy magiczną książkę przepisów (plik FASTA):**
    Na początku programu, jeśli nie mamy jeszcze naszej "książki przepisów"
    (czyli pliku o nazwie `sekwencje.txt`), program sam ją dla nas stworzy!
    Ta książka będzie zawierała od 31 do 50 losowych przepisów (sekwencji DNA).
    Każdy przepis ma swoją nazwę (np. "Sekwencja_1", "Sekwencja_2")
    i różną długość – od krótkich na 50 liter do dłuższych na 199 liter.
    
    Ważne jest to, że w tej książce specjalnie umieścimy kilka "sztuczek":
    * **Kilka takich samych przepisów:** To tak, jakbyśmy mieli kilka kopii
        tego samego przepisu w książce. Nazywamy je "duplikatami".
        Program doda ich 6, np. "Duplikat_A_1", "Duplikat_A_2" itd.,
        które będą kopiami jednej z normalnych sekwencji.
    * **Kilka przepisów z "literówkami":** Czasem zdarza się, że w DNA
        pojawią się błędy, tak jak literówki w przepisie. Nasz program
        celowo umieści 6 takich przepisów z błędnymi literami.
        Nazywamy je "niepoprawnymi" sekwencjami.
    * **Losowe literówki:** Nawet w "normalnych" przepisach (tych bez
        specjalnych nazw duplikatów czy błędnych), czasem przypadkowo
        pojawią się drobne "literówki" (nieprawidłowe litery).
        To sprawia, że nasza książka jest bardziej realistyczna,
        bo w prawdziwym DNA też zdarzają się niespodzianki.

2.  **Pyta, czy chcesz zacząć od nowa:**
    Jeśli już masz taką "książkę przepisów" (`sekwencje.txt`), program
    zapyta Cię, czy chcesz ją zresetować i stworzyć zupełnie nową,
    czy wolisz użyć tej, którą już masz. To tak, jakbyś mógł wybrać,
    czy zaczynasz nową historię, czy kontynuujesz czytanie tej, którą już znasz.

3.  **Czyta Twoją książkę:**
    Program dokładnie przeczyta całą "książkę przepisów" i zapamięta
    każdy przepis (sekwencję DNA) w swojej "pamięci" (w komputerze).
    Jeśli znajdzie jakieś dziwne rzeczy w książce (np. pustą stronę
    po tytule przepisu), to Ci o tym powie, ale spróbuje mimo to
    czytać dalej.

4.  **Dodajesz swój własny przepis:**
    Możesz być "kucharzem"! Program zapyta Cię, czy chcesz dodać
    swój własny, wymyślony przepis (sekwencję DNA). Podajesz jego nazwę
    i same literki DNA (A, T, C, G). Program dopisze Twój przepis do
    książki i zapamięta go. Jeśli Twój przepis będzie miał "literówki",
    program Cię o tym uprzedzi.

5.  **Pokazuje informacje o przepisach:**
    Program wyświetli Ci na ekranie najważniejsze informacje o kilku
    pierwszych przepisach, a także o tym, który Ty dodałeś. Będziesz
    wiedział:
    * **Nazwę przepisu**
    * **Jak długi jest przepis** (ile ma liter/nukleotydów)
    * **Ile w nim jest literek 'G' i 'C' razem** (to ważna cecha DNA,
        nazywamy ją "zawartością GC")
    * **Czy przepis jest "poprawny"** (czy nie ma literówek)
    * **Jaki jest "typ" przepisu** (czy ma dużo G i C, mało, czy tak
        w sam raz).

6.  **"Sprząta" i porządkuje przepisy:**
    To bardzo ważna część! Program przeanalizuje wszystkie przepisy i:
    * **Wyrzuci przepisy z literówkami:** Te "niepoprawne" sekwencje
        zostaną usunięte, bo są wadliwe.
    * **Zostawi tylko jeden egzemplarz każdego przepisu:** Jeśli znajdzie
        kilka identycznych przepisów (duplikatów), zostawi tylko jeden,
        a resztę usunie. Dzięki temu w naszej książce będzie porządek
        i nie będziemy mieli zbędnych kopii.
    Program powie Ci, ile przepisów usunął i ile zostało po sprzątaniu.

7.  **Tworzy "Tabelę Wyników":**
    Wszystkie czyste, uporządkowane przepisy trafiają do specjalnej
    "Tabeli Wyników" (nazywamy ją "DataFrame"). To taka tabela,
    w której każda linijka to jeden przepis, a kolumny to jego cechy
    (nazwa, długość, zawartość GC, typ). Ta tabela jest bardzo łatwa do
    analizowania. Program wyświetli ją całą na ekranie, abyś mógł ją zobaczyć.

8.  **Rysuje "Obrazki" z danych:**
    Żeby łatwiej było zrozumieć nasze przepisy, program narysuje cztery
    "obrazki" (wykresy):
    * **Wykres 1 (Długości):** Pokazuje, ile jest przepisów krótkich,
        ile średnich, a ile długich.
    * **Wykres 2 (Zawartości GC):** Pokazuje, ile przepisów ma dużo G i C,
        ile mało, a ile tak w sam raz.
    * **Wykres 3 (Długość a GC):** To taki "tajny" wykres, który pokazuje,
        czy dłuższe przepisy mają więcej G i C, czy może nie ma na to reguły.
        Każda kropka to jeden przepis.
    * **Wykres 4 (Typy przepisów):** Pokazuje, ile mamy przepisów typu
        "GC-rich" (dużo G i C), "AT-rich" (mało G i C) i "Standard".
