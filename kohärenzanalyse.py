# Folgende Befehle in Kommandozeile ausführen, falls 'gensim' und 'Levenshtein' noch nicht installiert sind
# pip install gensim
# pip install Levenshtein

import csv
from Levenshtein import distance
import re
from math import comb
from tqdm import tqdm
import heapq
import multiprocessing
from tqdm.contrib.concurrent import process_map

## Bite anpassen: Pfad zu Ordner, in der sich der durchzumessende Text befindet (muss im txt-Format vorliegen) ##
filepath =  'C:/Users/Kohärenz/hss/H-000.txt'

#### Zu modifizierende Parameter ####

#  Wie soll die Hs. in der fertig erstellten Tabelle benannt werden?
Name_fuer_Hs_in_Tabelle = 'H-000'

# Auf Grundlage welches Distanzmaßes sollen die Berechnungen durchgeführt werden?
# Nicht normalisierte Levenshtein Distanz (lev_nn) = die absolute Anzahl gezählter Edit-Operation zwischen zwei Zeichenketten
# Normalisierte Levenshtein Distanz (lev_n) = die absolute Anzahl gezählter Edit-Operation zwischen zwei Zeichenketten geteilt durch die Länge der längeren Zeichenkette
# Gewichtete Levenshtein Distanz (lev_nn + dictionary 'weights' mit selbstdefinierten Gewichtungen)
# Jaccard-Distanz (jac + selbstdefinierter Wert der Variable 'laenge_ngrams') =
distance_type = 'lev_nn'

# wieviele der ähnlichsten Verse sollen für jeden Vers gefunden  werden?
n = 1

grenzwerte_versch_distance_types = {'lev_n': [0.4, 0.7], 'lev_nn': [5, 10], 'jac': [0.4, 0.7]}

# Modifizierbares Parameter speziell für Jaccard-Distanz:
laenge_ngrams = 2

# Modifizierbares Parameter speziell für nicht normalisierte Levenshteindistanz:
# hier können spezifische Gewichtungen festgelegt werden in der Form  {('ʒ', 'z'): 0, ('z', 'ʒ'): 0, etc.}
weights = {}


#### Ab hier beginnt das Programm - hier  keine Änderungen vornehmen ####
schwellenwert_orange = None
schwellenwert_rot = None
def main(path=filepath):
    # Benennung der Zieldatei
    benennung_ergebnisdatei = 'koh_' + Name_fuer_Hs_in_Tabelle + '_' + distance_type
    global schwellenwert_orange, schwellenwert_rot
    if schwellenwert_orange is None:
        schwellenwert_orange = grenzwerte_versch_distance_types[distance_type][0]
    if schwellenwert_rot is None:
        schwellenwert_rot = grenzwerte_versch_distance_types[distance_type][1]

    # Einlesen der txt-Datei
    hs = preparing_text_verse_level(open_text_as_list_of_lines(filepath))

   # number_of_posbl_comb = comb(len(hs), 2)

    indices_hs = [x for x in range(0, len(hs))]

    # csv-Datei
    csv_file = open("levg_260624.csv", 'w', encoding='utf-8', newline='')

    # Benennung der Spalten
    column_1 = 'Versnummer'
    column_2 = 'Vers'
    columns = [column_1, column_2]
    n_columns_versnummer = [str(x) + '. Versnummer' for x in range(1, n + 1)]
    n_columns_vers = [str(x) + '. Vers' for x in range(1, n + 1)]
    n_columns_dis = [str(x) + f'. {distance_type}' for x in range(1, n + 1)]
    for r in range(0, len(n_columns_vers)):
        columns.append(n_columns_versnummer[r])
        columns.append(n_columns_vers[r])
        columns.append(n_columns_dis[r])
    columns.append('Ende')  # Endmarkierung für die Funktion 'find_n_ms' -> wird in Ergebnisdatei entfernt


    # parallelize
    print("Number of cpu : ", multiprocessing.cpu_count())
    verse_list = []
    print("\nMapping indices")
    for v_index in indices_hs:
        verse_list.append(map_indices_to_verses(v_index, hs, columns, n, distance_type))

    p = multiprocessing.Pool(multiprocessing.cpu_count())

    #results = p.map(calculate_distance_for_verse, verse_list)

    # testing with first 60 indices
    # do this to test your output structure, check for runtime errors or quickly estimate runtime for all indices
    # verse_list = verse_list[:60]
    results = process_map(calculate_distance_for_verse, verse_list, max_workers=multiprocessing.cpu_count())
    sorted_results = sorted(results, key=lambda x: x[0])

    print("\nWriting results to " + csv_file.name + ":")
    with csv_file:
        header = columns
        writer = csv.DictWriter(csv_file, fieldnames=header)
        writer.writeheader()
        for result in tqdm(sorted_results):
            v_index = result[0]
            first_two_columns = {columns[0]: v_index, columns[1]: hs[v_index]}
            final_columns = merge_two_dicts(first_two_columns, result[1])
            writer.writerow(final_columns)

def map_indices_to_verses(v_index, hs, columns, n, distance_type):
    return [v_index, hs, columns, n, distance_type]

def calculate_distance_for_verse(verse):
    return [verse[0], compare_verses(verse[0], verse[1], verse[2], verse[3], verse[4])]

# Berechnung Ähnlichkeit
def calculate_distance(vers1, vers2, distance_type, weights=weights):
    if distance_type == "lev_n":
        return distance(vers1, vers2) / max(len(vers1), len(vers2))
    elif distance_type == "lev_nn":
        if len(list(weights.keys())) < 1:
             return distance(vers1, vers2)
        else:
             return weighted_levenshtein(vers1, vers2, weights)
    elif distance_type == "jac":
        list1, list2 = word2ngrams(vers1, laenge_ngrams), word2ngrams(vers2, laenge_ngrams)
        intersection = len(list(set(list1).intersection(list2)))
        anzahl_verschiedener_vorkommender_woerter = len(set(list1 + list2))
        return 1 - (intersection / anzahl_verschiedener_vorkommender_woerter)

# notwendig für Jaccard Distanz
def word2ngrams(text, n):
    return ["".join(j) for j in zip(*[text[i:] for i in range(n)])]


# mithilfe von ChatGpt optimierte Funktion, ca. 10s/it
def weighted_levenshtein(s1, s2, weights):
    len1, len2 = len(s1), len(s2)

    # Edge cases for empty strings
    if len1 == 0:
        return sum(weights.get(('', s2[j]), 1) for j in range(len2))
    if len2 == 0:
        return sum(weights.get((s1[i], ''), 1) for i in range(len1))

    # Initialize the row for dynamic programming
    current_row = [0] * (len2 + 1)

    # Initialize the first row
    for j in range(1, len2 + 1):
        current_row[j] = current_row[j - 1] + weights.get(('', s2[j - 1]), 1)

    for i in range(1, len1 + 1):
        prev_val = current_row[0]
        current_row[0] += weights.get((s1[i - 1], ''), 1)
        for j in range(1, len2 + 1):
            insert_cost = current_row[j - 1] + weights.get(('', s2[j - 1]), 1)
            delete_cost = current_row[j] + weights.get((s1[i - 1], ''), 1)
            if s1[i - 1] == s2[j - 1]:
                substitute_cost = prev_val
            else:
                substitute_cost = prev_val + weights.get((s1[i - 1], s2[j - 1]), 1)

            prev_val = current_row[j]
            current_row[j] = min(insert_cost, delete_cost, substitute_cost)

    return current_row[len2]

# Texte in Liste einzelner Verse verwandeln
def preparing_text_verse_level(text: list):
    prep_text = list(filter(None, [t.replace('\n', '') for t in text]))
    return prep_text


# Text als Liste einlesen (notwendig für Funktion 'preparing_text_verse_level')
def open_text_as_list_of_lines(filename_hs):
    with open(filename_hs, mode='r', encoding='utf-8-sig') as f:
        liste_verse = f.readlines()
        regex = re.compile(r'^### ?\d{4}')
        filtered = [i for i in liste_verse if not regex.match(i)]
        return filtered


def color_rule(val, lim_max_rot=schwellenwert_rot, lim_max_orange=schwellenwert_orange):
    result = []
    for x in val:
        if x == 'n.a.':
            result.append('background-color: grey')
        elif float(x) >= lim_max_rot:
            result.append('background-color: red')
        elif lim_max_rot > float(x) > lim_max_orange:
            result.append('background-color: orange')
        elif float(x) <= lim_max_orange:
            result.append('background-color: green')
        else:
            result.append(None)
    return result


def merge_two_dicts(x, y):
    z = x.copy()  # start with keys and values of x
    z.update(y)  # modifies z with keys and values of y
    return z


def compare_verses(matching_vers_id: int, hs: list, columns: list, n: int,
                   distance_type: str):
    heap = []
    push = heapq.heappush
    pop = heapq.heappushpop

    for v_index in range(len(hs)):
        if v_index != matching_vers_id:
            calc_distance = calculate_distance(hs[v_index], hs[matching_vers_id], distance_type)
            if len(heap) < n:
                push(heap, (-calc_distance, v_index))
            else:
                pop(heap, (-calc_distance, v_index))

    n_aehnlichste = {}
    columns_copy = columns[2:]  # die ersten beiden columns ignorieren

    for i, (neg_lev, v_id) in enumerate(sorted(heap, reverse=True)):
        idx = i * 3
        n_aehnlichste[columns_copy[idx]] = v_id
        n_aehnlichste[columns_copy[idx + 1]] = hs[v_id]
        n_aehnlichste[columns_copy[idx + 2]] = -neg_lev

    return n_aehnlichste

if __name__ == "__main__":
    main()
