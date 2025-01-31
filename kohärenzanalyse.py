import csv
from Levenshtein import distance
import re
from math import comb
from tqdm import tqdm
import heapq
import multiprocessing
from tqdm.contrib.concurrent import process_map

laenge_ngrams = 2
schwellenwert_orange = None
schwellenwert_rot = None
distance_type = 'lev_gew'

grenzwerte_versch_distance_types = {'lev_n': [0.4, 0.7], 'lev_nn': [5, 10],
                                        'lev_gew': [5, 10], 'jac': [0.4, 0.7],
                                        'wmd': [0.25, 0.4]}

# Ao: Pfad zur Text-Datei (muss im txt-Format vorliegen)
dateiname = 'C:/Users/kiara/PycharmProjects/Kohärenz/hss/H-000.txt'
# Ao: Wie soll die Hs. in der finalen Tabelle genannt werden?
Name_fuer_Hs_in_Tabelle = 'H-000'
# Benennung der Zieldatei
benennung_ergebnisdatei = 'H-000_kohä_levg' + Name_fuer_Hs_in_Tabelle

weights = {
    ('ʒ', 'z'): 0,
    ('z', 'ʒ'): 0,
    ('uͦ', 'uo'): 0.5,
    ('uo', 'uͦ'): 0.5,
    ('vͦ', 'vu'): 0.5,
    ('vu', 'vͦ'): 0.5,
    ('wͦ', 'wo'): 0.5,
    ('wo', 'wͦ'): 0.5,
    ('Wͦ', 'Wo'): 0.5,
    ('Wo', 'Wͦ'): 0.5,
    ('Vͦ', 'Vo'): 0.5,
    ('Vo', 'Vͦ'): 0.5,
    ('ͤoͤ', 'oe'): 0.5,
    ('oe', 'ͤoͤ'): 0.5,
    ('vͤ', 've'): 0.5,
    ('ve', 'vͤ'): 0.5,
    ('iͤ', 'ie'): 0.5,
    ('ie', 'iͤ'): 0.5,
    ('uͤ', 'ue'): 0.5,
    ('ue', 'uͤ'): 0.5,
    ('ſ', 's'): 0,
    ('s', 'ſ'): 0,
    ('ſ', 'z'): 0,
    ('z', 'ſ'): 0,
    ('ä', 'ae'): 0,
    ('ae', 'ä'): 0,
    ('æ', 'ae'): 0,
    ('ae', 'æ'): 0,
    ('æ', 'ä'): 0,
    ('ä', 'æ'): 0,
    ('ö', 'oe'): 0.5,
    ('oe', 'ö'): 0.5,
    ('ü', 'ue'): 0.5,
    ('ue', 'ü'): 0.5,
    ('ß', 'ss'): 0.5,
    ('ss', 'ß'): 0.5,
    ('a', 'ä'): 0.75,
    ('ä', 'a'): 0.75,
    ('o', 'ö'): 0.75,
    ('ö', 'o'): 0.75,
    ('u', 'ü'): 0.75,
    ('ü', 'u'): 0.75,
    ('n̄', 'nn'): 0.5,
    ('nn', 'n̄'): 0.5,
    ('n̄', 'nt'): 0.5,
    ('nt', 'n̄'): 0.5,
    ('n̄', 'nd'): 0.5,
    ('nd', 'n̄'): 0.5,
    ('n̄', 'nde'): 0.5,
    ('nde', 'n̄'): 0.5,
    ('ē', 'en'): 0.5,
    ('en', 'ē'): 0.5,
    ('ē', 'em'): 0.5,
    ('em', 'ē'): 0.5,
    ('ā', 'an'): 0.5,
    ('an', 'ā'): 0.5,
    ('ā', 'am'): 0.5,
    ('am', 'ā'): 0.5,
    ('ū', 'un'): 0.5,
    ('un', 'ū'): 0.5,
    ('ū', 'um'): 0.5,
    ('um', 'ū'): 0.5,
    ('ẏ', 'i'): 0.5,
    ('i', 'ẏ'): 0.5,
    ('ẏ', 'j'): 0.5,
    ('j', 'ẏ'): 0.5,
    ('d͛', 'der'): 0,
    ('der', 'd͛'): 0,
    ('m͛', 'mer'): 0,
    ('mer', 'm͛'): 0,
    ('v͛', 'ver'): 0,
    ('ver', 'v͛'): 0,
    ('g͛', 'ger'): 0,
    ('ger', 'g͛'): 0,
    ('n͛', 'ner'): 0,
    ('ner', 'n͛'): 0,
    ('w͛', 'wer'): 0,
    ('wer', 'w͛'): 0,
    ('b͛', 'ber'): 0,
    ('ber', 'b͛'): 0,
    ('D͛', 'Der'): 0,
    ('Der', 'D͛'): 0,
    ('t͛', 'ter'): 0,
    ('ter', 't͛'): 0,
    ('oͮ', 'ov'): 0.5,
    ('ov', 'oͮ'): 0.5,
    ('oͮ', 'ou'): 0.5,
    ('ou', 'oͮ'): 0.5,
    ('vͥ', 'iv'): 0.5,
    ('iv', 'vͥ'): 0.5,
    ('vͥ', 'iu'): 0.5,
    ('iu', 'vͥ'): 0.5,
    ('ꝰ', 'us'): 0.5,
    ('us', 'ꝰ'): 0.5,
    ('v', 'u'): 0.5,
    ('u', 'v'): 0.5,
    ('v', 'f'): 0.5,
    ('f', 'v'): 0.5,
}
filepath = "/"


def main(path=filepath):
    zeilenumbrueche_hs = True
    zeilenbeginn_abgesetzt = True

    n = 1

    display_LevD = True

    global schwellenwert_orange, schwellenwert_rot
    if schwellenwert_orange is None:
        schwellenwert_orange = grenzwerte_versch_distance_types[distance_type][0]
    if schwellenwert_rot is None:
        schwellenwert_rot = grenzwerte_versch_distance_types[distance_type][1]

    # Fehlermeldung: 'UnpicklingError: pickle data was truncated'
    # path_to_model = "drive/MyDrive/corpora/mhdbdb-24-1.p"
    # with open(path_to_model, 'rb') as pickle_file:
    #    model = pickle.load(pickle_file)

    if zeilenumbrueche_hs is False:
        if zeilenbeginn_abgesetzt is False:
            # Einlesen der txt-Datei + Einfügen der Zeilenumbrüchen
            hs = preparing_text_if_no_linebreaks(dateiname)
        else:
            hs = entferne_absetzung(preparing_text_if_no_linebreaks(dateiname))
    else:
        if zeilenbeginn_abgesetzt is False:
            # Einlesen der txt-Datei
            hs = preparing_text_verse_level(open_text_as_list_of_lines(dateiname))
        else:
            hs = entferne_absetzung(preparing_text_verse_level(
                open_text_as_list_of_lines(dateiname)))

    number_of_posbl_comb = comb(len(hs), 2)
    #print(hs)
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
    #print(indices_hs)
    verse_list = []
    print("\nMapping indices")
    for v_index in indices_hs:
        verse_list.append(map_indices_to_verses(v_index, hs, columns, n, distance_type))

    #p = multiprocessing.Pool(multiprocessing.cpu_count-4)


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


def sort_verses(verses):
    return sorted(verses, key=lambda x: x[0])

def calculate_distance_for_verse(verse):
    return [verse[0], compare_verses(verse[0], verse[1], verse[2], verse[3], verse[4])]


# Berechnung Ähnlichkeit
def calculate_distance(vers1, vers2, distance_type, we_model=None, weights=weights):
    if distance_type == "lev_n":
        return distance(vers1, vers2) / max(len(vers1), len(vers2))
    elif distance_type == "lev_nn":
        return distance(vers1, vers2)
    elif distance_type == "lev_gew":
        return weighted_levenshtein(vers1, vers2, weights)
    elif distance_type == "jac":
        list1, list2 = word2ngrams(vers1, laenge_ngrams), word2ngrams(vers2, laenge_ngrams)
        intersection = len(list(set(list1).intersection(list2)))
        anzahl_verschiedener_vorkommender_woerter = len(set(list1 + list2))
        return 1 - (intersection / anzahl_verschiedener_vorkommender_woerter)


#    elif distance_type == "wmd":
#        return model.wv.wmdistance(vers1, vers2)

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


# Preprocessing Sonderfall: Funktion zum Entfernen der Lücken
# nach abgesetzten Buchstaben am Zeilenbeginn
def entferne_absetzung(liste_verse: list):
    result = []
    for l in liste_verse:
        l_chars = [a for a in l]
        if l[1] == ' ':
            del l_chars[1]
            result.append(''.join(l_chars))
        else:
            result.append(l)
    return result


# Einfügen von Versumbrüchen
def preparing_text_if_no_linebreaks(filename_hs):
    with open(filename_hs, mode='r', encoding='utf-8-sig') as f:
        text_as_string_with_no_lb = f.read()
        delimiter_replaced_by_lb = text_as_string_with_no_lb.translate({ord('.'): '\n', ord('/'): '\n', ord('·'): '\n'})
        whitespace_removed = ''.join(line.lstrip(' \t') for line in delimiter_replaced_by_lb.splitlines(True))
        prep_text = whitespace_removed.splitlines()
        return prep_text


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


def color_rule(val, lim_max_rot=schwellenwert_rot,
               lim_max_orange=schwellenwert_orange):
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
                   distance_type: str, display_LevD: bool = False):
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
