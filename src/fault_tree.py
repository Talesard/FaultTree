from copy import copy
import numpy as np
import itertools
from collections import Counter

import graphviz

def find_paths_algorithm(a, b, G, paths=[], q=[], visited=set()): # рекурсивный алгоритм поиска путей
    if a == b:
        paths.append(copy(q))
        return paths
    for u in G[a]:
        if u not in visited:
            q.append(u)
            visited.add(u)
            paths = find_paths_algorithm(u, b, G, paths, q, visited)
            q.pop()
            visited.remove(u)
    return paths

def find_all_paths(v1, v2, G) -> list: # обертка
    return  [x[::-1]+[v1] for x in find_paths_algorithm(v1, v2, G)]

def find_largest_path(paths):
    largest_path = paths[0]
    for path in paths:
        if len(path) > len(largest_path):
            largest_path = path
    return largest_path

def is_can_be_merged(path1, path2, vertex_index):
    for level in range(vertex_index, 0, -1):
        if not (path1[level] == path2[vertex_index] and path1[level-1] == path2[level-1]):
            return False
    return True

def get_all_pairs_on_level(paths, level):
    level_vertexes = []
    for i in range(len(paths)):
        try:
            level_vertexes.append({'vertex':paths[i][level], 'path_id': i})
        except:
            pass
    pairs = list(itertools.combinations(level_vertexes, r=2))
    return pairs

def merge_paths_to_tree(paths) -> dict: # построение дерева всех путей
    G = {paths[0][0]: [],}

    largest_path = find_largest_path(paths)
    for level in range(1, len(largest_path)):
        all_pairs = get_all_pairs_on_level(paths, level)
        to_merge_paths_ids = set()
        to_merge_then_add = []
        for pair in all_pairs:
            pair_can_be_merged = is_can_be_merged(paths[pair[0]["path_id"]], paths[pair[1]["path_id"]], level)
            if pair_can_be_merged:
                to_merge_then_add.append(pair)
                to_merge_paths_ids.add(pair[0]['path_id'])
                to_merge_paths_ids.add(pair[1]['path_id'])
        not_to_merge_paths_ids = list(set([i for i in range(len(paths))]).difference(to_merge_paths_ids))
        to_merge_paths_ids = list(to_merge_paths_ids)

        for path_id in not_to_merge_paths_ids:
            if paths[path_id][level-1] in G.keys():
                G[paths[path_id][level-1]].append(paths[path_id][level])
            else:
                G[paths[path_id][level-1]] = [paths[path_id][level]]
        
        for path_id in to_merge_paths_ids:
            if paths[path_id][level-1] in G.keys():
                if paths[path_id][level] not in G[paths[path_id][level-1]]:
                    G[paths[path_id][level-1]].append(paths[path_id][level])
            else:
                G[paths[path_id][level-1]] = [paths[path_id][level]]

    changed_parents_keys = []
    for key in list(G):
        counts = Counter(G[key])
        for counts_key in counts.keys():
            if counts[counts_key] > 1:
                changed_parents_keys.append(key)
                G[key] = [x for x in G[key] if x != counts_key]
                G[key].append(counts_key)
                for i in range(1, counts[counts_key]):
                    G[key + '_'+str(i)] = [counts_key + '_'+str(i)]

    counts_changed = [0 for i in range(len(changed_parents_keys))]
    for key in G.keys():
        for i_changed in range(len(changed_parents_keys)):
            try:
                index = G[key].index(changed_parents_keys[i_changed])
                counts_changed[i_changed] += 1
                if counts_changed[i_changed] > 1:
                    G[key][index] = G[key][index] + '_' + str(index+1)
            except:
                pass

    child_set = set()
    counts = {}
    for key in G.keys():
        for child_id in range(len(G[key])):
            if G[key][child_id] not in child_set:
                child_set.add(G[key][child_id])
                counts[G[key][child_id]] = 1
            else:
                counts[G[key][child_id]] += 1
                G[key][child_id] = G[key][child_id] + '_' + str(counts[G[key][child_id]])
    return G

# +названия, +причины, визуализация
def plot_graph(G, paths, names, reasons, debug=False) -> None:
    dot = graphviz.Digraph(comment='Fault Tree')
    for reasons_block in reasons.keys():
        counter = 0
        for reason in reasons[reasons_block]:
            dot.node(f'{reasons_block}_{counter}', reason)
            counter += 1
    for key in G.keys():
        if debug:
            dot.node(key, key, shape='box')
            for child in G[key]:
                dot.node(child, child)
                dot.edge(child, key)
                # dot.edge(key, child)
        else:
            if key == paths[0][0]:
                dot.node(key, 'Отказ блока ' + names[key] + ' <=> Отказ системы', shape='box')
            else:
                dot.node(key, 'Отказ блока ' + names[key.split('_')[0]], shape='box')
            for child in G[key]:
                dot.node(child, 'Отказ блока ' + names[child.split('_')[0]], shape='box')
                # dot.edge(key, child)
                dot.edge(child, key)

    for key in G.keys():
        for child in G[key]:
            child_clear = child.split('_')[0]
            for i in range(len(reasons[f'r_{child_clear}'])):
                v1 = child
                v2 = f'r_{child_clear}_{i}'
                dot.edge(v2, v1)
                # dot.edge(v1, v2)
    last_block = paths[0][0]
    for i in range(len(reasons[f'r_{last_block}'])):
        dot.edge(f'r_{last_block}_{i}', last_block)
        # dot.edge(last_block, f'r_{last_block}_{i}')
    dot.view('out/tree.gz')

def generate_fault_tree(G, names, reasons):
    paths = find_all_paths('1', '8', G)
    tree = merge_paths_to_tree(paths)
    plot_graph(tree, paths, names, reasons)