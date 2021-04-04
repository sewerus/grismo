import time
import bisect

from simple_graph import SimpleGraph


class IsomorfismTest:
    def __init__(self, matrix_1, matrix_2, method, wl_dim=None):
        self.graph_1 = SimpleGraph(matrix_1)
        self.graph_2 = SimpleGraph(matrix_2)
        self.n = len(matrix_1)
        self.method = method
        self.time = 0
        self.result = False
        self.wl_dim = wl_dim

    def make_test(self):
        time_start = time.time()

        if self.method == "Brute":
            self.result = self.brute_method()
        elif self.method == "BFS":
            self.result = self.bfs_method()
        elif self.method == "DFS":
            self.result = self.dfs_method()
        elif self.method == "WL":
            self.result = self.wl_method()
        elif self.method == "Spectrum":
            self.result = self.spectrum_method()

        self.time = time.time() - time_start + 0.000000000001

    def brute_method(self):
        # return self.brute_try(0, self.graph_1.matrix)
        return self.brute_try_fix([])

    def brute_try_fix(self, bijection_list):
        if not (SimpleGraph(self.graph_2.subgraph_by_bijection(bijection_list)).matrix_compare_to_another(
                SimpleGraph(self.graph_1.subgraph_by_bijection(range(len(bijection_list)))))):
            return False

        if len(bijection_list) == self.n:
            return SimpleGraph(self.graph_2.subgraph_by_bijection(bijection_list)).matrix_compare_to_another(
                self.graph_1)
        else:
            unused_vertices = list(range(self.n))
            for used_vertex in bijection_list:
                unused_vertices.remove(used_vertex)
            for vertex in unused_vertices:
                if self.brute_try_fix(bijection_list + [vertex]):
                    return True
            return False

    def brute_try(self, n, matrix):
        if n == self.n - 1:
            return SimpleGraph(matrix).matrix_compare_to_another(self.graph_2)
        else:
            # without swap
            result = self.brute_try(n + 1, matrix)
            if result:
                return True
            # with any swap
            for i in range(n + 1, self.n):
                swaped_graph = SimpleGraph(matrix)
                swaped_graph.swap_vertices(n, i)
                return self.brute_try(n + 1, swaped_graph.matrix)

    def bfs_method(self):
        neighbour_list_1 = self.graph_1.neighbour_list()
        neighbour_list_2 = self.graph_2.neighbour_list()

        # find the rarest vertices

        min_degree_1 = self.n - 1
        min_degree_2 = self.n - 1
        min_degree_vertices_1 = []
        min_degree_vertices_2 = []
        for i in range(self.n):
            if min_degree_1 == len(neighbour_list_1[i]):
                min_degree_vertices_1.append(i)
            elif min_degree_1 > len(neighbour_list_1[i]):
                min_degree_1 = len(neighbour_list_1[i])
                min_degree_vertices_1.clear()
                min_degree_vertices_1.append(i)
            if min_degree_2 == len(neighbour_list_2[i]):
                min_degree_vertices_2.append(i)
            elif min_degree_2 > len(neighbour_list_2[i]):
                min_degree_2 = len(neighbour_list_2[i])
                min_degree_vertices_2.clear()
                min_degree_vertices_2.append(i)

        if len(min_degree_vertices_1) != len(min_degree_vertices_2):
            return False

        # try by bfs by stages

        start_vertex_1 = min_degree_vertices_1[0]
        for start_vertex_2 in min_degree_vertices_2:
            stage_of_search = 0
            bijections_for_graph_1 = [start_vertex_1]
            bijections_for_graph_2 = [[start_vertex_2]]

            while len(bijections_for_graph_2) * len(
                    self.graph_1.bfs_stage_subgraph(start_vertex_1, stage_of_search + 1)):
                stage_of_search += 1

                next_vertices_for_graph_1 = self.graph_1.bfs_stage_subgraph(start_vertex_1, stage_of_search)
                next_vertices_for_graph_2 = self.graph_2.bfs_stage_subgraph(start_vertex_2, stage_of_search)

                if len(next_vertices_for_graph_1) != len(next_vertices_for_graph_2):
                    break

                bijections_for_graph_1 += next_vertices_for_graph_1
                subgraph_1_matrix = self.graph_1.subgraph_by_bijection(bijections_for_graph_1)
                subgraph_1 = SimpleGraph(subgraph_1_matrix)

                next_probable_bijections_2 = permutation(next_vertices_for_graph_2)
                new_bijections_for_graphs_2 = []

                for old_bijection in bijections_for_graph_2:
                    for next_bijection in next_probable_bijections_2:
                        probable_bijection = old_bijection + next_bijection
                        subgraph_2_matrix = self.graph_2.subgraph_by_bijection(probable_bijection)
                        if SimpleGraph(subgraph_2_matrix).matrix_compare_to_another(subgraph_1):
                            if len(probable_bijection) == self.n:
                                return True
                            new_bijections_for_graphs_2 += [probable_bijection]

                bijections_for_graph_2 = new_bijections_for_graphs_2.copy()
        return False

    def old_dfs_method(self):
        lengths_branches_1, branches_groups_1 = self.graph_1.all_dfs_branches_by_length(0)
        # print("graph 1")
        # print("lengths_branches_1: " + str(lengths_branches_1))
        # print("branches_groups_1: " + str(branches_groups_1))
        # print("----------------------")
        for start_vertex_2 in range(self.n):
            lengths_branches_2, branches_groups_2 = self.graph_2.all_dfs_branches_by_length(start_vertex_2)

            result = True

            # compare branches by lengths groups
            groups_amount = len(lengths_branches_1)
            if groups_amount != len(lengths_branches_2):
                result = False
                continue

            for i in range(groups_amount):
                if lengths_branches_1[i] != lengths_branches_2[i]:
                    result = False
                    break
                if len(branches_groups_1[i]) != len(branches_groups_2[i]):
                    result = False
                    break
            if not result:
                continue

            # print("PASS")
            # print("start_vertex2: " + str(start_vertex_2))
            # print("lengths_branches_2: " + str(lengths_branches_2))
            # print("branches_groups_2: " + str(branches_groups_2))

            # find bijection of vertices from graph_2 into graph_1
            main_possible_bijections = []

            # iterate through each group
            for group_index in range(groups_amount):
                group_1 = branches_groups_1[group_index]
                group_2 = branches_groups_2[group_index]
                # print("group_1: " + str(group_1))
                # print("group_2: " + str(group_2))

                branches_permutations = permutation(list(range(len(group_2))))
                # print("branches_permutations: " + str(branches_permutations))

                until_group_bijections = main_possible_bijections.copy()
                after_group_bijections = []

                # iterate through all permutations of this groups and create
                for single_permutation in branches_permutations:
                    # print("single_permutation: " + str(single_permutation))

                    permutation_bijections = until_group_bijections.copy()

                    # iterate through all branches in this group
                    for branch_index in range(len(single_permutation)):
                        # print("branch_index: " + str(branch_index))
                        branch_1 = group_1[branch_index]
                        branch_2 = group_2[single_permutation[branch_index]]

                        # print("possibles_for_this_permutation: " + str(permutation_bijections))

                        permutation_compatible = True
                        permutations_to_remove = []

                        # print("######################")
                        # print(branch_1)
                        # print(branch_2)
                        # print("######################")

                        # if it's first group and first branch - can add new main possibles
                        if group_index == 0 and branch_index == 0:
                            # print("HEAD")
                            permutation_bijections.append([])
                            for vertex_index in range(lengths_branches_1[group_index]):
                                permutation_bijections[-1].append([branch_1[vertex_index], branch_2[vertex_index]])
                            # print("updated main_possible_bijections:" + str(main_possible_bijections))
                        # if it's next branch add vertices if not present or destroy possible bijection
                        else:
                            # print("CHECK")
                            for previous_bijection_list in permutation_bijections:
                                # print("previous_bijection_list " + str(previous_bijection_list))
                                # check if branch is compatible
                                is_compatible = True
                                new_bijection_pairs = []
                                for vertex_index in range(lengths_branches_1[group_index]):
                                    # print("vertex_index " + str(vertex_index))
                                    single_new_bijection = [branch_1[vertex_index], branch_2[vertex_index]]
                                    # print("single_new_bijection: " + str(single_new_bijection))
                                    existed_before = False
                                    for single_old_bijection in previous_bijection_list:
                                        # print("single_old_bijection " + str(single_old_bijection))
                                        if single_new_bijection[0] == single_old_bijection[0] and single_new_bijection[1] != single_old_bijection[1]:
                                            is_compatible = False
                                            break
                                        if single_new_bijection[1] == single_old_bijection[1] and single_new_bijection[0] != single_old_bijection[0]:
                                            is_compatible = False
                                            break
                                        if single_new_bijection[0] == single_old_bijection[0] and single_new_bijection[                                            1] == single_old_bijection[1]:
                                            existed_before = True
                                            # print("pair already existed")
                                            break
                                    if not existed_before:
                                        new_bijection_pairs.append(single_new_bijection)
                                        # print("new pair: " + str(single_new_bijection))
                                        # print("updated new_bijection_pairs: " + str(new_bijection_pairs))
                                    if not is_compatible:
                                        permutation_compatible = False
                                        break
                                if is_compatible:
                                    # print("is compatible")
                                    for new_bijection in new_bijection_pairs:
                                        previous_bijection_list.append(new_bijection)
                                    new_bijection_pairs.clear()
                                    # print("updated previous_bijection_list: " + str(previous_bijection_list))
                                else:
                                    # print("not compatible!!!!")
                                    permutations_to_remove.append(previous_bijection_list)
                                    continue
                            for i in permutations_to_remove:
                                permutation_bijections.remove(i)

                    # check all bijection before start dfs_calc
                    elements_to_remove_from_permutation_bijections = []
                    for possible_to_check in [x for x in permutation_bijections if x not in main_possible_bijections]:
                        subgraph_1 = SimpleGraph(
                            self.graph_1.subgraph_by_bijection([row[0] for row in possible_to_check]))
                        subgraph_2 = SimpleGraph(
                            self.graph_2.subgraph_by_bijection([row[1] for row in possible_to_check]))
                        stupid_result = subgraph_1.matrix_compare_to_another(subgraph_2)
                        if len(possible_to_check) and stupid_result:
                            return True
                        elif not stupid_result:
                            elements_to_remove_from_permutation_bijections.append(possible_to_check)
                    for to_remove in elements_to_remove_from_permutation_bijections:
                        permutation_bijections.remove(to_remove)

                    for new_bijections in permutation_bijections:
                        after_group_bijections.append(new_bijections)
                        # print("updated after_group_bijections: " + str(after_group_bijections))

                main_possible_bijections = after_group_bijections.copy()
                # print("updated main_possible_bijections:" + str(main_possible_bijections))

                if len(main_possible_bijections) == 0:
                    # print("result FALSE")
                    result = False
                    break
            if result:
                return True

        # not returned True earlier
        return False

    def dfs_method(self):
        lengths_branches_1, branches_groups_1 = self.graph_1.all_dfs_branches_by_length(0)
        # print("graph 1")
        # print("lengths_branches_1: " + str(lengths_branches_1))
        # print("branches_groups_1: " + str(branches_groups_1))
        # print("----------------------")
        for start_vertex_2 in range(self.n):
            lengths_branches_2, branches_groups_2 = self.graph_2.all_dfs_branches_by_length(start_vertex_2)

            result = True

            # compare branches by lengths groups
            groups_amount = len(lengths_branches_1)
            if groups_amount != len(lengths_branches_2):
                result = False
                continue

            for i in range(groups_amount):
                if lengths_branches_1[i] != lengths_branches_2[i]:
                    result = False
                    break
                if len(branches_groups_1[i]) != len(branches_groups_2[i]):
                    result = False
                    break
            if not result:
                continue

            # print("PASS")
            # print("start_vertex2: " + str(start_vertex_2))
            # print("lengths_branches_2: " + str(lengths_branches_2))
            # print("branches_groups_2: " + str(branches_groups_2))

            # find bijection of vertices from graph_2 into graph_1
            main_possible_bijections = []

            # iterate through each group
            for group_index in range(groups_amount):
                group_1 = branches_groups_1[group_index]
                group_2 = branches_groups_2[group_index]
                # print("group_1: " + str(group_1))
                # print("group_2: " + str(group_2))

                until_group_bijections = main_possible_bijections.copy()
                after_group_bijections = []

                # print("START TEST")
                # iterate through all permutations of this groups and modify after_group_bijections
                if self.dfs_branches_permutation_test([], len(group_2), until_group_bijections, group_1, group_2, group_index, lengths_branches_1, after_group_bijections, main_possible_bijections):
                    return True
                # print("END TEST")
                # print("updated after_group_bijections: " + str(after_group_bijections))

                main_possible_bijections = after_group_bijections.copy()
                # print("updated main_possible_bijections:" + str(main_possible_bijections))

                if len(main_possible_bijections) == 0:
                    # print("result FALSE")
                    result = False
                    break
            if result:
                return True

        # not returned True earlier
        return False

    def dfs_branches_permutation_test(self, single_permutation, full_permutation_length, until_group_bijections, group_1, group_2, group_index, lengths_branches_1, after_group_bijections, main_possible_bijections):
        if len(single_permutation) < full_permutation_length:
            not_used = []
            for element in range(full_permutation_length):
                if not(element in single_permutation):
                    not_used.append(element)
            for element in not_used:
                if self.dfs_branches_permutation_test(single_permutation + [element], full_permutation_length, until_group_bijections, group_1, group_2, group_index, lengths_branches_1, after_group_bijections, main_possible_bijections):
                    return True
        else :
            # print("single_permutation: " + str(single_permutation))

            permutation_bijections = until_group_bijections.copy()

            # iterate through all branches in this group
            for branch_index in range(len(single_permutation)):
                # print("branch_index: " + str(branch_index))
                branch_1 = group_1[branch_index]
                branch_2 = group_2[single_permutation[branch_index]]

                # print("possibles_for_this_permutation: " + str(permutation_bijections))

                permutation_compatible = True
                permutations_to_remove = []

                # print("######################")
                # print(branch_1)
                # print(branch_2)
                # print("######################")

                # if it's first group and first branch - can add new main possibles
                if group_index == 0 and branch_index == 0:
                    # print("HEAD")
                    permutation_bijections.append([])
                    for vertex_index in range(lengths_branches_1[group_index]):
                        permutation_bijections[-1].append([branch_1[vertex_index], branch_2[vertex_index]])
                    # print("updated main_possible_bijections:" + str(main_possible_bijections))
                # if it's next branch add vertices if not present or destroy possible bijection
                else:
                    # print("CHECK")
                    for previous_bijection_list in permutation_bijections:
                        # print("previous_bijection_list " + str(previous_bijection_list))
                        # check if branch is compatible
                        is_compatible = True
                        new_bijection_pairs = []
                        for vertex_index in range(lengths_branches_1[group_index]):
                            # print("vertex_index " + str(vertex_index))
                            single_new_bijection = [branch_1[vertex_index], branch_2[vertex_index]]
                            # print("single_new_bijection: " + str(single_new_bijection))
                            existed_before = False
                            for single_old_bijection in previous_bijection_list:
                                # print("single_old_bijection " + str(single_old_bijection))
                                if single_new_bijection[0] == single_old_bijection[0] and single_new_bijection[1] != \
                                        single_old_bijection[1]:
                                    is_compatible = False
                                    break
                                if single_new_bijection[1] == single_old_bijection[1] and single_new_bijection[0] != \
                                        single_old_bijection[0]:
                                    is_compatible = False
                                    break
                                if single_new_bijection[0] == single_old_bijection[0] and single_new_bijection[1] == \
                                        single_old_bijection[1]:
                                    existed_before = True
                                    # print("pair already existed")
                                    break
                            if not existed_before:
                                new_bijection_pairs.append(single_new_bijection)
                                # print("new pair: " + str(single_new_bijection))
                                # print("updated new_bijection_pairs: " + str(new_bijection_pairs))
                            if not is_compatible:
                                permutation_compatible = False
                                break
                        if is_compatible:
                            # print("is compatible")
                            for new_bijection in new_bijection_pairs:
                                previous_bijection_list.append(new_bijection)
                            new_bijection_pairs.clear()
                            # print("updated previous_bijection_list: " + str(previous_bijection_list))
                        else:
                            # print("not compatible!!!!")
                            permutations_to_remove.append(previous_bijection_list)
                            continue
                    for i in permutations_to_remove:
                        permutation_bijections.remove(i)

            # check all bijection before start dfs_calc
            elements_to_remove_from_permutation_bijections = []
            for possible_to_check in [x for x in permutation_bijections if x not in main_possible_bijections]:
                subgraph_1 = SimpleGraph(
                    self.graph_1.subgraph_by_bijection([row[0] for row in possible_to_check]))
                subgraph_2 = SimpleGraph(
                    self.graph_2.subgraph_by_bijection([row[1] for row in possible_to_check]))
                stupid_result = subgraph_1.matrix_compare_to_another(subgraph_2)
                if len(possible_to_check) and stupid_result:
                    return True
                elif not stupid_result:
                    elements_to_remove_from_permutation_bijections.append(possible_to_check)
            for to_remove in elements_to_remove_from_permutation_bijections:
                permutation_bijections.remove(to_remove)

            for new_bijections in permutation_bijections:
                after_group_bijections.append(new_bijections)
            # print("updated after_group_bijections: " + str(after_group_bijections))
            return False

    def wl_method(self):
        graph_1_neighbour_list = self.graph_1.neighbour_list()
        graph_2_neighbour_list = self.graph_2.neighbour_list()

        # vertices' colors
        graph_1_colors = [0] * self.n
        graph_2_colors = [0] * self.n

        # repeat method method_dim times
        for i in range(self.wl_dim):
            # calc collections of neighbours' colors for each vertex
            graph_1_collection = []
            graph_2_collection = []
            for vertex in range(self.n):
                neighbours_colors = []
                for neighbour in graph_1_neighbour_list[vertex]:
                    bisect.insort(neighbours_colors, graph_1_colors[neighbour])
                graph_1_collection.append(neighbours_colors)
                neighbours_colors = []
                for neighbour in graph_2_neighbour_list[vertex]:
                    bisect.insort(neighbours_colors, graph_2_colors[neighbour])
                graph_2_collection.append(neighbours_colors)

            # prepare color - collection pairs
            pairs = []
            color_index = 0
            for vertex in range(self.n):
                collection = graph_1_collection[vertex]
                if not collection in [row[1] for row in pairs]:
                    pairs.append([color_index, collection])
                    color_index += 1

            # check if all collections from graph_2 are in prepared pairs
            for vertex in range(self.n):
                if not graph_2_collection[vertex] in [row[1] for row in pairs]:
                    return False

            # assign new colors
            for vertex in range(self.n):
                graph_1_colors[vertex] = pairs[[row[1] for row in pairs].index(graph_1_collection[vertex])][0]
                graph_2_colors[vertex] = pairs[[row[1] for row in pairs].index(graph_2_collection[vertex])][0]

        # on the end compare vertices' colors
        for vertex in range(self.n):
            if graph_1_colors[vertex] != graph_1_colors[vertex]:
                return False

        # not returned False before
        return True

    def spectrum_method(self):
        # calc spectrum
        spectrum_1 = self.graph_1.spectrum()
        spectrum_2 = self.graph_2.spectrum()

        spectrum_len = len(spectrum_1[0])

        # compare
        if spectrum_len != len(spectrum_2[0]):
            return False
        else:
            for i in range(spectrum_len):
                if spectrum_1[0][i] != spectrum_2[0][i]:
                    return False
                if spectrum_1[1][i] != spectrum_2[1][i]:
                    return False
        return True


def permutation(lst):
    if len(lst) == 0:
        return []
    if len(lst) == 1:
        return [lst]
    result = []
    for i in range(len(lst)):
        m = lst[i]
        rem_list = lst[:i] + lst[i + 1:]
        for p in permutation(rem_list):
            result.append([m] + p)
    return result
