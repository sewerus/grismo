import numpy as np


class SimpleGraph:
    def __init__(self, matrix):
        self.all_dfs_branches = []
        self.matrix = []
        for row in matrix:
            self.matrix.append([])
            for x in row:
                self.matrix[-1].append(x)
        self.n = len(matrix)

    def matrix_compare_to_another(self, another_graph):
        if self.n == len(another_graph.matrix):
            end_test = False
            result = True
            for i in range(self.n - 1):
                for j in range(i + 1, self.n):
                    if self.matrix[i][j] != another_graph.matrix[i][j]:
                        result = False
                        end_test = True
                        break
                if end_test:
                    break
        else:
            result = False
        return result

    def swap_vertices(self, a: int, b: int):
        old_a_neighbours = []
        old_b_neighbours = []
        for i in range(self.n):
            if self.matrix[i][a]:
                old_a_neighbours.append(i)
            if self.matrix[i][b]:
                old_b_neighbours.append(i)
        for i in range(self.n):
            self.matrix[i][a] = False
            self.matrix[a][i] = False
            self.matrix[i][b] = False
            self.matrix[b][i] = False
        for x in old_a_neighbours:
            if x == b:
                self.matrix[b][a] = True
                self.matrix[a][b] = True
            else:
                self.matrix[b][x] = True
                self.matrix[x][b] = True
        for x in old_b_neighbours:
            if x == a:
                self.matrix[a][b] = True
                self.matrix[b][a] = True
            else:
                self.matrix[a][x] = True
                self.matrix[x][a] = True

    def spectrum(self):
        # calc eigenvalues of matrix
        eigenvalues = np.linalg.eigvals(self.matrix)

        # count values and sort spectrum
        spectrum = [[], []]
        for value in eigenvalues:
            value = round(value, 2)
            if not(value in spectrum[0]):
                spectrum[0].append(value)
                spectrum[0].sort()
                spectrum[1].insert(spectrum[0].index(value), 1)
            else:
                spectrum[1][spectrum[0].index(value)] += 1

        return spectrum

    def neighbour_list(self):
        result = []
        for i in range(self.n):
            result.append([])
        for i in range(self.n - 1):
            for j in range(i, self.n):
                if self.matrix[i][j]:
                    result[i].append(j)
                    result[j].append(i)
        return result

    def stage_subgraph(self, start_vertex:int, stage:int):
        if stage == 0:
            return start_vertex

        neighbour_list = self.neighbour_list()

        # mark all the vertices as not visited
        visited = [False] * self.n

        # create a queue for BFS
        queue = []

        # modify BFS to know about stage of search
        stage_of_search = 0
        queue_of_stage = []

        # mark the source node as visited and enqueue it
        queue.append(start_vertex)
        visited[start_vertex] = True

        while queue:
            # dequeue a vertex from queue and print it
            s = queue.pop(0)

            # get all adjacent vertices of the dequeued vertex s.
            # if a adjacent has not been visited, then mark it visited and enqueue it
            for i in neighbour_list[s]:
                if not(visited[i]):
                    queue.append(i)
                    visited[i] = True

            # if queue_of_step is empty - it's a new step
            if s in queue_of_stage:
                queue_of_stage.remove(s)
            if len(queue_of_stage) == 0:
                queue_of_stage = queue.copy()
                stage_of_search += 1
            if stage_of_search == stage:
                return queue_of_stage

    def bfs_stage_subgraph(self, start_vertex:int, stage:int):
        if stage == 0:
            return [start_vertex]

        neighbour_list = self.neighbour_list()

        # mark all the vertices as not visited
        visited = [False] * self.n

        # create a queue for BFS
        queue = []

        # modify BFS to know about stage of search
        stage_of_search = 0
        queue_of_stage = []

        # mark the source node as visited and enqueue it
        queue.append(start_vertex)
        visited[start_vertex] = True

        while queue:
            # dequeue a vertex from queue and print it
            s = queue.pop(0)

            # get all adjacent vertices of the dequeued vertex s.
            # if a adjacent has not been visited, then mark it visited and enqueue it
            for i in neighbour_list[s]:
                if not(visited[i]):
                    queue.append(i)
                    visited[i] = True

            # if queue_of_step is empty - it's a new step
            if s in queue_of_stage:
                queue_of_stage.remove(s)
            if len(queue_of_stage) == 0:
                queue_of_stage = queue.copy()
                stage_of_search += 1
            if stage_of_search == stage:
                return queue_of_stage

    def all_dfs_branches_by_length(self, start_vertex:int):
        neighbour_list = self.neighbour_list()
        visited = [start_vertex]

        # call the recursive helper function to find DFS traversal
        self.all_dfs_branches = []
        self.find_dfs_branches(start_vertex, visited, neighbour_list)

        # group branches by length and order it
        lengths = []
        branches_groups = []
        for branch in self.all_dfs_branches:
            branch_len = len(branch)
            if branch_len in lengths:
                branches_groups[lengths.index(branch_len)].append(branch)
            else:
                lengths.append(branch_len)
                branches_groups.append([branch])

        sorted(branches_groups, key=len)
        branches_groups.reverse()
        lengths = []
        for group in branches_groups:
            lengths.append(len(group[0]))

        return [lengths, branches_groups]

    def find_dfs_branches(self, vertex, visited, neighbour_list):
        neighbours_to_visit = []

        for neighbour in neighbour_list[vertex]:
            if not(neighbour in visited):
                neighbours_to_visit.append(neighbour)

        if len(neighbours_to_visit) == 0:
            self.all_dfs_branches.append(visited)

        # recur for all the vertices adjacent to this vertex
        for neighbour in neighbours_to_visit:
                self.find_dfs_branches(neighbour, visited + [neighbour], neighbour_list)

    def subgraph_by_bijection(self, bijection_list):
        new_matrix = []
        new_n = len(bijection_list)
        for i in range(new_n):
            new_matrix.append([])
            for j in range(new_n):
                new_matrix[-1].append(False)
        for i in range(new_n-1):
            for j in range(i, new_n):
                if self.matrix[bijection_list[i]][bijection_list[j]]:
                    new_matrix[i][j] = True
                    new_matrix[j][i] = True
        return new_matrix

    def is_consistent(self):
        neighbour_list = self.neighbour_list()
        start_vertex = 0
        # mark all the vertices as not visited
        visited = [False] * len(neighbour_list)
        # create a queue for BFS
        queue = []
        # mark the source node as visited and enqueue it
        queue.append(start_vertex)
        visited[start_vertex] = True
        while queue:
            # dequeue a vertex from queue and print it
            s = queue.pop(0)
            # get all adjacent vertices of the dequeued vertex s.
            # if a adjacent has not been visited, then mark it visited and enqueue it
            for i in neighbour_list[s]:
                if not(visited[i]):
                    queue.append(i)
                    visited[i] = True
        result = True
        for vertex_visited in visited:
            if not vertex_visited:
                result = False
                break
        return result
