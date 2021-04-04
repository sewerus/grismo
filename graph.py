from random import randrange

from PyQt5.QtGui import QPen, QColor
from PyQt5.QtWidgets import *

import sip
import math

from simple_graph import SimpleGraph


class Graph:
    def __init__(self, n: int, draw_box: QVBoxLayout, list_box: QVBoxLayout, matrix_box: QVBoxLayout):
        self.n = n
        self.draw_box = draw_box
        self.list_box = list_box
        self.matrix_box = matrix_box
        self.scene = QGraphicsScene()
        self.neighbour_list = []
        self.neighbour_matrix = []
        self.checkbox_matrix = []
        for i in range(self.n):
            self.neighbour_list.append([])
            self.neighbour_matrix.append([])
            self.checkbox_matrix.append([])
            for j in range(self.n):
                self.neighbour_matrix[-1].append(False)
                self.checkbox_matrix[-1].append(QCheckBox())
                self.checkbox_matrix[-1][-1].setChecked(False)
                self.checkbox_matrix[-1][-1].stateChanged.connect(self.neighbour_matrix_change)

    def add_edge(self, a: int, b: int):
        self.neighbour_list[a].append(b)
        self.neighbour_list[b].append(a)
        a_list = self.list_box.children()[a].itemAt(1).widget()
        a_list.setText(a_list.text() + str(b+1) + ", ")
        b_list = self.list_box.children()[b].itemAt(1).widget()
        b_list.setText(b_list.text() + str(a+1) + ", ")
        self.neighbour_matrix[a][b] = True
        self.neighbour_matrix[b][a] = True
        self.checkbox_matrix[a][b].setChecked(True)
        self.checkbox_matrix[b][a].setChecked(True)

    def destroy_edge(self, a: int, b: int):
        self.neighbour_list[a].remove(b)
        self.neighbour_list[b].remove(a)
        a_list = self.list_box.children()[a].itemAt(1).widget()
        neighbours_string = ""
        for i in self.neighbour_list[a]:
            neighbours_string += str(i + 1) + ", "
        a_list.setText(neighbours_string)
        b_list = self.list_box.children()[b].itemAt(1).widget()
        neighbours_string = ""
        for i in self.neighbour_list[b]:
            neighbours_string += str(i + 1) + ", "
        b_list.setText(neighbours_string)
        self.neighbour_matrix[a][b] = False
        self.neighbour_matrix[b][a] = False
        self.checkbox_matrix[a][b].setChecked(False)
        self.checkbox_matrix[b][a].setChecked(False)

    def clear_all_edges(self):
        self.neighbour_list = []
        for i in range(self.n):
            self.neighbour_list.append([])
            i_list = self.list_box.children()[i].itemAt(1).widget()
            i_list.setText("")
            for j in range(self.n):
                self.neighbour_matrix[i][j] = False
                self.checkbox_matrix[i][j].setChecked(False)

    def add_random_edges(self, amount: int):
        for i in range(int(amount)):
            a = randrange(self.n)
            b = randrange(self.n)
            while a == b or self.neighbour_matrix[a][b]:
                a = randrange(self.n)
                b = randrange(self.n)
            # if not (a == b or self.neighbour_matrix[a][b]):
            self.add_edge(a, b)

    def random_graph(self):
        self.clear_all_edges()
        amount = self.n * (self.n - 1) / 4
        self.add_random_edges(int(amount))


    def random_graph_by_with_const_edges(self, edges_amount: int):
        self.clear_all_edges()
        self.add_random_edges(edges_amount)

        # check if graph is consistent
        while not(SimpleGraph(self.neighbour_matrix).is_consistent()):
            self.clear_all_edges()
            self.add_random_edges(edges_amount)

    def random_graph_by_with_var_edges(self, edges_ratio: int):
        self.clear_all_edges()
        edges_amount = edges_ratio * self.n * (self.n - 1) / 2
        self.add_random_edges(int(edges_amount))

        # check if graph is consistent
        while not(SimpleGraph(self.neighbour_matrix).is_consistent()):
            self.clear_all_edges()
            self.add_random_edges(int(edges_amount))

    def add_vertices(self, amount: int):
        matrix_grid_box = self.matrix_box.children()[0]
        for i in range(amount):
            self.neighbour_list.append([])
            vertex_label = QLabel(str(self.n + i + 1) + ": ")
            vertex_edit = QLineEdit("")
            vertex_edit.setReadOnly(True)
            new_h_box = QHBoxLayout()
            self.list_box.addLayout(new_h_box)
            new_h_box.addWidget(vertex_label)
            new_h_box.addWidget(vertex_edit)

            for row in self.neighbour_matrix:
                row.append(False)
            label_h = QLabel(str(self.n + i + 1))
            matrix_grid_box.addWidget(label_h, 0, self.n + i + 1)
            for j in range(self.n):
                row = self.checkbox_matrix[j]
                row.append(QCheckBox())
                row[-1].setChecked(False)
                row[-1].stateChanged.connect(self.neighbour_matrix_change)
                matrix_grid_box.addWidget(row[-1], j+1, self.n + i + 1)

        for i in range(amount):
            self.neighbour_matrix.append([])
            self.checkbox_matrix.append([])
            label_v = QLabel(str(self.n + i + 1))
            matrix_grid_box.addWidget(label_v, self.n + i + 1, 0)
            for j in range(self.n + amount):
                self.neighbour_matrix[-1].append(False)
                self.checkbox_matrix[-1].append(QCheckBox())
                self.checkbox_matrix[-1][j].setChecked(False)
                self.checkbox_matrix[-1][j].stateChanged.connect(self.neighbour_matrix_change)
                if self.n + i != j:
                    matrix_grid_box.addWidget(self.checkbox_matrix[-1][j], self.n + i + 1, j + 1)

        self.n = self.n + amount

    def destroy_vertices(self, amount: int):
        matrix_grid_box = self.matrix_box.children()[0]
        for i in range(amount):
            sip.delete(matrix_grid_box.itemAtPosition(0, self.n + i).widget())
            sip.delete(matrix_grid_box.itemAtPosition(self.n + i, 0).widget())
            sip.delete(self.list_box.children()[-1].itemAt(1).widget())
            sip.delete(self.list_box.children()[-1].itemAt(0).widget())
            sip.delete(self.list_box.children()[-1])

            for j in range(self.n):
                self.neighbour_matrix[j].pop()
                self.matrix_box.removeWidget(self.checkbox_matrix[j][-1])
                sip.delete(self.checkbox_matrix[j][-1])
                self.checkbox_matrix[j].pop()
                try:
                    self.neighbour_list[j].remove(self.n - 1)
                    neighbours_string = ""
                    for k in self.neighbour_list[j]:
                        neighbours_string += str(k + 1) + ", "
                    display_list = self.list_box.children()[j].itemAt(1).widget()
                    display_list.setText(neighbours_string)

                except ValueError:
                    pass
            self.neighbour_list.pop()
            self.neighbour_matrix.pop()
            for j in range(self.n - 1):
                self.matrix_box.removeWidget(self.checkbox_matrix[-1][j])
                sip.delete(self.checkbox_matrix[-1][j])
            self.checkbox_matrix.pop()
        self.n = self.n - amount

    def update_by_new_vertices_number(self, new_n: int):
        if self.n > new_n:
            for i in range(0, self.n - new_n):
                self.destroy_vertices(1)
        elif self.n < new_n:
            self.add_vertices(int(new_n - self.n))
        self.draw()

    def display_neighbour_list(self):
        for i in range(self.n):
            vertex_label = QLabel(str(i + 1) + ": ")
            neighbours_string = ""
            for j in self.neighbour_list[i]:
                neighbours_string += str(j + 1) + ", "
            vertex_edit = QLineEdit(neighbours_string[:-2])
            vertex_edit.setReadOnly(True)

            new_h_box = QHBoxLayout()
            self.list_box.addLayout(new_h_box)
            new_h_box.addWidget(vertex_label)
            new_h_box.addWidget(vertex_edit)

    def display_neighbour_matrix(self):
        new_grid_box = QGridLayout()
        new_grid_box.setVerticalSpacing(2)
        new_grid_box.setHorizontalSpacing(2)
        self.matrix_box.addLayout(new_grid_box)
        for i in range(self.n):
            label_h = QLabel(str(i + 1))
            label_v = QLabel(str(i + 1))
            new_grid_box.addWidget(label_h, 0, i + 1)
            new_grid_box.addWidget(label_v, i + 1, 0)

        for i in range(self.n):
            for j in range(self.n):
                if i != j:
                    new_grid_box.addWidget(self.checkbox_matrix[i][j], i + 1, j + 1)
                    self.checkbox_matrix[i][j].stateChanged.connect(self.neighbour_matrix_change)

    def display(self):
        self.display_neighbour_list()
        self.display_neighbour_matrix()
        graph_1_view = QGraphicsView(self.scene)
        graph_1_view.show()
        self.draw_box.addWidget(graph_1_view)
        self.draw()

    def draw(self):
        self.scene.clear()
        r = 120
        a = 2 * math.pi / self.n
        vertex_pen = QPen()
        vertex_pen.setWidth(9)
        edge_pen = QPen()
        edge_pen.setWidth(3)
        edge_pen.setColor(QColor(139, 69, 19))
        for i in range(self.n - 1):
            for j in range(i+1, self.n):
                if self.neighbour_matrix[i][j]:
                    self.scene.addLine(r * math.cos(i*a)+5, r * math.sin(i*a)+5, r * math.cos(j*a)+5,
                                       r * math.sin(j*a)+5, edge_pen)
        for i in range(self.n):
            self.scene.addEllipse(r * math.cos(i*a), r * math.sin(i*a), 10, 10, vertex_pen)
            title = self.scene.addText(str(i+1))
            title.setPos((r+20) * math.cos(i*a)-3, (r+20) * math.sin(i*a)-5)

    def neighbour_matrix_change(self):
        founded = False
        for i in range(self.n):
            for j in range(self.n):
                if i != j:
                    if self.checkbox_matrix[i][j].isChecked() and not (self.neighbour_matrix[i][j]):
                        self.add_edge(i, j)
                        founded = True
                    elif not (self.checkbox_matrix[i][j].isChecked()) and self.neighbour_matrix[i][j]:
                        self.destroy_edge(i, j)
                        founded = True
            if founded:
                self.draw()
                break

    def entry_iso_by(self, another_graph_matrix):
        self.clear_all_edges()
        for i in range(1, self.n):
            for j in range(i+1, self.n):
                if another_graph_matrix[i][j]:
                    self.add_edge(i-1, j-1)
            if another_graph_matrix[0][i]:
                self.add_edge(i-1, self.n-1)
