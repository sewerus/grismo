import os, sys
import time

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtWidgets import *
from graph import Graph
from isomorfism_test import IsomorfismTest
from result_dialogs import *


def resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath('.'), relative_path)


def random_graph(graph: Graph):
    graph.random_graph()


def entry_iso_graph(graph_1: Graph, graph_2: Graph):
    graph_1.entry_iso_by(graph_2.neighbour_matrix)


class Grismo(QDialog):
    def __init__(self, parent=None):
        super(Grismo, self).__init__(parent)

        layout = QHBoxLayout()

        # TOP LABELS

        graph_1_group_box = QGroupBox("Graf 1")
        graph_1_box = QVBoxLayout()
        graph_1_group_box.setLayout(graph_1_box)

        graph_2_group_box = QGroupBox("Graf 2")
        graph_2_box = QVBoxLayout()
        graph_2_group_box.setLayout(graph_2_box)

        menu_group_box = QGroupBox("MENU")
        menu_box = QVBoxLayout()
        menu_group_box.setLayout(menu_box)

        # GROUP BY SPLITTERS

        splitter_1 = QSplitter(Qt.Horizontal)
        splitter_1.addWidget(graph_1_group_box)
        splitter_1.addWidget(graph_2_group_box)
        splitter_1.setSizes([700, 700])
        splitter_2 = QSplitter(Qt.Horizontal)
        splitter_2.addWidget(splitter_1)
        splitter_2.addWidget(menu_group_box)
        splitter_2.setSizes([1400, 450])
        layout.addWidget(splitter_2)

        # GRAPHS

        graph_1_draw_group_box = QGroupBox("Wizualizacja:")
        self.graph_1_draw_box = QVBoxLayout()
        graph_1_draw_group_box.setLayout(self.graph_1_draw_box)
        graph_1_box.addWidget(graph_1_draw_group_box)

        graph_2_draw_group_box = QGroupBox("Wizualizacja:")
        self.graph_2_draw_box = QVBoxLayout()
        graph_2_draw_group_box.setLayout(self.graph_2_draw_box)
        graph_2_box.addWidget(graph_2_draw_group_box)

        graph_1_scroll = QScrollArea(self)
        graph_1_box.addWidget(graph_1_scroll)
        graph_1_scroll.setWidgetResizable(True)
        graph_1_scroll_content = QWidget(graph_1_scroll)
        graph_1_scroll_layout = QVBoxLayout(graph_1_scroll_content)
        graph_1_scroll_content.setLayout(graph_1_scroll_layout)
        graph_1_scroll.setWidget(graph_1_scroll_content)

        graph_2_scroll = QScrollArea(self)
        graph_2_box.addWidget(graph_2_scroll)
        graph_2_scroll.setWidgetResizable(True)
        graph_2_scroll_content = QWidget(graph_2_scroll)
        graph_2_scroll_layout = QVBoxLayout(graph_2_scroll_content)
        graph_2_scroll_content.setLayout(graph_2_scroll_layout)
        graph_2_scroll.setWidget(graph_2_scroll_content)

        graph_1_neighbours_list_group_box = QGroupBox("Lista sąsiedztwa:")
        self.graph_1_neighbours_list_box = QVBoxLayout()
        graph_1_neighbours_list_group_box.setLayout(self.graph_1_neighbours_list_box)
        graph_1_scroll_layout.addWidget(graph_1_neighbours_list_group_box)

        graph_2_neighbours_list_group_box = QGroupBox("Lista sąsiedztwa:")
        self.graph_2_neighbours_list_box = QVBoxLayout()
        graph_2_neighbours_list_group_box.setLayout(self.graph_2_neighbours_list_box)
        graph_2_scroll_layout.addWidget(graph_2_neighbours_list_group_box)

        graph_1_neighbours_matrix_group_box = QGroupBox("Macierz sąsiedztwa:")
        self.graph_1_neighbours_matrix_box = QVBoxLayout()
        graph_1_neighbours_matrix_group_box.setLayout(self.graph_1_neighbours_matrix_box)
        graph_1_scroll_layout.addWidget(graph_1_neighbours_matrix_group_box)

        graph_2_neighbours_matrix_group_box = QGroupBox("Macierz sąsiedztwa:")
        self.graph_2_neighbours_matrix_box = QVBoxLayout()
        graph_2_neighbours_matrix_group_box.setLayout(self.graph_2_neighbours_matrix_box)
        graph_2_scroll_layout.addWidget(graph_2_neighbours_matrix_group_box)

        # MENU

        algorithm_select_group_box = QGroupBox("Wybór porównywanych algorytmów:")
        algorithm_select = QVBoxLayout()
        algorithm_select_group_box.setLayout(algorithm_select)
        menu_box.addWidget(algorithm_select_group_box)

        self.brute_method_checkbox = QCheckBox("Metoda brutalna")
        self.brute_method_checkbox.setChecked(False)
        algorithm_select.addWidget(self.brute_method_checkbox)

        self.bfs_method_checkbox = QCheckBox("Metoda przeszukiwania BFS")
        self.bfs_method_checkbox.setChecked(False)
        algorithm_select.addWidget(self.bfs_method_checkbox)

        self.dfs_method_checkbox = QCheckBox("Metoda przeszukiwania DFS")
        self.dfs_method_checkbox.setChecked(False)
        algorithm_select.addWidget(self.dfs_method_checkbox)

        wl_method_widget = QWidget()
        wl_method_box = QHBoxLayout()
        wl_method_box.setContentsMargins(0,0,0,0)
        algorithm_select.addWidget(wl_method_widget)
        wl_method_widget.setLayout(wl_method_box)
        self.wl_method_checkbox = QCheckBox("Algorytm Weisfeilera-Lehmana stopnia")
        self.wl_method_checkbox.setChecked(False)
        wl_method_box.addWidget(self.wl_method_checkbox)
        self.wl_method_dim = QSpinBox()
        self.wl_method_dim.setMinimum(1)
        self.wl_method_dim.setValue(3)
        wl_method_box.addWidget(self.wl_method_dim)

        # self.wl_method_checkbox = QCheckBox("Algorytm Weisfeilera-Lehmana stopnia")
        # self.wl_method_checkbox.setChecked(False)
        # algorithm_select.addWidget(self.wl_method_checkbox)

        self.spectrum_method_checkbox = QCheckBox("Metoda obliczania spektrum")
        self.spectrum_method_checkbox.setChecked(False)
        algorithm_select.addWidget(self.spectrum_method_checkbox)

        # TEST TYPE TABS

        test_type_tab_box = QTabWidget()

        single_test_widget = QWidget()
        vertices_test_widget = QWidget()
        edges_test_widget = QWidget()

        single_test_tab = QVBoxLayout()
        vertices_test_tab = QVBoxLayout()
        edges_test_tab = QVBoxLayout()

        single_test_widget.setLayout(single_test_tab)
        vertices_test_widget.setLayout(vertices_test_tab)
        edges_test_widget.setLayout(edges_test_tab)

        test_type_tab_box.addTab(single_test_widget, "Prosty test")
        test_type_tab_box.addTab(vertices_test_widget, "Wierzchołki")
        test_type_tab_box.addTab(edges_test_widget, "Krawędzie")

        menu_box.addWidget(test_type_tab_box)

        # -- SINGLE TEST

        single_test_vertices_amount_group_box = QGroupBox("Liczba wierzchołków grafów:")
        single_test_vertices_amount_box = QVBoxLayout()
        single_test_vertices_amount_group_box.setLayout(single_test_vertices_amount_box)
        single_test_vertices_amount_group_box.setMaximumHeight(80)
        single_test_tab.addWidget(single_test_vertices_amount_group_box)

        self.single_test_vertices_amount = QSpinBox()
        self.single_test_vertices_amount.setMinimum(3)
        self.single_test_vertices_amount.setMaximum(100)
        self.single_test_vertices_amount.setValue(5)
        self.single_test_vertices_amount.lineEdit().setReadOnly(True)
        self.single_test_vertices_amount.valueChanged.connect(self.update_vertices_amount)
        single_test_vertices_amount_box.addWidget(self.single_test_vertices_amount)

        random_group_box = QGroupBox("Losuj grafy:")
        random_box = QVBoxLayout()
        random_group_box.setLayout(random_box)
        random_group_box.setMaximumHeight(120)
        single_test_tab.addWidget(random_group_box)

        random_graph_1_button = QPushButton("Losuj graf 1")
        random_graph_1_button.setCheckable(False)
        random_graph_1_button.clicked.connect(self.random_graph_1)
        random_box.addWidget(random_graph_1_button)

        random_graph_2_button = QPushButton("Losuj graf 2")
        random_graph_2_button.setCheckable(False)
        random_graph_2_button.clicked.connect(self.random_graph_2)
        random_box.addWidget(random_graph_2_button)

        entry_iso_group_box = QGroupBox("Wprowadź izomorfizm:")
        entry_iso_box = QVBoxLayout()
        entry_iso_group_box.setLayout(entry_iso_box)
        entry_iso_group_box.setMaximumHeight(120)
        single_test_tab.addWidget(entry_iso_group_box)

        entry_iso_graph_1_button = QPushButton("Generuj graf 2 na podstawie grafu 1")
        entry_iso_graph_1_button.setCheckable(False)
        entry_iso_graph_1_button.clicked.connect(self.entry_iso_graph_2)
        entry_iso_box.addWidget(entry_iso_graph_1_button)

        entry_iso_graph_2_button = QPushButton("Generuj graf 1 na podstawie grafu 2")
        entry_iso_graph_2_button.setCheckable(False)
        entry_iso_graph_2_button.clicked.connect(self.entry_iso_graph_1)
        entry_iso_box.addWidget(entry_iso_graph_2_button)

        make_single_test_group_box = QGroupBox("Porównaj wybrane metody:")
        single_test_box = QVBoxLayout()
        make_single_test_group_box.setLayout(single_test_box)
        make_single_test_group_box.setMaximumHeight(80)
        single_test_tab.addWidget(make_single_test_group_box)

        make_single_test_button = QPushButton("Przeprowadź badanie")
        make_single_test_button.setCheckable(False)
        make_single_test_button.clicked.connect(self.make_single_test)
        single_test_box.addWidget(make_single_test_button)

        # -- VERTICES TEST

        vertices_test_group_box = QGroupBox("Parametry testu:")
        vertices_test_tab.addWidget(vertices_test_group_box)

        vertices_test_form = QFormLayout()
        vertices_test_group_box.setLayout(vertices_test_form)

        min_vertices_label = QLabel("Minimalna liczba wierzchołków")
        self.min_vertices_amount = QSpinBox()
        self.min_vertices_amount.setMinimum(3)
        self.min_vertices_amount.setMaximum(100)
        self.min_vertices_amount.setValue(5)
        self.min_vertices_amount.lineEdit().setReadOnly(True)
        self.min_vertices_amount.valueChanged.connect(self.update_max_vertices_amount)
        vertices_test_form.addRow(min_vertices_label, self.min_vertices_amount)

        max_vertices_label = QLabel("Maksymalna liczba wierzchołków")
        self.max_vertices_amount = QSpinBox()
        self.max_vertices_amount.setMinimum(3)
        self.max_vertices_amount.setMaximum(100)
        self.max_vertices_amount.setValue(5)
        self.max_vertices_amount.lineEdit().setReadOnly(True)
        self.max_vertices_amount.valueChanged.connect(self.update_min_vertices_amount)
        vertices_test_form.addRow(max_vertices_label, self.max_vertices_amount)

        self.const_edges_radio = QRadioButton("Stała liczba krawędzi")
        self.const_edges_radio.setToolTip("W każdym teście ta sama liczba krawędzi o ile to możliwe")
        self.const_edges_radio.setChecked(True)
        self.const_edges_radio.toggled.connect(self.update_const_var_edges)
        vertices_test_form.addRow(self.const_edges_radio)

        self.var_edges_radio = QRadioButton("Zmienna liczba krawędzi")
        self.var_edges_radio.setToolTip("Stosunek w % do liczby wszystkim możliwych krawędzi")
        self.var_edges_radio.toggled.connect(self.update_const_var_edges)
        vertices_test_form.addRow(self.var_edges_radio)

        const_edges_frame = QFrame()
        const_edges_layout = QHBoxLayout()
        const_edges_frame.setLayout(const_edges_layout)
        vertices_test_form.addRow(const_edges_frame)

        var_edges_frame = QFrame()
        var_edges_layout = QHBoxLayout()
        var_edges_frame.setLayout(var_edges_layout)
        vertices_test_form.addRow(var_edges_frame)
        var_edges_frame.hide()

        const_edges_label = QLabel("Liczba krawędzi")
        const_edges_label.setToolTip("W każdym teście ta sama liczba krawędzi o ile to możliwe")
        self.const_edges_amount = QSpinBox()
        self.const_edges_amount.setMinimum(1)
        self.const_edges_amount.setValue(10)
        self.const_edges_amount.setToolTip("W każdym teście ta sama liczba krawędzi o ile to możliwe")
        const_edges_layout.addWidget(const_edges_label)
        const_edges_layout.addWidget(self.const_edges_amount)

        var_edges_label = QLabel("Procent krawędzi [%]")
        var_edges_label.setToolTip("Stosunek w % do liczby wszystkim możliwych krawędzi")
        self.var_edges_amount = QSpinBox()
        self.var_edges_amount.setMinimum(1)
        self.var_edges_amount.setMaximum(100)
        self.var_edges_amount.setValue(50)
        self.var_edges_amount.setToolTip("Stosunek w % do liczby wszystkim możliwych krawędzi")
        var_edges_layout.addWidget(var_edges_label)
        var_edges_layout.addWidget(self.var_edges_amount)

        vertices_test_entry_iso_group_box = QGroupBox("Czy losować grafy izomorficzne:")
        vertices_test_entry_iso_box = QVBoxLayout()
        vertices_test_entry_iso_group_box.setLayout(vertices_test_entry_iso_box)
        vertices_test_tab.addWidget(vertices_test_entry_iso_group_box)

        self.vertices_test_iso_radio = QRadioButton("Wszystkie grafy izomorficzne")
        self.vertices_test_iso_radio.setChecked(True)
        vertices_test_entry_iso_box.addWidget(self.vertices_test_iso_radio)

        self.vertices_test_no_iso_radio = QRadioButton("Grafy losowe")
        vertices_test_entry_iso_box.addWidget(self.vertices_test_no_iso_radio)

        vertices_test_samples_group_box = QGroupBox("Liczba testów dla każdej liczby wierchołków:")
        vertices_test_box = QHBoxLayout()
        vertices_test_samples_group_box.setLayout(vertices_test_box)
        vertices_test_tab.addWidget(vertices_test_samples_group_box)

        vertices_test_samples_label = QLabel("Liczba testów")
        self.vertices_test_samples_amount = QSpinBox()
        self.vertices_test_samples_amount.setMinimum(1)
        self.vertices_test_samples_amount.setMaximum(1000)
        self.vertices_test_samples_amount.setValue(50)
        vertices_test_box.addWidget(vertices_test_samples_label)
        vertices_test_box.addWidget(self.vertices_test_samples_amount)

        make_vertices_test_group_box = QGroupBox("Porównaj wybrane metody:")
        vertices_test_box = QVBoxLayout()
        make_vertices_test_group_box.setLayout(vertices_test_box)
        vertices_test_tab.addWidget(make_vertices_test_group_box)

        make_vertices_test_button = QPushButton("Przeprowadź badanie")
        make_vertices_test_button.setCheckable(False)
        make_vertices_test_button.clicked.connect(self.make_vertices_test)
        vertices_test_box.addWidget(make_vertices_test_button)

        # -- EDGES TEST

        edges_test_group_box = QGroupBox("Parametry testu:")
        edges_test_tab.addWidget(edges_test_group_box)

        edges_test_form = QFormLayout()
        edges_test_group_box.setLayout(edges_test_form)

        edges_test_vertices_label = QLabel("Liczba wierzchołków")
        self.edges_test_vertices_amount = QSpinBox()
        self.edges_test_vertices_amount.setMinimum(3)
        self.edges_test_vertices_amount.setMaximum(100)
        self.edges_test_vertices_amount.setValue(5)
        self.edges_test_vertices_amount.lineEdit().setReadOnly(True)
        self.edges_test_vertices_amount.valueChanged.connect(self.update_max_edges_by_vertices)
        edges_test_form.addRow(edges_test_vertices_label, self.edges_test_vertices_amount)

        min_edges_label = QLabel("Minimalna liczba krawędzi")
        self.min_edges_amount = QSpinBox()
        self.min_edges_amount.setMinimum(3)
        self.min_edges_amount.setMaximum(15)
        self.min_edges_amount.setValue(5)
        self.min_edges_amount.lineEdit().setReadOnly(True)
        self.min_edges_amount.valueChanged.connect(self.update_max_edges_amount)
        edges_test_form.addRow(min_edges_label, self.min_edges_amount)

        max_edges_label = QLabel("Maksymalna liczba krawędzi")
        self.max_edges_amount = QSpinBox()
        self.max_edges_amount.setMinimum(3)
        self.max_edges_amount.setMaximum(15)
        self.max_edges_amount.setValue(5)
        self.max_edges_amount.lineEdit().setReadOnly(True)
        self.max_edges_amount.valueChanged.connect(self.update_min_edges_amount)
        edges_test_form.addRow(max_edges_label, self.max_edges_amount)

        edges_test_entry_iso_group_box = QGroupBox("Czy losować grafy izomorficzne:")
        edges_test_entry_iso_box = QVBoxLayout()
        edges_test_entry_iso_group_box.setLayout(edges_test_entry_iso_box)
        edges_test_tab.addWidget(edges_test_entry_iso_group_box)

        self.edges_test_iso_radio = QRadioButton("Wszystkie grafy izomorficzne")
        self.edges_test_iso_radio.setChecked(True)
        edges_test_entry_iso_box.addWidget(self.edges_test_iso_radio)

        self.edges_test_no_iso_radio = QRadioButton("Grafy losowe")
        edges_test_entry_iso_box.addWidget(self.edges_test_no_iso_radio)

        edges_test_samples_group_box = QGroupBox("Liczba testów dla każdej liczby krawędzi:")
        edges_test_box = QHBoxLayout()
        edges_test_samples_group_box.setLayout(edges_test_box)
        edges_test_tab.addWidget(edges_test_samples_group_box)

        edges_test_samples_label = QLabel("Liczba testów")
        self.edges_test_samples_amount = QSpinBox()
        self.edges_test_samples_amount.setMinimum(1)
        self.edges_test_samples_amount.setMaximum(1000)
        self.edges_test_samples_amount.setValue(50)
        self.edges_test_samples_amount.valueChanged.connect(self.update_min_edges_amount)
        edges_test_box.addWidget(edges_test_samples_label)
        edges_test_box.addWidget(self.edges_test_samples_amount)

        make_edges_test_group_box = QGroupBox("Porównaj wybrane metody:")
        edges_test_box = QVBoxLayout()
        make_edges_test_group_box.setLayout(edges_test_box)
        edges_test_tab.addWidget(make_edges_test_group_box)

        make_edges_test_button = QPushButton("Przeprowadź badanie")
        make_edges_test_button.setCheckable(False)
        make_edges_test_button.clicked.connect(self.make_edges_test)
        edges_test_box.addWidget(make_edges_test_button)

        # prepare graph
        self.graph_1 = Graph(5, self.graph_1_draw_box, self.graph_1_neighbours_list_box,
                             self.graph_1_neighbours_matrix_box)
        self.graph_2 = Graph(5, self.graph_2_draw_box, self.graph_2_neighbours_list_box,
                             self.graph_2_neighbours_matrix_box)
        self.graph_1.display()
        self.graph_2.display()
        self.setLayout(layout)
        self.setWindowTitle("Grismo - Badanie izomorfizmu grafów")
        self.showMaximized()
        self.setFocus()
        systray_icon = QIcon(resource_path("icon.png"))
        systray = QSystemTrayIcon(systray_icon, self)
        menu = QMenu()
        close = QAction("Close", self)
        menu.addAction(close)
        systray.setContextMenu(menu)
        systray.show()
        systray.showMessage("Grismo", "Program powstał w czasie pandemii 2020.", QSystemTrayIcon.Information)
        close.triggered.connect(self.close)

    def random_graph_1(self):
        random_graph(self.graph_1)

    def random_graph_2(self):
        random_graph(self.graph_2)

    def entry_iso_graph_1(self):
        entry_iso_graph(self.graph_1, self.graph_2)

    def entry_iso_graph_2(self):
        entry_iso_graph(self.graph_2, self.graph_1)

    def update_vertices_amount(self):
        amount = self.single_test_vertices_amount.value()
        self.graph_1.update_by_new_vertices_number(amount)
        self.graph_2.update_by_new_vertices_number(amount)

    def update_min_vertices_amount(self):
        self.min_vertices_amount.setMaximum(self.max_vertices_amount.value())

    def update_max_vertices_amount(self):
        self.max_vertices_amount.setMinimum(self.min_vertices_amount.value())

    def update_min_edges_amount(self):
        self.min_edges_amount.setMaximum(self.max_edges_amount.value())

    def update_max_edges_amount(self):
        self.max_edges_amount.setMinimum(self.min_edges_amount.value())

    def update_max_edges_by_vertices(self):
        vertices = self.edges_test_vertices_amount.value()
        self.max_edges_amount.setMaximum(vertices * (vertices - 1) / 2)

    def update_const_var_edges(self):
        const_box = self.const_edges_radio.parent().children()[-2]
        var_box = self.const_edges_radio.parent().children()[-1]
        if self.const_edges_radio.isChecked():
            const_box.show()
            var_box.hide()
        if self.var_edges_radio.isChecked():
            var_box.show()
            const_box.hide()

    def chosen_methods(self):
        methods = []
        if self.brute_method_checkbox.isChecked():
            methods.append("Brute")
        if self.bfs_method_checkbox.isChecked():
            methods.append("BFS")
        if self.dfs_method_checkbox.isChecked():
            methods.append("DFS")
        if self.wl_method_checkbox.isChecked():
            methods.append("WL")
        if self.spectrum_method_checkbox.isChecked():
            methods.append("Spectrum")
        return methods

    def make_single_test(self):
        # results: [method, result, time]
        chosen_methods = self.chosen_methods()
        wl_dim_value = self.wl_method_dim.value()
        result_dialog = SimpleTestResultsDialog(self, chosen_methods)
        result_dialog.show()
        QApplication.processEvents()

        for method in chosen_methods:
            test = IsomorfismTest(self.graph_1.neighbour_matrix, self.graph_2.neighbour_matrix, method, wl_dim_value)
            test.make_test()
            result_dialog.add_result([method, test.result, test.time])
            QApplication.processEvents()

    def make_vertices_test(self):
        start_vertices = self.min_vertices_amount.value()
        end_vertices = self.max_vertices_amount.value()
        vertices_amounts = range(start_vertices, end_vertices + 1)
        is_const_edges = self.const_edges_radio.isChecked()
        const_edges_amount = self.const_edges_amount.value()
        is_var_edges = self.var_edges_radio.isChecked()
        var_edges_amount = self.var_edges_amount.value()
        iso_graphs = self.vertices_test_iso_radio.isChecked()
        no_iso_graphs = self.vertices_test_no_iso_radio.isChecked()
        samples_amount = self.vertices_test_samples_amount.value()
        chosen_methods = self.chosen_methods()
        wl_dim_value = self.wl_method_dim.value()
        self.stop_test = False
        result_dialog = VerticesTestResultsDialog(self, chosen_methods, vertices_amounts, samples_amount)
        result_dialog.show()

        for vertices_amount in vertices_amounts:
            self.graph_1.update_by_new_vertices_number(vertices_amount)
            self.graph_2.update_by_new_vertices_number(vertices_amount)
            for i in range(samples_amount):
                if is_const_edges:
                    self.graph_1.random_graph_by_with_const_edges(const_edges_amount)
                if is_var_edges:
                    self.graph_1.random_graph_by_with_var_edges(var_edges_amount / 100)
                if iso_graphs:
                    self.graph_2.entry_iso_by(self.graph_1.neighbour_matrix)
                if no_iso_graphs:
                    if is_const_edges:
                        self.graph_2.random_graph_by_with_const_edges(const_edges_amount)
                    if is_var_edges:
                        self.graph_2.random_graph_by_with_var_edges(var_edges_amount / 100)

                self.make_many_tests(chosen_methods, vertices_amount, result_dialog, i)

                if self.stop_test:
                    break

    def make_edges_test(self):
        vertices_amount = self.edges_test_vertices_amount.value()
        start_edges = self.min_edges_amount.value()
        end_edges = self.max_edges_amount.value()
        edges_amounts = range(start_edges, end_edges + 1)
        iso_graphs = self.edges_test_iso_radio.isChecked()
        no_iso_graphs = self.edges_test_no_iso_radio.isChecked()
        samples_amount = self.edges_test_samples_amount.value()
        chosen_methods = self.chosen_methods()
        self.stop_test = False
        result_dialog = EdgesTestResultsDialog(self, chosen_methods, edges_amounts, samples_amount)
        result_dialog.show()

        for edges_amount in edges_amounts:
            self.graph_1.update_by_new_vertices_number(vertices_amount)
            self.graph_2.update_by_new_vertices_number(vertices_amount)
            for i in range(samples_amount):
                self.graph_1.random_graph_by_with_const_edges(edges_amount)
                if iso_graphs:
                    self.graph_2.entry_iso_by(self.graph_1.neighbour_matrix)
                if no_iso_graphs:
                    self.graph_2.random_graph_by_with_const_edges(edges_amount)
                self.make_many_tests(chosen_methods, edges_amount, result_dialog, i)

            if self.stop_test:
                break

    def make_many_tests(self, chosen_methods, sth_amount, result_dialog, i):
        wl_dim_value = self.wl_method_dim.value()
        for method in chosen_methods:
            test = IsomorfismTest(self.graph_1.neighbour_matrix, self.graph_2.neighbour_matrix, method, wl_dim_value)
            test.make_test()
            result_dialog.add_result([method, test.result, test.time], sth_amount, i + 1)
            QApplication.processEvents()
            if self.stop_test:
                break


if __name__ == '__main__':
    app = QApplication(sys.argv)
    grismo_test = Grismo()

    splash_image = QPixmap(resource_path("icon.png"))
    splash = QSplashScreen(splash_image)
    splash.show()

    time.sleep(1)

    grismo_test.show()
sys.exit(app.exec_())
