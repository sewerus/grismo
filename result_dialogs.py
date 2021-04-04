from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import *
from pyqtgraph import PlotWidget, mkPen, AxisItem


def mean(array: list):
    return sum(array) / len(array)


class SimpleTestResultsDialog(QMainWindow):
    def __init__(self, parent=None, methods=None):
        super(SimpleTestResultsDialog, self).__init__(parent)
        self.setWindowTitle("Grismo - Wyniki dla prostego testu")
        self.setMinimumWidth(500)
        self.results_count = 0
        self.all_results = []
        self.methods = methods

        central = QWidget()
        self.setCentralWidget(central)

        result_dialog_layout = QVBoxLayout()
        central.setLayout(result_dialog_layout)

        self.results_table = QTableWidget()
        result_dialog_layout.addWidget(self.results_table)
        self.results_table.setRowCount(len(methods))
        self.results_table.setColumnCount(3)
        self.results_table.setColumnWidth(2, 200)
        self.results_table.setHorizontalHeaderLabels(["Metoda", "Wynik", "Czas [s]"])

        for i in range(len(methods)):
            self.results_table.setItem(i, 0, QTableWidgetItem(methods[i]))

        # plot box
        x_axis_dict = dict(enumerate(methods))
        x_axis = AxisItem(orientation='bottom')
        x_axis.setTicks([x_axis_dict.items()])
        self.results_plot = PlotWidget(axisItems={'bottom': x_axis})
        self.results_plot.setBackground('w')
        self.results_plot.setTitle("Porównanie metod dla wskazanych grafów")
        self.results_plot.setLabel('left', 'Czas obliczeń [s]', color='k', size=10)
        self.results_plot.setLabel('bottom', 'Metody ', color='k', size=10)
        self.results_plot.setMaximumWidth(600)
        self.results_plot.showGrid(y=True)
        result_dialog_layout.addWidget(self.results_plot)

        # prepare plot data
        pen = mkPen(color='k', width=2)
        self.plot_data = self.results_plot.plot([], [], pen=pen, symbol='+', symbolSize=10, symbolBrush='k')

    def add_result(self, result):
        for i in range(len(result)):
            self.results_table.setItem(self.results_count, i, QTableWidgetItem(str(result[i])))
        self.results_count = self.results_count + 1
        self.all_results.append(result[-1])
        self.plot_data.setData(range(len(self.all_results)), self.all_results)


class ComplicatedTestResultsDialog(QMainWindow):
    def __init__(self, parent=None, methods=None, vertices_amounts=None, samples_amount=None, name_1=None, name_2=None):
        super(ComplicatedTestResultsDialog, self).__init__(parent)
        self.parent_window = parent
        self.setWindowTitle("Grismo - Wyniki dla testu wpływu liczby " + name_1)
        self.setMinimumWidth(630)
        self.setMinimumHeight(800)
        self.setMaximumHeight(900)
        self.methods = methods
        self.results_count = 0
        self.start_vertices = vertices_amounts[0]
        self.all_vertices_amount = len(vertices_amounts)
        self.samples_amount = samples_amount
        # store all results to calculate means
        self.all_results = []

        central = QWidget()
        self.setCentralWidget(central)
        central_layout = QVBoxLayout()
        central.setLayout(central_layout)

        # plot box
        plot_label = QLabel("Wykres dla wszystkich testów:")
        central_layout.addWidget(plot_label)
        self.plot_log_mode = False
        self.results_plot = PlotWidget()
        self.results_plot.setLogMode(False, self.plot_log_mode)
        self.results_plot.setBackground('w')
        self.results_plot.setTitle("Badanie wpływu liczby " + name_1 + " na czas testu")
        self.results_plot.setLabel('left', 'Czas obliczeń [s]', color='k', size=10)
        self.results_plot.setLabel('bottom', 'Liczba ' + name_1, color='k', size=10)
        self.results_plot.setXRange(self.start_vertices, vertices_amounts[-1])
        self.results_plot.setMaximumWidth(600)
        self.results_plot.showGrid(y=True)
        central_layout.addWidget(self.results_plot)

        switch_plot_log_button = QPushButton("Zmień skalę osi Y (logarytmiczna/liniowa)")
        switch_plot_log_button.setCheckable(False)
        switch_plot_log_button.clicked.connect(self.switch_plot_log)
        central_layout.addWidget(switch_plot_log_button)

        # prepare plot lines
        self.plot_data = []
        method_colors = ['k', 'b', 'g', 'r', 'y']
        self.results_plot.plot([], [], name="")
        self.results_plot.addLegend()
        for method_index in range(len(self.methods)):
            method_name = self.methods[method_index]
            method_color = method_colors[method_index]
            pen = mkPen(color=method_color, width=2)
            self.plot_data.append(self.results_plot.plot([], [], name=method_name, pen=pen, symbol='+', symbolSize=10,
                                                         symbolBrush=method_color))
            self.results_plot.addLegend()

        # tables box
        tables_box = QScrollArea(self)
        tables_box.setWidgetResizable(True)
        tables_box_content = QWidget(tables_box)
        tables_box_layout = QVBoxLayout(tables_box_content)
        tables_box_content.setLayout(tables_box_layout)
        tables_box.setWidget(tables_box_content)
        central_layout.addWidget(tables_box)

        # for each vertices_amount prepare table: label -> table -> label -> progress_bar
        self.results_tables = []
        self.results_progress_bars = []
        bold_font = QFont()
        bold_font.setBold(True)
        for i in vertices_amounts:
            vertices_label = QLabel("Wyniki dla grafów o  " + str(i) + " " + name_2 + ":")
            vertices_label.setFont(bold_font)
            results_table = QTableWidget()
            results_table.setRowCount(len(methods))
            results_table.setColumnCount(4)
            results_table.setColumnWidth(1, 150)
            results_table.setColumnWidth(2, 150)
            results_table.setMinimumHeight((len(methods) + 1) * 30)
            results_table.setHorizontalHeaderLabels(["Metoda", "Wyniki pozytywne", "Wyniki negatywne",
                                                     "Średni czas [s]"])
            progress_label = QLabel("Postęp:")
            results_progress_bar = QProgressBar()
            results_progress_bar.setValue(0)

            self.results_tables.append(results_table)
            self.results_progress_bars.append(results_progress_bar)
            tables_box_layout.addWidget(vertices_label)
            tables_box_layout.addWidget(results_table)
            tables_box_layout.addWidget(progress_label)
            tables_box_layout.addWidget(results_progress_bar)

            self.all_results.append([])

            for method_index in range(len(methods)):
                method_title = methods[method_index]
                results_table.setItem(method_index, 0, QTableWidgetItem(method_title))
                self.all_results[-1].append([])

    def add_result(self, result, vertices_amount, sample_number):
        # result: [method, decision, time]

        # vertices table index
        table_index = vertices_amount - self.start_vertices

        # method index
        method_index = self.methods.index(result[0])

        # add result to all stored results
        self.all_results[table_index][method_index].append(result)

        # positive and negatives
        # firstly extract 2nd column from results matrix
        new_positives = sum([row[1] for row in self.all_results[table_index][method_index]])
        new_negatives = len(self.all_results[table_index][method_index]) - new_positives

        self.results_tables[table_index].setItem(method_index, 1, QTableWidgetItem(str(new_positives)))
        self.results_tables[table_index].setItem(method_index, 2, QTableWidgetItem(str(new_negatives)))

        # mean
        new_mean = mean([row[2] for row in self.all_results[table_index][method_index]])
        self.results_tables[table_index].setItem(method_index, 3, QTableWidgetItem(str(new_mean)))

        # progress_bar
        self.results_progress_bars[table_index].setValue(sample_number / self.samples_amount * 100)
        self.results_count = self.results_count + 1

        # update plot
        self.update_plot()

    def update_plot(self):
        for method_index in range(len(self.methods)):
            # for this method find all mean values
            x = []
            y = []
            for vertices_index in range(self.all_vertices_amount):
                vertices_amount = vertices_index + self.start_vertices
                if len(self.all_results[vertices_index][method_index]):
                    x.append(vertices_amount)
                    y.append(mean([row[2] for row in self.all_results[vertices_index][method_index]]))
            self.plot_data[method_index].setData(x, y)

    def switch_plot_log(self):
        self.plot_log_mode = not self.plot_log_mode
        self.results_plot.setLogMode(False, self.plot_log_mode)

    def closeEvent(self, event):
        self.parent_window.stop_test = True


class EdgesTestResultsDialog(ComplicatedTestResultsDialog):
    def __init__(self, parent=None, methods=None, edges_amounts=None, samples_amount=None):
        super(EdgesTestResultsDialog, self).__init__(parent, methods, edges_amounts, samples_amount,
                                                     "krawędzi", "krawędziach")


class VerticesTestResultsDialog(ComplicatedTestResultsDialog):
    def __init__(self, parent=None, methods=None, edges_amounts=None, samples_amount=None):
        super(VerticesTestResultsDialog, self).__init__(parent, methods, edges_amounts, samples_amount,
                                                        "wierzchołków", "wierzchołkach")
