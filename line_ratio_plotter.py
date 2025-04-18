import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QHBoxLayout, QWidget,
    QFileDialog, QLabel, QComboBox, QTableWidget, QTableWidgetItem, QCheckBox,
    QSizePolicy, QFrame, QGroupBox, QGridLayout, QMessageBox, QHeaderView
)
from PyQt5.QtCore import Qt

class LineRatioPlotter(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Line Ratio Plotter")
        self.setGeometry(100, 100, 1200, 800)
        
        # Store loaded data
        self.objects_df = None
        
        # Available lines (column names in the CSV)
        self.available_lines = ["Lya", "OVI", "CIV", "NV"]
        
        self.create_ui()
    
    def create_ui(self):
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        
        # Left panel for controls
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Load data section
        load_group = QGroupBox("Load Data")
        load_layout = QVBoxLayout()
        
        self.load_button = QPushButton("Load Objects CSV")
        self.load_button.clicked.connect(self.load_objects_csv)
        load_layout.addWidget(self.load_button)
        
        self.file_label = QLabel("No file loaded")
        load_layout.addWidget(self.file_label)
        
        # Object list
        self.object_list = QTableWidget()
        self.object_list.setColumnCount(self.available_lines.__len__() + 2)  # Object Name, Include, Lya, OVI, CIV, NV
        headers = ["Object Name", "Include"] + self.available_lines
        self.object_list.setHorizontalHeaderLabels(headers)
        self.object_list.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.object_list.horizontalHeader().setStretchLastSection(True)
        load_layout.addWidget(self.object_list)
        
        load_group.setLayout(load_layout)
        left_layout.addWidget(load_group)
        
        # Line ratio selection
        ratio_group = QGroupBox("Line Ratio Selection")
        ratio_layout = QGridLayout()

        # X-axis controls
        x_axis_label = QLabel("X-Axis Ratio:")
        ratio_layout.addWidget(x_axis_label, 0, 0)

        x_numerator_label = QLabel("Numerator:")
        self.numerator_combo = QComboBox()
        self.numerator_combo.addItems(self.available_lines)
        
        x_denominator_label = QLabel("Denominator:")
        self.denominator_combo = QComboBox()
        self.denominator_combo.addItems(self.available_lines)

        ratio_layout.addWidget(x_numerator_label, 1, 0)
        ratio_layout.addWidget(self.numerator_combo, 1, 1)
        ratio_layout.addWidget(x_denominator_label, 2, 0)
        ratio_layout.addWidget(self.denominator_combo, 2, 1)
        y_axis_label = QLabel("Y-Axis Ratio:")
        ratio_layout.addWidget(y_axis_label, 3, 0)

        y_numerator_label = QLabel("Numerator:")
        self.y_numerator_combo = QComboBox()
        self.y_numerator_combo.addItems(self.available_lines)

        y_denominator_label = QLabel("Denominator:")
        self.y_denominator_combo = QComboBox()
        self.y_denominator_combo.addItems(self.available_lines)

        ratio_layout.addWidget(y_numerator_label, 4, 0)
        ratio_layout.addWidget(self.y_numerator_combo, 4, 1)
        ratio_layout.addWidget(y_denominator_label, 5, 0)
        ratio_layout.addWidget(self.y_denominator_combo, 5, 1)
        self.plot_button = QPushButton("Plot Ratios")
        self.plot_button.clicked.connect(self.plot_line_ratio)
        ratio_layout.addWidget(self.plot_button, 6, 0, 1, 2)
        ratio_group.setLayout(ratio_layout)
        left_layout.addWidget(ratio_group)
        
        # Options
        options_group = QGroupBox("Plot Options")
        options_layout = QVBoxLayout()
        
        self.show_labels_cb = QCheckBox("Show Object Labels")
        self.show_labels_cb.setChecked(True)
        options_layout.addWidget(self.show_labels_cb)
        
        self.log_scale_cb = QCheckBox("Use Log Scale")
        self.log_scale_cb.setChecked(False)
        options_layout.addWidget(self.log_scale_cb)
        
        options_group.setLayout(options_layout)
        left_layout.addWidget(options_group)
        
        # Add export buttons
        export_group = QGroupBox("Export")
        export_layout = QVBoxLayout()
        
        self.export_plot_button = QPushButton("Export Plot")
        self.export_plot_button.clicked.connect(self.export_plot)
        export_layout.addWidget(self.export_plot_button)
        
        self.export_data_button = QPushButton("Export Data")
        self.export_data_button.clicked.connect(self.export_data)
        export_layout.addWidget(self.export_data_button)
        
        export_group.setLayout(export_layout)
        left_layout.addWidget(export_group)
        
        # Right panel for plot
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        # Matplotlib figure
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)
        
        # Add panels to main layout
        main_layout.addWidget(left_panel, 1)
        main_layout.addWidget(right_panel, 2)
        
        # Disable buttons until data is loaded
        self.plot_button.setEnabled(False)
        self.export_plot_button.setEnabled(False)
        self.export_data_button.setEnabled(False)
    
    def load_objects_csv(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Objects CSV", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        
        if not file_path:
            return
            
        try:
            # Load the objects CSV
            self.objects_df = pd.read_csv(file_path)
            self.file_label.setText(os.path.basename(file_path))
            
            # Update object list
            self.populate_object_list()
            
            # Enable buttons
            self.plot_button.setEnabled(True)
            self.export_plot_button.setEnabled(True)
            self.export_data_button.setEnabled(True)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading objects CSV: {str(e)}")
    
    def populate_object_list(self):
        """Populate the object list table with data from the loaded CSV"""
        if self.objects_df is None:
            return
            
        # Clear the table
        self.object_list.setRowCount(0)
        
        # Add each object as a row
        for i, row in self.objects_df.iterrows():
            row_position = self.object_list.rowCount()
            self.object_list.insertRow(row_position)
            
            # Object name
            object_name = row.get('Object Name', f"Object_{i}")
            self.object_list.setItem(row_position, 0, QTableWidgetItem(str(object_name)))
            
            # Checkbox for including in plot
            checkbox = QTableWidgetItem()
            checkbox.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            checkbox.setCheckState(Qt.Checked)
            self.object_list.setItem(row_position, 1, checkbox)
            
            # Add flux values for each line
            for col_idx, line in enumerate(self.available_lines):
                if line in self.objects_df.columns:
                    flux_value = row.get(line, "")
                    self.object_list.setItem(row_position, col_idx + 2, 
                                             QTableWidgetItem(str(flux_value)))
    
    def get_selected_objects(self):
        """Get list of selected objects for plotting"""
        selected_objects = []
        for row in range(self.object_list.rowCount()):
            checkbox_item = self.object_list.item(row, 1)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                object_data = {
                    'name': self.object_list.item(row, 0).text()
                }
                
                # Get flux values for each line
                for col_idx, line in enumerate(self.available_lines):
                    item = self.object_list.item(row, col_idx + 2)
                    if item and item.text():
                        try:
                            object_data[line] = float(item.text())
                        except ValueError:
                            object_data[line] = None
                    else:
                        object_data[line] = None
                
                selected_objects.append(object_data)
        
        return selected_objects
    
    def plot_line_ratio(self):
        # Get selected objects
        selected_objects = self.get_selected_objects()
        if not selected_objects:
            QMessageBox.warning(self, "Warning", "No objects selected for plotting")
            return
        
        # Get selected lines for X axis
        x_numerator = self.numerator_combo.currentText()
        x_denominator = self.denominator_combo.currentText()
        
        if x_numerator == x_denominator:
            QMessageBox.warning(self, "Warning", "X-axis numerator and denominator must be different")
            return
        
        # Get selected lines for Y axis
        y_numerator = self.y_numerator_combo.currentText()
        y_denominator = self.y_denominator_combo.currentText()
        
        if y_numerator == y_denominator:
            QMessageBox.warning(self, "Warning", "Y-axis numerator and denominator must be different")
            return
        
        # Clear the figure
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        # Extract data for plotting
        labels = []
        x_ratios = []
        y_ratios = []
        
        for obj in selected_objects:
            if (x_numerator in obj and x_denominator in obj and 
                y_numerator in obj and y_denominator in obj and
                obj[x_numerator] is not None and obj[x_denominator] is not None and
                obj[y_numerator] is not None and obj[y_denominator] is not None and
                obj[x_denominator] > 0 and obj[y_denominator] > 0):
                
                x_ratio = np.log10(obj[x_numerator] / obj[x_denominator])
                y_ratio = np.log10(obj[y_numerator] / obj[y_denominator])
                
                labels.append(obj['name'])
                x_ratios.append(x_ratio)
                y_ratios.append(y_ratio)
        
        # Plot scatter plot
        if not x_ratios or not y_ratios:
            QMessageBox.warning(self, "Warning", "No valid ratio data to plot")
            return
        
        scatter = ax.scatter(x_ratios, y_ratios)
        
        # Add labels if requested
        if self.show_labels_cb.isChecked():
            for i, label in enumerate(labels):
                ax.annotate(label, (x_ratios[i], y_ratios[i]), 
                        textcoords="offset points", 
                        xytext=(0,10), 
                        ha='center')
        
        if self.log_scale_cb.isChecked():
            ax.set_xscale('log')
            ax.set_yscale('log')
        
        # Set labels and title
        ax.set_xlabel(f'log{x_numerator}/{x_denominator} Ratio')
        ax.set_ylabel(f'log{y_numerator}/{y_denominator} Ratio')
        ax.set_title(f'Line Ratio Comparison')
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Adjust layout and draw
        self.figure.tight_layout()
        self.canvas.draw()
    
    def export_plot(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "", "PNG Files (*.png);;PDF Files (*.pdf);;All Files (*)", options=options
        )
        
        if file_path:
            try:
                self.figure.savefig(file_path, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Success", f"Plot saved to {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error saving plot: {str(e)}")
    
    def export_data(self):
        # Get selected objects
        selected_objects = self.get_selected_objects()
        if not selected_objects:
            QMessageBox.warning(self, "Warning", "No objects selected for export")
            return
        
        # Create data for all possible ratios
        ratio_data = {}
        for obj in selected_objects:
            obj_name = obj['name']
            ratio_data[obj_name] = {}
            
            for num_line in self.available_lines:
                for denom_line in self.available_lines:
                    if num_line == denom_line:
                        continue
                        
                    ratio_name = f"{num_line}/{denom_line}"
                    
                    if (num_line in obj and denom_line in obj and 
                        obj[num_line] is not None and obj[denom_line] is not None and 
                        obj[denom_line] > 0):
                        ratio = obj[num_line] / obj[denom_line]
                        ratio_data[obj_name][ratio_name] = ratio
                    else:
                        ratio_data[obj_name][ratio_name] = None
        
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(ratio_data, orient='index')
        
        # Save to CSV
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Data", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        
        if file_path:
            try:
                df.to_csv(file_path)
                QMessageBox.information(self, "Success", f"Data saved to {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error saving data: {str(e)}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = LineRatioPlotter()
    window.show()
    sys.exit(app.exec_())