import sys
import csv
from functools import partial
import numpy as np
from scipy.integrate import simpson
from astropy.table import vstack
from astropy.io import fits
import lmfit
import os
import math
import matplotlib.pyplot as plt
import logging
import logging.handlers as handlers

from PyQt5.QtWidgets import *
from gaussian_fitfunctions import SpectralLine, DoubletInfo
import csv
from astropy.stats import poisson_conf_interval as pcf
from PyQt5.QtCore import QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from gaussian_fitfunctions import SpectralRegion
from matplotlib.figure import Figure
from prepare_spec import coadd, prepare, prepare_other_grating, combine_tables
import pickle

class GaussianWindow(QWidget):
    """
    This "window" is a QWidget. Since it has no parent, it
    will appear as a free-floating window for comparison w/ other flux
    results. (THEORETICALLY. in reality it does not work.)
    """
    def __init__(self, fitter, mcmc_result, spectral_lines):
        super().__init__()
        print(f'GaussianWindow init')
        # layout = QVBoxLayout()
        label = QLabel("Gaussian Fit Results")
        # layout.addWidget(self.label)
        # self.setLayout(layout)
        
        # # Create tabs
        # tabs = QTabWidget()
        
        # # Flux results tab
        # flux_tab = QWidget()
        # flux_layout = QVBoxLayout(flux_tab)
        
        # # Create table for flux results
        # flux_table = QTableWidget()
        # flux_table.setColumnCount(6)
        # flux_table.setHorizontalHeaderLabels([
        #     "Line", "Component", "Flux", "Error (Low)", "Error (High)", "EW (Å)"
        # ])
        
        # # Analyze and add results for each line
        # row_count = 0
        # for line in spectral_lines:
        #     # Calculate line flux
        #     line_result = fitter.get_line_flux(mcmc_result, line.name)
            
        #     if not line_result:
        #         continue
            
        #     # Add total flux row
        #     flux_table.insertRow(row_count)
        #     flux_table.setItem(row_count, 0, QTableWidgetItem(line.name))
        #     flux_table.setItem(row_count, 1, QTableWidgetItem("Total"))
        #     flux_table.setItem(row_count, 2, QTableWidgetItem(f"{line_result['flux']:.3e}"))
        #     flux_table.setItem(row_count, 3, QTableWidgetItem(f"{line_result['error_down']:.3e}"))
        #     flux_table.setItem(row_count, 4, QTableWidgetItem(f"{line_result['error_up']:.3e}"))
        #     flux_table.setItem(row_count, 5, QTableWidgetItem("N/A"))
        #     row_count += 1
            
        #     # Add component rows
        #     for comp, comp_result in line_result['components'].items():
        #         flux_table.insertRow(row_count)
        #         flux_table.setItem(row_count, 0, QTableWidgetItem(""))
        #         flux_table.setItem(row_count, 1, QTableWidgetItem(comp))
        #         flux_table.setItem(row_count, 2, QTableWidgetItem(f"{comp_result['flux']:.3e}"))
        #         flux_table.setItem(row_count, 3, QTableWidgetItem(f"{comp_result['error_down']:.3e}"))
        #         flux_table.setItem(row_count, 4, QTableWidgetItem(f"{comp_result['error_up']:.3e}"))
        #         flux_table.setItem(row_count, 5, QTableWidgetItem("N/A"))
        #         row_count += 1
        
        # # Adjust table layout
        # flux_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # flux_table.horizontalHeader().setStretchLastSection(True)
        
        # flux_layout.addWidget(flux_table)
        
        # # Add export button
        # export_btn_layout = QHBoxLayout()
        # export_flux_btn = QPushButton("Export Flux Results")
        # export_btn_layout.addWidget(export_flux_btn)
        # flux_layout.addLayout(export_btn_layout)
        
        # # Component parameters tab
        # params_tab = QWidget()
        # params_layout = QVBoxLayout(params_tab)
        
        # # Create table for parameter results
        # params_table = QTableWidget()
        # params_table.setColumnCount(4)
        # params_table.setHorizontalHeaderLabels([
        #     "Component", "Parameter", "Value", "Error"
        # ])
        
        # # Get unique components
        # all_components = set()
        # for line in spectral_lines:
        #     all_components.update(line.components)
        
        # # Add rows for each component parameter
        # row_count = 0
        # for comp in sorted(all_components):
        #     # Get redshift
        #     z_key = f"z_{comp}"
        #     if z_key in mcmc_result.var_names:
        #         z_percentiles = np.percentile(mcmc_result.flatchain[z_key], [16, 50, 84])
        #         z_value = z_percentiles[1]
        #         z_error = (z_percentiles[2] - z_percentiles[0]) / 2
                
        #         params_table.insertRow(row_count)
        #         params_table.setItem(row_count, 0, QTableWidgetItem(comp))
        #         params_table.setItem(row_count, 1, QTableWidgetItem("Redshift"))
        #         params_table.setItem(row_count, 2, QTableWidgetItem(f"{z_value:.5f}"))
        #         params_table.setItem(row_count, 3, QTableWidgetItem(f"±{z_error:.5f}"))
        #         row_count += 1
            
        #     # Get sigma
        #     sigma_key = f"sigma_{comp}"
        #     if sigma_key in mcmc_result.var_names:
        #         sigma_percentiles = np.percentile(mcmc_result.flatchain[sigma_key], [16, 50, 84])
        #         sigma_value = sigma_percentiles[1]
        #         sigma_error = (sigma_percentiles[2] - sigma_percentiles[0]) / 2
                
        #         params_table.insertRow(row_count)
        #         params_table.setItem(row_count, 0, QTableWidgetItem(comp))
        #         params_table.setItem(row_count, 1, QTableWidgetItem("Sigma (km/s)"))
        #         params_table.setItem(row_count, 2, QTableWidgetItem(f"{sigma_value:.1f}"))
        #         params_table.setItem(row_count, 3, QTableWidgetItem(f"±{sigma_error:.1f}"))
        #         row_count += 1
        
        # # Adjust table layout
        # params_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # params_table.horizontalHeader().setStretchLastSection(True)
        
        # params_layout.addWidget(params_table)
        
        # # Export parameters button
        # export_params_btn = QPushButton("Export Parameter Results")
        # params_layout.addWidget(export_params_btn)
        
        # # Add tabs to tab widget
        # tabs.addTab(flux_tab, "Flux Results")
        # tabs.addTab(params_tab, "Component Parameters")
        
        # # Add tab widget to dialog
        # layout.addWidget(tabs)
        
        # Add close button
        # close_btn = QPushButton("Close")
        # close_btn.clicked.connect(dialog.accept)
        # layout.addWidget(close_btn)
        
        # Define export functions
        def export_flux_results():
            print(f'buggy for now')
            # file_path, _ = QFileDialog.getSaveFileName(
            #     dialog, "Save Flux Results", "", "CSV Files (*.csv);;All Files (*)"
            # )
            
            # if file_path:
            #     import csv
            #     with open(file_path, 'w', newline='') as f:
            #         writer = csv.writer(f)
            #         writer.writerow(["Line", "Component", "Flux", "Error (Low)", "Error (High)"])
                    
            #         for row in range(flux_table.rowCount()):
            #             line_name = flux_table.item(row, 0).text()
            #             component = flux_table.item(row, 1).text()
            #             flux_val = flux_table.item(row, 2).text()
            #             error_low = flux_table.item(row, 3).text()
            #             error_high = flux_table.item(row, 4).text()
                        
            #             writer.writerow([line_name, component, flux_val, error_low, error_high])
                
            #     self.statusBar().showMessage(f"Flux results exported to {file_path}")
        


class SpectralFluxApp(QMainWindow):
    # Main app with all the doohickeys and functionalities
    def __init__(self, csv_file_path=None):
        super().__init__()

        self.init_variables()
        self.create_main_widget()
        self.create_plot_area()
        self.create_controls()
        self.create_table()
        self.create_right_panel()  # Add this line

        self.setWindowTitle('Spectral Flux Measurement')
        self.setGeometry(100, 100, 1200, 800)  # Make the window a bit wider

        self.setMinimumSize(800, 600)  # Reasonable minimum size
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # self.setWindowTitle('Spectral Flux Measurement')
        # self.setGeometry(100, 100, 1000, 800)

    def init_variables(self):
        self.spectrum_data = []
        self.grating_data = []
        self.coadded_spectrum = None
        self.fits_file_paths = []
        self.all_plotted_models = []
        self.redshift = 0.0
        self.selection_step = 0
        self.cid = None  # Connection ID for event handler
        self.scale_value = 1e15 #value to scale the flux by
        self.flux_method = "Direct Integration"

        self.line_wavelengths = np.array([0.0,
            1215.67, 1031.92, 1037.61, 1238.82, 1242.8,
            1393.75, 1402.77, 1548.19, 1550.77
        ])
        self.line_labels = [
            'Full Spectrum', 
            'Lya', 'O VI 1031', 'O VI 1037', 'N V 1238',
            'N V 1242', 'Si 1393', 'Si 1402', 'C IV 1548', 'C IV 1550'
        ]

        # Variables for plotting selection overlays
        self.drawn_spectral_lines = []
        self.left_continuum_lines = []
        self.left_continuum_patch = None
        self.right_continuum_lines = []
        self.right_continuum_patch = None
        self.integration_lines = []
        self.integration_patch = None

    def create_main_widget(self):
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        
        # Create a horizontal layout for the entire application
        self.app_layout = QHBoxLayout(self.main_widget)
        
        # Create a vertical layout for the existing components (left side)
        self.main_layout = QVBoxLayout()
        self.app_layout.addLayout(self.main_layout)
        
        # Create a vertical layout for the new components (right side)
        self.right_layout = QVBoxLayout()
        self.app_layout.addLayout(self.right_layout)
        
        # Set the size ratio between left and right panels (2:1)
        self.app_layout.setStretch(0, 2)
        self.app_layout.setStretch(1, 1)

    def create_plot_area(self):
        self.figure = Figure(figsize=(10, 8))
        self.canvas = FigureCanvas(self.figure)
        self.main_layout.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.main_layout.addWidget(self.toolbar)

    def create_controls(self):
        self.controls_layout = QHBoxLayout()
        self.main_layout.addLayout(self.controls_layout)

        # # Load data button
        # self.load_button = QPushButton('Load CSV with FITS Paths')
        # self.load_button.clicked.connect(self.load_csv_with_fits_paths)
        # self.controls_layout.addWidget(self.load_button)

        # # Co-add and plot button
        # self.coadd_button = QPushButton('Co-add and Plot')
        # self.coadd_button.clicked.connect(self.coadd_and_plot)
        # self.controls_layout.addWidget(self.coadd_button)


        # Save Spectrum button
        self.save_button = QPushButton('Save Spectrum')
        self.save_button.clicked.connect(self.save_spectrum)
        self.controls_layout.addWidget(self.save_button)

        # Load Spectrum button
        self.load_spectrum_button = QPushButton('Load Spectrum')
        self.load_spectrum_button.clicked.connect(self.load_spectrum)
        self.controls_layout.addWidget(self.load_spectrum_button)

        # Redshift label and input
        self.redshift_label = QLabel('Redshift:')
        self.redshift_label.setWordWrap(True)
        self.controls_layout.addWidget(self.redshift_label)

        self.redshift_input = QLineEdit(self)
        # self.redshift_input.setPlaceholderText("Enter redshift value")
        self.redshift_input.setText(str(self.redshift))
        self.redshift_input.returnPressed.connect(self.update_redshift)
        self.redshift_input.setMinimumWidth(100)
        self.controls_layout.addWidget(self.redshift_input)

        self.scale_label = QLabel('Scale by:')
        self.scale_label.setWordWrap(True)
        self.controls_layout.addWidget(self.scale_label)

        self.scale_input = QLineEdit(self)
        # self.scale_input.setPlaceholderText(str(self.scale_value))
        self.scale_input.setText(f'{self.scale_value:.2e}')
        self.scale_input.returnPressed.connect(self.update_scale)
        self.scale_input.setMinimumWidth(100)
        self.controls_layout.addWidget(self.scale_input)

        # # Plot expected line locations button
        # self.plot_lines_button = QPushButton('Plot Expected Line Locations')
        # self.plot_lines_button.setEnabled(False)
        # self.controls_layout.addWidget(self.plot_lines_button)

        # Line selection dropdown
        self.surrounding_label = QLabel('Plot area surrounding:')
        self.surrounding_label.setWordWrap(True)
        self.controls_layout.addWidget(self.surrounding_label)

        self.line_dropdown = QComboBox(self)
        self.line_dropdown.addItems(self.line_labels)
        self.line_dropdown.currentIndexChanged.connect(self.zoom_to_line)
        self.controls_layout.addWidget(self.line_dropdown)

        # Progress bar
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setValue(0)
        self.main_layout.addWidget(self.progress_bar)

        # Status bar
        self.statusBar().showMessage('Ready')

    def create_table(self):
        # Table to display flux and flux error for each line
        self.table = QTableWidget(self)
        self.table.setColumnCount(10)  # Reduced from 10 to 7 columns
        self.table.setHorizontalHeaderLabels([
            'Line', 'Flux', 'Flux Error', 
            'Left Lower', 'Left Upper', 
            'Right Lower', 'Right Upper',
            'Continuum Slope', 'Continuum Intercept',
            'Gaussian Fit File'
        ])
        self.table.setRowCount(len(self.line_labels) - 1)  # Skip 'Full Spectrum'

        # Populate the table
        for i, line_label in enumerate(self.line_labels[1:]):
            self.table.setItem(i, 0, QTableWidgetItem(line_label))
            self.table.setItem(i, 1, QTableWidgetItem(''))
            self.table.setItem(i, 2, QTableWidgetItem(''))
            self.table.setItem(i, 3, QTableWidgetItem(''))  # Left lower bound
            self.table.setItem(i, 4, QTableWidgetItem(''))  # Left upper bound
            self.table.setItem(i, 5, QTableWidgetItem(''))  # Right lower bound
            self.table.setItem(i, 6, QTableWidgetItem(''))  # Right upper bound
            self.table.setItem(i, 7, QTableWidgetItem(''))  # Continuum Slope
            self.table.setItem(i, 8, QTableWidgetItem(''))  # Continuum Intercept
            self.table.setItem(i, 9, QTableWidgetItem(''))  # Gaussian fit file (if applicable)

        self.main_layout.addWidget(self.table)
        self.table.resizeColumnsToContents()
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Connect the row selection signal
        self.table.itemSelectionChanged.connect(self.on_table_selection_change)

    def on_table_selection_change(self):
        # Get the selected row(s)
        selected_rows = self.table.selectionModel().selectedRows()
        if selected_rows:
            row_index = selected_rows[0].row()
            # Get the line label from the first column
            line_item = self.table.item(row_index, 0)
            if line_item:
                line_label = line_item.text()
                # Update the display label
                self.selected_line_display.setText(line_label)
                # Find the matching index in line_labels (add 1 because we skip 'Full Spectrum')
                self.current_line = self.line_labels.index(line_label)
                # Enable the action buttons
                self.update_button_states(True)

                # self.action_button.setEnabled(True)
                # self.subtract_continuum_button.setEnabled(True)
                # self.undo_button.setEnabled(True)
    # Modified on_object_selection_change method to handle the object CSV files
    def on_object_selection_change(self):
        """Handle when a row in the objects table is selected"""
        # Get the selected row
        selected_rows = self.objects_table.selectionModel().selectedRows()
        if not selected_rows:
            return
                
        row_index = selected_rows[0].row()
        
        # Get the object name from column 0
        object_name_item = self.objects_table.item(row_index, 0)
        if not object_name_item or not object_name_item.text():
            self.statusBar().showMessage("No object name found for this row")
            return
            
        object_name = object_name_item.text()
        
        # Get or create the object CSV file path
        csv_dir = "/Users/eliseke/Research/UV_spectral_fitting/data/Object_csvs"
        csv_filename = f"{object_name}_lines.csv"
        csv_path = os.path.join(csv_dir, csv_filename)
        
        # Ensure the directory exists
        os.makedirs(csv_dir, exist_ok=True)
        
        # Check if the CSV file exists, and create it if not
        if not os.path.exists(csv_path):
            self.create_object_csv(csv_path)
        else:
            # Load data from existing CSV to left table
            self.update_line_values_from_object_csv(row_index)
        
        # Update the object_csv cell in the right table
        self.objects_table.setItem(row_index, 9, QTableWidgetItem(csv_path))
        
        # Get the fits file path from column 10 (Object_fits)
        fits_item = self.objects_table.item(row_index, 10)
        if hasattr(self, 'spectral_lines') and self.spectral_lines:
            self.save_spectral_lines()
        if hasattr(self, 'component_parameters') and self.component_parameters:
            self.save_component_parameters()
        
        # Clear existing spectral lines and parameters
        #note: this does not work
        self.spectral_lines = []
        self.component_parameters = {}
        self.clear_line_entries()
        # try now?
        self.clear_all_lines()
        
        # Then load the ones for the newly selected object
        self.load_spectral_lines()
        self.load_component_parameters()
        if fits_item and fits_item.text():
            fits_path = fits_item.text()
            
            # Check if the path exists - if not, just use the path as reference without loading
            if os.path.exists(fits_path):
                # Load the spectrum from the fits file
                try:
                    self.load_spectrum_from_path(fits_path)
                    self.statusBar().showMessage(f"Loaded spectrum from {fits_path}")
                except Exception as e:
                    self.statusBar().showMessage(f"Error loading spectrum: {str(e)}")
        
        self.statusBar().showMessage(f"Working with object: {object_name}, CSV: {csv_path}")

    def initialize_line_with_redshift(self, line, update_existing=False):
        """Initialize line components with the current redshift"""
        for comp in line.components:
            param_key = f"{line.name}_{comp}"
            
            # Skip if the parameter already exists and we're not updating
            if param_key in self.component_parameters and not update_existing:
                continue
                
            # Create or update the parameter entry
            if param_key not in self.component_parameters:
                self.component_parameters[param_key] = {}
            
            # Set redshift parameters
            self.component_parameters[param_key]['z_value'] = self.redshift
            self.component_parameters[param_key]['z_min'] = max(0, self.redshift - 0.001)
            self.component_parameters[param_key]['z_max'] = self.redshift + 0.001
            self.component_parameters[param_key]['z_vary'] = True
            
            # Set default values for other parameters if they don't exist
            if 'sigma_value' not in self.component_parameters[param_key]:
                self.component_parameters[param_key]['sigma_value'] = 200
                self.component_parameters[param_key]['sigma_min'] = 10
                self.component_parameters[param_key]['sigma_max'] = 1000
                self.component_parameters[param_key]['flux_value'] = 100
        
        # Save the updated parameters
        self.save_component_parameters()
    def create_object_csv(self, csv_path):
        """Create a new object CSV file with headers matching the left table"""
        try:
            # Create headers based on the left table structure
            headers = ["Redshift"]
            
            # Create line data structure
            line_data = [str(self.redshift)]  # First row is redshift
            
            # Write the CSV file
            with open(csv_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(headers)
                writer.writerow(line_data)
                
                # Add placeholder rows for fitting data (can be populated later)
                for i in range(len(self.fits_file_paths)):
                    writer.writerow(["placeholder_path", "placeholder_grating"])
                    
            self.statusBar().showMessage(f"Created new object CSV file: {csv_path}")
        except Exception as e:
            self.statusBar().showMessage(f"Error creating object CSV: {str(e)}")

    def load_object_csv_to_left_table(self, csv_path):
        """Load data from an object CSV file to the left table"""
        try:
            with open(csv_path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                rows = list(reader)
                
                if len(rows) > 0:
                    # First row should be headers
                    # Second row should have redshift
                    if len(rows) > 1:
                        try:
                            self.redshift = float(rows[1][0])
                            self.redshift_input.setText(str(self.redshift))
                        except (ValueError, IndexError):
                            self.statusBar().showMessage("Could not read redshift from CSV")
                    
                    # Read line measurements if they exist in the CSV
                    self.read_line_measurements_from_csv(csv_path)
        except Exception as e:
            self.statusBar().showMessage(f"Error loading object CSV: {str(e)}")

    # New method to read line measurements from CSV
    def read_line_measurements_from_csv(self, csv_path):
        """Read line measurements from CSV and update the left table"""
        try:
            # Define the mapping between CSV rows and table rows
            # This assumes your CSV has a specific structure for line measurements
            line_row_mapping = {
                "Lya": 0,
                "OVI": 1,
                "NV": 2,
                "CIV": 3,
            }
            
            # Try to read the CSV as a table with line names as rows
            import pandas as pd
            try:
                df = pd.read_csv(csv_path, index_col=0)
                
                # Update the left table with values from the DataFrame
                for line_name, table_row in line_row_mapping.items():
                    if line_name in df.index:
                        line_data = df.loc[line_name]
                        
                        # Update flux value if available
                        if 'Flux' in line_data:
                            self.table.setItem(table_row, 1, QTableWidgetItem(str(line_data['Flux'])))
                        
                        # Update error value if available
                        if 'Flux_Error' in line_data:
                            self.table.setItem(table_row, 2, QTableWidgetItem(str(line_data['Flux_Error'])))
                        
                        # Update continuum bounds if available
                        if 'Left_Lower' in line_data:
                            self.table.setItem(table_row, 3, QTableWidgetItem(str(line_data['Left_Lower'])))
                        if 'Left_Upper' in line_data:
                            self.table.setItem(table_row, 4, QTableWidgetItem(str(line_data['Left_Upper'])))
                        if 'Right_Lower' in line_data:
                            self.table.setItem(table_row, 5, QTableWidgetItem(str(line_data['Right_Lower'])))
                        if 'Right_Upper' in line_data:
                            self.table.setItem(table_row, 6, QTableWidgetItem(str(line_data['Right_Upper'])))
                        if 'Slope' in line_data:
                            self.table.setItem(table_row, 7, QTableWidgetItem(str(line_data['Slope'])))
                        if 'Intercept' in line_data:
                            self.table.setItem(table_row, 8, QTableWidgetItem(str(line_data['Intercept'])))
            except Exception as e:
                self.statusBar().showMessage(f"Could not read line measurements: {str(e)}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error reading line measurements from CSV: {str(e)}")

    def update_object_csv_from_left_table(self, row_index):
        """Update the individual object CSV file with data from the left table"""
        # Get the object CSV path
        csv_item = self.objects_table.item(row_index, 9)
        if not csv_item or not csv_item.text():
            return
            
        csv_path = csv_item.text()
        
        try:
            # First get the column headers from the table
            headers = []
            for col in range(self.table.columnCount()):
                headers.append(self.table.horizontalHeaderItem(col).text())
            
            # Create a list to store the rows
            rows = [headers]  # First row is the header
            
            # Extract all data from the left table
            for row in range(self.table.rowCount()):
                row_data = []
                for col in range(self.table.columnCount()):
                    item = self.table.item(row, col)
                    row_data.append(item.text() if item and item.text() else "")
                rows.append(row_data)
            
            # Save to CSV with redshift as the first row
            with open(csv_path, 'w', newline='') as csvfile:
                # First write the redshift
                writer = csv.writer(csvfile)
                writer.writerow(["Redshift"])
                writer.writerow([str(self.redshift)])
                
                # Then write the table data
                for row in rows:
                    writer.writerow(row)
                    
            self.statusBar().showMessage(f"Updated object CSV: {csv_path}")
        except Exception as e:
            self.statusBar().showMessage(f"Error updating object CSV: {str(e)}")

    def update_line_values_from_object_csv(self, row_index):
        """Update the line values in the left table from the selected object's CSV file"""
        # Get the object CSV path
        csv_item = self.objects_table.item(row_index, 9)  # Object_csv column
        
        if not csv_item or not csv_item.text():
            self.statusBar().showMessage("No CSV file found for this object")
            return
        
        csv_path = csv_item.text()
        
        try:
            # Clear the left table values first
            for row in range(self.table.rowCount()):
                for col in range(1, self.table.columnCount()):  # Skip the first column (line name)
                    self.table.setItem(row, col, QTableWidgetItem(""))
            
            # Read the CSV file
            with open(csv_path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                
                # First two rows should be "Redshift" and the redshift value
                header = next(reader, None)
                if header and header[0] == "Redshift":
                    redshift_row = next(reader, None)
                    if redshift_row and redshift_row[0]:
                        self.redshift = float(redshift_row[0])
                        self.redshift_input.setText(str(self.redshift))
                
                # Next row should be column headers
                headers = next(reader, None)
                if not headers:
                    self.statusBar().showMessage(f"CSV file {csv_path} has invalid format")
                    return
                
                # Map column names to their indices
                col_indices = {name: idx for idx, name in enumerate(headers) if name}
                
                # Process data rows
                for row_data in reader:
                    if not row_data or len(row_data) < len(headers):
                        continue
                    
                    # Get the line name from the first column
                    line_name = row_data[0]
                    
                    # Map line names to their row indices in the left table
                    line_row_mapping = {
                        "Lya": 0,
                        "O VI 1031": 1,  
                        "O VI 1037": 2,
                        "N V 1238": 3,
                        "N V 1242": 4,
                        "Si 1393": 5,
                        "Si 1402": 6,
                        "C IV 1548": 7,
                        "C IV 1550": 8,
                    }
                    
                    # If this line exists in our mapping
                    if line_name in line_row_mapping:
                        left_row = line_row_mapping[line_name]
                        
                        # Update each column value
                        for col_name, col_idx in col_indices.items():
                            if col_idx < len(row_data) and col_name != "Line":  # Skip the line name column
                                # Map column names to their column indices in the left table
                                col_mapping = {
                                    "Flux": 1,
                                    "Flux Error": 2,
                                    "Left Lower": 3,
                                    "Left Upper": 4,
                                    "Right Lower": 5,
                                    "Right Upper": 6,
                                    "Continuum Slope": 7,
                                    "Continuum Intercept": 8,
                                    "Gaussian Fit File": 9
                                }
                                
                                if col_name in col_mapping:
                                    left_col = col_mapping[col_name]
                                    self.table.setItem(left_row, left_col, QTableWidgetItem(row_data[col_idx]))
            
            self.statusBar().showMessage(f"Updated left table from CSV: {csv_path}")
        except Exception as e:
            self.statusBar().showMessage(f"Error reading object CSV: {str(e)}")

    def init_gaussian_fit_mode(self):
            """Initialize variables needed for Gaussian fit mode"""
            # Create internal storage directory if it doesn't exist
            os.makedirs("internal", exist_ok=True)
            
            # Initialize spectral lines list if not exists
            if not hasattr(self, 'spectral_lines'):
                self.spectral_lines = []
            
            # Initialize component parameters dictionary if not exists
            if not hasattr(self, 'component_parameters'):
                self.component_parameters = {}
                self.load_component_parameters()
            
            # Load any existing spectral lines
            self.load_spectral_lines()

    def create_new_object_csv(self, csv_path, row_index):
        """Create a new CSV file with the default structure for this object"""
        try:
            # First get the column headers from the table
            headers = []
            for col in range(self.table.columnCount()):
                headers.append(self.table.horizontalHeaderItem(col).text())
            
            # Get values from the right table for this object
            line_mappings = {
                "Lya": (1, 0),       # (Column in objects table, Row in left table)
                "Lya_err": (2, 0),
                "OVI": (3, 1),
                "OVI_err": (4, 1),
                "CIV": (5, 3),
                "CIV_err": (6, 3),
                "NV": (7, 2),
                "NV_err": (8, 2)
            }
            
            # Map line names as they appear in the left table
            line_names = ["Lya", "O VI 1031", "N V 1238", "C IV 1548"]
            
            # Create rows with default values
            rows = [headers]  # First row is the header
            
            # Initialize with values from the objects table where available
            for i, line_name in enumerate(line_names):
                row_data = [line_name]  # First column is the line name
                
                # Fill the rest with empty strings by default
                for j in range(1, len(headers)):
                    row_data.append("")
                
                # Update Flux and Flux Error if available in objects table
                for obj_key, (obj_col, left_row) in line_mappings.items():
                    if left_row == i:  # If this mapping is for the current line
                        value_item = self.objects_table.item(row_index, obj_col)
                        if value_item and value_item.text():
                            # Determine which column to update (Flux or Flux Error)
                            update_col = 1 if "_err" not in obj_key else 2
                            row_data[update_col] = value_item.text()
                
                rows.append(row_data)
            
            # Get redshift value
            redshift = 0.0
            try:
                redshift = float(self.redshift_input.text())
            except:
                pass
            
            # Save to CSV with redshift as the first row
            with open(csv_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["Redshift"])
                writer.writerow([str(redshift)])
                
                # Then write the table data
                for row in rows:
                    writer.writerow(row)
                    
            self.statusBar().showMessage(f"Created new object CSV: {csv_path}")
        except Exception as e:
            self.statusBar().showMessage(f"Error creating object CSV: {str(e)}")
    def action_on_selected_line(self):
        # Get the flux method
        flux_method = self.flux_method_dropdown.currentText()
        
        # Depending on the method, call the appropriate function
        if flux_method == "Direct Integration":
            self.find_flux_for_line(self.current_line)
        else:  # "Gaussian Fit"
            # This would be a new method you'd need to implement
            self.fit_gaussian_for_line(self.current_line)

    def subtract_continuum_on_selected_line(self):
        self.subtract_continuum_for_line(self.current_line)

    def undo_continuum_on_selected_line(self):
        self.undo_continuum_for_line(self.current_line)


    def create_right_panel(self):
        # Create a table for object data
        self.objects_table = QTableWidget()
        self.objects_table.setColumnCount(11)
        self.objects_table.setHorizontalHeaderLabels([
            "Object Name", "Lya", "Lya_err", "OVI", "OVI_err", 
            "CIV", "CIV_err", "NV", "NV_err", "Object_csv", "Object_fits"
        ])
        
        # Make the table take up most of the right panel
        self.right_layout.addWidget(self.objects_table)
        
        # Connect the selection signal
        self.objects_table.itemSelectionChanged.connect(self.on_object_selection_change)
        
        # Create a horizontal layout for the buttons
        button_layout = QHBoxLayout()
        
        # Add the "load csv" button
        self.load_csv_button = QPushButton("Load CSV")
        self.load_csv_button.clicked.connect(self.load_objects_csv)
        button_layout.addWidget(self.load_csv_button)
        
        # Add the "save csv" button
        self.save_csv_button = QPushButton("Save CSV")
        self.save_csv_button.clicked.connect(self.save_objects_csv)
        button_layout.addWidget(self.save_csv_button)
        
        # Add the button layout to the right panel
        self.right_layout.addLayout(button_layout)
        
        self.create_detailed_controls()
        
        # Set the stretch factors to make the table take up most of the space
        self.right_layout.setStretch(0, 4)  # Table gets most space
        self.right_layout.setStretch(1, 1)  # Buttons get less space
        self.right_layout.setStretch(2, 2)  # Detailed controls get more space

    def load_spectrum_from_path(self, fits_path):
        """Load a spectrum directly from a fits file path"""
        try:
            if fits_path.endswith('.fits'):
                # Load the spectrum from the FITS file
                hdul = fits.open(fits_path)
                data = hdul[1].data
                self.coadded_spectrum = {
                    'wave': data['Wavelength'],
                    'flux': data['Flux']*self.scale_value,
                    'error_up': data['Error_Up']*self.scale_value,
                    'error_down': data['Error_Down']*self.scale_value,
                    'counts': data['gcounts'] # Fallback if not present
                }
                self.plot_spectrum()
                self.statusBar().showMessage(f'Spectrum loaded from {fits_path}')
                
                # Try to get redshift from header
                try:
                    self.redshift = hdul[1].header.get('REDSHIFT', self.redshift)
                    self.redshift_input.setText(str(self.redshift))
                except:
                    pass
                
            elif fits_path.endswith('.npz'):
                # Load the spectrum data from the NPZ file
                data = np.load(fits_path)

                # Restore the co-added spectrum and redshift
                self.coadded_spectrum = {
                    'wave': data['wave'],
                    'flux': data['flux']*self.scale_value,
                    'error_up': data['error_up']*self.scale_value,
                    'error_down': data['error_down']*self.scale_value,
                    'counts': data.get('counts', np.zeros_like(data['wave']))  # Fallback if not present
                }
                self.redshift = data['redshift']
                self.redshift_input.setText(str(self.redshift))

                self.plot_spectrum()
                self.statusBar().showMessage(f'Spectrum loaded from {fits_path}')
            else:
                self.statusBar().showMessage('Unsupported file format.')
        except Exception as e:
            self.statusBar().showMessage(f'Error loading spectrum: {e}')
    def load_csv_with_fits_paths(self):
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getOpenFileName(
            self, "Open CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        if csv_file:
            self.load_csv(csv_file)

    def load_objects_csv(self):
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getOpenFileName(
            self, "Open Objects CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        if csv_file:
            try:
                with open(csv_file, 'r') as file:
                    import csv
                    reader = csv.reader(file)
                    header = next(reader)  # Skip header row
                    
                    # Clear the table
                    self.objects_table.setRowCount(0)
                    
                    # Add rows to the table
                    for row_num, row_data in enumerate(reader):
                        self.objects_table.insertRow(row_num)
                        for col_num, data in enumerate(row_data):
                            item = QTableWidgetItem(data)
                            self.objects_table.setItem(row_num, col_num, item)
                    
                self.statusBar().showMessage(f'Loaded objects CSV: {csv_file}')
            except Exception as e:
                self.statusBar().showMessage(f'Error loading objects CSV: {e}')

    def create_detailed_controls(self):
        # Create a frame with a border for the detailed controls
        self.detailed_controls_frame = QFrame()
        self.detailed_controls_frame.setFrameShape(QFrame.StyledPanel)
        self.detailed_controls_frame.setFrameShadow(QFrame.Raised)
        self.detailed_controls_frame.setMinimumHeight(200)
    
        # Create a layout for the detailed controls
        detailed_layout = QVBoxLayout(self.detailed_controls_frame)
        
        # Add "Selected line:" label and value display
        line_selection_layout = QHBoxLayout()
        line_selection_label = QLabel("Selected line:")
        line_selection_layout.addWidget(line_selection_label)
        
        self.selected_line_display = QLabel("None")
        self.selected_line_display.setStyleSheet("font-weight: bold;")
        line_selection_layout.addWidget(self.selected_line_display)
        line_selection_layout.addStretch()
        
        detailed_layout.addLayout(line_selection_layout)
        
        # Add "Flux method:" label and dropdown
        method_layout = QHBoxLayout()
        method_label = QLabel("Flux method:")
        method_layout.addWidget(method_label)
        
        self.flux_method_dropdown = QComboBox()
        self.flux_method_dropdown.addItems(["Direct Integration", "Gaussian Fit"])
        self.flux_method_dropdown.currentIndexChanged.connect(self.on_flux_method_changed)
        method_layout.addWidget(self.flux_method_dropdown)
        method_layout.addStretch()
    
        detailed_layout.addLayout(method_layout)
        
        self.method_controls_container = QWidget()
        self.method_layout = QVBoxLayout(self.method_controls_container)

        self.create_direct_integration_controls()

        detailed_layout.addWidget(self.method_controls_container)
        detailed_layout.addStretch()

        self.right_layout.addWidget(self.detailed_controls_frame)

        self.action_button.setEnabled(False)
        self.subtract_continuum_button.setEnabled(False)
        self.undo_button.setEnabled(False)

    def create_direct_integration_controls(self):
        # Clear existing controls
        self.clear_method_controls()
        
        # Create the buttons for Direct Integration
        buttons_layout = QHBoxLayout()
        
        self.action_button = QPushButton("Find Flux and Error")
        self.action_button.clicked.connect(self.action_on_selected_line)
        buttons_layout.addWidget(self.action_button)
        self.subtract_continuum_button = QPushButton("Subtract Continuum")
        self.subtract_continuum_button.clicked.connect(self.subtract_continuum_on_selected_line)
        buttons_layout.addWidget(self.subtract_continuum_button)
        
        self.undo_button = QPushButton("Undo Continuum")
        self.undo_button.clicked.connect(self.undo_continuum_on_selected_line)
        buttons_layout.addWidget(self.undo_button)
        
        self.method_layout.addLayout(buttons_layout)
        
        # Initially disable the buttons until a line is selected
        self.action_button.setEnabled(False)
        self.subtract_continuum_button.setEnabled(False)
        self.undo_button.setEnabled(False)
    def create_gaussian_fit_controls(self):
        # Clear existing controls
        self.clear_method_controls()
        
        # Initialize storage for spectral lines and components
        if not hasattr(self, 'spectral_lines'):
            self.spectral_lines = []
        
        # Create top row of buttons
        buttons_layout = QHBoxLayout()
        
        self.fit_model_button = QPushButton("Fit Model")
        self.fit_model_button.clicked.connect(self.fit_gaussian_model)
        buttons_layout.addWidget(self.fit_model_button)
        
        self.add_line_button = QPushButton("Add Line")
        self.add_line_button.clicked.connect(self.add_gaussian_line)
        buttons_layout.addWidget(self.add_line_button)
        
        self.view_results_button = QPushButton("View Fit Results")
        self.view_results_button.clicked.connect(self.view_fit_results)
        buttons_layout.addWidget(self.view_results_button)
        
        self.method_layout.addLayout(buttons_layout)
        
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setMinimumHeight(150)
        self.scroll_area.setFrameShape(QFrame.StyledPanel)
        
        # Create a container widget for the scroll area
        self.line_entries_container = QWidget()
        self.line_entries_layout = QVBoxLayout(self.line_entries_container)
        self.line_entries_layout.setSpacing(5)
        self.line_entries_layout.addStretch()
        
        # Set the container as the scroll area's widget
        self.scroll_area.setWidget(self.line_entries_container)
        
        # Add the scroll area to the method layout
        self.method_layout.addWidget(self.scroll_area)
        
        # Load existing spectral lines if any
        self.load_spectral_lines()
        
        # Initially disable buttons until a line is selected
        self.fit_model_button.setEnabled(False)
        self.add_line_button.setEnabled(True)  # We should allow adding lines at any time
        self.view_results_button.setEnabled(False)

    def clear_method_controls(self):
        # Remove all widgets from the method layout
        while self.method_layout.count():
            item = self.method_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                # Recursively clear layouts
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
    # Modified save_objects_csv method to also save individual object CSVs
    def save_objects_csv(self):
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getSaveFileName(
            self, "Save Objects CSV", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        if csv_file:
            try:
                with open(csv_file, 'w', newline='') as file:
                    import csv
                    writer = csv.writer(file)
                    
                    # Write header
                    headers = []
                    for col in range(self.objects_table.columnCount()):
                        headers.append(self.objects_table.horizontalHeaderItem(col).text())
                    writer.writerow(headers)
                    
                    # Write data and update individual CSVs for each row
                    for row in range(self.objects_table.rowCount()):
                        row_data = []
                        for col in range(self.objects_table.columnCount()):
                            item = self.objects_table.item(row, col)
                            if item is not None:
                                row_data.append(item.text())
                            else:
                                row_data.append('')
                        writer.writerow(row_data)
                        
                        # Update individual object CSV for this row
                        self.update_object_csv_from_left_table(row)
                        
                    self.statusBar().showMessage(f'Saved objects CSV: {csv_file} and updated individual object CSVs')
                
            except Exception as e:
                self.statusBar().showMessage(f'Error saving objects CSV: {e}')
    def load_csv(self, csv_file_path):
        try:
            # Store the current CSV path
            self.current_csv_path = csv_file_path
            
            with open(csv_file_path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                lines = list(reader)

                # The first line should contain the redshift
                if lines:
                    self.redshift = float(lines[0][0])
                    self.confirm_redshift()

                # The remaining lines should contain FITS file paths and grating information
                self.fits_file_paths = []
                for row in lines[1:]:
                    print(f"Row: {row}")
                    if len(row) == 2:  # Ensure that each row contains exactly two items
                        self.fits_file_paths.append((row[0], row[1]))
                    else:
                        self.statusBar().showMessage(f'Invalid CSV format on line: {row}')
                        return

            self.statusBar().showMessage(f'Loaded {len(self.fits_file_paths)} FITS file paths from CSV')
            self.load_fits_files()
        except Exception as e:
            self.statusBar().showMessage(f'Error loading CSV: {e}')

    def update_redshift(self):
        try:
            new_redshift = float(self.redshift_input.text())
            old_redshift = self.redshift
            self.redshift = new_redshift

            if old_redshift != new_redshift and hasattr(self, 'spectral_lines'):
                for line in self.spectral_lines:
                    self.initialize_line_with_redshift(line, update_existing=True)

            self.statusBar().showMessage(f'Redshift updated to {self.redshift}')
            if self.coadded_spectrum:
                self.plot_spectrum()
        except ValueError:
            self.statusBar().showMessage('Invalid redshift value')
            self.redshift_input.setText(str(self.redshift))
    # Add method to handle flux method changes
    def on_flux_method_changed(self, index):
        selected_method = self.flux_method_dropdown.currentText()
        
        if selected_method == "Direct Integration":
            self.create_direct_integration_controls()
        elif selected_method == "Gaussian Fit":
            self.create_gaussian_fit_controls()
            self.init_gaussian_fit_mode()  # Initialize Gaussian fit mode
        
        # Update button states based on line selection
        if self.selected_line_display.text() != "None":
            self.update_button_states(True)
        else:
            self.update_button_states(False)
    # Add method to update button states
    def update_button_states(self, enabled):
        selected_method = self.flux_method_dropdown.currentText()
            
        if selected_method == "Direct Integration":
            if hasattr(self, 'action_button'):
                self.action_button.setEnabled(enabled)
            if hasattr(self, 'subtract_continuum_button'):
                self.subtract_continuum_button.setEnabled(enabled)
            if hasattr(self, 'undo_button'):
                self.undo_button.setEnabled(enabled)
        elif selected_method == "Gaussian Fit":
            if hasattr(self, 'fit_model_button'):
                self.fit_model_button.setEnabled(enabled)
            if hasattr(self, 'add_line_button'):
                self.add_line_button.setEnabled(enabled)
            if hasattr(self, 'view_results_button'):
                self.view_results_button.setEnabled(enabled)
    def update_scale(self):
        try:
            new_scale = float(self.scale_input.text())
            self.scale_value = new_scale
            self.statusBar().showMessage(f'Scale value updated to {self.scale_value:.2e}')
            self.scale_input.setText(f"{self.scale_value:.2e}")  # Update with scientific notation
            if self.coadded_spectrum:
                self.plot_spectrum()
        except ValueError:
            self.statusBar().showMessage('Invalid scale value')
            self.scale_input.setText(f"{self.scale_value:.2e}")  # Revert with scientific notation

    def create_spectral_line_entry(self, name, rest_wavelength, doublet_info=None, geocoronal=False):
        """Create a visual entry for a spectral line and add to the UI"""
        print(f'CREATING SPECTRAL LINE ENTRY: {name}, {rest_wavelength}, {doublet_info}, {geocoronal}')
        # Create a SpectralLine object
        components = ["narrow"]  # Default component
        
        # Create DoubletInfo object if needed
        doublet = None
        if doublet_info:
            doublet = DoubletInfo(
                secondary_wavelength=doublet_info["secondary_wavelength"],
                ratio=doublet_info["ratio"]
            )
        
        # Create the SpectralLine object
        spectral_line = SpectralLine(
            rest_wavelength=rest_wavelength,
            name=name,
            components=components,
            doublet=doublet,
            geocoronal=geocoronal
        )
        
        # Store the spectral line
        self.spectral_lines.append(spectral_line)
        self.initialize_line_with_redshift(spectral_line)
        for comp in spectral_line.components:
            param_key = f"{name}_{comp}"
            if param_key not in self.component_parameters:
                self.component_parameters[param_key] = {}
            self.component_parameters[param_key]['z_value'] = self.redshift
            self.component_parameters[param_key]['z_min'] = self.redshift - 0.001
            self.component_parameters[param_key]['z_max'] = self.redshift + 0.001
            self.component_parameters[param_key]['sigma_value'] = 200
            self.component_parameters[param_key]['sigma_min'] = 10
            self.component_parameters[param_key]['sigma_max'] = 1000
            self.component_parameters[param_key]['flux_value'] = 100
        self.save_component_parameters()
        
        # Create visual frame
        line_entry = QFrame()
        line_entry.setFrameShape(QFrame.StyledPanel)
        line_entry.setFrameShadow(QFrame.Raised)
        line_entry.setStyleSheet("background-color: #f0f0f0;")
        line_entry.setMinimumHeight(40)
        line_entry.setMaximumHeight(40)
        
        entry_layout = QHBoxLayout(line_entry)
        entry_layout.setContentsMargins(10, 5, 10, 5)
        
        # Display format: Name - Wavelength
        if doublet:
            display_text = f"{name} ({rest_wavelength:.2f}, {doublet.secondary_wavelength:.2f} Å)"
        else:
            display_text = f"{name} ({rest_wavelength:.2f} Å)"
            
        line_label = QLabel(display_text)
        entry_layout.addWidget(line_label)
        
        entry_layout.addStretch()
        
        # Store a reference to the spectral line object in the widget
        line_entry.spectral_line = spectral_line
        
        view_params_button = QPushButton("Parameters")
        view_params_button.clicked.connect(lambda: self.view_line_parameters(spectral_line))
        entry_layout.addWidget(view_params_button)
        
        # Remove the stretch from the end of line_entries_layout
        if self.line_entries_layout.count() > 0:
            item = self.line_entries_layout.itemAt(self.line_entries_layout.count() - 1)
            if item and item.spacerItem():
                self.line_entries_layout.removeItem(item)
        
        # Add the line entry
        self.line_entries_layout.addWidget(line_entry)
        self.line_entries_layout.addStretch()
        
        # Enable the fit button if we have at least one line
        self.fit_model_button.setEnabled(True)
        
        return spectral_line

    def save_spectral_lines(self):
        """Save the spectral lines to a CSV file specific to the current object"""
        print(f'SAVING SPECTRAL LINES: {self.spectral_lines}')
        
        unique_lines = []
        seen_names = set()
        for line in self.spectral_lines:
            if line.name not in seen_names:
                unique_lines.append(line)
                seen_names.add(line.name)
        self.spectral_lines = unique_lines
        print(f'SAVING SPECTRAL LINES: {self.spectral_lines}')

        os.makedirs("internal", exist_ok=True)
        
        # Get the current object name
        object_name = self.get_current_object_name()
        if not object_name:
            # Fall back to default if no object is selected
            csv_path = os.path.join("internal", "spectral_lines.csv")
        else:
            csv_path = os.path.join("internal", f"{object_name}_spectral_lines.csv")
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow(["name", "rest_wavelength", "components", "doublet_wavelength", "doublet_ratio", "geocoronal"])
            
            # Write each line
            for line in self.spectral_lines:
                # Components are stored as a comma-separated list
                components_str = ",".join(line.components)
                
                # Doublet info
                doublet_wavelength = ""
                doublet_ratio = ""
                if line.doublet:
                    doublet_wavelength = line.doublet.secondary_wavelength
                    doublet_ratio = line.doublet.ratio
                    
                writer.writerow([
                    line.name,
                    line.rest_wavelength,
                    components_str,
                    doublet_wavelength,
                    doublet_ratio,
                    line.geocoronal
                ])
        
        self.statusBar().showMessage(f"Saved {len(self.spectral_lines)} spectral lines to {csv_path}")

    def save_component_parameters(self):
        """Save component parameters to a CSV file specific to the current object"""

        
        os.makedirs("internal", exist_ok=True)
        
        # Get the current object name
        object_name = self.get_current_object_name()
        if not object_name:
            # Fall back to default if no object is selected
            csv_path = os.path.join("internal", "component_parameters.csv")
        else:
            csv_path = os.path.join("internal", f"{object_name}_component_parameters.csv")
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow([
                "key", "z_value", "z_min", "z_max", 
                "sigma_value", "sigma_min", "sigma_max", 
                "flux_value", "z_vary"  # Added z_vary
            ])
            
            # Write each component's parameters
            for key, params in self.component_parameters.items():
                writer.writerow([
                    key,
                    params.get('z_value', self.redshift),
                    params.get('z_min', max(0, self.redshift - 0.05)),
                    params.get('z_max', self.redshift + 0.05),
                    params.get('sigma_value', 200),
                    params.get('sigma_min', 10),
                    params.get('sigma_max', 1000),
                    params.get('flux_value', 100),
                    params.get('z_vary', True)  # Added z_vary
                ])
        
        self.statusBar().showMessage(f"Saved parameters for {len(self.component_parameters)} components to {csv_path}")
    def load_spectral_lines(self):
        """Load spectral lines from object-specific CSV file"""
        
        print(f'LOADING SPECTRAL LINES')
        # Get the current object name
        object_name = self.get_current_object_name()
        if not object_name:
            # Fall back to default if no object is selected
            csv_path = os.path.join("internal", "spectral_lines.csv")
        else:
            csv_path = os.path.join("internal", f"{object_name}_spectral_lines.csv")
        
        self.spectral_lines = []
        self.clear_line_entries()
        if not os.path.exists(csv_path):
            # No saved lines yet for this object
            return
        
        # Load from CSV
        try:
            with open(csv_path, 'r') as f:
                reader = csv.reader(f)
                
                # Skip header
                next(reader, None)
                
                for row in reader:
                    if len(row) >= 6:
                        name = row[0]
                        rest_wavelength = float(row[1])
                        components = row[2].split(",") if row[2] else ["narrow"]
                        
                        # Parse doublet info
                        doublet_info = None
                        if row[3] and row[4]:
                            doublet_info = {
                                "secondary_wavelength": float(row[3]),
                                "ratio": float(row[4])
                            }
                        
                        # Parse geocoronal flag
                        geocoronal = row[5].lower() == "true"
                    
                        
                        # Create DoubletInfo object if needed
                        doublet = None
                        if doublet_info:
                            doublet = DoubletInfo(
                                secondary_wavelength=doublet_info["secondary_wavelength"],
                                ratio=doublet_info["ratio"]
                            )
                        
                        # Create the SpectralLine object
                        spectral_line = SpectralLine(
                            rest_wavelength=rest_wavelength,
                            name=name,
                            components=components,
                            doublet=doublet,
                            geocoronal=geocoronal
                        )
                        
                        # Add to the UI
                        self.create_line_entry_from_spectral_line(spectral_line)
            
            self.statusBar().showMessage(f"Loaded {len(self.spectral_lines)} spectral lines from {csv_path}")
        except Exception as e:
            print(f"Error loading spectral lines: {e}")
            import traceback
            print(traceback.format_exc())

    def load_component_parameters(self):
        """Load component parameters from the object-specific CSV file"""
        
        # Get the current object name
        object_name = self.get_current_object_name()
        if not object_name:
            # Fall back to default if no object is selected
            csv_path = os.path.join("internal", "component_parameters.csv")
        else:
            csv_path = os.path.join("internal", f"{object_name}_component_parameters.csv")
        
        if not os.path.exists(csv_path):
            # No saved parameters yet for this object
            self.component_parameters = {}
            return
        
        # Initialize parameters dictionary
        self.component_parameters = {}
        
        # Load from CSV
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            
            # Skip header
            next(reader, None)
            
            for row in reader:
                if len(row) >= 8:  # Check we have at least 8 columns (more if z_vary is included)
                    key = row[0]
                    print(f'HERE IS ROW REDSHIFT: {row[1]}, HERE IS FLOAT {float(row[1])}')
                    print(f'AND THIS IS OBJECT REDSHIFT {self.redshift}')
                    self.component_parameters[key] = {
                        'z_value': float(row[1]),
                        'z_min': float(row[2]),
                        'z_max': float(row[3]),
                        'sigma_value': float(row[4]),
                        'sigma_min': float(row[5]),
                        'sigma_max': float(row[6]),
                        'flux_value': float(row[7])
                    }
                    #print(f'PARAMETER SAVED: {self.component_parameters['z_value']}')
                    # Add z_vary if available (in newer files)
                    if len(row) >= 9:
                        self.component_parameters[key]['z_vary'] = row[8].lower() == 'true'
        
        self.statusBar().showMessage(f"Loaded parameters for {len(self.component_parameters)} components from {csv_path}")
    def view_line_parameters(self, spectral_line):
        """Open a dialog to view and edit a line's parameters and components"""
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Parameters for {spectral_line.name}")
        dialog.setMinimumWidth(500)
        main_layout = QVBoxLayout(dialog)
        
        # Create tabs for Line Properties and Components
        tab_widget = QTabWidget()
        
        # Line Properties Tab
        line_tab = QWidget()
        line_layout = QFormLayout(line_tab)
        
        # Basic line properties
        name_input = QLineEdit(spectral_line.name)
        line_layout.addRow("Name:", name_input)
        
        wavelength_input = QDoubleSpinBox()
        wavelength_input.setRange(1000, 10000)
        wavelength_input.setValue(spectral_line.rest_wavelength)
        wavelength_input.setDecimals(2)
        line_layout.addRow("Rest Wavelength (Å):", wavelength_input)
        
        # Doublet properties
        is_doublet = QCheckBox("Is a doublet")
        is_doublet.setChecked(spectral_line.doublet is not None)
        line_layout.addRow(is_doublet)
        
        secondary_wavelength = QDoubleSpinBox()
        secondary_wavelength.setRange(1000, 10000)
        secondary_wavelength.setValue(spectral_line.doublet.secondary_wavelength if spectral_line.doublet else spectral_line.rest_wavelength + 5)
        secondary_wavelength.setDecimals(2)
        secondary_wavelength.setEnabled(spectral_line.doublet is not None)
        line_layout.addRow("Secondary Wavelength (Å):", secondary_wavelength)
        
        ratio_input = QDoubleSpinBox()
        ratio_input.setRange(0.01, 10)
        ratio_input.setValue(spectral_line.doublet.ratio if spectral_line.doublet else 0.5)
        ratio_input.setDecimals(3)
        ratio_input.setEnabled(spectral_line.doublet is not None)
        line_layout.addRow("Flux Ratio:", ratio_input)
        
        # Connect doublet checkbox
        is_doublet.stateChanged.connect(lambda state: secondary_wavelength.setEnabled(state))
        is_doublet.stateChanged.connect(lambda state: ratio_input.setEnabled(state))
        
        # Geocoronal option
        is_geocoronal = QCheckBox("Is a geocoronal line")
        is_geocoronal.setChecked(spectral_line.geocoronal)
        line_layout.addRow(is_geocoronal)
        
        tab_widget.addTab(line_tab, "Line Properties")
        
        # Components Tab
        components_tab = QWidget()
        components_layout = QVBoxLayout(components_tab)
        
        # Component list
        component_list = QListWidget()
        for comp in spectral_line.components:
            component_list.addItem(comp)
        components_layout.addWidget(component_list)
        
        # Buttons for add/remove components
        components_btn_layout = QHBoxLayout()
        
        add_comp_btn = QPushButton("Add Component")
        components_btn_layout.addWidget(add_comp_btn)
        
        remove_comp_btn = QPushButton("Remove Selected")
        components_btn_layout.addWidget(remove_comp_btn)
        
        components_layout.addLayout(components_btn_layout)
        
        # Component parameter group
        comp_params_group = QGroupBox("Component Parameters")
        comp_params_layout = QFormLayout(comp_params_group)
        
        comp_name_input = QLineEdit()
        comp_params_layout.addRow("Name:", comp_name_input)
        
        # Initial redshift from main GUI
        redshift_input = QDoubleSpinBox()
        # redshift_input.setRange(0, 10)
        redshift_input.setDecimals(7)
        redshift_input.setValue(self.redshift)
        
        comp_params_layout.addRow("Initial Redshift:", redshift_input)
        
        # Redshift constraints
        redshift_min = QDoubleSpinBox()
        # redshift_min.setRange(0, 10)
        redshift_min.setDecimals(7)
        redshift_min.setValue(max(0, self.redshift - 0.05))
        comp_params_layout.addRow("Min Redshift:", redshift_min)
        
        redshift_max = QDoubleSpinBox()
        redshift_max.setDecimals(7)
        redshift_max.setValue(self.redshift + 0.05)
        comp_params_layout.addRow("Max Redshift:", redshift_max)
        
        # Sigma (velocity dispersion)
        sigma_input = QDoubleSpinBox()
        sigma_input.setRange(0, 10000)
        sigma_input.setDecimals(10)
        sigma_input.setValue(300) 
        comp_params_layout.addRow("Initial Sigma (km/s):", sigma_input)
        
        # Sigma constraints
        sigma_min = QDoubleSpinBox()
        sigma_min.setRange(0, 10000)
        sigma_min.setDecimals(10)
        sigma_min.setValue(10)
        comp_params_layout.addRow("Min Sigma:", sigma_min)
        
        sigma_max = QDoubleSpinBox()
        sigma_max.setRange(0, 10000)
        sigma_max.setDecimals(8)
        sigma_max.setValue(1000)
        comp_params_layout.addRow("Max Sigma:", sigma_max)
        
        # Initial flux guess
        flux_input = QDoubleSpinBox()
        flux_input.setRange(0, 1e15)
        flux_input.setDecimals(10)
        flux_input.setValue(10)
        comp_params_layout.addRow("Initial Flux:", flux_input)
        
        # Add to layout
        components_layout.addWidget(comp_params_group)
        
        # Add/Save component button
        save_comp_btn = QPushButton("Add/Update Component")
        components_layout.addWidget(save_comp_btn)
        
        # Initialize component parameters dictionary if not exists
        if not hasattr(self, 'component_parameters'):
            self.component_parameters = {}
        
        # Handler for component selection
        def on_component_selected():
            if not component_list.currentItem():
                return
                
            comp_name = component_list.currentItem().text()
            comp_name_input.setText(comp_name)
            
            # Load saved parameters if they exist
            param_key = f"{spectral_line.name}_{comp_name}"
            if param_key in self.component_parameters:
                params = self.component_parameters[param_key]
                redshift_input.setValue(params.get('z_value', self.redshift))
                redshift_min.setValue(params.get('z_min', max(0, self.redshift - 0.05)))
                redshift_max.setValue(params.get('z_max', self.redshift + 0.05))
                sigma_input.setValue(params.get('sigma_value', 200))
                sigma_min.setValue(params.get('sigma_min', 10))
                sigma_max.setValue(params.get('sigma_max', 1000))
                flux_input.setValue(params.get('flux_value', 100))
            
        
        # Connect component selection
        component_list.itemClicked.connect(lambda: on_component_selected())
        
        # Handler for adding a new component
        def add_component():
            comp_name = f"component_{component_list.count() + 1}"
            component_list.addItem(comp_name)
            component_list.setCurrentRow(component_list.count() - 1)
            comp_name_input.setText(comp_name)
            # Add to spectral line components
            spectral_line.components.append(comp_name)
        
        # Handler for removing a component
        def remove_component():
            if not component_list.currentItem():
                return
                
            row = component_list.currentRow()
            comp_name = component_list.currentItem().text()
            
            # Remove from components list and UI
            component_list.takeItem(row)
            if comp_name in spectral_line.components:
                spectral_line.components.remove(comp_name)
                
            # Also remove parameters
            param_key = f"{spectral_line.name}_{comp_name}"
            if param_key in self.component_parameters:
                del self.component_parameters[param_key]
        
        # Handler for saving component parameters
        def save_component():
            if not component_list.currentItem():
                return
                
            old_comp_name = component_list.currentItem().text()
            new_comp_name = comp_name_input.text()
            
            # Update component name if changed
            if old_comp_name != new_comp_name:
                # Update in spectral_line.components
                idx = spectral_line.components.index(old_comp_name)
                spectral_line.components[idx] = new_comp_name
                
                # Update in list widget
                component_list.currentItem().setText(new_comp_name)
                
                # Transfer parameters to new name if they exist
                old_param_key = f"{spectral_line.name}_{old_comp_name}"
                if old_param_key in self.component_parameters:
                    self.component_parameters[f"{spectral_line.name}_{new_comp_name}"] = \
                        self.component_parameters[old_param_key]
                    del self.component_parameters[old_param_key]
            
            # Save component parameters
            param_key = f"{spectral_line.name}_{new_comp_name}"
            self.component_parameters[param_key] = {
                'z_value': redshift_input.value(),
                'z_min': redshift_min.value(),
                'z_max': redshift_max.value(),
                'sigma_value': sigma_input.value(),
                'sigma_min': sigma_min.value(),
                'sigma_max': sigma_max.value(),
                'flux_value': flux_input.value()
            }
        
        # Connect buttons
        add_comp_btn.clicked.connect(add_component)
        remove_comp_btn.clicked.connect(remove_component)
        save_comp_btn.clicked.connect(save_component)
        
        tab_widget.addTab(components_tab, "Components")
        
        # Add tabs to main layout
        main_layout.addWidget(tab_widget)
        
        # Dialog buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        main_layout.addWidget(button_box)
        
        # Show the dialog
        if dialog.exec_() == QDialog.Accepted:
            # Update line properties
            spectral_line.name = name_input.text()
            spectral_line.rest_wavelength = wavelength_input.value()
            
            # Update doublet info
            if is_doublet.isChecked():
                from gaussian_fitfunctions import DoubletInfo
                spectral_line.doublet = DoubletInfo(
                    secondary_wavelength=secondary_wavelength.value(),
                    ratio=ratio_input.value()
                )
            else:
                spectral_line.doublet = None
                
            # Update geocoronal flag
            spectral_line.geocoronal = is_geocoronal.isChecked()
            
            # Save updated spectral lines
            self.save_spectral_lines()
            
            # Update the UI to reflect changes
            self.refresh_line_entries()
            
            # Save component parameters to file
            self.save_component_parameters()


    def view_fit_results(self):
        """Display the fitting results with flux measurements"""
        if not hasattr(self, 'fit_results'):
            self.statusBar().showMessage("No fit results available.")
            return
        
        # Get the fit results
        fitter = self.fit_results['fitter']
        mcmc_result = self.fit_results['mcmc_result']

        print(f'In this function, now attempting to open Gaussian Window')
        
        new_window = GaussianWindow(fitter, mcmc_result, self.spectral_lines)
        new_window.show()
        # Create dialog
        # dialog = QDialog(self)
        # dialog.setWindowTitle("Fit Results")
        # dialog.setMinimumSize(600, 500)
        # layout = QVBoxLayout(dialog)
        
        # # Create tabs
        # tabs = QTabWidget()
        
        # # Flux results tab
        # flux_tab = QWidget()
        # flux_layout = QVBoxLayout(flux_tab)
        
        # # Create table for flux results
        # flux_table = QTableWidget()
        # flux_table.setColumnCount(6)
        # flux_table.setHorizontalHeaderLabels([
        #     "Line", "Component", "Flux", "Error (Low)", "Error (High)", "EW (Å)"
        # ])
        
        # # Analyze and add results for each line
        # row_count = 0
        # for line in self.spectral_lines:
        #     # Calculate line flux
        #     line_result = fitter.get_line_flux(mcmc_result, line.name)
            
        #     if not line_result:
        #         continue
            
        #     # Add total flux row
        #     flux_table.insertRow(row_count)
        #     flux_table.setItem(row_count, 0, QTableWidgetItem(line.name))
        #     flux_table.setItem(row_count, 1, QTableWidgetItem("Total"))
        #     flux_table.setItem(row_count, 2, QTableWidgetItem(f"{line_result['flux']:.3e}"))
        #     flux_table.setItem(row_count, 3, QTableWidgetItem(f"{line_result['error_down']:.3e}"))
        #     flux_table.setItem(row_count, 4, QTableWidgetItem(f"{line_result['error_up']:.3e}"))
        #     flux_table.setItem(row_count, 5, QTableWidgetItem("N/A"))
        #     row_count += 1
            
        #     # Add component rows
        #     for comp, comp_result in line_result['components'].items():
        #         flux_table.insertRow(row_count)
        #         flux_table.setItem(row_count, 0, QTableWidgetItem(""))
        #         flux_table.setItem(row_count, 1, QTableWidgetItem(comp))
        #         flux_table.setItem(row_count, 2, QTableWidgetItem(f"{comp_result['flux']:.3e}"))
        #         flux_table.setItem(row_count, 3, QTableWidgetItem(f"{comp_result['error_down']:.3e}"))
        #         flux_table.setItem(row_count, 4, QTableWidgetItem(f"{comp_result['error_up']:.3e}"))
        #         flux_table.setItem(row_count, 5, QTableWidgetItem("N/A"))
        #         row_count += 1
        
        # # Adjust table layout
        # flux_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # flux_table.horizontalHeader().setStretchLastSection(True)
        
        # flux_layout.addWidget(flux_table)
        
        # # Add export button
        # export_btn_layout = QHBoxLayout()
        # export_flux_btn = QPushButton("Export Flux Results")
        # export_btn_layout.addWidget(export_flux_btn)
        # flux_layout.addLayout(export_btn_layout)
        
        # # Component parameters tab
        # params_tab = QWidget()
        # params_layout = QVBoxLayout(params_tab)
        
        # # Create table for parameter results
        # params_table = QTableWidget()
        # params_table.setColumnCount(4)
        # params_table.setHorizontalHeaderLabels([
        #     "Component", "Parameter", "Value", "Error"
        # ])
        
        # # Get unique components
        # all_components = set()
        # for line in self.spectral_lines:
        #     all_components.update(line.components)
        
        # # Add rows for each component parameter
        # row_count = 0
        # for comp in sorted(all_components):
        #     # Get redshift
        #     z_key = f"z_{comp}"
        #     if z_key in mcmc_result.var_names:
        #         z_percentiles = np.percentile(mcmc_result.flatchain[z_key], [16, 50, 84])
        #         z_value = z_percentiles[1]
        #         z_error = (z_percentiles[2] - z_percentiles[0]) / 2
                
        #         params_table.insertRow(row_count)
        #         params_table.setItem(row_count, 0, QTableWidgetItem(comp))
        #         params_table.setItem(row_count, 1, QTableWidgetItem("Redshift"))
        #         params_table.setItem(row_count, 2, QTableWidgetItem(f"{z_value:.5f}"))
        #         params_table.setItem(row_count, 3, QTableWidgetItem(f"±{z_error:.5f}"))
        #         row_count += 1
            
        #     # Get sigma
        #     sigma_key = f"sigma_{comp}"
        #     if sigma_key in mcmc_result.var_names:
        #         sigma_percentiles = np.percentile(mcmc_result.flatchain[sigma_key], [16, 50, 84])
        #         sigma_value = sigma_percentiles[1]
        #         sigma_error = (sigma_percentiles[2] - sigma_percentiles[0]) / 2
                
        #         params_table.insertRow(row_count)
        #         params_table.setItem(row_count, 0, QTableWidgetItem(comp))
        #         params_table.setItem(row_count, 1, QTableWidgetItem("Sigma (km/s)"))
        #         params_table.setItem(row_count, 2, QTableWidgetItem(f"{sigma_value:.1f}"))
        #         params_table.setItem(row_count, 3, QTableWidgetItem(f"±{sigma_error:.1f}"))
        #         row_count += 1
        
        # # Adjust table layout
        # params_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # params_table.horizontalHeader().setStretchLastSection(True)
        
        # params_layout.addWidget(params_table)
        
        # # Export parameters button
        # export_params_btn = QPushButton("Export Parameter Results")
        # params_layout.addWidget(export_params_btn)
        
        # # Add tabs to tab widget
        # tabs.addTab(flux_tab, "Flux Results")
        # tabs.addTab(params_tab, "Component Parameters")
        
        # # Add tab widget to dialog
        # layout.addWidget(tabs)
        
        # # Add close button
        # close_btn = QPushButton("Close")
        # close_btn.clicked.connect(dialog.accept)
        # layout.addWidget(close_btn)
        
        # # Define export functions
        # def export_flux_results():
        #     file_path, _ = QFileDialog.getSaveFileName(
        #         dialog, "Save Flux Results", "", "CSV Files (*.csv);;All Files (*)"
        #     )
            
        #     if file_path:
        #         import csv
        #         with open(file_path, 'w', newline='') as f:
        #             writer = csv.writer(f)
        #             writer.writerow(["Line", "Component", "Flux", "Error (Low)", "Error (High)"])
                    
        #             for row in range(flux_table.rowCount()):
        #                 line_name = flux_table.item(row, 0).text()
        #                 component = flux_table.item(row, 1).text()
        #                 flux_val = flux_table.item(row, 2).text()
        #                 error_low = flux_table.item(row, 3).text()
        #                 error_high = flux_table.item(row, 4).text()
                        
        #                 writer.writerow([line_name, component, flux_val, error_low, error_high])
                
        #         self.statusBar().showMessage(f"Flux results exported to {file_path}")
        
        # def export_parameter_results():
        #     file_path, _ = QFileDialog.getSaveFileName(
        #         dialog, "Save Parameter Results", "", "CSV Files (*.csv);;All Files (*)"
        #     )
            
        #     if file_path:
        #         import csv
        #         with open(file_path, 'w', newline='') as f:
        #             writer = csv.writer(f)
        #             writer.writerow(["Component", "Parameter", "Value", "Error"])
                    
        #             for row in range(params_table.rowCount()):
        #                 component = params_table.item(row, 0).text()
        #                 parameter = params_table.item(row, 1).text()
        #                 value = params_table.item(row, 2).text()
        #                 error = params_table.item(row, 3).text()
                        
        #                 writer.writerow([component, parameter, value, error])

    def update_flux_results_in_table(self, mcmc_result):
        """Update the main table with flux results from the fit"""
        if not hasattr(self, 'fit_results'):
            return
        
        fitter = self.fit_results['fitter']
        
        # Process each line in the table
        for row in range(self.table.rowCount()):
            line_item = self.table.item(row, 0)
            if not line_item:
                continue
                
            line_name = line_item.text()
            
            # Find matching spectral line
            matching_line = None
            for line in self.spectral_lines:
                if line.name == line_name:
                    matching_line = line
                    break
            
            if not matching_line:
                continue
            
            # Get flux results for this line
            line_result = fitter.get_line_flux(mcmc_result, line_name)
            
            if not line_result:
                continue
            
            # Update flux and error in the table
            flux_value = line_result['flux']
            flux_error = max(line_result['error_up'], line_result['error_down'])
            
            self.table.setItem(row, 1, QTableWidgetItem(f"{flux_value:.3e}"))
            self.table.setItem(row, 2, QTableWidgetItem(f"{flux_error:.3e}"))
            
            # If we have parameters for this line, update them too
            # For simplicity, assume we've saved component parameters
            for comp in matching_line.components:
                z_key = f"z_{comp}"
                sigma_key = f"sigma_{comp}"
                
                if z_key in mcmc_result.var_names and sigma_key in mcmc_result.var_names:
                    z_value = np.median(mcmc_result.flatchain[z_key])
                    sigma_value = np.median(mcmc_result.flatchain[sigma_key])
                    
                    param_key = f"{line_name}_{comp}"
                    if not hasattr(self, 'component_parameters'):
                        self.component_parameters = {}
                    
                    self.component_parameters[param_key] = {
                        'z_value': z_value,
                        'z_min': max(0, z_value - 0.05),
                        'z_max': z_value + 0.05,
                        'sigma_value': sigma_value,
                        'sigma_min': 10,
                        'sigma_max': sigma_value * 2,
                        'flux_value': line_result['components'].get(comp, {}).get('flux', 100)
                    }
                    
            # Save the updated parameters
            self.save_component_parameters()
        
        # Update the objects table if a row is selected
        selected_rows = self.objects_table.selectionModel().selectedRows()
        if selected_rows:
            row_index = selected_rows[0].row()
            self.update_object_row_from_line_values(row_index)
    def plot_fit_results(self):
        """Plot the fitting results on the existing plot"""
        if not hasattr(self, 'fit_results'):
            return
        
        # clear other plots
        self.clear_fit_result_plot()

        legend_labels = []
        legend_handles = []
        # Get the fit results
        fitter = self.fit_results['fitter']
        mcmc_result = self.fit_results['mcmc_result']
        region = self.fit_results['region']
        
        # Get wavelength, flux, and error arrays
        wave = self.coadded_spectrum['wave']
        flux = self.coadded_spectrum['flux']
        error = (self.coadded_spectrum['error_up'] + self.coadded_spectrum['error_down']) / 2
        
        # Get Gaussian models
        models_data = fitter.get_gaussian_models(
            mcmc_result,
            wave=wave,
            use_mcmc_median=True
        )
        
        # Plot data on existing plot
        # Make sure we're using the current axis (create if needed)
        if not hasattr(self, 'ax') or self.ax is None:
            self.ax = self.figure.add_subplot(111)
        
        # Clear the plot if we want to start fresh
        # self.ax.clear()
        
        # Plot the original data if it's not already plotted
        # self.ax.plot(wave, flux, color='black', drawstyle='steps-mid', label='Data')
        
        # Process each region
        for region_key, region_data in models_data.items():
            region_wave = region_data['wave']
            # Plot total model
            total_model = self.ax.plot(region_wave, region_data['total_model'], 'r-', 
                    label='Model Fit', lw=2)
            self.all_plotted_models.append(total_model)
            legend_labels.append("Total Model")
            legend_handles.append(total_model[0])
            # Plot individual components
            
            colors = plt.cm.tab10(np.linspace(0, 1, len(region_data['components'])))
            for (comp_name, comp_data), color in zip(region_data['components'].items(), colors):
                plotted = self.ax.plot(region_wave, comp_data['model'] + region_data['continuum'],
                        '--', color=color, alpha=0.7)
                self.all_plotted_models.append(plotted)
                legend_labels.append(f"{comp_name}")
                legend_handles.append(plotted[0])
            continuum_model = self.ax.plot(region_wave, region_data['continuum'], 'k--', color='gray', alpha = 0.7)
            self.all_plotted_models.append(continuum_model)
            legend_labels.append("Continuum")
            legend_handles.append(continuum_model[0])
        # Set appropriate axis limits
        # Get the region limits from the first (and likely only) region
        first_region = list(models_data.values())[0]
        region_wave = first_region['wave']
        
        # Set x-limits to the region wavelength range (with some padding)
        x_min, x_max = np.min(region_wave), np.max(region_wave)
        x_padding = (x_max - x_min) * 0.05  # 5% padding
        self.ax.set_xlim([x_min - x_padding, x_max + x_padding])
        
        # Set appropriate y-limits
        total_model = first_region['total_model']
        y_max = np.max(total_model[~np.isnan(total_model)]) * 1.2  # 20% headroom
        y_min = np.min([0, np.min(total_model[~np.isnan(total_model)]) * 1.2])
        self.ax.set_ylim([y_min, y_max])
        
        # Add labels and legend
        self.ax.set_xlabel('Wavelength (Å)')
        self.ax.set_ylabel(f'Flux (scaled by {self.scale_value:.2e})')
        self.ax.set_title("Spectral Fit")
        
        # Update legend - get existing handles and labels, then add new ones
        handles, labels = self.ax.get_legend_handles_labels()
        self.ax.legend(labels = legend_labels, handles = legend_handles, loc='best')
        
        # Redraw the canvas
        self.canvas.draw()

    def clear_fit_result_plot(self):
        for line in self.ax.lines:
            if line in self.all_plotted_models:
                print(f"Removing line: {line}")
                line.remove()
        self.all_plotted_models = []
            
    def fit_gaussian_model(self):
        """Perform fitting using the JointSpectralFitter class"""
        print("\n=== Starting Gaussian Fitting Process ===")
        
        if not self.spectral_lines:
            self.statusBar().showMessage("No spectral lines to fit. Please add lines first.")
            print("Error: No spectral lines found")
            return
            
        if not self.coadded_spectrum:
            self.statusBar().showMessage("No spectrum loaded. Please load a spectrum first.")
            print("Error: No spectrum loaded")
            return
        
        # Load component parameters if not already loaded
        if not hasattr(self, 'component_parameters'):
            print("Loading component parameters")
            self.load_component_parameters()
        
        print(f"Current redshift: {self.redshift}")
        print(f"Number of spectral lines: {len(self.spectral_lines)}")
        for i, line in enumerate(self.spectral_lines):
            print(f"Line {i+1}: {line.name} at rest wavelength {line.rest_wavelength} Å")
            print(f"  Components: {line.components}")
            if line.doublet:
                print(f"  Doublet with secondary wavelength: {line.doublet.secondary_wavelength} Å, ratio: {line.doublet.ratio}")
        
        x_min, x_max = self.ax.get_xlim()
        print(f'Plot wavelength range: {x_min:.2f} - {x_max:.2f} Å')
        
        # Initialize min_wave and max_wave with the plot limits
        min_wave = x_min
        max_wave = x_max
        
        print("Observed wavelengths of lines:")
        for line in self.spectral_lines:
            # Calculate the observed wavelength based on redshift
            obs_wavelength = line.rest_wavelength * (1 + self.redshift)
            print(f"  {line.name}: {obs_wavelength:.2f} Å")
            
            # Also handle doublet if present
            if line.doublet:
                obs_doublet = line.doublet.secondary_wavelength * (1 + self.redshift)
                print(f"  {line.name} doublet: {obs_doublet:.2f} Å")
        
        # Check if any lines are within the plot range
        lines_in_range = []
        for line in self.spectral_lines:
            obs_wavelength = line.rest_wavelength * (1 + self.redshift)
            if x_min <= obs_wavelength <= x_max:
                lines_in_range.append(line.name)
                
            # Check doublet as well
            if line.doublet:
                obs_doublet = line.doublet.secondary_wavelength * (1 + self.redshift)
                if x_min <= obs_doublet <= x_max and line.name not in lines_in_range:
                    lines_in_range.append(f"{line.name} doublet")
        
        if not lines_in_range:
            print("Warning: No spectral lines appear to be within the current plot range!")
            print(f"Plot range: {x_min:.2f} - {x_max:.2f} Å")
        else:
            print(f"Lines in plot range: {', '.join(lines_in_range)}")
        
        print("\n=== Creating Spectral Region ===")
        # Create a single region with all lines
        try:
            from gaussian_fitfunctions import SpectralRegion
            region = SpectralRegion(
                wave_min=x_min,
                wave_max=x_max,
                lines=self.spectral_lines
            )
            print(f"Created SpectralRegion with range: {region.wave_min:.2f} - {region.wave_max:.2f} Å")
        except Exception as e:
            print(f"Error creating SpectralRegion: {str(e)}")
            import traceback
            print(traceback.format_exc())
            self.statusBar().showMessage(f"Error creating region: {str(e)}")
            return
        
        print("\n=== Creating Fitter ===")
        # Create the JointSpectralFitter
        try:
            from gaussian_fitfunctions import JointSpectralFitter
            fitter = JointSpectralFitter(
                z=self.redshift,
                regions=[region]
            )
            print(f"Created JointSpectralFitter with redshift: {self.redshift}")
        except Exception as e:
            print(f"Error creating JointSpectralFitter: {str(e)}")
            import traceback
            print(traceback.format_exc())
            self.statusBar().showMessage(f"Error creating fitter: {str(e)}")
            return
        
        print("\n=== Preparing Parameter Constraints ===")
        # Prepare parameter constraints based on saved component parameters
        parameter_constraints = {}
        
        # Process each line and component
        print("Component parameters being used:")
        for line in self.spectral_lines:
            for comp in line.components:
                # Get parameter key
                param_key = f"{line.name}_{comp}"
                
                # Set constraints for this component if saved
                if param_key in self.component_parameters:
                    params = self.component_parameters[param_key]
                    print(f"Parameters for {param_key}: {params}")
                    
                    # Redshift constraints
                    z_name = f"z_{comp}"
                    if 'z_value' in params:
                        print(f"  Using saved z_value: {params['z_value']}")
                    else:
                        print(f"  No saved z_value, using redshift: {self.redshift}")
                        
                    parameter_constraints[z_name] = {
                        'value': params['z_value'] if 'z_value' in params else self.redshift,
                        'min': params['z_min'] if 'z_min' in params else max(0, self.redshift - 0.05),
                        'max': params['z_max'] if 'z_max' in params else self.redshift + 0.05,
                        'vary': params.get('z_vary', True)  # Use get with default
                    }
                    print(f"  Set {z_name} constraints: {parameter_constraints[z_name]}")
                    
                    # Sigma constraints
                    sigma_name = f"sigma_{comp}"
                    parameter_constraints[sigma_name] = {
                        'value': params['sigma_value'] if 'sigma_value' in params else 200,
                        'min': params['sigma_min'] if 'sigma_min' in params else 10,
                        'max': params['sigma_max'] if 'sigma_max' in params else 1000,
                        'vary': True
                    }
                    print(f"  Set {sigma_name} constraints: {parameter_constraints[sigma_name]}")
                    
                    # Flux constraints
                    flux_name = f"flux_{line.name}_{comp}"
                    parameter_constraints[flux_name] = {
                        'value': params['flux_value'] if 'flux_value' in params else 100,
                        'min': 0,
                        'vary': True
                    }
                    print(f"  Set {flux_name} constraints: {parameter_constraints[flux_name]}")
                else:
                    print(f"No saved parameters for {param_key}, using defaults")
                    
                    # Create default constraints
                    z_name = f"z_{comp}"
                    sigma_name = f"sigma_{comp}"
                    flux_name = f"flux_{line.name}_{comp}"
                    
                    parameter_constraints[z_name] = {
                        'value': self.redshift,
                        'min': max(0, self.redshift - 0.05),
                        'max': self.redshift + 0.05,
                        'vary': True
                    }
                    
                    parameter_constraints[sigma_name] = {
                        'value': 200,
                        'min': 10,
                        'max': 1000,
                        'vary': True
                    }
                    
                    parameter_constraints[flux_name] = {
                        'value': 100,
                        'min': 0,
                        'vary': True
                    }
                    print(f"  Created default constraints for {param_key}")
        
        # Set up continuum parameters for the region
        parameter_constraints[f'm_0'] = {'value': 0, 'vary': True}
        parameter_constraints[f'c_0'] = {'value': 0, 'vary': True}
        print(f"Set continuum parameters: m_0={parameter_constraints['m_0']}, c_0={parameter_constraints['c_0']}")
        
        # Get wavelength, flux, and error arrays from the loaded spectrum
        print("\n=== Preparing Spectrum Data ===")
        wave = self.coadded_spectrum['wave']
        flux = self.coadded_spectrum['flux']
        error = (self.coadded_spectrum['error_up'] + self.coadded_spectrum['error_down']) / 2
        
        # Check for NaN or inf values
        wave_nans = np.isnan(wave).sum()
        flux_nans = np.isnan(flux).sum()
        error_nans = np.isnan(error).sum()
        
        if wave_nans > 0 or flux_nans > 0 or error_nans > 0:
            print(f"Warning: Found NaN values in data: wave={wave_nans}, flux={flux_nans}, error={error_nans}")
        
        # Check for negative or zero error values
        zero_errors = (error <= 0).sum()
        if zero_errors > 0:
            print(f"Warning: Found {zero_errors} non-positive error values")
            # placeholder -- something is wrong with errors
            error[error <= 0] = abs(0.1*flux[error <= 0])
            # error[error < = 0] = 
        # Calculate data statistics to verify reasonable values
        print(f"Wavelength range: {np.min(wave):.2f} - {np.max(wave):.2f} Å")
        print(f"Flux range: {np.min(flux):.2e} - {np.max(flux):.2e}")
        print(f"Error range: {np.min(error[error > 0]):.2e} - {np.max(error):.2e}")
        
        # Check that wavelength range includes the fit region
        if np.min(wave) > x_min or np.max(wave) < x_max:
            print(f"Warning: Fit region ({x_min:.2f} - {x_max:.2f} Å) extends beyond available data range ({np.min(wave):.2f} - {np.max(wave):.2f} Å)")
        
        # Show a status message
        self.statusBar().showMessage("Performing initial fit...")
        print("\n=== Starting Initial Fit ===")
        
        # Perform initial fit
        try:
            print("Calling fit_regions with:")
            print(f"  wave shape: {wave.shape}")
            print(f"  flux shape: {flux.shape}")
            print(f"  error shape: {error.shape}")
            print(f"  Number of parameter constraints: {len(parameter_constraints)}")
            
            initial_fit = fitter.fit_regions(
                wave=wave,
                flux=flux,
                error=error,
                parameter_constraints=parameter_constraints
            )
            
            print("\nInitial fit results:")
            print(f"Success: {initial_fit.success}")
            print(f"Number of function evaluations: {initial_fit.nfev}")
            print(f"Number of variables: {initial_fit.nvarys}")
            print(f"Message: {initial_fit.message}")
            print(f"Reduced chi-squared: {initial_fit.redchi:.2f}")
            
            print("\nBest fit parameter values:")
            for name, param in initial_fit.params.items():
                print(f"  {name}: {param.value:.5g} ± {param.stderr if param.stderr is not None else 'N/A'}")
            
            self.statusBar().showMessage("Initial fit complete. Running MCMC...")
            print("\n=== Starting MCMC Fit ===")
            
            # Perform MCMC fit
            mcmc_result = fitter.mcmc_joint_fit(
                wave=wave,
                flux=flux,
                error=error,
                initial_fit=initial_fit,
                steps=1000,  # Adjust as needed
                burn=100,
                thin=15,
                parameter_constraints=parameter_constraints
            )
            
            print("\nMCMC completed.")
            if hasattr(mcmc_result, 'flatchain'):
                print(f"MCMC chain shape: {mcmc_result.flatchain.shape}")
                print(f"Number of parameters: {len(mcmc_result.var_names)}")
                print(f"Parameter names: {mcmc_result.var_names}")
            else:
                print("Warning: MCMC result does not have a flatchain attribute")
            
            # Store the results
            self.fit_results = {
                'initial_fit': initial_fit,
                'mcmc_result': mcmc_result,
                'fitter': fitter,
                'region': region
            }
            
            print("\n=== Plotting Results ===")
            # Plot the results
            self.plot_fit_results()
            
            # Enable the View Results button
            self.view_results_button.setEnabled(True)
            
            self.statusBar().showMessage("Fitting complete!")
            
        except Exception as e:
            import traceback
            traceback_str = traceback.format_exc()
            print("\n=== ERROR DURING FITTING ===")
            print(traceback_str)
            print(f"Error message: {str(e)}")
            
            # Additional debugging for common errors
            if "Cannot handle mix of" in str(e):
                print("\nData type error detected. Checking data types:")
                print(f"wave dtype: {wave.dtype}")
                print(f"flux dtype: {flux.dtype}")
                print(f"error dtype: {error.dtype}")
                
            elif "out of bounds" in str(e):
                print("\nOut of bounds error detected. This might be related to parameter constraints.")
                print("Check if your initial parameter values are within the min/max bounds.")
                
            elif "division by zero" in str(e) or "divide by zero" in str(e):
                print("\nDivision by zero detected. Check for zero values in error array.")
                print(f"Number of zero values in error array: {(error == 0).sum()}")
                
            self.statusBar().showMessage(f"Error during fitting: {str(e)}")

    def refresh_line_entries(self):
        """Refresh all line entries in the UI"""
        # Save current lines
        current_lines = self.spectral_lines.copy()
        
        # Clear existing entries
        self.clear_line_entries()
        
        # Re-create entries for each line
        for line in current_lines:
            self.create_line_entry_from_spectral_line(line)

    def get_current_object_name(self):
        """Get the name of the currently selected object"""
        selected_rows = self.objects_table.selectionModel().selectedRows()
        if selected_rows:
            row_index = selected_rows[0].row()
            object_item = self.objects_table.item(row_index, 0)
            if object_item and object_item.text():
                return object_item.text()
        return None
    def create_line_entry_from_spectral_line(self, spectral_line):
        """Create a UI entry from an existing SpectralLine object"""
        
        # Check if the spectral line already exists
        for existing_line in self.spectral_lines:
            if existing_line.name == spectral_line.name:
                print(f'spectral line already exists, got caught here')
                # return
            
        # But... if the line does not exist yet, add it
        self.spectral_lines.append(spectral_line)
        
        # Skip creating UI elements if we're not in Gaussian Fit mode
        if not hasattr(self, 'line_entries_layout'):
            print(f"Note: Added spectral line {spectral_line.name} but skipping UI update (not in Gaussian Fit mode)")
            return
        
        # Create visual frame
        line_entry = QFrame()
        line_entry.setFrameShape(QFrame.StyledPanel)
        line_entry.setFrameShadow(QFrame.Raised)
        line_entry.setStyleSheet("background-color: #f0f0f0;")
        line_entry.setMinimumHeight(40)
        line_entry.setMaximumHeight(40)
        
        entry_layout = QHBoxLayout(line_entry)
        entry_layout.setContentsMargins(10, 5, 10, 5)
        
        # Display format: Name - Wavelength
        if spectral_line.doublet:
            display_text = f"{spectral_line.name} ({spectral_line.rest_wavelength:.2f}, {spectral_line.doublet.secondary_wavelength:.2f} Å)"
        else:
            display_text = f"{spectral_line.name} ({spectral_line.rest_wavelength:.2f} Å)"
            
        line_label = QLabel(display_text)
        entry_layout.addWidget(line_label)
        
        entry_layout.addStretch()
        
        # Store a reference to the spectral line object in the widget
        line_entry.spectral_line = spectral_line
        
        view_params_button = QPushButton("Parameters")
        view_params_button.clicked.connect(lambda: self.view_line_parameters(spectral_line))
        entry_layout.addWidget(view_params_button)

        include_line_button = QCheckBox("Include in Fit?")
        include_line_button.setChecked(True)
        include_line_button.stateChanged.connect(lambda state: spectral_line.set_include_in_fit(state))
        entry_layout.addWidget(include_line_button)
        
        
        # Remove the stretch from the end of line_entries_layout
        if self.line_entries_layout.count() > 0:
            item = self.line_entries_layout.itemAt(self.line_entries_layout.count() - 1)
            if item and item.spacerItem():
                self.line_entries_layout.removeItem(item)
        
        # Add the line entry
        self.line_entries_layout.addWidget(line_entry)
        self.line_entries_layout.addStretch()
        
        # Enable the fit button if we have at least one line
        if hasattr(self, 'fit_model_button'):
            self.fit_model_button.setEnabled(True)
    def add_gaussian_line(self):
        """Open a dialog to add a new spectral line for fitting with predefined options"""
        # Create directory for internal storage if it doesn't exist
        os.makedirs("internal", exist_ok=True)
        print(f'ADD GAUSSIAN LINE FUNCTION')
        
        # Define common spectral lines with their properties
        predefined_lines = {
            "Lya": {
                "wavelength": 1215.67,
                "is_doublet": False
            },
            "NV": {
                "wavelength": 1238.82,
                "is_doublet": True,
                "secondary_wavelength": 1242.80,
                "ratio": 2.95
            },
            "SiIV": {
                "wavelength": 1393.75,
                "is_doublet": True,
                "secondary_wavelength": 1402.77,
                "ratio": 2.95
            },
            "CIV": {
                "wavelength": 1548.19,
                "is_doublet": True,
                "secondary_wavelength": 1550.77,
                "ratio": 2.95
            },
            "OVI": {
                "wavelength": 1031.92,
                "is_doublet": True,
                "secondary_wavelength": 1037.61,
                "ratio": 2.95
            },
            "HeII": {
                "wavelength": 1640.42,
                "is_doublet": False
            },
            "geocoronal_Lya": {
                "wavelength": 1215.67,
                "is_doublet": False,
                "geocoronal": True
            }
        }
        
        # Create the dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Spectral Line")
        main_layout = QVBoxLayout(dialog)
        
        # Create a group for line selection method
        selection_group = QGroupBox("Line Selection Method")
        selection_layout = QVBoxLayout(selection_group)
        
        # Radio buttons for selection method
        predefined_radio = QRadioButton("Choose from common lines")
        predefined_radio.setChecked(True)  # Default to predefined
        custom_radio = QRadioButton("Enter custom line")
        
        selection_layout.addWidget(predefined_radio)
        selection_layout.addWidget(custom_radio)
        
        main_layout.addWidget(selection_group)
        
        # Create form layouts for both options
        predefined_form = QFormLayout()
        custom_form = QFormLayout()
        
        # Predefined lines form
        line_dropdown = QComboBox()
        line_dropdown.addItems(list(predefined_lines.keys()))
        line_dropdown.currentIndexChanged.connect(lambda: update_line_details(line_dropdown.currentText()))
        predefined_form.addRow("Select Line:", line_dropdown)
        
        # Custom line form
        name_input = QLineEdit()
        custom_form.addRow("Line Name:", name_input)
        
        wavelength_input = QDoubleSpinBox()
        wavelength_input.setRange(1000, 10000)
        wavelength_input.setValue(1215.67)  # Default to Lyman-alpha
        wavelength_input.setDecimals(2)
        custom_form.addRow("Rest Wavelength (Å):", wavelength_input)
        
        # Common elements for both forms
        is_doublet = QCheckBox("Is this a doublet?")
        
        secondary_wavelength = QDoubleSpinBox()
        secondary_wavelength.setRange(1000, 10000)
        secondary_wavelength.setValue(1215.67 + 5)  # Default slightly higher
        secondary_wavelength.setDecimals(2)
        secondary_wavelength.setEnabled(False)
        
        ratio_input = QDoubleSpinBox()
        ratio_input.setRange(0.01, 10)
        ratio_input.setValue(0.5)  # Default ratio
        ratio_input.setDecimals(3)
        ratio_input.setEnabled(False)
        
        is_doublet.stateChanged.connect(lambda state: secondary_wavelength.setEnabled(state))
        is_doublet.stateChanged.connect(lambda state: ratio_input.setEnabled(state))
        
        is_geocoronal = QCheckBox("Is this a geocoronal line?")
        
        # Add common elements to both forms
        for form in [predefined_form, custom_form]:
            form.addRow(is_doublet)
            form.addRow("Secondary Wavelength (Å):", secondary_wavelength)
            form.addRow("Flux Ratio (secondary/primary):", ratio_input)
            form.addRow(is_geocoronal)
            form.addRow(QLabel("Default component will be added automatically"))
        
        # Create containers for each form
        predefined_container = QGroupBox("Predefined Line Options")
        predefined_container.setLayout(predefined_form)
        
        custom_container = QGroupBox("Custom Line Options")
        custom_container.setLayout(custom_form)
        
        # Initially show only predefined
        main_layout.addWidget(predefined_container)
        main_layout.addWidget(custom_container)
        custom_container.setVisible(False)
        
        # Connect radio buttons to show/hide options
        predefined_radio.toggled.connect(lambda checked: predefined_container.setVisible(checked))
        custom_radio.toggled.connect(lambda checked: custom_container.setVisible(checked))
        
        # Function to update form when selecting predefined line
        def update_line_details(line_name):
            if line_name in predefined_lines:
                line_data = predefined_lines[line_name]
                # Update doublet status
                is_doublet.setChecked(line_data["is_doublet"])
                
                if line_data["is_doublet"]:
                    secondary_wavelength.setValue(line_data["secondary_wavelength"])
                    ratio_input.setValue(line_data["ratio"])
                    
        # Initialize with the first selection
        update_line_details(line_dropdown.currentText())
        
        # Dialog buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        main_layout.addWidget(button_box)
        
        # Show the dialog and process result
        if dialog.exec_() == QDialog.Accepted:
            # Determine if using predefined or custom
            if predefined_radio.isChecked():
                line_name = line_dropdown.currentText()
                line_data = predefined_lines[line_name]
                rest_wavelength = line_data["wavelength"]
            else:
                line_name = name_input.text()
                rest_wavelength = wavelength_input.value()
            
            # Create doublet info if needed
            doublet_info = None
            if is_doublet.isChecked():
                doublet_info = {
                    "secondary_wavelength": secondary_wavelength.value(),
                    "ratio": ratio_input.value()
                }
            
            # Create the line entry
            self.create_spectral_line_entry(
                line_name, 
                rest_wavelength, 
                doublet_info, 
                is_geocoronal.isChecked()
            )
            print(f'ADDED SPECTRAL LINE {line_name}, NOW SAVING')
            
            # Save to CSV
            self.save_spectral_lines()
    def clear_line_entries(self):
        if not hasattr(self, 'line_entries_layout'):
            return
        if hasattr(self, 'line_entries_layout'):
            while self.line_entries_layout.count():
                item = self.line_entries_layout.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
                elif item.spacerItem():
                    self.line_entries_layout.removeItem(item)
            
            # Add back the stretch
            self.line_entries_layout.addStretch()

    def clear_method_controls(self):
        # Clear the line entries if we're switching from Gaussian fit mode
        if hasattr(self, 'line_entries_layout'):
            self.clear_line_entries()
        
        # Remove all widgets from the method layout
        while self.method_layout.count():
            item = self.method_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                # Recursively clear layouts
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
    # def view_line_parameters(self, line_name):
    #     QMessageBox.information(self, "Line Parameters", f"Parameters for {line_name} - placeholder")

    def save_spectrum(self):
        if self.coadded_spectrum is None:
            self.statusBar().showMessage('No co-added spectrum to reference.')
            return

        # Get the save file path
        options = QFileDialog.Options()
        save_file, _ = QFileDialog.getSaveFileName(
            self, "Select Reference Path for Spectrum", "", "FITS Files (*.fits);;NPZ Files (*.npz);;All Files (*)", options=options
        )

        if save_file:
            try:
                # Update the objects table with the saved spectrum reference path
                self.update_objects_table_with_spectrum_path(save_file)
                self.statusBar().showMessage(f'Spectrum path recorded: {save_file}')
            except Exception as e:
                self.statusBar().showMessage(f'Error recording spectrum path: {e}')

    def update_objects_table_with_spectrum_path(self, fits_path):
        """Update the objects table with the path to the spectrum"""
        # Check if there's a selected row
        selected_rows = self.objects_table.selectionModel().selectedRows()
        
        if selected_rows:
            # Update the existing row
            row_index = selected_rows[0].row()
            self.objects_table.setItem(row_index, 10, QTableWidgetItem(fits_path))
            
            # Update the flux values from the left table
            self.update_object_row_from_line_values(row_index)
            
            # Save the object name from column 0 if it exists
            object_name_item = self.objects_table.item(row_index, 0)
            if object_name_item and object_name_item.text():
                object_name = object_name_item.text()
            else:
                # Generate a default name based on the file path
                object_name = os.path.basename(fits_path).split('.')[0]
                self.objects_table.setItem(row_index, 0, QTableWidgetItem(object_name))
        else:
            # Add a new row - we need to create a new row with the spectrum path
            self.add_new_object_row(fits_path)

    def add_new_object_row(self, fits_path):
        """Add a new row to the objects table with the spectrum path"""
        # Get current row count
        row_count = self.objects_table.rowCount()
        # Insert a new row
        self.objects_table.insertRow(row_count)
        
        # Generate a default object name from the file path
        object_name = os.path.basename(fits_path).split('.')[0]
        self.objects_table.setItem(row_count, 0, QTableWidgetItem(object_name))
        
        # Set the spectrum path in the Object_fits column (index 10)
        self.objects_table.setItem(row_count, 10, QTableWidgetItem(fits_path))
        
        # Set the current CSV path in the Object_csv column (index 9) if we have one loaded
        if hasattr(self, 'current_csv_path') and self.current_csv_path:
            self.objects_table.setItem(row_count, 9, QTableWidgetItem(self.current_csv_path))
        
        # Update the flux values from the left table
        self.update_object_row_from_line_values(row_count)
        
        # Select the newly added row
        self.objects_table.selectRow(row_count)

    # Add this method to the SpectralFluxApp class in fit_UV_spectrum.py

    def update_object_row_from_line_values(self, row_index):
        """Update the flux values in the objects table from the left table with doublet summing"""
        # Dictionary to track doublet components for summing
        doublet_sums = {
            "Lya": {"flux": 0.0, "error": 0.0, "count": 0},
            "OVI": {"flux": 0.0, "error": 0.0, "count": 0},
            "NV": {"flux": 0.0, "error": 0.0, "count": 0},
            "Si": {"flux": 0.0, "error": 0.0, "count": 0},
            "CIV": {"flux": 0.0, "error": 0.0, "count": 0}
        }
        
        # Line name mapping to doublet category
        line_to_doublet = {
            "Lya": "Lya",
            "O VI 1031": "OVI",
            "O VI 1037": "OVI",
            "N V 1238": "NV",
            "N V 1242": "NV",
            "Si 1393": "Si",
            "Si 1402": "Si",
            "C IV 1548": "CIV",
            "C IV 1550": "CIV"
        }
    
        # Mapping from doublet categories to object table columns
        doublet_to_columns = {
            "Lya": (1, 2),    # (Flux column, Error column)
            "OVI": (3, 4),
            "NV": (7, 8),
            "CIV": (5, 6),
            "Si": None        # Not included in right table, can be added if needed
        }
    
        # Process each row in the left table
        for left_row in range(self.table.rowCount()):
            # Get the line name from the first column
            line_item = self.table.item(left_row, 0)
            if not line_item or not line_item.text():
                continue
                
            line_name = line_item.text()
            # Check if this line belongs to a doublet category
            if line_name in line_to_doublet:
                doublet_category = line_to_doublet[line_name]
                
                # Get flux and error values if they exist
                flux_item = self.table.item(left_row, 1)
                error_item = self.table.item(left_row, 2)
                
                if flux_item and flux_item.text() and error_item and error_item.text():
                    try:
                        flux_value = float(flux_item.text())
                        error_value = float(error_item.text())
                        
                        # Add to the appropriate doublet sum
                        doublet_sums[doublet_category]["flux"] += flux_value
                        # Adding errors in quadrature for better statistical handling
                        doublet_sums[doublet_category]["error"] += error_value**2
                        doublet_sums[doublet_category]["count"] += 1
                    except ValueError:
                        # Skip if we can't convert to float
                        pass
    
        # Now update the object table with the summed values
        for doublet, sum_data in doublet_sums.items():
            if sum_data["count"] > 0 and doublet_to_columns.get(doublet):
                flux_col, error_col = doublet_to_columns[doublet]
                
                # Handle single line case (like Lya)
                if doublet == "Lya":
                    flux_value = sum_data["flux"]
                    error_value = sum_data["error"] if sum_data["count"] == 1 else math.sqrt(sum_data["error"])
                else:
                    # For doublets, use the sum of flux measurements
                    flux_value = sum_data["flux"]
                    # Take the square root of summed squared errors (error propagation)
                    error_value = math.sqrt(sum_data["error"])
                
                # Update the values in the objects table
                self.objects_table.setItem(row_index, flux_col, QTableWidgetItem(f"{flux_value:.2e}"))
                self.objects_table.setItem(row_index, error_col, QTableWidgetItem(f"{error_value:.2e}"))
        

    def load_spectrum(self):
        # Get the load file path
        options = QFileDialog.Options()
        load_file, _ = QFileDialog.getOpenFileName(
            self, "Load Spectrum", "", "FITS Files (*.fits);;NPZ Files (*.npz);;All Files (*)", options=options
        )

        if load_file:
            try:
                if load_file.endswith('.npz'):
                    # Load the spectrum data from the NPZ file
                    data = np.load(load_file)

                    # Restore the co-added spectrum and redshift
                    self.coadded_spectrum = {
                        'wave': data['wave'],
                        'flux': data['flux']*self.scale_value,
                        'error_up': data['error_up']*self.scale_value,
                        'error_down': data['error_down']*self.scale_value
                    }
                    self.redshift = data['redshift']
                    self.redshift_input.setText(str(self.redshift))
                    

                    self.plot_spectrum()
                    self.statusBar().showMessage(f'Spectrum loaded from {load_file}')
                elif load_file.endswith('.fits'):
                    # Load the spectrum from the FITS file
                    hdul = fits.open(load_file)
                    data = hdul[1].data
                    self.coadded_spectrum = {
                        'wave': data['Wavelength'],
                        'flux': data['Flux']*self.scale_value,
                        'error_up': data['Error_Up']*self.scale_value,
                        'error_down': data['Error_Down']*self.scale_value,
                        'counts': data['gcounts']
                    }
                    self.redshift = hdul[1].header['REDSHIFT']
                    self.redshift_input.setText(str(self.redshift))

                    self.plot_spectrum()
                    self.statusBar().showMessage(f'Spectrum loaded from {load_file}')
                else:
                    self.statusBar().showMessage('Unsupported file format.')
            except Exception as e:
                self.statusBar().showMessage(f'Error loading spectrum: {e}')


    def confirm_redshift(self):
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle('Redshift Confirmation')
        msg_box.setText(f'Is this your redshift? {self.redshift}')
        msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        result = msg_box.exec_()

        if result == QMessageBox.No:
            redshift, ok = QInputDialog.getDouble(
                self, "Enter Redshift", "Redshift:", self.redshift, -10.0, 10.0, 3
            )
            if ok:
                self.redshift = redshift
                self.redshift_input.setText(str(self.redshift))
                if self.coadded_spectrum:
                    self.plot_spectrum()
        else:
            self.redshift_input.setText(str(self.redshift))

    def load_fits_files(self):
        self.spectrum_data = [self.load_and_prepare_fits(file_info) for file_info in self.fits_file_paths]
        if self.spectrum_data:
            self.plot_lines_button.setEnabled(True)

    def load_and_prepare_fits(self, file_info):
        try:
            file_name, grating = file_info  # Unpack file path and grating from the tuple
            self.grating_data.append(grating)
            if grating == 'G140L':
                # Use the standard prepare function for G140L grating
                return prepare(file_name)
            else:
                # Use prepare_other_grating for non-G140L gratings
                # try:
                if grating not in self.grating_data:
                    self.grating_data.append(grating)
                return prepare_other_grating(file_name, grating=grating)
                # except Exception as e:
                #     self.statusBar().showMessage(f'Error prepping non-G140L FITS file {file_name}: {e}')
                #     print((f'Error prepping non-G140L FITS file {file_name}: {e}'))
        except Exception as e:
            self.statusBar().showMessage(f'Error loading FITS file {file_name}: {e}')
            return None

    def subtract_continuum_for_line(self, line_index):
        """Subtract the selected continuum for the given line and update the plot."""
        # Retrieve the continuum bounds from the table
        table_row = line_index - 1
        left_lower = self.table.item(table_row, 3).text()
        left_upper = self.table.item(table_row, 4).text()
        right_lower = self.table.item(table_row, 5).text()
        right_upper = self.table.item(table_row, 6).text()
        slope = self.table.item(table_row, 7).text()
        intercept = self.table.item(table_row, 8).text()

        if not (left_lower and left_upper and right_lower and right_upper):
            self.statusBar().showMessage('Continuum bounds are not set for this line.')
            return

        try:
            slope = float(slope)
            intercept = float(intercept)

            # Subtract the continuum using the bounds from the table and plot the new spectrum
            self.apply_continuum_subtraction(slope, intercept)
            self.statusBar().showMessage(f'Continuum subtracted for line: {self.line_labels[line_index]}')
        except ValueError:
            self.statusBar().showMessage('Invalid continuum bounds input.')


    def apply_continuum_subtraction(self, slope, intercept):
        """Perform the continuum subtraction and update the plot."""
        if self.coadded_spectrum is None:
            self.statusBar().showMessage('No co-added spectrum to subtract continuum.')
            return

        wave = self.coadded_spectrum['wave']
        flux = self.coadded_spectrum['flux']

        # # Create masks for the continuum regions
        # left_mask = (wave >= left_lower) & (wave <= left_upper)
        # right_mask = (wave >= right_lower) & (wave <= right_upper)

        # # Calculate the continuum levels
        # left_flux = np.mean(flux[left_mask])
        # right_flux = np.mean(flux[right_mask])

        # # Linear continuum between the two regions
        # continuum = np.interp(wave, [left_lower, right_upper], [left_flux, right_flux])
        continuum = slope * wave + intercept
        # Subtract the continuum from the flux
        self.flux_after_subtraction = flux - continuum

        # Replot the spectrum with the continuum subtracted
        # self.figure.clear()
        # ax = self.figure.add_subplot(111)
        self.ax.plot(wave, self.flux_after_subtraction, color='blue')
        self.ax.set_xlabel('Wavelength')
        self.ax.set_ylabel(f'Flux (erg / s / cm^2 / A ) (scaled by {self.scale_value})')
        self.canvas.draw()

    def undo_continuum_for_line(self, line_index):
        """Undo the continuum subtraction and restore the original spectrum."""
        if self.coadded_spectrum is None:
            self.statusBar().showMessage('No co-added spectrum to undo.')
            return

        #self.statusBar().showMessage(f'Continuum subtraction undone for line: {self.line_labels[line_index]}')
        self.statusBar().showMessage(f'This is buggy... cut for now. Sorry!!')


    def on_coadd_complete(self, coadded_spectrum):
        if coadded_spectrum is None:
            self.statusBar().showMessage('Co-addition failed')
            self.progress_bar.setValue(0)
            return

        self.coadded_spectrum = coadded_spectrum
        self.plot_spectrum()
        self.statusBar().showMessage('Co-addition and plotting complete')
        self.progress_bar.setValue(0)

    def plot_spectrum(self):
        if self.coadded_spectrum:
            wave = self.coadded_spectrum['wave']
            flux = self.coadded_spectrum['flux']
            error_up = self.coadded_spectrum['error_up']
            error_down = self.coadded_spectrum['error_down']

            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            ax = self.ax
            ax.plot(wave, flux,  color='black', drawstyle='steps-mid')
            ax.plot(wave, (error_up + error_down)/2, color='grey', drawstyle='steps-mid')
            ax.set_xlabel('Wavelength')
            ax.set_xlim(1100, 1880)
            # print(np.max(flux[~np.isnan(flux)]))
            # ax.set_ylim(-0.5, np.max(flux[~np.isnan(flux)]))
            ax.set_ylabel(f'Flux (erg / s / cm^2 / A ) (scaled by {self.scale_value:0.2e})')
            self.plot_expected_lines(ax)
            self.canvas.draw()
        else:
            self.statusBar().showMessage('No co-added spectrum to plot')

    def clear_all_lines(self):
        try:
            if len(self.ax.lines) > 0:
                for line in self.ax.lines:
                    line.remove()
        except Exception:
            print(f'no axis defined yet/no lines to remove')
            
            
    def plot_expected_lines(self, ax=None):
        try:
            redshift = float(self.redshift_input.text()) if self.redshift_input.text() else self.redshift
            observed_lines = self.line_wavelengths * (1 + redshift)
            colors = ['r', 'g', 'b', 'orange', 'purple', 'pink', 'yellow', 'darkblue', 'black']

            if ax is None:
                ax = self.figure.gca()

            for i, line in enumerate(observed_lines[1:]):  # Skip 'Full Spectrum' placeholder
                thisline = ax.axvline(x=line, linestyle='--', color=colors[i % len(colors)])
                self.drawn_spectral_lines.append(thisline)
                ax.text(line, ax.get_ylim()[1] * 0.5, self.line_labels[i+1],
                        color=colors[i % len(colors)], rotation=90, verticalalignment='bottom')
        except ValueError:
            self.statusBar().showMessage('Invalid redshift input. Please enter a valid number.')

    def zoom_to_line(self):
        if self.coadded_spectrum:
            selected_line = self.line_dropdown.currentIndex()
            rest_frame_wavelength = self.line_wavelengths[selected_line] * (1 + self.redshift)
            wave = self.coadded_spectrum['wave']
            flux = self.coadded_spectrum['flux']
            error_up = self.coadded_spectrum['error_up']
            error_down = self.coadded_spectrum['error_down']

            if self.line_labels[selected_line] == 'Full Spectrum':
                # Plot the full spectrum
                # self.figure.clear()
                # self.ax = self.figure.add_subplot(111)
                ax = self.ax
                ax.plot(wave, flux, color='black', drawstyle='steps-mid')
                ax.plot(wave, (error_up + error_down)/2, color='grey', drawstyle='steps-mid')
                ax.set_xlabel('Wavelength')
                ax.set_ylabel(f'Flux (erg / s / cm^2 / A ) (scaled by {self.scale_value:0.2e})')
                self.plot_expected_lines(ax)
                self.canvas.draw()
            else:
                zoom_range = (rest_frame_wavelength - 50, rest_frame_wavelength + 50)
                mask = (wave >= zoom_range[0]) & (wave <= zoom_range[1])
                wave_zoom = wave[mask]
                flux_zoom = flux[mask]

                # self.figure.clear()
                # self.ax = self.figure.add_subplot(111)
                ax = self.ax
                ax.plot(wave, flux, color='black', drawstyle='steps-mid')
                ax.plot(wave, (error_up + error_down)/2, color='grey', drawstyle='steps-mid')
                ax.set_xlim(zoom_range)
                ax.set_ylim(0, np.max(flux_zoom[~np.isnan(flux_zoom)]) * 1.1)
                ax.set_xlabel('Wavelength')
                ax.set_ylabel(f'Flux (erg / s / cm^2 / A ) (scaled by {self.scale_value:0.2e})')
                self.plot_expected_lines(ax)
                self.canvas.draw()
        else:
            self.statusBar().showMessage('No co-added spectrum to zoom into')

    def find_flux_for_line(self, line_index):
        self.current_line = line_index
        self.selection_step = 0
        self.statusBar().showMessage(f'Select left continuum start for {self.line_labels[line_index]}')

        # Clear previous selections from the plot
        self.clear_selection_overlays()

        # Connect event handler
        if self.cid is None:
            self.cid = self.canvas.mpl_connect('button_press_event', self.on_click)

    def clear_selection_overlays(self):
        # Remove existing lines and patches
        try:
            for line in self.left_continuum_lines + self.right_continuum_lines + self.integration_lines:
                line.remove()
            self.left_continuum_lines = []
            self.right_continuum_lines = []
            self.integration_lines = []
        except Exception:
            pass
        try:
            if self.left_continuum_patch:
                self.left_continuum_patch.remove()
                self.left_continuum_patch = None
            if self.right_continuum_patch:
                self.right_continuum_patch.remove()
                self.right_continuum_patch = None
            if self.integration_patch:
                self.integration_patch.remove()
                self.integration_patch = None
        except Exception:
            pass
        self.canvas.draw()

    def on_click(self, event):
        if event.inaxes is not None:
            x = event.xdata
            if self.selection_step == 0:
                self.left_continuum_start = x
                self.statusBar().showMessage(f'Select left continuum end for {self.line_labels[self.current_line]}')
            elif self.selection_step == 1:
                self.left_continuum_end = x
                self.plot_continuum_region('left')
                self.statusBar().showMessage(f'Select right continuum start for {self.line_labels[self.current_line]}')
            elif self.selection_step == 2:
                self.right_continuum_start = x
                self.statusBar().showMessage(f'Select right continuum end for {self.line_labels[self.current_line]}')
            elif self.selection_step == 3:
                self.right_continuum_end = x
                self.plot_continuum_region('right')
                self.statusBar().showMessage(f'Select integration start for {self.line_labels[self.current_line]}')
            elif self.selection_step == 4:
                self.integration_start = x
                self.statusBar().showMessage(f'Select integration end for {self.line_labels[self.current_line]}')
            elif self.selection_step == 5:
                self.integration_end = x
                self.plot_integration_region()
                self.canvas.mpl_disconnect(self.cid)
                self.cid = None
                self.calculate_flux_for_line()
                return
            self.selection_step += 1

    def plot_continuum_region(self, side):
        if side == 'left':
            start = self.left_continuum_start
            end = self.left_continuum_end
            color = 'blue'
            lines_attr = 'left_continuum_lines'
            patch_attr = 'left_continuum_patch'
        elif side == 'right':
            start = self.right_continuum_start
            end = self.right_continuum_end
            color = 'green'
            lines_attr = 'right_continuum_lines'
            patch_attr = 'right_continuum_patch'

        # Remove previous lines and patches if they exist
        for line in getattr(self, lines_attr):
            try:
                line.remove()
            except Exception as e:
                print(f'error removing line: {e}')
        setattr(self, lines_attr, [])
        try:
            if getattr(self, patch_attr):
                getattr(self, patch_attr).remove()
                setattr(self, patch_attr, None)
        except Exception as e:
            print(f'error removing patch: {e}')

        # Plot vertical lines
        line1 = self.ax.axvline(x=start, color=color, alpha=0.5)
        line2 = self.ax.axvline(x=end, color=color,alpha = 0.5)
        setattr(self, lines_attr, [line1, line2])

        # Shade region
        setattr(self, patch_attr, self.ax.axvspan(start, end, color=color, alpha=0.3))
        self.canvas.draw()

    def plot_integration_region(self):
        # Remove previous lines and patches if they exist
        # for line in self.integration_lines:
        #     line.remove()
        self.integration_lines = []

        if self.integration_patch:
            try:
                self.integration_patch.remove()
                self.integration_patch = None
            except Exception as e:
                print(f'weird integration patch error')

        # Plot vertical lines
        line1 = self.ax.axvline(x=self.integration_start, color='red',alpha=0.5)
        line2 = self.ax.axvline(x=self.integration_end, color='red', alpha=0.5)
        self.integration_lines = [line1, line2]

        # Shade region
        self.integration_patch = self.ax.axvspan(self.integration_start, self.integration_end, color='red', alpha=0.3)
        self.canvas.draw()

    def calculate_flux_for_line(self):
        try:
            table_row = self.current_line - 1
            
            # Get the flux and error values
            print("About to calculate flux...")
            flux_value, flux_error, cont_slope, cont_intercept = self.get_flux(
                self.integration_start, self.integration_end,
                self.left_continuum_start, self.left_continuum_end,
                self.right_continuum_start, self.right_continuum_end
            )
            print(f"Flux calculation returned: {flux_value}, {flux_error}")
            # Update the objects table if a row is selected

            if flux_value is None:
                self.statusBar().showMessage(f'Flux calculation failed for {self.line_labels[self.current_line]}')
                return

            # Update the table with calculated values
            flux_item = QTableWidgetItem(f'{flux_value:.2e}')
            self.table.setItem(table_row, 1, flux_item)
            flux_error_item = QTableWidgetItem(f'{flux_error:.2e}')
            self.table.setItem(table_row, 2, flux_error_item)

            # Update the table with continuum bounds
            self.table.setItem(table_row, 3, QTableWidgetItem(f'{self.left_continuum_start:.2f}'))
            self.table.setItem(table_row, 4, QTableWidgetItem(f'{self.left_continuum_end:.2f}'))
            self.table.setItem(table_row, 5, QTableWidgetItem(f'{self.right_continuum_start:.2f}'))
            self.table.setItem(table_row, 6, QTableWidgetItem(f'{self.right_continuum_end:.2f}'))
            self.table.setItem(table_row, 7, QTableWidgetItem(f'{cont_slope:.2e}'))
            self.table.setItem(table_row, 8, QTableWidgetItem(f'{cont_intercept:.2e}'))

            self.table.resizeColumnsToContents()

            selected_rows = self.objects_table.selectionModel().selectedRows()
            if selected_rows:
                row_index = selected_rows[0].row()
                self.update_object_row_from_line_values(row_index)
                self.update_object_csv_from_left_table(row_index)
            # Update the objects table with the calculated flux
            #self.update_object_row_from_line_values(table_row)
            self.statusBar().showMessage(f'Calculated flux for {self.line_labels[self.current_line]}')
        except Exception as e:
            import traceback
            print(f"Error in calculate_flux_for_line: {str(e)}")
            print(traceback.format_exc())
            self.statusBar().showMessage(f'Error calculating flux: {str(e)}')


    def get_continuum(self, continuum_wave, continuum_flux, continuum_error, weights = None):
        
        print(f"Start of get_continuum")
        
        wavelength_array = continuum_wave

        flux_array = continuum_flux
        print(f'Successfully masked for continuum')
        
        if (weights == None):
            weights = np.nan_to_num(1/(continuum_error))
        
        print(f'Weights worked')
        def linear(x, m, b):
            return x*m + b
        

        linear_model = lmfit.Model(linear)
        print(f'linear model worked')
        params_continuum = lmfit.Parameters()
        params_continuum.add('m', value = (np.nanmean(flux_array)/2))
        params_continuum.add('b', value = 0 )
        nan_fluxes = 0.
        total_fluxes = 0.
        for flux in flux_array:
            total_fluxes += 1
            if np.isnan(flux):
                nan_fluxes += 1
        print(f'total fluxes: {total_fluxes}, nan fluxes: {nan_fluxes}')
        try:
            result = linear_model.fit(flux_array, x=wavelength_array, params=params_continuum, weights = weights)
        except Exception as e:
            print(f'Issue with result... error is: {e}')
            result = None
        print(f'Got result')
        print(result.fit_report())
        values = {}
        values['slope'] = result.best_values['m']
        values['intercept'] = result.best_values['b']
        m = values['slope']
        b = values['intercept'] 

        # print(f'slope is {m}, intercept is {b}')
        return values
    def get_flux(self, integ_start, integ_end, cont1_start, cont1_end, cont2_start, cont2_end):
        """
        Calculates the flux and flux error for a given spectral line.

        Parameters
        ----------
        integ_start : float
            Start of the integration range.
        integ_end : float
            End of the integration range.
        cont1_start : float
            Start of the first continuum region.
        cont1_end : float
            End of the first continuum region.
        cont2_start : float
            Start of the second continuum region.
        cont2_end : float
            End of the second continuum region.

        Returns
        -------
        flux : float
            Calculated flux.
        flux_error : float
            Estimated flux error.
        """
        try:
            wavelength_array = self.coadded_spectrum['wave']
            flux_array = self.coadded_spectrum['flux']
            error_array = (self.coadded_spectrum['error_up'] + self.coadded_spectrum['error_down']) / 2
            counts_array = self.coadded_spectrum['counts']
            print(f'Ranges successfully found')
            # Masks for the integration and continuum ranges
            integ_mask = (wavelength_array >= integ_start) & (wavelength_array <= integ_end)
            cont1_mask = (wavelength_array >= cont1_start) & (wavelength_array <= cont1_end)
            cont2_mask = (wavelength_array >= cont2_start) & (wavelength_array <= cont2_end)

            cont_mask = ((wavelength_array >= cont2_start) & (wavelength_array <= cont2_end) | (wavelength_array >= cont2_start) & (wavelength_array <= cont2_end))
            # Mask the flux and error arrays
            flux_array_continuum = flux_array[cont_mask]
            error_array_continuum = error_array[cont_mask]
            wavelength_array_continuum = wavelength_array[cont_mask]
            print(f'mask done')
            # Check for empty masks
            if not np.any(integ_mask):
                self.statusBar().showMessage('No data in the selected integration range.')
                return None, None
            if not np.any(cont1_mask) or not np.any(cont2_mask):
                self.statusBar().showMessage('No data in one or both of the selected continuum ranges.')
                return None, None

            # Calculate continuum levels
            print(f'Finding continuum...')
            continuum_values = self.get_continuum(wavelength_array_continuum, flux_array_continuum, error_array_continuum)
            def linear(x, m, b):
                return x*m + b
            
            slope = continuum_values['slope']
            intercept = continuum_values['intercept']
            print(f'slope is {slope:0.3f}, intercept is {intercept:0.3f}')
            continuum = linear(wavelength_array, slope, intercept)

            # Subtract continuum from flux
            flux_corrected = flux_array[integ_mask] - continuum[integ_mask]
            flux = simpson(flux_corrected, wavelength_array[integ_mask])
            print(f"Calculated flux: {flux}, scaled by {self.scale_value}")
            # Estimate flux error
            counts_array = counts_array[integ_mask]
            total_counts = np.sum(counts_array)
            print(f'total counts worked')
            poisson_conf = pcf(total_counts, interval = 'frequentist-confidence')
            print(f'pcf worked')
            gcounts_error_up = poisson_conf[1] - total_counts
            gcounts_error_down = total_counts - poisson_conf[0] 
            error_poisson_up = (gcounts_error_up/total_counts)*flux
            error_poisson_down = (gcounts_error_down/total_counts)*flux
            # error_integ = error_array[integ_mask]
            flux_error = np.max([error_poisson_up, error_poisson_down])
            print(f'found flux error in get_flux: {flux_error}')

            return flux, flux_error, slope, intercept
        except Exception as e:
            self.statusBar().showMessage(f'Error calculating flux: {e}')
            return None, None


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_win = SpectralFluxApp()
    main_win.show()
    sys.exit(app.exec_())
