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
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget,
    QFileDialog, QLabel, QHBoxLayout, QLineEdit, QProgressBar,
    QInputDialog, QMessageBox, QComboBox, QTableWidget,
    QTableWidgetItem, QSizePolicy, QFrame
)
import csv
from astropy.stats import poisson_conf_interval as pcf
from PyQt5.QtCore import QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from matplotlib.figure import Figure
from prepare_spec import coadd, prepare, prepare_other_grating, combine_tables
import pickle

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

    def fit_gaussian_for_line(self, line_index):
        # This is a placeholder for the Gaussian fitting functionality
        self.statusBar().showMessage(f"Gaussian fitting not yet implemented for {self.line_labels[line_index]}")
        # jj implement the Gaussian fitting here
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

    # def on_object_selection_change(self):
    #     # Get the selected row
    #     selected_rows = self.objects_table.selectionModel().selectedRows()
    #     if not selected_rows:
    #         return
            
    #     row_index = selected_rows[0].row()
        
    #     # Get the csv and fits file paths
    #     csv_item = self.objects_table.item(row_index, 9)  # Object_csv column
    #     fits_item = self.objects_table.item(row_index, 10)  # Object_fits column
        
    #     if not csv_item or not fits_item:
    #         return
            
    #     csv_path = csv_item.text()
    #     fits_path = fits_item.text()
        
    #     # Check if the paths exist
    #     if os.path.exists(csv_path) and os.path.exists(fits_path):
    #         # Load the CSV file (which should contain redshift info and possibly other metadata)
    #         self.load_csv(csv_path)
            
    #         # Load the spectrum from the fits file
    #         self.load_spectrum_from_path(fits_path)
            
    #         # Update the line values in the table
    #         self.update_line_values_from_object_row(row_index)

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
        
        # Create a scroll area for line entries
        from PyQt5.QtWidgets import QScrollArea
        
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
        
        # Initially disable buttons until a line is selected
        self.fit_model_button.setEnabled(False)
        self.add_line_button.setEnabled(False)
        self.view_results_button.setEnabled(False)
    # def create_gaussian_fit_controls(self):
    #     # Clear existing controls
    #     self.clear_method_controls()
        
    #     # Create top row of buttons
    #     buttons_layout = QHBoxLayout()
        
    #     self.fit_model_button = QPushButton("Fit Model")
    #     self.fit_model_button.clicked.connect(self.fit_gaussian_model)
    #     buttons_layout.addWidget(self.fit_model_button)
        
    #     self.add_line_button = QPushButton("Add Line")
    #     self.add_line_button.clicked.connect(self.add_gaussian_line)
    #     buttons_layout.addWidget(self.add_line_button)
        
    #     self.view_results_button = QPushButton("View Fit Results")
        # self.view_results_button.clicked.connect(self.view_fit_results)
        # buttons_layout.addWidget(self.view_results_button)
        # self.method_layout.addLayout(buttons_layout)
        
        # # Create a container for line entries
        # self.line_entries_container = QWidget()
        # self.line_entries_layout = QVBoxLayout(self.line_entries_container)
        # self.method_layout.addWidget(self.line_entries_container)
        
        # # Initially disable buttons until a line is selected
        # self.fit_model_button.setEnabled(False)
        # self.add_line_button.setEnabled(False)
        # self.view_results_button.setEnabled(False)
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
            self.redshift = new_redshift
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
    # Add placeholder functions for Gaussian Fit mode
    def fit_gaussian_model(self):
        self.statusBar().showMessage("Gaussian model fitting - placeholder")

    def view_fit_results(self):
        self.statusBar().showMessage("View fit results - placeholder")
    def add_gaussian_line(self):
        self.statusBar().showMessage("Adding Gaussian line - placeholder")
        
        # Create a new line entry row (black box with a button)
        line_entry = QFrame()
        line_entry.setFrameShape(QFrame.StyledPanel)
        line_entry.setFrameShadow(QFrame.Raised)
        line_entry.setStyleSheet("background-color: black; color: white;")
        line_entry.setMinimumHeight(40)
        line_entry.setMaximumHeight(40)
        
        entry_layout = QHBoxLayout(line_entry)
        entry_layout.setContentsMargins(10, 5, 10, 5)
        
        line_label = QLabel("Line Name")
        line_label.setStyleSheet("color: white;")
        entry_layout.addWidget(line_label)
        
        entry_layout.addStretch()
        
        view_params_button = QPushButton("View Parameters")
        view_params_button.clicked.connect(lambda: self.view_line_parameters(line_label.text()))
        entry_layout.addWidget(view_params_button)
        
        # Remove the stretch from the end to allow new widgets to appear
        if self.line_entries_layout.count() > 0:
            # Check if the last item is a spacer
            item = self.line_entries_layout.itemAt(self.line_entries_layout.count() - 1)
            if item and item.spacerItem():
                self.line_entries_layout.removeItem(item)
        
        # Add the line entry to the container
        self.line_entries_layout.addWidget(line_entry)
        
        # Add the stretch back
        self.line_entries_layout.addStretch()
        
        # Update the scroll area to ensure the new widget is visible
        self.scroll_area.ensureWidgetVisible(line_entry)

    def clear_line_entries(self):
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
    def view_line_parameters(self, line_name):
        QMessageBox.information(self, "Line Parameters", f"Parameters for {line_name} - placeholder")

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
        self.ax.plot(wave, self.flux_after_subtraction, label='Spectrum with Continuum Subtracted', color='blue')
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
            ax.plot(wave, flux, label='Co-added Spectrum', color='black', drawstyle='steps-mid')
            ax.plot(wave, (error_up + error_down)/2, label='Error', color='grey', drawstyle='steps-mid')
            ax.set_xlabel('Wavelength')
            ax.set_xlim(1100, 1880)
            # print(np.max(flux[~np.isnan(flux)]))
            # ax.set_ylim(-0.5, np.max(flux[~np.isnan(flux)]))
            ax.set_ylabel(f'Flux (erg / s / cm^2 / A ) (scaled by {self.scale_value:0.2e})')
            self.plot_expected_lines(ax)
            self.canvas.draw()
        else:
            self.statusBar().showMessage('No co-added spectrum to plot')

    def plot_expected_lines(self, ax=None):
        try:
            redshift = float(self.redshift_input.text()) if self.redshift_input.text() else self.redshift
            observed_lines = self.line_wavelengths * (1 + redshift)
            colors = ['r', 'g', 'b', 'orange', 'purple', 'pink', 'yellow', 'darkblue', 'black']

            if ax is None:
                ax = self.figure.gca()

            for i, line in enumerate(observed_lines[1:]):  # Skip 'Full Spectrum' placeholder
                ax.axvline(x=line, linestyle='--', color=colors[i % len(colors)], label=f'{self.line_labels[i+1]}')
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
                ax.plot(wave, flux, label='Co-added Spectrum', color='black', drawstyle='steps-mid')
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
                ax.plot(wave, flux, label='Co-added Spectrum', color='black', drawstyle='steps-mid')
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
        for line in self.integration_lines:
            line.remove()
        self.integration_lines = []

        if self.integration_patch:
            self.integration_patch.remove()
            self.integration_patch = None

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
