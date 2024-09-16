import sys
import argparse
import csv
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QFileDialog, QLabel, QHBoxLayout, QLineEdit, QProgressBar, QInputDialog, QMessageBox, QComboBox
)
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem, QPushButton
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QDoubleValidator
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from astropy.table import vstack
from prepare_spec import coadd, prepare, prepare_other_grating
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from UV_spec_class import spec
from scipy.integrate import simps
import matplotlib.pyplot as plt
from astropy.stats import mad_std


class CoaddWorker(QThread):
    update_progress = pyqtSignal(int)
    coadd_complete = pyqtSignal(object)

    def __init__(self, spectrum_data):
        super().__init__()
        self.spectrum_data = spectrum_data

    def run(self):
        combined_table = vstack(self.spectrum_data)
        delta = 1 # Example value, adjust as necessary

        # Simulate progress updates (since actual progress is hard to track with vstack and coadd)
        for i in range(1, 101):
            self.update_progress.emit(i)
            self.msleep(50)  # Simulate some delay for each step

        coadded_spectrum = coadd(combined_table, delta)
        self.coadd_complete.emit(coadded_spectrum)


class SpectralFluxApp(QMainWindow):
    def __init__(self, csv_file_path=None):
        super().__init__()

        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)

        # Create the main layout
        self.main_layout = QVBoxLayout(self.main_widget)

        # Plot area with a larger figure
        self.figure = Figure(figsize=(10, 8))  # Set larger figure size
        self.canvas = FigureCanvas(self.figure)
        self.main_layout.addWidget(self.canvas)

        # Add toolbar directly below the canvas
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.main_layout.addWidget(self.toolbar)

        # Controls layout for buttons, labels, etc.
        self.controls_layout = QHBoxLayout()
        self.main_layout.addLayout(self.controls_layout)

        # Load data button
        self.load_button = QPushButton('Load CSV with FITS Paths')
        self.load_button.clicked.connect(self.load_csv_with_fits_paths)
        self.controls_layout.addWidget(self.load_button)

        self.coadd_button = QPushButton('Co-add and Plot')
        self.coadd_button.clicked.connect(self.coadd_and_plot)
        self.controls_layout.addWidget(self.coadd_button)

        # Add Redshift label and input
        self.redshift_label = QLabel('Redshift:')
        self.controls_layout.addWidget(self.redshift_label)

        self.redshift_input = QLineEdit(self)
        self.redshift_input.setMaxLength(20)
        self.redshift_input.setPlaceholderText("Enter redshift value")
        self.redshift_input.setClearButtonEnabled(True)
        self.redshift_input.setDragEnabled(True)
        self.controls_layout.addWidget(self.redshift_input)

        # Add button to plot expected line locations
        self.plot_lines_button = QPushButton('Plot Expected Line Locations')
        self.plot_lines_button.setEnabled(False)
        self.controls_layout.addWidget(self.plot_lines_button)

        # Add drop down for "plot area surrounding"
        self.surrounding_label = QLabel('Plot area surrounding:')
        self.controls_layout.addWidget(self.surrounding_label)

        self.line_dropdown = QComboBox(self)
        self.line_dropdown.addItems(['Lya', 'O VI 1031', 'O VI 1037', 'N V 1238', 'N V 1242', 'Si 1393', 'Si 1402', 'C IV 1548', 'C IV 1550'])
        self.line_dropdown.currentIndexChanged.connect(self.zoom_to_line)
        self.controls_layout.addWidget(self.line_dropdown)

        # Progress bar
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setValue(0)
        self.main_layout.addWidget(self.progress_bar)

        # Status bar
        self.statusBar().showMessage('Ready')

        self.line_wavelengths = np.array([1215.67, 1031.92, 1037.61, 1238.82, 1242.8, 1393.75, 1402.77, 1548.19, 1550.77])
        self.line_labels = ['Lya', 'O VI 1031', 'O VI 1037', 'N V 1238', 'N V 1242', 'Si 1393', 'Si 1402', 'C IV 1548', 'C IV 1550']

        
        # Add table to display flux and flux error for each line
        self.table = QTableWidget(self)
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(['Line', 'Flux', 'Flux Error'])
        self.table.setRowCount(len(self.line_labels))

        # Populate the table with line names and buttons to find flux
        for i, line_label in enumerate(self.line_labels):
            self.table.setItem(i, 0, QTableWidgetItem(line_label))  # Line name

            # Create layout for Flux column
            flux_layout = QHBoxLayout()
            flux_item = QTableWidgetItem('')  # Placeholder for calculated flux
            self.table.setItem(i, 1, flux_item)
            
            find_flux_button = QPushButton(f'Find Flux')
            find_flux_button.clicked.connect(lambda _, l=i: self.find_flux_for_line(l))  # Connect button to find_flux_for_line
            flux_layout.addWidget(find_flux_button)

            # You need to use QWidget to wrap the layout and set it as cell widget
            flux_widget = QWidget()
            flux_widget.setLayout(flux_layout)
            self.table.setCellWidget(i, 1, flux_widget)

            # Create layout for Flux Error column
            flux_error_layout = QHBoxLayout()
            flux_error_item = QTableWidgetItem('')  # Placeholder for flux error
            self.table.setItem(i, 2, flux_error_item)

            find_flux_error_button = QPushButton(f'Find Flux Error')
            find_flux_error_button.clicked.connect(lambda _, l=i: self.find_flux_for_line(l))  # Modify if you want separate logic
            flux_error_layout.addWidget(find_flux_error_button)

            # Again, use QWidget to wrap the layout
            flux_error_widget = QWidget()
            flux_error_widget.setLayout(flux_error_layout)
            self.table.setCellWidget(i, 2, flux_error_widget)


        self.main_layout.addWidget(self.table)  # Add table to the layout

        self.setWindowTitle('Spectral Flux Measurement')
        self.setGeometry(100, 100, 1000, 800)  # Optionally adjust the window size

        # Initialize spectrum data
        self.spectrum_data = []
        self.coadded_spectrum = None
        self.fits_file_paths = []
        self.redshift = 0.0
    def load_csv_with_fits_paths(self):
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getOpenFileName(self, "Open CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if csv_file:
            self.load_csv(csv_file)

    def load_csv(self, csv_file_path):
        try:
            with open(csv_file_path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                lines = list(reader)

                # The first line is the redshift
                if lines:
                    self.redshift = float(lines[0][0])
                    self.confirm_redshift()

                # The remaining lines are FITS file paths
                self.fits_file_paths = [row[0] for row in lines[1:]]

            self.statusBar().showMessage(f'Loaded {len(self.fits_file_paths)} FITS file paths from CSV')
            self.load_fits_files()

        except Exception as e:
            self.statusBar().showMessage(f'Error loading CSV: {e}')

    def confirm_redshift(self):
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle('Redshift Confirmation')
        msg_box.setText(f'Is this your redshift? {self.redshift}')
        msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        result = msg_box.exec_()

        if result == QMessageBox.No:
            # If the user clicks No, prompt them to enter a new redshift
            redshift, ok = QInputDialog.getDouble(self, "Enter Redshift", "Redshift:", self.redshift, -10.0, 10.0, 3)
            if ok:
                self.redshift = redshift

        self.redshift_input.setText(str(self.redshift))

    def load_fits_files(self):
        self.spectrum_data = [self.load_and_prepare_fits(file) for file in self.fits_file_paths]
        if self.spectrum_data:  # If FITS files are loaded successfully
            self.plot_lines_button.setEnabled(True)

    def load_and_prepare_fits(self, file_name):
        if 'other_grating' in file_name:
            return prepare_other_grating(file_name)
        else:
            return prepare(file_name)

    def coadd_and_plot(self):
        if len(self.spectrum_data) > 0:
            # Create a worker thread for coadding and plotting
            self.worker = CoaddWorker(self.spectrum_data)
            self.worker.update_progress.connect(self.progress_bar.setValue)
            self.worker.coadd_complete.connect(self.on_coadd_complete)
            self.worker.start()
            self.statusBar().showMessage('Co-adding and plotting...')
        else:
            self.statusBar().showMessage('No spectra loaded')

    def on_coadd_complete(self, coadded_spectrum):
        self.coadded_spectrum = coadded_spectrum
        self.plot_spectrum()
        self.statusBar().showMessage('Co-addition and plotting complete')
        self.progress_bar.setValue(0)  # Reset the progress bar

    def plot_spectrum(self):
        if self.coadded_spectrum:
            wave = self.coadded_spectrum['wave']
            flux = self.coadded_spectrum['flux']

            # Plotting the spectrum
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.plot(wave, flux, label='Co-added Spectrum')
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')

            # Automatically plot expected lines
            self.plot_expected_lines(ax)

            # ax.legend()
            self.canvas.draw()
        else:
            self.statusBar().showMessage('No co-added spectrum to plot')

    def plot_expected_lines(self, ax=None):
        try:
            redshift = self.redshift

            # Apply redshift to get observed line locations
            observed_lines = self.line_wavelengths * (1 + redshift)
            colors = ['r', 'g', 'b', 'orange', 'purple', 'pink', 'yellow', 'dark blue', 'black']

            if ax is None:
                ax = self.figure.gca()

            for i, line in enumerate(observed_lines):
                ax.axvline(x=line, linestyle='--', color=colors[i], label=f'{self.line_labels[i]}')
                ax.text(line, ax.get_ylim()[1] * 0.9, self.line_labels[i], color=colors[i], rotation=90, verticalalignment='bottom')

            # ax.legend()

        except ValueError:
            self.statusBar().showMessage('Invalid redshift input. Please enter a valid number.')

    def zoom_to_line(self):
        if self.coadded_spectrum:
            selected_line = self.line_dropdown.currentIndex()
            rest_frame_wavelength = self.line_wavelengths[selected_line] * (1 + self.redshift)

            wave = self.coadded_spectrum['wave']
            flux = self.coadded_spectrum['flux']

            # Zoom into the region surrounding the selected line (+/- 50 Ångströms)
            zoom_range = (rest_frame_wavelength - 50, rest_frame_wavelength + 50)

            # Get the portion of the spectrum within this range
            mask = (wave >= zoom_range[0]) & (wave <= zoom_range[1])
            wave_zoom = wave[mask]
            flux_zoom = flux[mask]

            # Plot the zoomed-in region
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.plot(wave, flux, label='Co-added Spectrum')
            ax.set_xlim(zoom_range)
            ax.set_ylim(0, np.max(flux_zoom) * 1.1)  # Set y-axis limit to 10% more than the max flux in the region
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')

            # Replot expected lines in the zoomed region
            self.plot_expected_lines(ax)

            self.canvas.draw()
    def find_flux_for_line(self, line_index):
        self.current_line = line_index
        self.statusBar().showMessage(f'Select left continuum for {self.line_labels[line_index]}')

        # Listen for the next two clicks to set left continuum
        self.canvas.mpl_connect('button_press_event', self.on_select_left_continuum)

    def on_select_left_continuum(self, event):
        if event.inaxes is not None:
            self.left_continuum_start = event.xdata
            self.statusBar().showMessage(f'Select right bound of left continuum for {self.line_labels[self.current_line]}')
            self.canvas.mpl_connect('button_press_event', self.on_select_right_continuum)

    def on_select_right_continuum(self, event):
        if event.inaxes is not None:
            self.left_continuum_end = event.xdata
            self.statusBar().showMessage(f'Select left bound of right continuum for {self.line_labels[self.current_line]}')
            self.canvas.mpl_connect('button_press_event', self.on_select_left_right_continuum)

    def on_select_left_right_continuum(self, event):
        if event.inaxes is not None:
            self.right_continuum_start = event.xdata
            self.statusBar().showMessage(f'Select right bound of right continuum for {self.line_labels[self.current_line]}')
            self.canvas.mpl_connect('button_press_event', self.on_select_right_right_continuum)

    def on_select_right_right_continuum(self, event):
        if event.inaxes is not None:
            self.right_continuum_end = event.xdata
            self.statusBar().showMessage(f'Select integration range for {self.line_labels[self.current_line]}')
            self.canvas.mpl_connect('button_press_event', self.on_select_integration_range)

    def on_select_integration_range(self, event):
        if event.inaxes is not None:
            self.integration_start = event.xdata
            self.statusBar().showMessage(f'Select right bound of integration range for {self.line_labels[self.current_line]}')
            self.canvas.mpl_connect('button_press_event', self.on_select_right_integration)

    def on_select_right_integration(self, event):
        if event.inaxes is not None:
            self.integration_end = event.xdata

            # Now that all bounds are selected, calculate flux
            self.calculate_flux_for_line()

    def calculate_flux_for_line(self):
        # Call the get_flux method to calculate flux for the selected line
        flux_value = self.get_flux(self.integration_start, self.integration_end,
                                self.left_continuum_start, self.left_continuum_end,
                                self.right_continuum_start, self.right_continuum_end)
        
        # Display flux in the table
        flux_item = QTableWidgetItem(f'{flux_value:.2e}')
        self.table.setItem(self.current_line, 1, flux_item)

        # You can also calculate flux error and display it (as needed)
        flux_error = 0  # Replace with actual error calculation logic
        flux_error_item = QTableWidgetItem(f'{flux_error:.2e}')
        self.table.setItem(self.current_line, 2, flux_error_item)

        # Clear the status bar
        self.statusBar().showMessage(f'Calculated flux for {self.line_labels[self.current_line]}')


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_win = SpectralFluxApp()
    main_win.show()
    sys.exit(app.exec_())
