import sys
import argparse
import csv
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QFileDialog, QLabel, QHBoxLayout, QLineEdit
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from astropy.table import Table, vstack
from prepare_spec import coadd, prepare, prepare_other_grating
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


class SpectralFluxApp(QMainWindow):
    def __init__(self, csv_file_path=None):
        super().__init__()

        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)

        # Create the main layout
        self.main_layout = QVBoxLayout(self.main_widget)

        # Plot area
        self.figure = Figure()
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
        self.controls_layout.addWidget(self.redshift_input)

        # Add button to plot expected line locations
        self.plot_lines_button = QPushButton('Plot Expected Line Locations')
        self.plot_lines_button.clicked.connect(self.plot_expected_lines)
        self.controls_layout.addWidget(self.plot_lines_button)

        # Status bar
        self.statusBar().showMessage('Ready')

        self.setWindowTitle('Spectral Flux Measurement')
        self.setGeometry(100, 100, 800, 600)

        # Initialize spectrum data
        self.spectrum_data = []
        self.coadded_spectrum = None
        self.fits_file_paths = []

        # If CSV file path was provided, load it on startup
        if csv_file_path:
            self.load_csv(csv_file_path)

    def load_csv_with_fits_paths(self):
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getOpenFileName(self, "Open CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if csv_file:
            self.load_csv(csv_file)

    def load_csv(self, csv_file_path):
        try:
            # Load the paths of FITS files from the CSV
            with open(csv_file_path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                self.fits_file_paths = [row[0] for row in reader]  # Assuming each row contains one file path

            self.statusBar().showMessage(f'Loaded {len(self.fits_file_paths)} FITS file paths from CSV')
            self.load_fits_files()
        except Exception as e:
            self.statusBar().showMessage(f'Error loading CSV: {e}')

    def load_fits_files(self):
        # Load the data from the FITS files listed in the CSV, using prepare functions
        self.spectrum_data = [self.load_and_prepare_fits(file) for file in self.fits_file_paths]

    def load_and_prepare_fits(self, file_name):
        # Determine which preparation function to use based on file content or metadata
        # Example: if file contains a specific grating type, use prepare_other_grating

        # For simplicity, we'll assume that files with 'other_grating' in their name use prepare_other_grating
        if 'other_grating' in file_name:
            return prepare_other_grating(file_name)
        else:
            return prepare(file_name)

    def coadd_and_plot(self):
        if len(self.spectrum_data) > 0:
            # Combine all tables using vstack
            combined_table = vstack(self.spectrum_data)

            # Perform co-addition (you can change delta as needed)
            delta = 0.1  # Example value, adjust as necessary
            self.coadded_spectrum = coadd(combined_table, delta)

            # Plot the co-added spectrum
            self.plot_spectrum()

        else:
            self.statusBar().showMessage('No spectra loaded')

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
            ax.legend()
            self.canvas.draw()
        else:
            self.statusBar().showMessage('No co-added spectrum to plot')

    def plot_expected_lines(self):
        try:
            redshift = float(self.redshift_input.text())
            line_wavelengths = np.array([1215.67, 1031.92, 1037.61, 1238.82, 1242.8, 1393.75, 1402.77, 1548.19, 1550.77])  # Example line locations in rest frame
            line_labels = ['Lya', 'O VI 1031', 'O VI 1037', 'N V 1238', 'N V 1242', 'Si 1393', 'Si 1402', 'C IV 1548', 'C IV 1550']
            observed_lines = line_wavelengths * (1 + redshift)

            ax = self.figure.gca()
            for i in range(len(line_wavelengths)):
                ax.axvline(x=line_wavelengths[i], linestyle='--', label=f'{line_labels[i]}')
            ax.legend()
            self.canvas.draw()
        except ValueError as e:
            print(e)
            self.statusBar().showMessage('Invalid redshift input. Please enter a number.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Spectral Flux App')
    parser.add_argument('--csv', type=str, help='Path to the CSV file containing FITS file paths', default=None)
    args = parser.parse_args()

    app = QApplication(sys.argv)
    main_win = SpectralFluxApp(csv_file_path=args.csv)
    main_win.show()
    sys.exit(app.exec_())
