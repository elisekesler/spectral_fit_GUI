import sys
import argparse
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QFileDialog, QLabel, QHBoxLayout, QLineEdit, QComboBox
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from astropy.table import Table, vstack

class SpectralFluxApp(QMainWindow):
    
    def __init__(self, file_path = None):
        super().__init__()

        self.file_path = file_path
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)

        self.layout = QVBoxLayout(self.main_widget)

        # plotting area 
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)

        # controls
        self.controls_layout = QHBoxLayout()
        self.layout.addLayout(self.controls_layout)

        # load data button
        self.load_button = QPushButton('Load Spectrum')
        self.load_button.clicked.connect(self.load_spectrum)
        self.controls_layout.addWidget(self.load_button)

        # labels and inputs for measuring fluxes
        self.line_center_dropdown = QComboBox()
        self.line_center_dropdown.addItems(self.get_spectral_lines())
        self.controls_layout.addWidget(self.line_center_dropdown)

        self.line_center_input = QLineEdit('Redshift')
        self.controls_layout.addWidget(self.line_center_input)

        self.line_button = QPushButton('Find Line Location')
        self.controls_layout.addWidget(self.line_button)

        self.load_button = QPushButton('Locate Spectrum')
        self.load_button.clicked.connect(self.load_spectrum)
        self.controls_layout.addWidget(self.load_button)

        self.measure_button = QPushButton('Measure Flux')
        self.measure_button.clicked.connect(self.measure_flux)
        self.controls_layout.addWidget(self.measure_button)

        # status bar
        self.statusBar().showMessage('Ready')

        self.setWindowTitle('Spectral Flux Measurement')
        self.setGeometry(100, 100, 800, 600)

        # initialize spectrum data
        self.spectrum_data = None
    
        if file_path:
                self.masked_spec = self.prepare(file = file_path)
                self.plot_spectrum(file_name = file_path)

    def get_spectral_lines(self):
        lines = [
            "H I Lya",
            "N V 1238.82/1242.80",
            "O VI 1031.92/1037.61",
            "C IV 1548.19/1550.77"
        ]
        return lines
    def load_spectrum(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Open Spectrum File", "", "All Files (*);;Text Files (*.txt)", options=options)
        if file_name:
            self.statusBar().showMessage(f'Loaded {file_name}')
            self.plot_spectrum(file_name)

    def prepare(self, file):
        """Nicely formats file into columns. Also does bitmasking"""
        
        spec = Table.read(file)

        spec_formatted = Table()
        
        spec_formatted['wave'] = np.array(spec['WAVELENGTH'][0])
        spec_formatted['flux'] = np.array(spec['FLUX'][0])
        spec_formatted['flag'] = np.array(spec['DQ'][0])
        
        if (spec['EXPTIME'].shape != (1,)):
            spec_formatted['exp_time'] = np.array(np.max(spec['EXPTIME']))
            
        else: 
            spec_formatted['exp_time'] = np.array(spec['EXPTIME'])
        spec_formatted['error_upper'] = np.array(spec['ERROR'][0])
        spec_formatted['error_lower'] = np.array(spec['ERROR_LOWER'][0])
        spec_formatted['G230L_error_up'] = np.zeros_like(spec_formatted['error_lower'])
        spec_formatted['G230L_error_down'] = np.zeros_like(spec_formatted['error_lower'])
        spec_formatted['gcounts'] = np.array(spec['GCOUNTS'][0])
        spec_formatted['background'] = np.array(spec['BACKGROUND'][0])
        spec_formatted['net'] = np.array(spec['NET'][0])

        net_exp_array = np.array(spec['NET'][0])*spec_formatted['exp_time'][0] 
        bkg_exp_array = np.array(spec['BACKGROUND'][0])*spec_formatted['exp_time'][0] 
        spec_formatted['bkg_times_exp'] = np.nan_to_num(bkg_exp_array)
        spec_formatted['net_times_exp_time'] = np.nan_to_num(net_exp_array)
        
        variance_flat = np.array(spec['VARIANCE_FLAT'][0])
        flags = spec_formatted['flag']
        mask = np.where((flags == 0) | (flags == 4))
        masked_spec = spec_formatted[(mask)]
        
        return masked_spec 
    

    def plot_spectrum(self, file_name=None):
        if file_name is None:
            file_name = self.file_path  # Use the file path stored in the instance

        if file_name:
            wavelength = self.spectrum_data[:, 0]
            flux = self.spectrum_data[:, 1]

            # Plotting the spectrum
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.plot(wavelength, flux, label='Spectrum')
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')
            ax.legend()
            self.canvas.draw()
        else:
            self.statusBar().showMessage('No file path provided.')


    def measure_flux(self):
        if self.spectrum_data is not None:
            try:
                line_center = float(self.line_center_input.text())
                flux = self.calculate_flux(line_center)
                self.statusBar().showMessage(f'Measured Flux: {flux:.2e}')
            except ValueError:
                self.statusBar().showMessage('Invalid input for line center')
        else:
            self.statusBar().showMessage('No spectrum loaded')

    def calculate_flux(self, line_center):
        # Placeholder for actual flux measurement logic
        # Here, you would integrate the flux around the line_center
        # For now, let's return a dummy value
        return np.random.random()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Spectral Flux App')
    parser.add_argument('--file', type=str, help='Path to the spectrum file to load at startup', default=None)
    args = parser.parse_args()

    app = QApplication(sys.argv)
    main_win = SpectralFluxApp(file_path=args.file)  # Pass file_path here
    main_win.show()
    sys.exit(app.exec_())
