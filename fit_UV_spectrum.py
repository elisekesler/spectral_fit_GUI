import sys
import csv
from functools import partial
import numpy as np
from scipy.integrate import simps
from astropy.table import vstack
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget,
    QFileDialog, QLabel, QHBoxLayout, QLineEdit, QProgressBar,
    QInputDialog, QMessageBox, QComboBox, QTableWidget,
    QTableWidgetItem
)
from PyQt5.QtCore import QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from matplotlib.figure import Figure
from prepare_spec import coadd, prepare, prepare_other_grating


class CoaddWorker(QThread):
    update_progress = pyqtSignal(int)
    coadd_complete = pyqtSignal(object)

    def __init__(self, spectrum_data):
        super().__init__()
        self.spectrum_data = spectrum_data

    def run(self):
        try:
            combined_table = vstack(self.spectrum_data)
            delta = 1  # Example value, adjust as necessary

            # Simulate progress updates
            for i in range(1, 101):
                self.update_progress.emit(i)
                self.msleep(50)  # Simulate some delay for each step

            coadded_spectrum = coadd(combined_table, delta)
            self.coadd_complete.emit(coadded_spectrum)
        except Exception as e:
            self.coadd_complete.emit(None)


class SpectralFluxApp(QMainWindow):
    def __init__(self, csv_file_path=None):
        super().__init__()

        self.init_variables()
        self.create_main_widget()
        self.create_plot_area()
        self.create_controls()
        self.create_table()

        self.setWindowTitle('Spectral Flux Measurement')
        self.setGeometry(100, 100, 1000, 800)

    def init_variables(self):
        # Initialize variables
        self.spectrum_data = []
        self.coadded_spectrum = None
        self.fits_file_paths = []
        self.redshift = 0.0
        self.selection_step = 0
        self.cid = None  # Connection ID for event handler

        self.line_wavelengths = np.array([
            1215.67, 1031.92, 1037.61, 1238.82, 1242.8,
            1393.75, 1402.77, 1548.19, 1550.77
        ])
        self.line_labels = [
            'Lya', 'O VI 1031', 'O VI 1037', 'N V 1238',
            'N V 1242', 'Si 1393', 'Si 1402', 'C IV 1548', 'C IV 1550'
        ]

    def create_main_widget(self):
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        self.main_layout = QVBoxLayout(self.main_widget)

    def create_plot_area(self):
        self.figure = Figure(figsize=(10, 8))
        self.canvas = FigureCanvas(self.figure)
        self.main_layout.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.main_layout.addWidget(self.toolbar)

    def create_controls(self):
        self.controls_layout = QHBoxLayout()
        self.main_layout.addLayout(self.controls_layout)

        # Load data button
        self.load_button = QPushButton('Load CSV with FITS Paths')
        self.load_button.clicked.connect(self.load_csv_with_fits_paths)
        self.controls_layout.addWidget(self.load_button)

        # Co-add and plot button
        self.coadd_button = QPushButton('Co-add and Plot')
        self.coadd_button.clicked.connect(self.coadd_and_plot)
        self.controls_layout.addWidget(self.coadd_button)

        # Redshift label and input
        self.redshift_label = QLabel('Redshift:')
        self.controls_layout.addWidget(self.redshift_label)

        self.redshift_input = QLineEdit(self)
        self.redshift_input.setPlaceholderText("Enter redshift value")
        self.controls_layout.addWidget(self.redshift_input)

        # Plot expected line locations button
        self.plot_lines_button = QPushButton('Plot Expected Line Locations')
        self.plot_lines_button.setEnabled(False)
        self.controls_layout.addWidget(self.plot_lines_button)

        # Line selection dropdown
        self.surrounding_label = QLabel('Plot area surrounding:')
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
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(['Line', 'Flux', 'Flux Error', 'Action'])
        self.table.setRowCount(len(self.line_labels))

        # Populate the table
        for i, line_label in enumerate(self.line_labels):
            self.table.setItem(i, 0, QTableWidgetItem(line_label))
            self.table.setItem(i, 1, QTableWidgetItem(''))
            self.table.setItem(i, 2, QTableWidgetItem(''))
            find_flux_error_button = QPushButton('Find Flux and Error')
            find_flux_error_button.clicked.connect(partial(self.find_flux_for_line, i))
            self.table.setCellWidget(i, 3, find_flux_error_button)

        self.main_layout.addWidget(self.table)

    def load_csv_with_fits_paths(self):
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getOpenFileName(
            self, "Open CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
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
            redshift, ok = QInputDialog.getDouble(
                self, "Enter Redshift", "Redshift:", self.redshift, -10.0, 10.0, 3
            )
            if ok:
                self.redshift = redshift

        self.redshift_input.setText(str(self.redshift))

    def load_fits_files(self):
        self.spectrum_data = [self.load_and_prepare_fits(file) for file in self.fits_file_paths]
        if self.spectrum_data:
            self.plot_lines_button.setEnabled(True)

    def load_and_prepare_fits(self, file_name):
        try:
            if 'other_grating' in file_name:
                return prepare_other_grating(file_name)
            else:
                return prepare(file_name)
        except Exception as e:
            self.statusBar().showMessage(f'Error loading FITS file {file_name}: {e}')
            return None

    def coadd_and_plot(self):
        if len(self.spectrum_data) > 0:
            self.worker = CoaddWorker(self.spectrum_data)
            self.worker.update_progress.connect(self.progress_bar.setValue)
            self.worker.coadd_complete.connect(self.on_coadd_complete)
            self.worker.start()
            self.statusBar().showMessage('Co-adding and plotting...')
        else:
            self.statusBar().showMessage('No spectra loaded')

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
            ax = self.figure.add_subplot(111)
            ax.plot(wave, flux, label='Co-added Spectrum', drawstyle='steps-mid')
            ax.fill_between(wave, flux, flux + error_up, color='blue', alpha=0.3, label='Upper Error')
            ax.fill_between(wave, flux - error_down, flux, color='red', alpha=0.2, label='Lower Error')
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')
            self.plot_expected_lines(ax)
            ax.legend()
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

            for i, line in enumerate(observed_lines):
                ax.axvline(x=line, linestyle='--', color=colors[i % len(colors)], label=f'{self.line_labels[i]}')
                ax.text(line, ax.get_ylim()[1] * 0.9, self.line_labels[i],
                        color=colors[i % len(colors)], rotation=90, verticalalignment='bottom')
        except ValueError:
            self.statusBar().showMessage('Invalid redshift input. Please enter a valid number.')

    def zoom_to_line(self):
        if self.coadded_spectrum:
            selected_line = self.line_dropdown.currentIndex()
            rest_frame_wavelength = self.line_wavelengths[selected_line] * (1 + self.redshift)

            wave = self.coadded_spectrum['wave']
            flux = self.coadded_spectrum['flux']

            zoom_range = (rest_frame_wavelength - 50, rest_frame_wavelength + 50)
            mask = (wave >= zoom_range[0]) & (wave <= zoom_range[1])
            wave_zoom = wave[mask]
            flux_zoom = flux[mask]

            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.plot(wave, flux, label='Co-added Spectrum')
            ax.set_xlim(zoom_range)
            ax.set_ylim(0, np.max(flux_zoom) * 1.1)
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')
            self.plot_expected_lines(ax)
            self.canvas.draw()
        else:
            self.statusBar().showMessage('No co-added spectrum to zoom into')

    def find_flux_for_line(self, line_index):
        self.current_line = line_index
        self.selection_step = 0
        self.statusBar().showMessage(f'Select left continuum start for {self.line_labels[line_index]}')

        # Connect event handler
        if self.cid is None:
            self.cid = self.canvas.mpl_connect('button_press_event', self.on_click)

    def on_click(self, event):
        if event.inaxes is not None:
            x = event.xdata
            if self.selection_step == 0:
                self.left_continuum_start = x
                self.statusBar().showMessage(f'Select left continuum end for {self.line_labels[self.current_line]}')
            elif self.selection_step == 1:
                self.left_continuum_end = x
                self.statusBar().showMessage(f'Select right continuum start for {self.line_labels[self.current_line]}')
            elif self.selection_step == 2:
                self.right_continuum_start = x
                self.statusBar().showMessage(f'Select right continuum end for {self.line_labels[self.current_line]}')
            elif self.selection_step == 3:
                self.right_continuum_end = x
                self.statusBar().showMessage(f'Select integration start for {self.line_labels[self.current_line]}')
            elif self.selection_step == 4:
                self.integration_start = x
                self.statusBar().showMessage(f'Select integration end for {self.line_labels[self.current_line]}')
            elif self.selection_step == 5:
                self.integration_end = x
                self.canvas.mpl_disconnect(self.cid)
                self.cid = None
                self.calculate_flux_for_line()
                return
            self.selection_step += 1

    def calculate_flux_for_line(self):
        flux_value, flux_error = self.get_flux(
            self.integration_start, self.integration_end,
            self.left_continuum_start, self.left_continuum_end,
            self.right_continuum_start, self.right_continuum_end
        )

        if flux_value is None:
            self.statusBar().showMessage(f'Flux calculation failed for {self.line_labels[self.current_line]}')
            return

        # Display the calculated flux and error in the table
        flux_item = QTableWidgetItem(f'{flux_value:.2e}')
        self.table.setItem(self.current_line, 1, flux_item)
        flux_error_item = QTableWidgetItem(f'{flux_error:.2e}')
        self.table.setItem(self.current_line, 2, flux_error_item)

        self.statusBar().showMessage(f'Calculated flux for {self.line_labels[self.current_line]}')

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

            # Masks for the integration and continuum ranges
            integ_mask = (wavelength_array >= integ_start) & (wavelength_array <= integ_end)
            cont1_mask = (wavelength_array >= cont1_start) & (wavelength_array <= cont1_end)
            cont2_mask = (wavelength_array >= cont2_start) & (wavelength_array <= cont2_end)

            # Check for empty masks
            if not np.any(integ_mask):
                self.statusBar().showMessage('No data in the selected integration range.')
                return None, None
            if not np.any(cont1_mask) or not np.any(cont2_mask):
                self.statusBar().showMessage('No data in one or both of the selected continuum ranges.')
                return None, None

            # Calculate continuum levels
            cont1_flux = np.mean(flux_array[cont1_mask])
            cont2_flux = np.mean(flux_array[cont2_mask])

            # Calculate the slope and intercept of the continuum
            cont1_center = np.mean(wavelength_array[cont1_mask])
            cont2_center = np.mean(wavelength_array[cont2_mask])
            slope = (cont2_flux - cont1_flux) / (cont2_center - cont1_center)
            intercept = cont1_flux - slope * cont1_center

            # Continuum model over integration range
            continuum = slope * wavelength_array[integ_mask] + intercept

            # Subtract continuum from flux
            flux_corrected = flux_array[integ_mask] - continuum

            # Integrate the flux using Simpson's rule
            flux = simps(flux_corrected, wavelength_array[integ_mask])

            # Estimate flux error
            error_integ = error_array[integ_mask]
            flux_error = np.sqrt(np.sum(error_integ**2)) * np.mean(np.diff(wavelength_array[integ_mask]))

            # Plot for visualization
            ax = self.figure.gca()
            ax.plot(wavelength_array[integ_mask], flux_array[integ_mask], label="Flux")
            ax.plot(wavelength_array[integ_mask], continuum, label="Continuum", linestyle='--')
            ax.fill_between(wavelength_array[integ_mask], continuum, flux_array[integ_mask],
                            alpha=0.3, label=f'Integrated Flux = {flux:.2e}')
            ax.legend()
            self.canvas.draw()

            return flux, flux_error
        except Exception as e:
            self.statusBar().showMessage(f'Error calculating flux: {e}')
            return None, None


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_win = SpectralFluxApp()
    main_win.show()
    sys.exit(app.exec_())
