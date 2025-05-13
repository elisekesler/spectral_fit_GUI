import sys
import os
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget,
    QFileDialog, QLabel, QHBoxLayout, QLineEdit, QProgressBar,
    QMessageBox, QComboBox, QListWidget, QGroupBox, QGridLayout,
    QSpinBox, QDoubleSpinBox
)
from PyQt5.QtCore import QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from matplotlib.figure import Figure
from astropy.table import Table, vstack
from astropy.io import fits
from prepare_spec import prepare, prepare_other_grating, combine_tables, coadd

class CoadderThread(QThread):
    progress_updated = pyqtSignal(int)
    finished_signal = pyqtSignal(object)
    error_signal = pyqtSignal(str)
    
    def __init__(self, file_data, delta, parent=None):
        super().__init__(parent)
        self.file_data = file_data
        self.delta = delta
        
    def prepare(self, file):
        """Nicely formats file into columns. Also does bitmasking"""
        print(f'preparing {file}')
        spec = Table.read(file)
        print(f'file read')
        spec_formatted = Table()
        try:
            if spec['DQ'][0] == 55:
                print(f'this is a test spectrum')
                # this is the flag for a test spectrum. maybe not the best way to do it, but it works
                spec['DQ'][0] = 0
                spec_formatted['wave'] = np.array(spec['WAVELENGTH'])
                spec_formatted['flux'] = np.array(spec['FLUX'])
                spec_formatted['flag'] = np.array(spec['DQ'])
                spec_formatted['exp_time'] = np.array(spec['EXPTIME'])

                spec_formatted['error_upper'] = np.array(spec['ERROR'])
                spec_formatted['error_lower'] = np.array(spec['ERROR_LOWER'])
                spec_formatted['G230L_error_up'] = np.zeros_like(spec_formatted['error_lower'])
                spec_formatted['G230L_error_down'] = np.zeros_like(spec_formatted['error_lower'])
                spec_formatted['gcounts'] = np.array(spec['GCOUNTS'])
                spec_formatted['background'] = np.array(spec['BACKGROUND'])
                spec_formatted['net'] = np.array(spec['NET'])
                net_exp_array = np.array(spec['NET'])*spec_formatted['exp_time'] 
                bkg_exp_array = np.array(spec['BACKGROUND'])*spec_formatted['exp_time']
                variance_flat = np.array(spec['VARIANCE_FLAT'])
                spec_formatted['bkg_times_exp'] = np.nan_to_num(bkg_exp_array)
                spec_formatted['net_times_exp_time'] = np.nan_to_num(net_exp_array)
                
                flags = spec_formatted['flag']
                mask = np.where((flags == 0) | (flags == 4))
                masked_spec = spec_formatted[(mask)]

                return masked_spec

        except Exception as e:
            print(f'not a test')

        print(f'wave')
        spec_formatted['wave'] = np.array(spec['WAVELENGTH'][0])
        print(f'flux')
        spec_formatted['flux'] = np.array(spec['FLUX'][0])
        print('flag')
        spec_formatted['flag'] = np.array(spec['DQ'][0])
        

        if (spec['EXPTIME'].shape != (1,)):
            print(f'exptime loop 1')
            spec_formatted['exp_time'] = np.array(np.max(spec['EXPTIME']))
            
        else: 
            print(f'exptime loop 2')
            spec_formatted['exp_time'] = np.array(spec['EXPTIME'])
        
        print('err up')
        spec_formatted['error_upper'] = np.array(spec['ERROR'][0])
        spec_formatted['error_lower'] = np.array(spec['ERROR_LOWER'][0])
        print('err low')
        spec_formatted['G230L_error_up'] = np.zeros_like(spec_formatted['error_lower'])
        spec_formatted['G230L_error_down'] = np.zeros_like(spec_formatted['error_lower'])
        print('gcounts')
        spec_formatted['gcounts'] = np.array(spec['GCOUNTS'][0])
        print('background')
        spec_formatted['background'] = np.array(spec['BACKGROUND'][0])
        print('net')
        spec_formatted['net'] = np.array(spec['NET'][0])
        print('net_exp_array')
        net_exp_array = np.array(spec['NET'][0])*spec_formatted['exp_time'][0] 
        print('bkg_exp_array')
        bkg_exp_array = np.array(spec['BACKGROUND'][0])*spec_formatted['exp_time'][0] 
        print('variance flat')
        variance_flat = np.array(spec['VARIANCE_FLAT'][0])

        print('bkg times exp')
        spec_formatted['bkg_times_exp'] = np.nan_to_num(bkg_exp_array)
        print('net times exp')
        spec_formatted['net_times_exp_time'] = np.nan_to_num(net_exp_array)
        
        flags = spec_formatted['flag']
        mask = np.where((flags == 0) | (flags == 4))
        print(f'mask shape {mask}')
        masked_spec = spec_formatted[(mask)]
        print(f'made it through prep')
        
        return masked_spec 
    def run(self):
        try:
            # Process FITS files
            processed_data = []
            grating_types = []
            
            total_files = len(self.file_data)
            for i, (file_path, grating) in enumerate(self.file_data):
                print(f'have {i} of {total_files} files, grating {grating}')
                try:
                    if grating == 'G140L':
                        print(f'made it to grating ')
                        table = self.prepare(file_path)
                    else:
                        table = prepare_other_grating(file_path, grating=grating)
                    
                    print(f'now appending')
                    processed_data.append(table)
                    grating_types.append(grating)
                    print('appended successfully')
                    # Update progress (50% for loading and processing)
                    self.progress_updated.emit(int(50 * (i + 1) / total_files))
                except Exception as e:
                    self.error_signal.emit(f"Error processing file {file_path}: {str(e)}")
                    return
            
            # Separate by grating type
            G140L_tables = [data for data, grating in zip(processed_data, grating_types) 
                           if grating == 'G140L']
            G130M_tables = [data for data, grating in zip(processed_data, grating_types) 
                           if grating == 'G130M']
            
            # Handle case with mixed grating types or single grating type
            if G140L_tables and G130M_tables:
                # For mixed gratings, would need ranges for combining
                # This would typically be handled through UI input
                # For now, assuming default ranges
                G140L_table = vstack(G140L_tables)
                G130M_table = vstack(G130M_tables)
                
                # Default ranges for combining tables - these should ideally come from UI
                low1, high1 = 1100, 1150
                low2, high2 = 1150, 1450
                low3, high3 = 1450, 1800
                
                combined_table = combine_tables(
                    G140L_table, G130M_table,
                    low1, high1, low2, high2, low3, high3
                )
            elif G140L_tables:
                combined_table = vstack(G140L_tables)
            elif G130M_tables:
                combined_table = vstack(G130M_tables)
            else:
                self.error_signal.emit("No valid spectra found to coadd.")
                return
            
            # Update progress (75% after stacking)
            self.progress_updated.emit(75)
            
            # Coadd the spectra
            coadded_spectrum = coadd(combined_table, self.delta)
            
            # Update progress (100% after coadding)
            self.progress_updated.emit(100)
            
            # Emit the completed result
            self.finished_signal.emit(coadded_spectrum)
            
        except Exception as e:
            self.error_signal.emit(f"Error during coadding: {str(e)}")


class SpectraCoadderApp(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.file_data = []  # Store (file_path, grating_type) tuples
        self.coadded_spectrum = None
        
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('UV Spectra Coadder')
        self.setGeometry(100, 100, 1000, 700)
        
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)
        
        # Create sections
        self.createFileSection(main_layout)
        self.createPlotSection(main_layout)
        self.createParametersSection(main_layout)
        self.createActionsSection(main_layout)
        
        # Status bar and progress bar
        self.statusBar().showMessage('Ready')
        self.progress_bar = QProgressBar()
        main_layout.addWidget(self.progress_bar)
        
    def createFileSection(self, parent_layout):
        file_group = QGroupBox("Spectrum Files")
        file_layout = QVBoxLayout()
        
        # File list
        self.file_list = QListWidget()
        file_layout.addWidget(self.file_list)
        
        # Control buttons
        btn_layout = QHBoxLayout()
        
        self.add_file_btn = QPushButton("Add Files")
        self.add_file_btn.clicked.connect(self.addFiles)
        btn_layout.addWidget(self.add_file_btn)
        
        self.remove_file_btn = QPushButton("Remove Selected")
        self.remove_file_btn.clicked.connect(self.removeSelectedFile)
        btn_layout.addWidget(self.remove_file_btn)
        
        self.clear_files_btn = QPushButton("Clear All")
        self.clear_files_btn.clicked.connect(self.clearFiles)
        btn_layout.addWidget(self.clear_files_btn)
        
        file_layout.addLayout(btn_layout)
        file_group.setLayout(file_layout)
        parent_layout.addWidget(file_group)
        
    def createPlotSection(self, parent_layout):
        plot_group = QGroupBox("Spectrum Preview")
        plot_layout = QVBoxLayout()
        
        # Matplotlib Figure
        self.figure = Figure(figsize=(10, 4))
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel('Wavelength (Å)')
        self.ax.set_ylabel('Flux')
        self.ax.set_title('Coadded Spectrum Preview')
        
        plot_layout.addWidget(self.canvas)
        
        # Navigation toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        
        plot_group.setLayout(plot_layout)
        parent_layout.addWidget(plot_group)
        
    def createParametersSection(self, parent_layout):
        param_group = QGroupBox("Coadding Parameters")
        param_layout = QGridLayout()
        
        # Delta parameter
        param_layout.addWidget(QLabel("Delta (Å):"), 0, 0)
        self.delta_spin = QDoubleSpinBox()
        self.delta_spin.setRange(0.01, 10.0)
        self.delta_spin.setValue(0.5)
        self.delta_spin.setSingleStep(0.1)
        param_layout.addWidget(self.delta_spin, 0, 1)
        
        # Redshift
        param_layout.addWidget(QLabel("Redshift:"), 1, 0)
        self.redshift_input = QLineEdit("0.0")
        param_layout.addWidget(self.redshift_input, 1, 1)
        
        param_group.setLayout(param_layout)
        parent_layout.addWidget(param_group)
        
    def createActionsSection(self, parent_layout):
        actions_layout = QHBoxLayout()
        
        self.coadd_btn = QPushButton("Coadd Spectra")
        self.coadd_btn.clicked.connect(self.startCoadding)
        actions_layout.addWidget(self.coadd_btn)
        
        self.save_fits_btn = QPushButton("Save as FITS")
        self.save_fits_btn.clicked.connect(lambda: self.saveSpectrum('fits'))
        self.save_fits_btn.setEnabled(False)
        actions_layout.addWidget(self.save_fits_btn)
        
        self.save_npz_btn = QPushButton("Save as NPZ")
        self.save_npz_btn.clicked.connect(lambda: self.saveSpectrum('npz'))
        self.save_npz_btn.setEnabled(False)
        actions_layout.addWidget(self.save_npz_btn)
        
        self.export_csv_btn = QPushButton("Export File List")
        self.export_csv_btn.clicked.connect(self.exportFileListToCSV)
        actions_layout.addWidget(self.export_csv_btn)
        
        parent_layout.addLayout(actions_layout)

    def exportFileListToCSV(self):
        """Export the current file list to a CSV file"""
        if not self.file_data:
            QMessageBox.warning(self, "No Files", "No files to export.")
            return
            
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getSaveFileName(
            self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        
        if not csv_file:
            return
            
        if not csv_file.lower().endswith('.csv'):
            csv_file += '.csv'
            
        try:
            import csv
            with open(csv_file, 'w', newline='') as f:
                writer = csv.writer(f)
                
                # First row is redshift
                writer.writerow([self.redshift_input.text()])
                
                # Remaining rows are file path and grating type
                for file_path, grating_type in self.file_data:
                    writer.writerow([file_path, grating_type])
                    
            self.statusBar().showMessage(f'Exported file list to {csv_file}')
            QMessageBox.information(self, "Export Successful", f"File list exported to {csv_file}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error exporting to CSV: {str(e)}")

    def createFileSection(self, parent_layout):
        file_group = QGroupBox("Spectrum Files")
        file_layout = QVBoxLayout()
        
        # File list
        self.file_list = QListWidget()
        file_layout.addWidget(self.file_list)
        
        # Control buttons
        btn_layout = QHBoxLayout()
        
        self.add_file_btn = QPushButton("Add Files")
        self.add_file_btn.clicked.connect(self.addFiles)
        btn_layout.addWidget(self.add_file_btn)
        
        self.import_csv_btn = QPushButton("Import from CSV")
        self.import_csv_btn.clicked.connect(self.importFromCSV)
        btn_layout.addWidget(self.import_csv_btn)
        
        self.remove_file_btn = QPushButton("Remove Selected")
        self.remove_file_btn.clicked.connect(self.removeSelectedFile)
        btn_layout.addWidget(self.remove_file_btn)
        
        self.clear_files_btn = QPushButton("Clear All")
        self.clear_files_btn.clicked.connect(self.clearFiles)
        btn_layout.addWidget(self.clear_files_btn)
        
        file_layout.addLayout(btn_layout)
        file_group.setLayout(file_layout)
        parent_layout.addWidget(file_group)

    def importFromCSV(self):
        """Import file paths and grating types from a CSV file"""
        options = QFileDialog.Options()
        csv_file, _ = QFileDialog.getOpenFileName(
            self, "Open CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options
        )
        
        if not csv_file:
            return
            
        try:
            import csv
            with open(csv_file, 'r') as f:
                reader = csv.reader(f)
                
                # First line should contain redshift (optional)
                try:
                    first_line = next(reader)
                    if len(first_line) == 1:
                        try:
                            redshift = float(first_line[0])
                            self.redshift_input.setText(str(redshift))
                            # Skip to next line which should have file info
                            first_file_line = next(reader)
                            rows = [first_file_line] + list(reader)
                        except ValueError:
                            # If first line isn't a redshift number, treat it as file info
                            rows = [first_line] + list(reader)
                    else:
                        #jj mplconnect
                        # First line has multiple columns, treat as file info
                        rows = [first_line] + list(reader)
                except StopIteration:
                    # Empty file
                    QMessageBox.warning(self, "Empty CSV", "The CSV file appears to be empty.")
                    return
                    
                # Process file info rows
                files_added = 0
                for row in rows:
                    if len(row) >= 2:  # Need at least path and grating
                        file_path = row[0].strip()
                        grating_type = row[1].strip()
                        
                        # Verify the file exists
                        if not os.path.isfile(file_path):
                            QMessageBox.warning(
                                self, "File Not Found", 
                                f"File not found: {file_path}\nSkipping this entry."
                            )
                            continue
                            
                        # Validate grating type
                        if grating_type not in ['G140L', 'G130M']:
                            reply = QMessageBox.question(
                                self, "Invalid Grating Type", 
                                f"Grating type '{grating_type}' for file {os.path.basename(file_path)} "
                                f"is not recognized. Use G140L instead?",
                                QMessageBox.Yes | QMessageBox.No
                            )
                            
                            grating_type = 'G140L' if reply == QMessageBox.Yes else 'G130M'
                        
                        self.file_data.append((file_path, grating_type))
                        self.file_list.addItem(f"{os.path.basename(file_path)} ({grating_type})")
                        files_added += 1
                        
                if files_added > 0:
                    self.statusBar().showMessage(f'Added {files_added} files from CSV')
                else:
                    QMessageBox.warning(self, "No Valid Files", "No valid file entries found in the CSV.")
                    
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error importing from CSV: {str(e)}")
    def addFiles(self):
        options = QFileDialog.Options()
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FITS Files", "", "FITS Files (*.fits);;All Files (*)", options=options
        )
        
        if file_paths:
            for file_path in file_paths:
                # Ask for grating type
                reply = QMessageBox.question(
                    self, "Grating Type", 
                    f"Is {os.path.basename(file_path)} a G140L grating?",
                    QMessageBox.Yes | QMessageBox.No
                )
                
                if reply == QMessageBox.Yes:
                    grating_type = 'G140L'
                else:
                    grating_type = 'G130M'  # Default to G130M if not G140L
                
                self.file_data.append((file_path, grating_type))
                self.file_list.addItem(f"{os.path.basename(file_path)} ({grating_type})")
                
            self.statusBar().showMessage(f'Added {len(file_paths)} files')
    def removeSelectedFile(self):
        selected_rows = self.file_list.selectedIndexes()
        if not selected_rows:
            return
            
        for index in sorted(selected_rows, reverse=True):
            del self.file_data[index.row()]
            self.file_list.takeItem(index.row())
    
    def clearFiles(self):
        self.file_data = []
        self.file_list.clear()
        
    def startCoadding(self):
        if not self.file_data:
            QMessageBox.warning(self, "No Files", "Please add FITS files first.")
            return
            
        # Get delta value
        delta = self.delta_spin.value()
        
        # Start coadding in a separate thread
        self.progress_bar.setValue(0)
        self.statusBar().showMessage('Coadding spectra...')
        
        self.coadder_thread = CoadderThread(self.file_data, delta)
        self.coadder_thread.progress_updated.connect(self.progress_bar.setValue)
        self.coadder_thread.finished_signal.connect(self.onCoaddingFinished)
        self.coadder_thread.error_signal.connect(self.onCoaddingError)
        self.coadder_thread.start()
        
    def onCoaddingFinished(self, coadded_spectrum):
        self.coadded_spectrum = coadded_spectrum
        self.statusBar().showMessage('Coadding complete')
        
        # Plot the result
        self.plotCoaddedSpectrum()
        
        # Enable save buttons
        self.save_fits_btn.setEnabled(True)
        self.save_npz_btn.setEnabled(True)
        
    def onCoaddingError(self, error_message):
        self.statusBar().showMessage('Error during coadding')
        QMessageBox.critical(self, "Error", error_message)
        self.progress_bar.setValue(0)
        
    def plotCoaddedSpectrum(self):
        if self.coadded_spectrum is None:
            return
            
        self.ax.clear()
        
        # Plot spectrum
        wave = self.coadded_spectrum['wave']
        flux = self.coadded_spectrum['flux']
        error_up = self.coadded_spectrum['error_up']
        error_down = self.coadded_spectrum['error_down']
        
        self.ax.plot(wave, flux, color='black', drawstyle='steps-mid', label='Flux')
        self.ax.plot(wave, (error_up + error_down)/2, color='gray', drawstyle='steps-mid', label='Error')
        
        # Add redshift markers if available
        try:
            redshift = float(self.redshift_input.text())
            if redshift > 0:
                # Common emission lines
                line_wavelengths = [
                    1215.67,  # Lya
                    1031.92, 1037.61,  # O VI
                    1238.82, 1242.8,  # N V
                    1393.75, 1402.77,  # Si IV
                    1548.19, 1550.77  # C IV
                ]
                line_labels = [
                    'Lyα', 'O VI 1031', 'O VI 1037', 'N V 1238',
                    'N V 1242', 'Si IV 1393', 'Si IV 1402', 'C IV 1548', 'C IV 1550'
                ]
                
                colors = ['r', 'g', 'b', 'orange', 'purple', 'pink', 'yellow', 'darkblue', 'teal']
                
                for i, (wavelength, label) in enumerate(zip(line_wavelengths, line_labels)):
                    observed = wavelength * (1 + redshift)
                    if observed > min(wave) and observed < max(wave):
                        self.ax.axvline(x=observed, linestyle='--', color=colors[i % len(colors)], 
                                       alpha=0.7, label=label)
                
        except ValueError:
            pass  # Invalid redshift input, skip line markers
            
        self.ax.set_xlabel('Wavelength (Å)')
        self.ax.set_ylabel('Flux')
        self.ax.set_title('Coadded Spectrum')
        self.ax.legend(loc='upper right', fontsize='small')
        
        self.canvas.draw()
        
    def saveSpectrum(self, format_type):
        if self.coadded_spectrum is None:
            QMessageBox.warning(self, "No Data", "No coadded spectrum to save.")
            return
            
        try:
            redshift = float(self.redshift_input.text())
        except ValueError:
            QMessageBox.warning(self, "Invalid Redshift", "Please enter a valid redshift value.")
            return
            
        options = QFileDialog.Options()
        if format_type == 'fits':
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save FITS File", "", "FITS Files (*.fits)", options=options
            )
            if file_path:
                if not file_path.lower().endswith('.fits'):
                    file_path += '.fits'
                self.saveFitsFile(file_path, redshift)
                
        elif format_type == 'npz':
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save NPZ File", "", "NPZ Files (*.npz)", options=options
            )
            if file_path:
                if not file_path.lower().endswith('.npz'):
                    file_path += '.npz'
                self.saveNpzFile(file_path, redshift)
    
    def saveFitsFile(self, file_path, redshift):
        try:
            # Create columns for the FITS file
            col1 = fits.Column(name='Wavelength', format='E', array=self.coadded_spectrum['wave'])
            col2 = fits.Column(name='Flux', format='E', array=self.coadded_spectrum['flux'])
            col3 = fits.Column(name='Error_Up', format='E', array=self.coadded_spectrum['error_up'])
            col4 = fits.Column(name='Error_Down', format='E', array=self.coadded_spectrum['error_down'])
            col5 = fits.Column(name='Gcounts', format='E', array=self.coadded_spectrum['gcounts'])
            
            # Create HDU
            hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
            hdu.header['REDSHIFT'] = redshift
            hdu.writeto(file_path, overwrite=True)
            
            self.statusBar().showMessage(f'Saved spectrum to {file_path}')
            QMessageBox.information(self, "File Saved", f"Spectrum saved to {file_path}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving FITS file: {str(e)}")
    
    def saveNpzFile(self, file_path, redshift):
        try:
            # Save as NPZ file
            np.savez(
                file_path,
                wave=self.coadded_spectrum['wave'],
                flux=self.coadded_spectrum['flux'],
                error_up=self.coadded_spectrum['error_up'],
                error_down=self.coadded_spectrum['error_down'],
                gcounts=self.coadded_spectrum['gcounts'],
                redshift=redshift
            )
            
            self.statusBar().showMessage(f'Saved spectrum to {file_path}')
            QMessageBox.information(self, "File Saved", f"Spectrum saved to {file_path}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving NPZ file: {str(e)}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = SpectraCoadderApp()
    main_window.show()
    sys.exit(app.exec_())