import numpy as np
from astropy.table import Table, vstack
import matplotlib as mpt
from matplotlib import pyplot as plt
import matplotlib as mp
from astropy.stats import mad_std
from scipy import stats
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import lmfit
from astropy.io.fits import getheader,getdata
from astropy.io import fits
from astropy.stats import poisson_conf_interval as pcf
from scipy.integrate import simpson as simps
from csv import writer
import extinction
from statsmodels import robust

class spec:
    def __init__(self, table, z, options='normal', dereddened = False):
        if options=='normal':
            self.table = table
            if dereddened == False:
                self.flux = table['flux']
            if dereddened == True:
                self.flux = table['dereddened']
            self.z = z
            self.wave = table['wave']
            self.exp_time = table['exp_time']
            self.gcounts = table['gcounts']
            
            for pixel in table:
                if pixel['G230L_error_up'] == 0:
                    pixel['error_up'] = pixel['error_poisson_up']
                elif pixel['G230L_error_up'] != 0:
                    pixel['error_up'] = pixel['G230L_error_up']    
                
                if pixel['G230L_error_down'] == 0:
                    pixel['error_down'] = pixel['error_poisson_down']
                elif pixel['G230L_error_down'] != 0:
                    pixel['error_down'] = pixel['G230L_error_down']
                    
                    
            self.error_up = table['error_up']
            self.error_down = table['error_down']
            
            
                
        
            self.lya = 1215.67*(1 + self.z)
            self.OVI1031 = 1031.92 *(1 + self.z)
            self.OVI1037 = 1037.61 *(1 + self.z)
        
            self.NV1238 = 1238.82*(1 + self.z)
            self.NV1242 = 1242.8*(1 + self.z)
            self.Si1393 = 1393.75*(1 + self.z)
            self.Si1402 = 1402.77*(1 + self.z)
            self.CIV1548 = 1548.19*(1 + self.z)
            self.CIV1550 = 1550.77*(1 + self.z)
            
        else:
            self.flux = table['flux']
            self.z = z
            self.wave = table['wave']
            self.error_up = table['errup']
            self.error_down = table['errdn']
            self.gcounts = table['gcounts']
            self.exp_time = 0.0
            
            self.lya = 1215.67*(1 + self.z)
            self.OVI1031 = 1031.92 *(1 + self.z)
            self.OVI1037 = 1037.61 *(1 + self.z)
        
            self.NV1238 = 1238.82*(1 + self.z)
            self.NV1242 = 1242.8*(1 + self.z)
            self.Si1393 = 1393.75*(1 + self.z)
            self.Si1402 = 1402.77*(1 + self.z)
            self.CIV1548 = 1548.19*(1 + self.z)
            self.CIV1550 = 1550.77*(1 + self.z)
    
        
        
    def plot(self,x_lower=1300, x_upper=1800, y_lower=-0.3e-15, y_upper=2e-15, style = 'filled'):
        fig, ax = plt.subplots(1, figsize=(10, 5))

        #data
        ax.plot(self.wave, self.flux, color='black', drawstyle='steps-mid')
        
        if style == ('filled'):
            ax.fill_between(self.wave, self.flux + self.error_up, self.flux - self.error_up, color = 'blue', alpha=0.3)
            ax.fill_between(self.wave, self.flux + self.error_down, self.flux - self.error_down, color = 'red', alpha=0.2)
        
        if style == ('lines'):
            
            ax.plot(self.wave, self.error_up - self.error_down, color='grey', drawstyle='steps-mid')
#             ax.plot(self.wave, self.flux + self.error_down, color='red', drawstyle='steps-mid')
#             ax.plot(self.wave, self.flux - self.error_up, color='blue', drawstyle='steps-mid')
#             ax.plot(self.wave, self.flux - self.error_down, color='red', drawstyle='steps-mid')
            
        # [Lyα] emission line
        ax.axvline(self.lya, color='red', linestyle='--', label="[Lya]")

        # [O VI] emission lines
        ax.axvline(self.OVI1031, color='maroon', linestyle='--', label="[O VI] 1031 ")
        ax.axvline(self.OVI1037, color='grey', linestyle=':', label="[O VI] 1037")

        # [N V] emission lines
        ax.axvline(self.NV1238, color='blue', linestyle=':', label="[N V] 1238")
        ax.axvline(self.NV1242, color='forestgreen', linestyle=':', label="[N V] 1242")


        # [Si IV]emission line
        ax.axvline(self.Si1393, color='orange', linestyle=':', label="[Si IV] 1393")
        ax.axvline(self.Si1402, color='darkblue', linestyle=':', label="[Si IV] 1402")


        # [C IV]emission lines
        ax.axvline(self.CIV1548, color='green', linestyle=':', label="[C IV] 1548")
        ax.axvline(self.CIV1550, color='teal', linestyle=':', label="[C IV] 1550")


        plt.legend(prop={'size': 10})

        ax.tick_params(axis='both', which='major', direction='in', top='on', bottom='on', left='on', right='on',labelsize=20, size=5)
        ax.tick_params(axis='both', which='minor', direction='in', top='on', bottom='on', left='on', right='on',size=3)
        ax.minorticks_on()

        ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$', size = 20)
        ax.set_ylabel(r'$\rm flux\ [1\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$', size = 20)
        ax.set_xlim([x_lower, x_upper])
        ax.set_ylim([y_lower, y_upper])
        
        
    def plot_lya(self, x_lower= 0, x_upper= 1600, y_lower=0, y_upper=2e-15, style='filled'):
        fig, ax = plt.subplots(1, figsize=(10, 5))

        #data
        ax.plot(self.wave, self.flux, color='black', drawstyle='steps-mid')
        
        if style == ('filled'):
            ax.fill_between(self.wave, self.flux + self.error_up, self.flux - self.error_up, color = 'blue', alpha=0.3)
            ax.fill_between(self.wave, self.flux + self.error_down, self.flux - self.error_down, color = 'red', alpha=0.2)
        
        if style == ('lines'):
            
            ax.plot(self.wave, self.flux + self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux + self.error_down, color='red', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_down, color='red', drawstyle='steps-mid')
        
        # [Lyα] emission line
        ax.axvline(self.lya, color='red', linestyle='--', label="[Lya]")
        plt.legend()


        ax.minorticks_on()
        ax.set_xlabel(r'$\rm{\rm observed\ wavelength\ [\AA]}$')
        ax.set_ylabel(r'$\rm {flux\ [1\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]}$')
        
        if x_lower == 0:
            x_lower = self.lya - 50
            x_upper = self.lya + 50
#             y_upper = np.max(self.flux) + (np.max(self.flux)*0.2)
            
        ax.set_xlim([x_lower, x_upper])
        ax.set_ylim([y_lower, y_upper])
        
    def plot_NV(self, x_lower= 0, x_upper= 1600, y_lower=0, y_upper=2e-15, style = 'filled'):
        fig, ax = plt.subplots(1, figsize=(10, 5))
        
        if style == ('filled'):
            ax.fill_between(self.wave, self.flux + self.error_up, self.flux - self.error_up, color = 'blue', alpha=0.3)
            ax.fill_between(self.wave, self.flux + self.error_down, self.flux - self.error_down, color = 'red', alpha=0.2)
        
        if style == ('lines'):
            
            ax.plot(self.wave, self.flux + self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux + self.error_down, color='red', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_down, color='red', drawstyle='steps-mid')
        

        #data
        ax.plot(self.wave, self.flux, color='black', drawstyle='steps-mid')

         # [N V] emission lines
        ax.axvline(self.NV1238, color='blue', linestyle=':', label="[N V] 1238")
        ax.axvline(self.NV1242, color='forestgreen', linestyle=':', label="[N V] 1242")

        ax.minorticks_on()
        ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
        ax.set_ylabel(r'$\rm flux\ [1\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
        
        if x_lower == 0:
            x_lower = ((self.NV1242 + self.NV1238)/2) - 50
            x_upper = ((self.NV1242 + self.NV1238)/2) + 50
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', direction='in', top='on', bottom='on', left='on', right='on',labelsize=20, size=5)
        ax.tick_params(axis='both', which='minor', direction='in', top='on', bottom='on', left='on', right='on',size=3)
        ax.set_xlim([x_lower, x_upper])
        ax.set_ylim([y_lower, y_upper])
        
    def plot_OVI(self, x_lower= 0, x_upper= 1600, y_lower=0, y_upper=1e-15, style = 'filled'):
        fig, ax = plt.subplots(1, figsize=(10, 5))
        
        if style == ('filled'):
            ax.fill_between(self.wave, self.flux + self.error_up, self.flux - self.error_up, color = 'blue', alpha=0.3)
            ax.fill_between(self.wave, self.flux + self.error_down, self.flux - self.error_down, color = 'red', alpha=0.2)
        
        if style == ('lines'):
            
            ax.plot(self.wave, self.flux + self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux + self.error_down, color='red', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_down, color='red', drawstyle='steps-mid')
        

        #data
        ax.plot(self.wave, self.flux, color='black', drawstyle='steps-mid')

        # [O VI] emission lines
        ax.axvline(self.OVI1031, color='maroon', linestyle='--', label="[O VI] 1031 ")
        ax.axvline(self.OVI1037, color='grey', linestyle=':', label="[O VI] 1037")

        ax.minorticks_on()
        ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
        ax.set_ylabel(r'$\rm flux\ [1\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
        
        if x_lower == 0:
            x_lower = ((self.OVI1031 + self.OVI1037)/2) - 50
            x_upper = ((self.OVI1031 + self.OVI1037)/2) + 50
            
        ax.set_xlim([x_lower, x_upper])
        ax.set_ylim([y_lower, y_upper])
        
    def plot_Si(self, x_lower= 0, x_upper= 1600, y_lower=0, y_upper=1e-15, style = 'filled'):
        fig, ax = plt.subplots(1, figsize=(10, 5))
        
        if style == ('filled'):
            ax.fill_between(self.wave, self.flux + self.error_up, self.flux - self.error_up, color = 'blue', alpha=0.3)
            ax.fill_between(self.wave, self.flux + self.error_down, self.flux - self.error_down, color = 'red', alpha=0.2)
        
        if style == ('lines'):
            
            ax.plot(self.wave, self.flux + self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux + self.error_down, color='red', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_down, color='red', drawstyle='steps-mid')
        
        #data
        ax.plot(self.wave, self.flux, color='black', drawstyle='steps-mid')

        ax.axvline(self.Si1393, color='orange', linestyle=':', label="[Si IV] 1393")
        ax.axvline(self.Si1402, color='darkblue', linestyle=':', label="[Si IV] 1402")

        ax.minorticks_on()
        ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
        ax.set_ylabel(r'$\rm flux\ [1\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
        
        if x_lower == 0:
            x_lower = ((self.Si1393 + self.Si1402)/2) - 50
            x_upper = ((self.Si1393 + self.Si1402)/2) + 50
            
        ax.set_xlim([x_lower, x_upper])
        ax.set_ylim([y_lower, y_upper])
        
    def plot_CIV(self, x_lower= 0, x_upper= 1600, y_lower=0, y_upper=1e-15, style='filled'):
        fig, ax = plt.subplots(1, figsize=(10, 5))
        
        if style == ('filled'):
            ax.fill_between(self.wave, self.flux + self.error_up, self.flux - self.error_up, color = 'blue', alpha=0.3)
            ax.fill_between(self.wave, self.flux + self.error_down, self.flux - self.error_down, color = 'red', alpha=0.2)
        
        if style == ('lines'):
            
            ax.plot(self.wave, self.flux + self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux + self.error_down, color='red', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_up, color='blue', drawstyle='steps-mid')
            ax.plot(self.wave, self.flux - self.error_down, color='red', drawstyle='steps-mid')
               
        
        #data
        ax.plot(self.wave, self.flux, color='black', drawstyle='steps-mid')

        # [C IV]emission lines
        ax.axvline(self.CIV1548, color='green', linestyle=':', label="[C IV] 1548")
        ax.axvline(self.CIV1550, color='teal', linestyle=':', label="[C IV] 1550")

        ax.minorticks_on()
        ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
        ax.set_ylabel(r'$\rm flux\ [1\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
        
        if x_lower == 0:
            x_lower = ((self.CIV1548 + self.CIV1550)/2) - 50
            x_upper = ((self.CIV1548 + self.CIV1550)/2) + 50
            
        ax.set_xlim([x_lower, x_upper])
        ax.set_ylim([y_lower, y_upper])
        
    def get_continuum(self, mask_low_1, mask_low_2, mask_high_1, mask_high_2, weights = None, option = 'normal'):
        
        values = {}
        if option == 'normal':
            mask = np.where(((self.wave >= mask_low_1) & (self.wave < mask_high_1)) | ((self.wave >= mask_low_2) & (self.wave < mask_high_2)))
        else:
            mask = np.where(((self.wave >= mask_low_1) & (self.wave < mask_high_1)))
        wavelength_array = self.wave[mask]
        flux_array =self.flux[mask]
        
        if (weights == None):
            weights = 1/self.error_up[mask]
            
        def linear(x, m, b):
            return x*m + b
        
        linear_model = lmfit.Model(linear)
        params_continuum = lmfit.Parameters()
        params_continuum.add('m', value = (np.nanmean(flux_array)/2))
        params_continuum.add('b', value = 0 )
        
        result = linear_model.fit(flux_array, x=wavelength_array, params=params_continuum, weights = weights)

        values['slope'] = result.best_values['m']
        values['intercept'] = result.best_values['b']
        m = values['slope']
        b = values['intercept'] 
        #print(f'slope is {m}, intercept is {b}')
        return values


    def get_flux(self, lower_wavelength_limit, upper_wavelength_limit, mask_low_1, mask_low_2, mask_high_1, mask_high_2, N=100, get_monte_errors = False, option = 'normal', print_out = True, error_array = None, array_errors = None):
        '''Gets flux for a given area. Intended to be used in tandem with a plot, so one can manually choose the zone to integrate over for a given spectral line. 
        Parameters
        ----------
        lower_wavelength_limit : integer
            value in angstroms of the estimated start of a spectral line; lower limit of the flux integral
        upper_wavelength_limit : integer
            value in angstroms of the estimated end of a spectral line; upper limit of the flux integral
        mask_low_1 : integer
            value in angstroms of lower part of one continuum portion to be subtracted. 
        mask_low_2 : integer
            value in angstroms of lower part of another continuum portion to be subtracted. 
        mask_high_1: integer
            value in angstroms of upper part of one continuum portion to be subtracted, corresponds with mask_low_1
        mask_high_2: integer
            value in angstroms of upper part of the other continuum portion to be subtracted, corresponds with mask_low_2

        Returns
        -------
        not done with this description
            
         '''
        if option == 'normal':
            continuum_values = self.get_continuum(mask_low_1, mask_low_2, mask_high_1, mask_high_2)
        else:
            continuum_values = self.get_continuum(mask_low_1, mask_low_2, mask_high_1, mask_high_2, option = 'not normal')
        wavelength_array = self.wave 
        mask = ((wavelength_array) >= (lower_wavelength_limit)) & ((wavelength_array) < (upper_wavelength_limit))
        
        #determine an estimation for sigma by getting a flux per pixel from continuum area
        mask_for_flux_per_pixel = (((wavelength_array) >= (mask_low_1)) & ((wavelength_array) < (mask_high_1)) | ((wavelength_array) >= (mask_low_2)) & ((wavelength_array) >= (mask_high_2)))

        
        
#         for pixel in flux_per_pixel:
#             deviation = (pixel-mean_flux_per_pixel)**2
#             deviation_from_mean.append(deviation)
            
#         sigma_estimate = np.sum(deviation_from_mean)/len(deviation_from_mean)
            
            
        
        
        gcounts_array = np.array(self.gcounts[mask])
        wave_array = np.array(self.wave[mask])
        flux_array = np.array(self.flux[mask])
        exp_array = np.array(self.exp_time[mask])
        otherflux = np.array(self.flux[mask] + self.error_up[mask])
        
        if error_array != None:
            flux_array = array_errors
        
        def linear(x, m, b):
            return x*m + b
        
        slope = continuum_values['slope']
        intercept = continuum_values['intercept']
        continuum = linear(wave_array, slope, intercept)

        flux_new = flux_array - continuum
        otherflux_new = otherflux - continuum

        
        other_index = ((wavelength_array) >= (lower_wavelength_limit - 20)) & ((wavelength_array) < (upper_wavelength_limit + 20))
        other_wave = self.wave[other_index]
        other_continuum = linear(other_wave, slope, intercept)
        flux_area_new = self.flux[other_index] - other_continuum
        errors = np.array(self.error_up[mask])
        poisson_table = Table()
        poisson_table['gcounts'] = gcounts_array
        poisson_table['wave'] = wave_array
        poisson_table['flux']=flux_new
        poisson_table['exp_time']=exp_array
        shape = np.shape(poisson_table['exp_time'])
        poisson_table['number']=np.zeros(shape)
#         poisson_table['gcounts_error_up']=np.zeros(shape)
#         poisson_table['gcounts_error_down']=np.zeros(shape)
#         poisson_table['error_poisson_up']=np.zeros(shape)
#         poisson_table['error_poisson_down']=np.zeros(shape)
        
        
        #coadd with desired lambda
        poisson_table.write('poisson_table.txt', overwrite=True, format='ascii.fixed_width')
#         plt.plot(poisson_table['wave'], poisson_table['gcounts'])
#         plt.savefig('photoofcounts.jpg')
#         plt.close()
        
        lim = np.max(flux_new)
        flux = simps(flux_new, wave_array)
        otherflux_value = simps(otherflux_new, wave_array)
        
        
        if print_out == True:
            fig, ax = plt.subplots(1, figsize=(10, 5))
            ax.plot(wave_array, (flux_new), color='red', drawstyle='steps-mid')
            ax.plot(other_wave, flux_area_new, color='black', alpha=0.8, drawstyle='steps-mid')
            ax.set_xlim((lower_wavelength_limit - 20), (upper_wavelength_limit + 20))
            ax.set_ylim((0, (lim + (lim*0.25))))
            print(f'Calculated flux is {flux} erg cm^{-2} s^{-1}')
            
        
        if get_monte_errors == True:
              fluxes = []
              for flux in range(N):  
                  self.random = np.random.randn(int(len(flux_array)),)
  
                  array = flux_array + flux_array*self.random                 
                  result = self.get_flux(lower_wavelength_limit, upper_wavelength_limit, mask_low_1, mask_low_2, mask_high_1, mask_high_2, get_errors = False, print_out = False, error_array = True, array_errors = array)
                  fluxes.append(result)
      
    
              average = np.sum(fluxes)/N
              print(f'Monte carlo average sum of fluxes is {average:0.5}')
              error = robust.mad(fluxes)/np.sqrt(2)
              print(f'Monte carlo standard deviation estimate is {error:0.5}')
            
        total_counts = np.sum(poisson_table['gcounts'])
        print(f'total counts is {total_counts}')
        total_s = np.sum(poisson_table['exp_time'])
        poisson_conf = pcf(total_counts, interval='frequentist-confidence')
        gcounts_error_up = poisson_conf[1] - total_counts
        gcounts_error_down = total_counts - poisson_conf[0] 
        error_poisson_up = (gcounts_error_up/total_counts)*flux
        error_poisson_down = (gcounts_error_down/total_counts)*flux
        
        jackknife_fluxes = []
        
        for pixel in poisson_table:
            flux_of = [x for i,x in enumerate(flux_new) if x!=pixel['flux']]
            wave_jackknife = [x for i,x in enumerate(wave_array) if x!=pixel['wave']]
            pixel_flux = simps(flux_of, wave_jackknife)
            jackknife_fluxes.append(pixel_flux)
            

        jackknife_sigma2 = robust.mad(jackknife_fluxes)
        jackknife_sigma = mad_std(jackknife_fluxes)

        
        
        string = f'Poisson estimates give lower error of {error_poisson_down} and upper error of {error_poisson_up}. Jackknife estimate is {jackknife_sigma}, robust.mad gives {jackknife_sigma2}'
        print(string)
        return flux, error_poisson_up
        