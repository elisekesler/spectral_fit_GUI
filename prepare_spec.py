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
from astropy.stats import poisson_conf_interval as pcf

def prepare(file):
    """Nicely formats file into columns. Also does bitmasking"""
    
    spec = Table.read(file)

    spec_formatted = Table()

    if spec['DQ'][0] == 55:
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
    else:

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
        variance_flat = np.array(spec['VARIANCE_FLAT'][0])

        spec_formatted['bkg_times_exp'] = np.nan_to_num(bkg_exp_array)
        spec_formatted['net_times_exp_time'] = np.nan_to_num(net_exp_array)
        
        flags = spec_formatted['flag']
        mask = np.where((flags == 0) | (flags == 4))
        masked_spec = spec_formatted[(mask)]
    
    return masked_spec 


# In[64]:


def prepare_other_grating(file, type = 0, grating = 'G230L'):
    """Nicely formats file into columns. Also does bitmasking"""
    
    spec = Table.read(file)
#     print(spec['FLUX'][0], spec['FLUX'][1] )
#     print(spec['WAVELENGTH'][0], spec['WAVELENGTH'][1])

    spec_formatted = Table()
    
    range_arrays = [0, 1]
    spec_array = ['WAVELENGTH', 'FLUX', 'DQ', 'ERROR', 'ERROR_LOWER', 'GCOUNTS', 'BACKGROUND', 'NET', 'VARIANCE_FLAT']
    formatted = []
#     print(spec.colnames)
    if type == 0:
        for column in spec_array: 
            for i in range_arrays:
#                 print(column)
                formatted.append(np.array(spec[column][i]))
#                 print(spec[column][i])
    if type == 1:
        for column in spec_array: 
             formatted.append(np.array(spec[column][0]))  

    spec_formatted['wave'] = np.concatenate((formatted[0], formatted[1]))
    

    spec_formatted['flux'] = np.concatenate((formatted[2], formatted[3]))
    spec_formatted['flag'] = np.concatenate((formatted[4], formatted[5]))
    
    if (spec['EXPTIME'].shape != (1,)):
        spec_formatted['exp_time'] = np.array(np.max(spec['EXPTIME']))
        
    else: 
         spec_formatted['exp_time'] = np.array(spec['EXPTIME'])
            
    spec_formatted['error_upper'] = np.concatenate((formatted[6], formatted[7]))
#     print(formatted[6], formatted[7])
    spec_formatted['error_lower'] = np.concatenate((formatted[8], formatted[9]))
    if grating == 'G230L':
        spec_formatted['G230L_error_up'] = np.concatenate((formatted[6], formatted[7]))
    #     print(formatted[6], formatted[7])
        spec_formatted['G230L_error_down'] = np.concatenate((formatted[8], formatted[9]))
    elif grating == 'G130M':
        spec_formatted['error_upper'] = np.concatenate((formatted[6], formatted[7]))
        spec_formatted['error_lower'] = np.concatenate((formatted[8], formatted[9]))
        spec_formatted['G230L_error_up'] = np.zeros_like(spec_formatted['error_lower'])
        spec_formatted['G230L_error_down'] = np.zeros_like(spec_formatted['error_lower'])
        
    spec_formatted['gcounts'] = np.concatenate((formatted[10], formatted[11]))
    spec_formatted['background'] = np.concatenate((formatted[12], formatted[13]))
    spec_formatted['net'] = np.concatenate((formatted[14], formatted[15]))
    variance_flat = np.concatenate((formatted[16], formatted[17]))

    net_exp_array = np.array(spec_formatted['net'])*spec_formatted['exp_time'][0] 
    spec_formatted['net_times_exp_time'] = np.nan_to_num(net_exp_array)

    bkg_exp_array = np.array(spec_formatted['background'])*spec_formatted['exp_time'][0] 
    spec_formatted['bkg_times_exp'] = np.nan_to_num(bkg_exp_array)
    
    
    flags = spec_formatted['flag']
    mask = np.where((flags == 0) | (flags == 4))
    masked_spec = spec_formatted[(mask)]
    
    return masked_spec 
    


# In[78]:

def combine_tables(table1, table2, low1, high1, low2, high2, low3, high3):
    flux1 = table1['flux']
    flux2 = table2['flux']
    
    wave1 = table1['wave']
    wave2 = table2['wave']
    
    
    mask1 = np.where((((wave1 >= low1) & (wave1 < high1)) | ((wave1 > low3) & (wave1 <= high3))))
    mask2 = np.where((wave2 >= low2) & (wave2 < high2))
    
    new_table1 = table1[mask1]
    new_table2 = table2[mask2]
    
    final_table = vstack([new_table1, new_table2])
    

    
    return final_table

def coadd(combined_table, delta, subtract_background=False):
    """
    Coadds spectra from a combined table of exposures.

    Parameters
    ----------
    combined_table : stacked table
        Combined data from exposures with the same grating.
    delta : float
        Delta value between points. Can be adjusted as desired/needed.

    Returns
    -------
    coadded_spectrum : table
        New table with coadded flux and wavelength.
    """
    # Create a wavelength array over which to coadd
    min_wave = np.min(combined_table['wave'])
    max_wave = np.max(combined_table['wave'])
    wave = np.arange(min_wave, max_wave, delta)
    
    # Initialize columns of the dataset
    coadded_spectrum = Table()
    coadded_spectrum['wave'] = wave
    coadded_spectrum['flux'] = 0.0
    coadded_spectrum['exp_time'] = 0.0
    coadded_spectrum['gcounts'] = 0.0
    coadded_spectrum['number'] = 0.0
    coadded_spectrum['error_poisson_up'] = 0.0
    coadded_spectrum['error_poisson_down'] = 0.0
    coadded_spectrum['error_up'] = 0.0
    coadded_spectrum['error_down'] = 0.0
    coadded_spectrum['G230L_error_up'] = 0.0
    coadded_spectrum['G230L_error_down'] = 0.0
    coadded_spectrum['error_down'] = 0.0
    coadded_spectrum['background'] = 0.0
    coadded_spectrum['pcf_n'] = 0.0
    coadded_spectrum['gcounts_error_up'] = 0.0
    coadded_spectrum['gcounts_error_down'] = 0.0
    coadded_spectrum['net'] = 0.0
    coadded_spectrum['net_times_exp'] = 0.0
    coadded_spectrum['background_times_exp'] = 0.0
    coadded_spectrum['wavediff_ratio'] = 0.0

    #coadd with desired lambda
    for pixel in coadded_spectrum:
        thispixel = combined_table[(combined_table['wave'] >= (pixel['wave'] - delta/2)) & (combined_table['wave'] < (pixel['wave'] + delta/2))]
        pixel['flux'] = np.nanmean(thispixel['flux'])
        pixel['G230L_error_down'] = np.nanmean(thispixel['G230L_error_down'])
        pixel['G230L_error_up'] = np.nanmean(thispixel['G230L_error_up'])
        
        pixel['exp_time'] = np.sum(thispixel['exp_time'])
        pixel['gcounts'] = np.sum(thispixel['gcounts'])
        pixel['background'] = np.sum(thispixel['background'])
        pixel['background_times_exp'] = np.nanmean(thispixel['bkg_times_exp'])
        pixel['number'] = len(thispixel)
        pixel['net'] = np.sum(thispixel['net'])
        pixel['net_times_exp'] = np.sum(thispixel['net_times_exp_time'])
        if subtract_background:
            print(pixel['flux'], pixel['background_times_exp'], pixel['net'])
        # if np.isnan(pixel['net_times_exp']):
        #     exp_net = thispixel['net_times_exp_time']
        #     print(f'this pixel had a nan net times exp time. it has net exptime {exp_net}')
        
        poisson_conf = pcf(pixel['gcounts'], interval='frequentist-confidence')
        pixel['gcounts_error_up'] = poisson_conf[1] - pixel['gcounts']
        pixel['gcounts_error_down'] = pixel['gcounts'] - poisson_conf[0] 
        pixel['error_poisson_up'] = pixel['gcounts_error_up']/pixel['gcounts']*pixel['flux']
        pixel['error_poisson_down'] = pixel['gcounts_error_down']/pixel['gcounts']*pixel['flux']

        if (pixel['gcounts'] < 1) | (pixel['flux'] < 0.0):
              pixel['error_poisson_up'] = np.nan
              pixel['error_poisson_down'] = 0 

    good_errors = coadded_spectrum[np.isfinite(coadded_spectrum['error_poisson_up'])]
    print(good_errors['error_poisson_up'])
    for pixel in coadded_spectrum:
        if np.isnan(pixel['error_poisson_up']):
      
            # Find the nearest pixel with a good error estimate
            minIndex = np.argmin(np.abs(good_errors['wave'] - pixel['wave']))
            pixel_nearest_good_error = good_errors[minIndex]
            # Calculate the ratio of flux to gcounts for the nearest pixel with a good error
            ratio = pixel_nearest_good_error['flux']/pixel_nearest_good_error['gcounts']
        
            # Multiple this ratio by the upward going error in gcounts
            pixel['error_poisson_up'] = pixel['gcounts_error_up']*ratio

          # Record the wavelenght mis-match between the nearest good error pixel and this pixel. Need to make sure this isn't too big (a few  to ~10 angstroms is ok.

            pixel['wavediff_ratio'] = np.abs(pixel['wave'] - pixel_nearest_good_error['wave'])

        if np.isnan(pixel['net_times_exp']):
            print(f'pixel had nan net*exp, setting that one to 0')
            pixel['net_times_exp'] = 0.0

    return coadded_spectrum
