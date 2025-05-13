import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy import stats
from scipy import integrate
from scipy.optimize import curve_fit
import lmfit
from statsmodels import robust
import emcee
import corner
import extinction
from lmfit import Model, Parameters
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Any, Sequence
import pandas as pd
import os


c_kms = 2.998e5

def fwhm_ang_to_velocity(fwhm_ang, central_wavelength):
    """
    Convert FWHM from Angstroms to km/s
    
    Parameters:
    fwhm_ang: FWHM in Angstroms
    central_wavelength: Central wavelength of the line in Angstroms
    
    Returns:
    fwhm_velocity: FWHM in km/s
    """
    c = 2.998e5  # speed of light in km/s
    fwhm_velocity = (fwhm_ang / central_wavelength) * c
    print(f'FWHM velocity of {fwhm_velocity} km/s')
    return fwhm_velocity
def find_FWHM(lower_bound, upper_bound, dereddened_table, plot=True):
    """
    Find the Full Width at Half Maximum (FWHM) of a spectral line
    using improved interpolation.
    """
    from scipy import interpolate
    
    # Create initial mask
    mask = np.where((dereddened_table['wave'] >= lower_bound) & 
                    (dereddened_table['wave'] <= upper_bound))
    masked_table = dereddened_table[mask]
    
    # Create a finer wavelength grid for smooth interpolation
    wave_interp = np.linspace(min(masked_table['wave']), 
                             max(masked_table['wave']), 
                             1000)
    
    # Use cubic spline interpolation for smoother curve
    spline = interpolate.splrep(masked_table['wave'], masked_table['flux'], k=3)
    flux_interp = interpolate.splev(wave_interp, spline)
    
    # Find maximum
    max_idx = np.argmax(flux_interp)
    max_val = flux_interp[max_idx]
    max_wave = wave_interp[max_idx]
    half_max = max_val/2
    
    # Find crossing points
    above_half = flux_interp >= half_max
    crossing_indices = np.where(np.diff(above_half))[0]
    
    if len(crossing_indices) >= 2:
        left_idx = crossing_indices[0]
        right_idx = crossing_indices[-1]
        
        # Get precise crossing points through linear interpolation at crossings
        left_x = [wave_interp[left_idx], wave_interp[left_idx + 1]]
        left_y = [flux_interp[left_idx], flux_interp[left_idx + 1]]
        right_x = [wave_interp[right_idx], wave_interp[right_idx + 1]]
        right_y = [flux_interp[right_idx], flux_interp[right_idx + 1]]
        
        half_max_left = np.interp(half_max, left_y, left_x)
        half_max_right = np.interp(half_max, right_y, right_x)
    else:
        print("Warning: Could not find clear half-maximum crossing points")
        return None
    
    fwhm = half_max_right - half_max_left
    print(f'FWHM of {fwhm:.2f} angstroms')
    
    if plot:
        plt.figure(figsize=(10, 6))
        # Original data points
        plt.plot(masked_table['wave'], masked_table['flux'], 
                'o', label='Data points', markersize=5)
        
        # Interpolated curve
        plt.plot(wave_interp, flux_interp, '-', 
                label='Interpolated', alpha=0.7)
        
        # Half maximum line
        plt.axhline(y=half_max, color='r', linestyle=':', 
                   label='Half maximum', alpha=0.7)
        
        # FWHM points
        plt.plot([half_max_left, half_max_right], 
                [half_max, half_max], 'go', 
                label='FWHM points', markersize=8)
        
        # Peak point
        plt.plot(max_wave, max_val, 'ro', 
                label='Peak', markersize=8)
        
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.title('FWHM Measurement')
        
    return fwhm

def continuum(wave, m, c):
  return m*wave + c


def genericgaussian(wave, rest_wavelength, redshift, sigma_kms, flux):
  wave_observed = rest_wavelength*(1+redshift)
  sigma_Ang = sigma_kms/c_kms*wave_observed
  peak = flux/np.sqrt(2*sigma_Ang**2*np.pi)
  gaussian = peak*np.exp(-(wave - wave_observed)**2/2/sigma_Ang**2)


@dataclass
class DoubletInfo:
    """Class to hold information about a spectral line doublet"""
    secondary_wavelength: float  # Rest wavelength of the secondary line
    ratio: float  # Ratio of secondary to primary line flux (e.g., 1/2.95 for [O III])

@dataclass
class FitComponent:
    """Class to hold fit component parameters"""
    redshift: float
    sigma: float
    flux: float
    name: str

@dataclass
class SpectralLine:
    """Class to hold spectral line parameters"""
    rest_wavelength: float
    name: str
    components: List[str]
    doublet: Optional[DoubletInfo] = None
    geocoronal: bool = False
    include_in_fit: bool = True  # Flag to include this line in the fit
    def set_include_in_fit(self, include: bool):
        """Set whether to include this line in the fit"""
        self.include_in_fit = include
    
@dataclass
class FluxResult:
    """Class to hold flux measurement results"""
    total: float
    err_low: float
    err_high: float
    components: Dict[str, Tuple[float, float, float]]  # component name -> (value, err_low, err_high)
    samples: np.ndarray  # For corner plots if needed

class SpectralFitter:
    """Class for fitting multiple Gaussian components to spectral lines"""
    
    def __init__(self, z: float, wavelength_range: Tuple[float, float]):
        """
        Initialize the fitter
        
        Args:
            z: Initial redshift guess
            wavelength_range: (min_wavelength, max_wavelength) for fitting
        """
        self.z = z
        self.wave_min, self.wave_max = wavelength_range
        
    def _gaussian(self, wave: np.ndarray, rest_wave: float, 
                redshift: float, sigma_kms: float, flux: float,
                doublet: Optional[DoubletInfo] = None, geocoronal: Optional[bool] = False, ratio: Optional[float] = None) -> np.ndarray:
        """
        Calculate Gaussian profile in flux units, handling doublets if specified
        
        Args:
            wave: Wavelength array
            rest_wave: Rest wavelength of primary line
            redshift: Redshift
            sigma_kms: Velocity dispersion in km/s
            flux: Flux of primary line
            doublet: Optional DoubletInfo for doublet lines
        """
        # Convert velocity dispersion to wavelength units
        center = rest_wave * (1 + redshift)


        if geocoronal:
            center = rest_wave
        sigma_wave = rest_wave * sigma_kms / 3e5
        
        # Primary line
        gaussian = (flux / (sigma_wave * np.sqrt(2*np.pi))) * \
                np.exp(-0.5 * ((wave - center) / sigma_wave)**2)
        
        # Add secondary line for doublets
        if doublet is not None:
            rest_wave_2 = doublet.secondary_wavelength
            center_2 = rest_wave_2 * (1 + redshift)
            sigma_wave_2 = rest_wave_2 * sigma_kms / 3e5
            
            # Add secondary line with specified ratio
            if ratio is not None:
                flux_2 = flux * ratio
            else:
                flux_2 = flux * doublet.ratio
            gaussian_2 = (flux_2 / (sigma_wave_2 * np.sqrt(2*np.pi))) * \
                        np.exp(-0.5 * ((wave - center_2) / sigma_wave_2)**2)
            gaussian += gaussian_2
        
        return gaussian
    
    def _continuum(self, wave: np.ndarray, slope: float, intercept: float) -> np.ndarray:
        """Linear continuum model"""
        return slope * wave + intercept

    def create_model(self, lines: List[SpectralLine], 
                    parameter_constraints: Dict = None, ratio: Optional[bool] = False) -> Tuple[lmfit.Model, lmfit.Parameters]:
        """
        Create lmfit model and parameters for multiple emission lines
        
        Args:
            lines: List of SpectralLine objects to fit
            parameter_constraints: Dictionary of parameter constraints
                Format: {
                    'parameter_name': {
                        'value': initial_value,
                        'min': min_value,
                        'max': max_value,
                        'vary': bool,
                        'expr': expression
                    }
                }
            
        Returns:
            model: lmfit Model object
            params: lmfit Parameters object with initial guesses
        """
        if parameter_constraints is None:
            parameter_constraints = {}
            
        def model_func(wave, m, c, **kwargs):
            """Combined model of continuum + emission lines"""
            total = self._continuum(wave, m, c)
            
            # Add each line component
            for line in lines:
                if not line.include_in_fit:
                    continue
                for comp in line.components:
                    z_param = f"z_{comp}"
                    sigma_param = f"sigma_{comp}"
                    flux_param = f"flux_{line.name}_{comp}"
                    
                    if ratio:
                        try:
                            ratio_param = f"ratio_{line.name}"
                            total += self._gaussian(
                                wave, 
                                line.rest_wavelength,
                                kwargs[z_param], 
                                kwargs[sigma_param],
                                kwargs[flux_param],
                                doublet=line.doublet,
                                geocoronal=line.geocoronal,
                                ratio = kwargs[ratio_param]
                            )
                        except Exception:
                            pass
                    else:
                        total += self._gaussian(
                            wave, 
                            line.rest_wavelength,
                            kwargs[z_param], 
                            kwargs[sigma_param],
                            kwargs[flux_param],
                            doublet=line.doublet,
                            geocoronal=line.geocoronal,
                        )
            return total
            
        self._current_model_func = model_func
        # Create parameters
        params = lmfit.Parameters()
        
        # Add continuum parameters
        params.add('m', value=parameter_constraints.get('m', {}).get('value', 0.0),
                  vary=parameter_constraints.get('m', {}).get('vary', True))
        params.add('c', value=parameter_constraints.get('c', {}).get('value', 0.0),
                  vary=parameter_constraints.get('c', {}).get('vary', True))

        # Add parameters for each line component
        for line in lines:
            try:
                ratio_name = f"ratio_{line.name}"
                parameter_constraints.get(ratio_name).get('value')
                params.add(ratio_name, value=parameter_constraints.get(ratio_name, {}).get('value', 0.0),
                        vary=parameter_constraints.get(ratio_name, {}).get('vary', True),
                        min=parameter_constraints.get(ratio_name, {}).get('min', 0.0),
                        max=parameter_constraints.get(ratio_name, {}).get('max', 1.0))
            except Exception:
                pass
            for comp in line.components:
                # Redshift parameter (shared between lines)
                z_name = f"z_{comp}"
                if z_name not in params:
                    z_constraints = parameter_constraints.get(z_name, {})
                    params.add(z_name, 
                             value=z_constraints.get('value', self.z),
                             min=z_constraints.get('min', self.z-0.01),
                             max=z_constraints.get('max', self.z+0.01),
                             vary=z_constraints.get('vary', True),
                             expr=z_constraints.get('expr', None))
                
                
                # Velocity dispersion (shared between lines)    
                sigma_name = f"sigma_{comp}"
                if sigma_name not in params:
                    sigma_constraints = parameter_constraints.get(sigma_name, {})
                    params.add(sigma_name,
                             value=sigma_constraints.get('value', 200.0),
                             min=sigma_constraints.get('min', 0),
                             max=sigma_constraints.get('max', 1000),
                             vary=sigma_constraints.get('vary', True),
                             expr=sigma_constraints.get('expr', None))
                    
                # Flux (separate for each line)
                flux_name = f"flux_{line.name}_{comp}"
                flux_constraints = parameter_constraints.get(flux_name, {})
                params.add(flux_name,
                         value=flux_constraints.get('value', 100.0),
                         min=flux_constraints.get('min', 0),
                         vary=flux_constraints.get('vary', True),
                         expr=flux_constraints.get('expr', None))
                
        return Model(model_func), params
    
    def fit_spectrum(self, wave: np.ndarray, flux: np.ndarray, 
                    error: np.ndarray, lines: List[SpectralLine],
                    parameter_constraints: Optional[Dict] = None, ratio: Optional[bool] = False) -> Any:
        """
        Fit the spectrum with multiple emission lines
        
        Args:
            wave: Wavelength array
            flux: Flux array
            error: Error array
            lines: List of SpectralLine objects to fit
            parameter_constraints: Dictionary of parameter constraints
            
        Returns:
            result: lmfit ModelResult object
        """
        # Create mask for fitting region
        mask = (wave > self.wave_min) & (wave < self.wave_max)
        mask &= np.isfinite(flux) & np.isfinite(error)
        
        # Create model and parameters
        model, params = self.create_model(lines, parameter_constraints, ratio=ratio)
        
        # Perform fit
        result = model.fit(flux[mask], wave=wave[mask], 
                         weights=1/error[mask], params=params)
        
        return result

    def mcmc_fit(self, wave: np.ndarray, flux: np.ndarray, 
                error: np.ndarray, lines: List[SpectralLine],
                initial_fit: Any = None,
                steps: int = 1000, 
                burn: int = 100, 
                thin: int = 15,
                is_weighted: bool = True,
                parameter_constraints: Optional[Dict] = None, 
                ratio: Optional[bool] = False) -> Any:
        """
        Perform MCMC fit to the spectrum
        
        Args:
            wave: Wavelength array
            flux: Flux array
            error: Error array
            lines: List of SpectralLine objects to fit
            initial_fit: Optional initial fit result to use as starting point
            steps: Number of MCMC steps
            burn: Number of burn-in steps to discard
            thin: Thinning factor for chain
            is_weighted: Whether to use error weights in the fit
            parameter_constraints: Optional parameter constraints
            
        Returns:
            MCMC fit result
        """
        # Create mask for fitting region
        mask = (wave > self.wave_min) & (wave < self.wave_max)
        mask &= np.isfinite(flux) & np.isfinite(error)
        
        # Get model and parameters
        if initial_fit is not None:
            # Use parameters from initial fit
            emcee_params = initial_fit.params.copy()
        else:
            # Create fresh model and parameters
            model, emcee_params = self.create_model(lines, parameter_constraints, ratio=ratio)
            
        # Optional: add log-sigma parameter for error estimation
        # emcee_params.add('__lnsigma', value=np.log(0.1), 
        #                  min=np.log(0.001), max=np.log(2.0))
        
        # Setup MCMC keywords
        emcee_kws = dict(
            steps=steps,
            burn=burn,
            thin=thin,
            is_weighted=is_weighted
        )
        
        # Create model and run MCMC
        if not hasattr(self, '_current_model_func'):
            _, _ = self.create_model(lines, parameter_constraints)

        model = Model(self._current_model_func)

        result_emcee = model.fit(
            flux[mask],
            wave=wave[mask],
            weights=1/error[mask] if is_weighted else None,
            params=emcee_params,
            method='emcee',
            nan_policy='omit',
            fit_kws=emcee_kws
        )
        
        return result_emcee

    def analyze_line_flux(self, 
                         line: SpectralLine,
                         mcmc_result: Any,
                         plot_corner: bool = False) -> FluxResult:
        """
        Analyze the total flux and uncertainties for a given line
        
        Args:
            line: SpectralLine object to analyze
            mcmc_result: MCMC fit result
            plot_corner: Whether to create corner plot of flux parameters
            
        Returns:
            FluxResult object containing flux measurements and uncertainties
        """
        # Get flux parameter names for this line
        flux_keys = [f"flux_{line.name}_{comp}" for comp in line.components]
        
        # Get samples for all flux parameters
        flux_samples = mcmc_result.flatchain[flux_keys]
        
        # Calculate total flux samples and percentiles
        total_flux_samples = flux_samples.sum(axis=1)
        percentiles = np.percentile(total_flux_samples, [16, 50, 84])
        
        # Calculate component-wise percentiles
        component_results = {}
        for comp, key in zip(line.components, flux_keys):
            comp_percentiles = np.percentile(mcmc_result.flatchain[key], [16, 50, 84])
            component_results[comp] = (
                comp_percentiles[1],  # median
                comp_percentiles[1] - comp_percentiles[0],  # lower error
                comp_percentiles[2] - comp_percentiles[1]   # upper error
            )
        
        # Create FluxResult object
        result = FluxResult(
            total=percentiles[1],
            err_low=percentiles[1] - percentiles[0],
            err_high=percentiles[2] - percentiles[1],
            components=component_results,
            samples=flux_samples.values
        )
        
        # Create corner plot if requested
        if plot_corner:
            fig = corner.corner(
                flux_samples,
                labels=flux_keys,
                quantiles=[0.16, 0.5, 0.84],
                show_titles=True,
                title_fmt='.3f',
                title_kwargs={"fontsize": 12}
            )
            plt.show()
            
        return result

    def print_flux_summary(self, line: SpectralLine, 
                          flux_result: FluxResult,
                          initial_fit: Any = None):
        """Print a summary of the flux measurements"""
        print(f"\nResults for {line.name}:")
        print(f"Total flux: {flux_result.total:.3f} (+{flux_result.err_high:.3f}, -{flux_result.err_low:.3f})")
        
        print("\nComponent breakdown:")
        for comp, (val, err_low, err_high) in flux_result.components.items():
            print(f"{comp}: {val:.3f} (+{err_high:.3f}, -{err_low:.3f})")
            
        if initial_fit is not None:
            total_initial = sum(
                initial_fit.best_values[f"flux_{line.name}_{comp}"]
                for comp in line.components
            )
            print(f"\nInitial fit total flux: {total_initial:.3f}")

    def access_flux_value(self, line: SpectralLine, 
                          flux_result: FluxResult,
                          initial_fit: Any = None):
        """Print a summary of the flux measurements"""
        print(f"\nReturning results for {line.name}:")
        print(f"Total flux: {flux_result.total:.3f} (+{flux_result.err_high:.3f}, -{flux_result.err_low:.3f})")
        
        return flux_result.total, flux_result.err_high, flux_result.err_low


    def plot_fit(self, result, wave: np.ndarray, 
                flux: np.ndarray, error: np.ndarray, 
                lines: List[SpectralLine], ax=None,
                use_mcmc_median: bool = False,
                subtract_continuum: bool = False,
                show_line_markers: bool = True,
                save_file = None, 
                redshift = None,
                ylim = None):
        """Plot the fit results"""
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))
        
        # Get continuum parameters
        if use_mcmc_median:
            slope = np.median(result.flatchain['m'])
            intercept = np.median(result.flatchain['c'])
        else:
            slope = result.best_values['m']
            intercept = result.best_values['c']
        
        # Calculate continuum
        continuum = self._continuum(wave, slope, intercept)
        
        # Plot data and error (with or without continuum)
        plot_flux = flux - continuum if subtract_continuum else flux
        ax.errorbar(wave, plot_flux, yerr=error, color='black', 
                drawstyle='steps-mid', label='Data')
        
        # Plot best-fit model
        if use_mcmc_median:
            params = result.params.copy()
            for name in result.var_names:
                params[name].value = np.median(result.flatchain[name])
            model_fit = result.model.eval(params, wave=wave)
            if subtract_continuum:
                model_fit -= continuum
            ax.plot(wave, model_fit, color='red', label='Best fit (MCMC median)')
        else:
            model_fit = result.best_fit
            if subtract_continuum:
                model_fit -= continuum
            ax.plot(wave, model_fit, color='red', label='Best fit')
        
        # Define some colors for different components
        all_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        color_dict = {
            f"{line.name}_{comp}": all_colors[i + j*len(lines[0].components)]
            for i, comp in enumerate(lines[0].components)
            for j, line in enumerate(lines)
        }
        
        # Plot components
        for line in lines:
            for i, comp in enumerate(line.components):
                try:
                    # Get parameter values based on fit type
                    if use_mcmc_median:
                        try:
                            z = np.median(result.flatchain[f"z_{comp}"])
                            sigma = np.median(result.flatchain[f"sigma_{comp}"])
                        except Exception:
                            z= result.best_values[f"z_{comp}"]
                            sigma = result.best_values[f"sigma_{comp}"]
                        flux_val = np.median(result.flatchain[f"flux_{line.name}_{comp}"])
                    else:
                        z = result.best_values[f"z_{comp}"]
                        sigma = result.best_values[f"sigma_{comp}"]
                        flux_val = result.best_values[f"flux_{line.name}_{comp}"]
                    
                    # Plot both primary and secondary components in same color and style
                    model_primary = self._gaussian(
                        wave, 
                        line.rest_wavelength,
                        z, sigma, flux_val,
                        doublet=None
                    )
                    
                    if line.doublet is not None:
                        # For doublets, add both components before plotting
                        model_secondary = self._gaussian(
                            wave,
                            line.doublet.secondary_wavelength,
                            z, sigma, flux_val * line.doublet.ratio,
                            doublet=None
                        )
                        # Plot sum of both components
                        ax.plot(wave, model_primary + model_secondary, '--', 
                            color=color_dict[f"{line.name}_{comp}"], label=f'{line.name} {comp}')
                        if save_file is not None:
                            fits_file = pd.DataFrame({'wavelength': wave, f'{line.name}': model_primary + model_secondary})
                            fits_file.to_csv(f'{save_file}_{line.name}.csv', index=False)
                            print(f'Saving to {save_file}_{line.name}.csv')
                    else:
                        # For single lines, plot as normal
                        ax.plot(wave, model_primary, '--', 
                            color=color_dict[f"{line.name}_{comp}"], label=f'{line.name} {comp}')
                        if save_file is not None:
                            fits_file = pd.DataFrame({'wavelength': wave, f'{line.name}': model_primary})
                            fits_file.to_csv(f'{save_file}_{line.name}.csv', index=False)
                            print(f'Saving to {save_file}_{line.name}.csv')
                except Exception as e:
                    print(f'May be a mismatch in components between lines; carry on')
        try:
            if redshift == None:
                redshift = z
        except Exception:
            print('issue with redshift')
        if ylim is not None:
            ax.set_ylim(ylim[0], ylim[1])
        if show_line_markers:
            ymin, ymax = ax.get_ylim()
            text_y = ymax * 0.9  # Position for text labels
            
            for line in lines:
                # Plot primary line marker
                if line.geocoronal:
                    obs_wavelength = line.rest_wavelength
                else:
                    obs_wavelength = line.rest_wavelength * (1 + redshift)
                ax.axvline(obs_wavelength, color='gray', linestyle=':', alpha=0.5)
                ax.text(obs_wavelength, text_y, 
                    f'{line.name}\n{line.rest_wavelength:.1f}Å', 
                    rotation=90, ha='right', va='bottom')
                
                # Plot secondary line marker if it's a doublet
                if line.doublet is not None:
                    obs_wavelength_2 = line.doublet.secondary_wavelength * (1 + redshift)
                    ax.axvline(obs_wavelength_2, color='gray', linestyle=':', alpha=0.5)
                    ax.text(obs_wavelength_2, text_y,
                        f'{line.name}\n{line.doublet.secondary_wavelength:.1f}Å',
                        rotation=90, ha='right', va='bottom')

        ax.minorticks_on()
        ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
        ax.set_ylabel(r'$\rm flux\ [10^{-17}\ erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
        ax.legend()
        fig.tight_layout()
        return ax
    
@dataclass
class SpectralRegion:
    """Class to hold information about a spectral fitting region"""
    wave_min: float
    wave_max: float
    lines: List[SpectralLine]
    
class JointSpectralFitter(SpectralFitter):
    """Extension of SpectralFitter to handle joint fitting across multiple regions"""
    
    def __init__(self, z: float, regions: List[SpectralRegion]):
        """
        Initialize the joint fitter
        
        Args:
            z: Initial redshift guess
            regions: List of SpectralRegion objects defining the fitting regions
        """
        self.z = z
        self.regions = regions
        
    def create_joint_model(self, parameter_constraints: Dict = None) -> Tuple[lmfit.Model, lmfit.Parameters]:
        """
        Create model for joint fitting with separate continua for each region
        
        Args:
            parameter_constraints: Dictionary of parameter constraints
        """
        if parameter_constraints is None:
            parameter_constraints = {}
            
        def joint_model(wave, **kwargs):
            """Model function handling multiple regions with separate continua"""
            # Initialize output array
            total = np.zeros_like(wave)
            
            # Process each region
            for i, region in enumerate(self.regions):
                # Create mask for this region
                region_mask = (wave >= region.wave_min) & (wave <= region.wave_max)
                
                if not np.any(region_mask):
                    continue
                
                # Add continuum for this region
                m_param = f"m_{i}"
                c_param = f"c_{i}"
                total[region_mask] += self._continuum(
                    wave[region_mask], 
                    kwargs[m_param], 
                    kwargs[c_param]
                )
                
                # Add emission lines for this region
                for line in region.lines:
                    for comp in line.components:
                        z_param = f"z_{comp}"
                        sigma_param = f"sigma_{comp}"
                        flux_param = f"flux_{line.name}_{comp}"
                        
                        total[region_mask] += self._gaussian(
                            wave[region_mask],
                            line.rest_wavelength,
                            kwargs[z_param],
                            kwargs[sigma_param],
                            kwargs[flux_param],
                            doublet=line.doublet,
                            geocoronal=line.geocoronal
                        )
            
            return total
            
        self._current_model_func = joint_model
        
        # Create parameters
        params = lmfit.Parameters()
        
        # Add continuum parameters for each region
        for i in range(len(self.regions)):
            params.add(f'm_{i}', value=parameter_constraints.get(f'm_{i}', {}).get('value', 0.0),
                      vary=parameter_constraints.get(f'm_{i}', {}).get('vary', True))
            params.add(f'c_{i}', value=parameter_constraints.get(f'c_{i}', {}).get('value', 0.0),
                      vary=parameter_constraints.get(f'c_{i}', {}).get('vary', True))
        
        # Get unique components across all regions
        all_components = set()
        all_lines = []
        for region in self.regions:
            for line in region.lines:
                all_components.update(line.components)
                all_lines.append(line)
        
        # Add shared parameters for all components
        for comp in all_components:
            # Redshift parameter (shared between lines)
            z_name = f"z_{comp}"
            z_constraints = parameter_constraints.get(z_name, {})
            params.add(z_name,
                      value=z_constraints.get('value', self.z),
                      min=z_constraints.get('min', self.z-0.01),
                      max=z_constraints.get('max', self.z+0.01),
                      vary=z_constraints.get('vary', True),
                      expr=z_constraints.get('expr', None))
            
            # Velocity dispersion (shared between lines)
            sigma_name = f"sigma_{comp}"
            sigma_constraints = parameter_constraints.get(sigma_name, {})
            params.add(sigma_name,
                      value=sigma_constraints.get('value', 200.0),
                      min=sigma_constraints.get('min', 0),
                      max=sigma_constraints.get('max', 1000),
                      vary=sigma_constraints.get('vary', True),
                      expr=sigma_constraints.get('expr', None))
        
        # Add flux parameters for each line and component
        for line in all_lines:
            for comp in line.components:
                flux_name = f"flux_{line.name}_{comp}"
                flux_constraints = parameter_constraints.get(flux_name, {})
                params.add(flux_name,
                          value=flux_constraints.get('value', 100.0),
                          min=flux_constraints.get('min', 0),
                          vary=flux_constraints.get('vary', True),
                          expr=flux_constraints.get('expr', None))
        
        return Model(joint_model), params
    
    def fit_regions(self, wave: np.ndarray, flux: np.ndarray, 
                   error: np.ndarray, parameter_constraints: Optional[Dict] = None):
        """
        Fit all regions simultaneously
        """
        # Create mask for all fitting regions
        mask = np.zeros_like(wave, dtype=bool)
        for region in self.regions:
            mask |= (wave >= region.wave_min) & (wave <= region.wave_max)
        
        mask &= np.isfinite(flux) & np.isfinite(error)
        
        # Create model and parameters
        model, params = self.create_joint_model(parameter_constraints)
        
        # Perform fit
        result = model.fit(flux[mask], wave=wave[mask], 
                         weights=1/error[mask], params=params)
        
        return result
    
    def mcmc_joint_fit(self, wave: np.ndarray, flux: np.ndarray, 
                      error: np.ndarray, initial_fit: Any = None,
                      steps: int = 1000, burn: int = 100, 
                      thin: int = 15, is_weighted: bool = True,
                      parameter_constraints: Optional[Dict] = None):
        """
        Perform MCMC fit across all regions
        """
        # Create mask for all fitting regions
        mask = np.zeros_like(wave, dtype=bool)
        for region in self.regions:
            mask |= (wave >= region.wave_min) & (wave <= region.wave_max)
        
        mask &= np.isfinite(flux) & np.isfinite(error)
        
        # Get model and parameters
        if initial_fit is not None:
            emcee_params = initial_fit.params.copy()
        else:
            model, emcee_params = self.create_joint_model(parameter_constraints)
        
        # Setup MCMC keywords
        emcee_kws = dict(steps=steps, burn=burn, thin=thin, is_weighted=is_weighted)
        
        # Create model and run MCMC
        if not hasattr(self, '_current_model_func'):
            _, _ = self.create_joint_model(parameter_constraints)
            
        model = Model(self._current_model_func)
        
        result_emcee = model.fit(
            flux[mask],
            wave=wave[mask],
            weights=1/error[mask] if is_weighted else None,
            params=emcee_params,
            method='emcee',
            nan_policy='omit',
            fit_kws=emcee_kws
        )
        
        return result_emcee
    def calculate_equivalent_width(self, wave: np.ndarray, line_model: np.ndarray, continuum: np.ndarray) -> float:
        """
        Calculate the equivalent width of an emission line.
        
        Args:
            wave: Wavelength array
            line_model: Line emission model without continuum
            continuum: Continuum level
            
        Returns:
            equivalent_width: The equivalent width in wavelength units
        """
        # Calculate normalized flux (line / continuum)
        normalized_flux = line_model / continuum
        
        eq_width = -np.trapz(normalized_flux, wave)
        
        return eq_width
    def plot_joint_fit(self, result, wave: np.ndarray, flux: np.ndarray, 
                    error: np.ndarray, use_mcmc_median: bool = True,
                    subtract_continuum: bool = False, show_components: bool = True, save_file=None,
                    figsize=(15, 10), labelit: bool = True, calculate_ew: bool = False):
        """
        Plot joint fit results showing each region
        
        Args:
            result: Fit result from fit_regions or mcmc_joint_fit
            wave: Wavelength array
            flux: Flux array
            error: Error array
            use_mcmc_median: Whether to use MCMC median parameters (if available)
            subtract_continuum: Whether to subtract continuum from data and model
            show_components: Whether to show individual line components
            figsize: Figure size tuple
        """
        # Create figure with a subplot for each region
        n_regions = len(self.regions)
        fig, axes = plt.subplots(n_regions, 1, figsize=figsize)
        if n_regions == 1:
            axes = [axes]
        equivalent_widths = {}
        # Plot each region
        for i, (region, ax) in enumerate(zip(self.regions, axes)):
            # Create mask for this region
            mask = (wave >= region.wave_min) & (wave <= region.wave_max)
            region_wave = wave[mask]
            region_flux = flux[mask]
            region_error = error[mask]
            
            # Get parameters (either MCMC median or best fit)
            if use_mcmc_median and hasattr(result, 'flatchain'):
                params = result.params.copy()
                for name in result.var_names:
                    params[name].value = np.median(result.flatchain[name])
            else:
                params = result.params
            
            # Get continuum for this region
            slope = params[f'm_{i}'].value
            intercept = params[f'c_{i}'].value
            continuum = slope * region_wave + intercept
            
            # Plot data
            plot_flux = region_flux - continuum if subtract_continuum else region_flux
            # ax.errorbar(region_wave, plot_flux, yerr=region_error, 
            #         color='black', drawstyle='steps-mid', alpha=0.5, label='Data')
            
            # Calculate and plot total model
            model_comp = {}
            total_model = np.zeros_like(region_wave)
            
            # Add continuum to model if not subtracting
            if not subtract_continuum:
                total_model += continuum
            
            # Add lines
            for line in region.lines:
                for comp in line.components:
                    z = params[f'z_{comp}'].value
                    sigma = params[f'sigma_{comp}'].value
                    flux_val = params[f'flux_{line.name}_{comp}'].value
                    
                    line_model = self._gaussian(
                        region_wave,
                        line.rest_wavelength,
                        z, sigma, flux_val,
                        doublet=line.doublet,
                        geocoronal=line.geocoronal
                    )
                    
                    total_model += line_model
                    
                    if show_components:
                        model_comp[f'{line.name}_{comp}'] = line_model
                        # if name == 'Ha_broad' or name == 'Ha_narrow':
                        #     print(f'sum of {line.name}_{comp} is {np.sum([x for x in line_model])}')
                        #     if np.sum([x for x in line_model])==0:
                        #         print(f'everything is 0 for {line.name}_{comp}')

                if save_file is not None:
                        fits_file = pd.DataFrame({'wavelength': region_wave, f'{line.name}': total_model - continuum})
                        fits_file.to_csv(f'{save_file}_{line.name}.csv', index=False)
                        print(f'Saving to {save_file}_{line.name}.csv')
            # Calculate equivalent width for each line component
            for line in region.lines:
                for comp in line.components:
                    z = params[f'z_{comp}'].value
                    sigma = params[f'sigma_{comp}'].value
                    flux_val = params[f'flux_{line.name}_{comp}'].value
        
                    line_model = self._gaussian(
                        region_wave,
                        line.rest_wavelength,
                        z, sigma, flux_val,
                        doublet=line.doublet,
                        geocoronal=line.geocoronal
                    )
                    
                    # Calculate equivalent width
                    if calculate_ew:
                        ew = self.calculate_equivalent_width(region_wave, line_model, continuum)
                        equivalent_widths[f'{line.name}_{comp}_EW'] = ew
                    
            # if save_file is not None:
            #     # Add equivalent widths to saved files
            #     fits_file = pd.DataFrame({
            #         'wavelength': region_wave, 
            #         f'{line.name}': line_model,
            #         f'{line.name}_EW': [ew] * len(region_wave)  # Constant value across all rows
            #     })
            #     fits_file.to_csv(f'{save_file}_{line.name}.csv', index=False)
            #     print(f'Saving to {save_file}_{line.name}.csv with EW = {ew:.2f} Å')
            # Plot total model
            # ax.plot(region_wave, total_model, 'r-', label='Model', lw=2)
            
            # Plot components if requested
            if show_components:
                colors = plt.cm.tab10(np.linspace(0, 1, len(model_comp)))
                for (name, comp_model), color in zip(model_comp.items(), colors):
                    ax.plot(region_wave, comp_model + (0 if subtract_continuum else continuum),
                        '--', color=color, label=name, alpha=0.7)
                    # if name == 'Ha_broad' or name == 'Ha_narrow':
                    #     print(f'plotting model component {name}')
                    #     print(f'max is {np.max(comp_model)}')
                    #     print(f'{region_wave[np.where(np.max(comp_model))]}, {comp_model[np.where(np.max(comp_model))]}')
            
            # Add wavelength markers for lines
            if not subtract_continuum:
                ymin, ymax = ax.get_ylim()
                text_y = ymax * 0.9
                
                for line in region.lines:
                    obs_wavelength = line.rest_wavelength * (1 + self.z)
                    if obs_wavelength >= region.wave_min and obs_wavelength <= region.wave_max:
                        ax.axvline(obs_wavelength, color='gray', linestyle=':', alpha=0.5)
                        if labelit:
                            ax.text(obs_wavelength, text_y,
                                f'{line.name}\n{line.rest_wavelength:.1f}Å',
                                rotation=90, ha='right', va='bottom')
            
            # Customize plot
            ax.set_xlim(region.wave_min, region.wave_max)
            ax.set_xlabel('Wavelength (Å)')
            ax.set_ylabel('Flux')
            ax.legend(ncol=2)
            ax.minorticks_on()
            
        plt.tight_layout()
        if calculate_ew:
            print(f'saving equivalent widths')
            return fig, axes, equivalent_widths
        
        else:
            print(f'not saving equivalent widths')
            return fig, axes
    def get_gaussian_models(self, result, wave: np.ndarray, use_mcmc_median: bool = True):
        """
        Calculate Gaussian models for plotting without creating a new figure
        
        Args:
            result: Fit result from fit_regions or mcmc_joint_fit
            wave: Wavelength array
            use_mcmc_median: Whether to use MCMC median parameters (if available)
            
        Returns:
            Dictionary containing model data for each region and component
        """
        models_data = {}
        
        # Process each region
        for i, region in enumerate(self.regions):
            # Create mask for this region
            mask = (wave >= region.wave_min) & (wave <= region.wave_max)
            region_wave = wave[mask]
            
            # Get parameters (either MCMC median or best fit)
            if use_mcmc_median and hasattr(result, 'flatchain'):
                params = result.params.copy()
                for name in result.var_names:
                    params[name].value = np.median(result.flatchain[name])
            else:
                params = result.params
            
            # Get continuum for this region
            slope = params[f'm_{i}'].value
            intercept = params[f'c_{i}'].value
            continuum = slope * region_wave + intercept
            
            # Calculate total model and components
            models_data[f'region_{i}'] = {
                'wave': region_wave,
                'continuum': continuum,
                'total_model': np.zeros_like(region_wave) + continuum,  # Start with continuum
                'components': {}
            }
            
            # Add lines
            for line in region.lines:
                # Initialize line model for this specific line
                line_model = np.zeros_like(region_wave)
                
                for comp in line.components:
                    z = params[f'z_{comp}'].value
                    sigma = params[f'sigma_{comp}'].value
                    flux_val = params[f'flux_{line.name}_{comp}'].value
                    
                    # Calculate component model
                    component_model = self._gaussian(
                        region_wave,
                        line.rest_wavelength,
                        z, sigma, flux_val,
                        doublet=line.doublet,
                        geocoronal=line.geocoronal
                    )
                    
                    # Add to total model
                    models_data[f'region_{i}']['total_model'] += component_model
                    
                    # Store component model
                    component_key = f'{line.name}_{comp}'
                    models_data[f'region_{i}']['components'][component_key] = {
                        'model': component_model,
                        'params': {
                            'z': z,
                            'sigma': sigma,
                            'flux': flux_val
                        }
                    }
                    
                    # Add to line-specific model
                    line_model += component_model
                
                # Store combined line model (all components for this line)
                models_data[f'region_{i}'][f'line_{line.name}'] = line_model
            
        return models_data
    def print_all_line_fluxes(self, mcmc_result):
        """
        Print fluxes and errors for all lines in all regions from MCMC results
        
        Args:
            mcmc_result: Result from mcmc_joint_fit
        """
        print("\nMCMC Flux Results:")
        print("-" * 60)
        
        # Process each region
        for region in self.regions:
            for line in region.lines:
                # Get flux parameter names for this line
                flux_keys = [f"flux_{line.name}_{comp}" for comp in line.components]
                
                # Get samples for all flux parameters
                flux_samples = mcmc_result.flatchain[flux_keys]
                
                # Calculate total flux samples and percentiles
                total_flux_samples = flux_samples.sum(axis=1)
                percentiles = np.percentile(total_flux_samples, [16, 50, 84])
                
                # Print total flux
                print(f"\n{line.name}:")
                print(f"Total flux: {percentiles[1]:.2e} (+{percentiles[2]-percentiles[1]:.2e}, -{percentiles[1]-percentiles[0]:.2e})")
                
                # Print components
                print("Components:")
                for comp in line.components:
                    flux_key = f"flux_{line.name}_{comp}"
                    comp_percentiles = np.percentile(mcmc_result.flatchain[flux_key], [16, 50, 84])
                    print(f"  {comp}: {comp_percentiles[1]:.2e} (+{comp_percentiles[2]-comp_percentiles[1]:.2e}, -{comp_percentiles[1]-comp_percentiles[0]:.2e})")
        
        print("\nLine component parameters:")
        print("-" * 60)
        
        # Get unique components
        all_components = set()
        for region in self.regions:
            for line in region.lines:
                all_components.update(line.components)
                
        # Print parameters for each component
        for comp in sorted(all_components):
            # Get redshift
            z_key = f"z_{comp}"
            z_percentiles = np.percentile(mcmc_result.flatchain[z_key], [16, 50, 84])
            
            # Get sigma
            sigma_key = f"sigma_{comp}"
            sigma_percentiles = np.percentile(mcmc_result.flatchain[sigma_key], [16, 50, 84])
            
            print(f"\n{comp} component:")
            print(f"  z = {z_percentiles[1]:.5f} (+{z_percentiles[2]-z_percentiles[1]:.5f}, -{z_percentiles[1]-z_percentiles[0]:.5f})")
            print(f"  σ = {sigma_percentiles[1]:.1f} km/s (+{sigma_percentiles[2]-sigma_percentiles[1]:.1f}, -{sigma_percentiles[1]-sigma_percentiles[0]:.1f})")

    def get_line_flux(self, mcmc_result, line_name):
        """
        Get flux and errors for a specific line by name
        
        Args:
            mcmc_result: Result from mcmc_joint_fit
            line_name: Name of the line (e.g., 'Ha', 'Hb', 'OIII')
            
        Returns:
            dict with flux and error values, or None if line not found
        """
        # Find the line in any region
        target_line = None
        for region in self.regions:
            for line in region.lines:
                if line.name == line_name:
                    target_line = line
                    break
            if target_line:
                break
        
        if not target_line:
            print(f"Warning: Line '{line_name}' not found in any region")
            return None
        
        # Get flux parameter names for this line
        flux_keys = [f"flux_{line_name}_{comp}" for comp in target_line.components]
        
        # Get samples for all flux parameters
        flux_samples = mcmc_result.flatchain[flux_keys]
        
        # Calculate total flux samples and percentiles
        total_flux_samples = flux_samples.sum(axis=1)
        percentiles = np.percentile(total_flux_samples, [16, 50, 84])
        
        result = {
            'flux': percentiles[1],
            'error_up': percentiles[2] - percentiles[1],
            'error_down': percentiles[1] - percentiles[0],
            'components': {}
        }
        
        # Add component-wise results
        for comp in target_line.components:
            flux_key = f"flux_{line_name}_{comp}"
            comp_percentiles = np.percentile(mcmc_result.flatchain[flux_key], [16, 50, 84])
            result['components'][comp] = {
                'flux': comp_percentiles[1],
                'error_up': comp_percentiles[2] - comp_percentiles[1],
                'error_down': comp_percentiles[1] - comp_percentiles[0]
            }
        
        return result
    
    def get_line_ew(self, ew_result, line_name):
        """
        Get equivalent width for a specific line by name, combining components
        
        Args:
            ew_result: Dictionary of equivalent widths from plot_joint_fit
            line_name: Name of the line (e.g., 'Ha', 'Hb', 'OIII')
            
        Returns:
            dict with total equivalent width and component-wise values
        """
        # Find the line in any region
        target_line = None
        for region in self.regions:
            for line in region.lines:
                if line.name == line_name:
                    target_line = line
                    break
            if target_line:
                break
        
        if not target_line:
            print(f"Warning: Line '{line_name}' not found in any region")
            return None
            
        # Get EW values for this line's components
        component_ews = {}
        total_ew = 0
        
        for comp in target_line.components:
            ew_key = f'{line_name}_{comp}_EW'
            if ew_key in ew_result:
                component_ews[comp] = ew_result[ew_key]
                total_ew += ew_result[ew_key]
        
        # Create result dictionary
        result = {
            'ew': total_ew,
            'components': component_ews
        }
        
        return result