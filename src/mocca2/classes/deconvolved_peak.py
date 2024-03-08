from typing import List, Dict
from numpy.typing import NDArray

import numpy as np

from mocca2.classes import Peak, Component
from copy import deepcopy

class DeconvolvedPeak(Peak):
    """Information about peak and its deconvolved components"""

    components : List[Component]
    """Deconvolved components of the peak"""

    residual_mse: float
    """Residual MSE after deconvolution"""

    r2: float
    """R2 after deconvolution"""

    resolved: bool
    """This specifies, whether the deconvolution sufficiently explains the peak"""

    def __init__(self, peak: Peak, residual_mse: float, r2: float, resolved: bool, concentrations: NDArray = None, spectra: NDArray = None, json_reconstruct: bool = False):
        """
        Initializes DeconvolvedPeak based on the Peak and deconvolved data

        Parameters
        ----------
        
        concentrations: NDArray = None
            Concentrations of the individual components [component, time]
        
        spectra: NDArray = None
            Spectra of the individual components [component, wavelength]

        residual_mse: float
            Residual MSE after deconvolution

        r2: float
            R2 after deconvolution

        resolved: bool
            This specifies, whether the deconvolution sufficiently explains the peak

        json_reconstruct: bool
            This specifies, if the class is being reconstructed from a json file, in
            which case some additional functionality in init should not be ran

        """

        self.__dict__ = peak.__dict__

        self.residual_mse = residual_mse
        self.r2 = r2
        self.components = []
        self.resolved = resolved

        if not json_reconstruct:
            assert concentrations.shape[0] == spectra.shape[0], "The number of components in concentrations and spectra does not match"

            total_conc = np.sum(concentrations)

            for conc, spec in zip(concentrations, spectra):
                peak_frac = np.sum(conc)/total_conc
                self.components.append(Component(conc, spec, peak.left, peak_frac))

    def merge_same_components(self) -> None:
        """
        Merges all the components with identical ID (not None).
        
        Concentrations are added, spectra are averaged, weighted by concentration integral.
        """

        # Assign all components to a dict by compound_id
        components : Dict[int,List[Component]] = dict()
        for component in self.components:
            id = component.compound_id
            if id in components:
                components[id].append(component)
            elif id is not None:
                components[id] = [component]

        # Merge components with identical compound_id
        for compound_id, comps in components.items():
            if len(comps) == 1 or compound_id is None:
                continue

            integrals = np.array([c.integral for c in comps])
            spectra = np.array([c.spectrum for c in comps])
            
            conc = np.sum([c.concentration for c in comps], axis=0)
            spectrum = np.sum(spectra.T*integrals, axis=1)/np.sum(integrals)

            peak_frac = np.sum(c.peak_fraction for c in comps)

            merged = Component(conc, spectrum, self.left, peak_frac, compound_id)

            components[compound_id] = [merged]

        # Flatten the dict into list and save
        self.components = [component for comps in components.values() for component in comps]
        
    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        if json_dict['components']:
            json_dict['components'] = [comp.to_json() for comp in json_dict['components']]
        return json_dict

    def from_json(json_dict):

        peak_vars_init = ['left', 'right', 'maximum', 'height', 'prominence', 'all_maxima']
        deconvolved_init_vars = ['peak', 'concentrations', 'spectra', 'residual_mse', 'r2', 'resolved']
        instance_vars = ['components']

        peak_vars_dict = {}
        for var in peak_vars_init:
            if var in json_dict.keys():
                peak_vars_dict[var] = json_dict[var]

        deconvolved_init_vars_dict = {}
        for var in deconvolved_init_vars:
            if var in json_dict.keys():
                deconvolved_init_vars_dict[var] = json_dict[var]

        # Convert components into json depending on peak type (Peak or DeconvolvedPeak)
        if len(json_dict['components']) > 0:
            comps = [Component.from_json(comp) for comp in json_dict['components']]

        peak = Peak.from_json(peak_vars_dict)
        deconvolvedpeak = DeconvolvedPeak(peak, **deconvolved_init_vars_dict, json_reconstruct=True)
        deconvolvedpeak.components = comps

        return deconvolvedpeak
