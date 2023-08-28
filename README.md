# BeamGasAnalysis
In the CERN accelerator complex the ion beams suffer from losses. In this repository we study the losses coming from interactions between the ion beams and residual gas molecules. These effects is modelled by the rest gas collision cross section, which consists of electron loss (EL) cross section, and electron capture (EC) cross section. Hence, we have $\sigma = \sigma_{EC} + \sigma_{EL}$. The EL cross section is calculated using a semi-emperical formula combining the studies from [^1] and [^2]. The EC is calculated using the Schlachter formula[^3]. 

The cross section is inversely proportional to the lifetime of an ion in the following sense, [^4]

$$
\tau = \frac{1}{\beta  c\sigma n},
$$

where $\beta$ is the relativistic beta ($\beta=v/c$ with $v$ being the velocity of the ion), $c$ is the speed of light, $\sigma$ is the combined EL and EC cross section, and $n$ is the molecular density of the gas in the beam pipe. Hence, it is relatively simple to go from cross section estimates to lifetime estimates and vica versa. From beam data it is also possible to estimate a lifetime for a given ion. This is done with the following beam intensity formula,

$$
I(t) = a\ \text{e}^{-bt},
$$

where $\tau = 1/b$ is the lifetime of the ion, and $a$ is a constant. By fitting this formula to the beam data, we can estimate the lifetime of the ion.

## Usage and Results
All code is run from the file `Executor.py` found in the `LifetimeExecutor`, where different plots are made from the `Tests` class. The current methods implemented in the `Tests` class is explained below.

### `Tests.get_lifetimes_from_sigma_estimates()`
This method calculates a lifetime for a given ring type (LEIR, PS, or SPS). The estimated cross sections are calculated using the methods referred to above [^1][^2][^3]. The lifetimes are calculated using the inversely proportional relation between lifetime and cross section. The plot that is obtained using this method is shown below (for PS).


### `Tests.get_sigma_estimates_with_varying_i_p()`
This method calculates the cross section for a given ring type (LEIR, PS, or SPS). The cross sections are calculated using the methods referred to above [^1][^2][^3]. The cross sections are calculated for different ionization potentials $I_p$ to see the impact of this paremeter.

### `Tests.get_sigma_estimates_with_varying_n_0()`
This method calculates the cross section for a given ring type (LEIR, PS, or SPS). The cross sections are calculated using the methods referred to above [^1][^2][^3]. The cross sections are calculated for different principle quantum numbers $n_0$ to see the impact of this paremeter.

### `Tests.get_lifetime_from_data()`
This method gets an estimate of the lifetime based on the beam data from `DataObject().elements`. The lifetime is estimated by fitting the beam intensity formula to the beam data. This lifetime is then compared to the cross section estimates in the cross section domain. Note that when the lifetime is found from the cross section then from the inverse proportionality we have an $n$ for each residual gas in the beam pipe. Hence we have a cross section comparison for each residual gas.

## Data
The data used for processing is found in the `DataHandler`, and the data is stored in the `DataObject` class. This class contains the following: 
- `gas_fractions`, which is percentage of the different gases in the beam pipe (also referred to as residual gasses). 
- `pressure_data`, which is the pressure in the beam pipe. It is assumed that the pressure is constant in the beam pipe.
- `projectile_data`, which constists of the atomic number $Z$, the charge $q$, the kinetic energy $E_{kin}$, the first ionization potential of the projectile $I_p$, the principle quantum number of the outermost projectile electron $n_0$, and the relativistic beta $\beta$. The projectile is the term used for the ion beam.
  - $I_p$ is found using the [NIST database](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) and converting the value for ionization energy to keV.
  - $n_0$ is found by sorting the electron configuration in ascending order according to the energy level of the atom in question. The energy level order can be found [here](https://chem.libretexts.org/Bookshelves/General_Chemistry/Map%3A_A_Molecular_Approach_(Tro)/08%3A_Periodic_Properties_of_the_Elements/8.03%3A_Electron_Configurations-_How_Electrons_Occupy_Orbitals) in Figure 8.3.1. Starting from the highest energies then we remove from the superscript. For instance the electron configuration for oxygen is $1s^2 2s^2 2p^4 \ \Rightarrow \ 1s^2 2s^2$, which gives $n_0=2$ (because this is the energy level $2s$). 
- `elements`, which is the beam data and is used to find the lifetime of an ion through fitting of the exponential function given above.

[^1]: [DuBois et al. Electron capture cross section](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.84.022702)
[^2]: [Shevelko et al. Electron loss cross section](https://www.sciencedirect.com/science/article/pii/S0168583X11003272)
[^3]: [Schlachter Formula. Electron loss cross section](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.27.3372)
[^4]: https://cds.cern.ch/record/902810/files/at-2005-013.pdf
