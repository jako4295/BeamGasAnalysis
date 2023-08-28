# BeamGasAnalysis
In the CERN accelerator complex the ion beams suffer from losses. In this repository we study the losses coming from interactions between the ion beams and residual gas molecules. These effects is modelled by the rest gas collision cross section, which consists of electron loss cross section, and electron capture cross section. Hence, we have 
$$
\sigma = \sigma_{EC} + \sigma_{EL}
$$ 

[^1]

## Usage
All code is run from the file `Executor.py` where different plots are made from the `Tests` class. The current methods implemented in the `Tests` class is the following.

### `Tests.get_lifetimes_from_sigma_estimates()`
This method calculates a lifetime for a given ring type (LEIR, PS, or SPS). This is done using the 



[^1]: [Schlachter Formula](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.27.3372)