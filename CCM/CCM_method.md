# SOME CONSIDERATIONS ABOUT THE USE OF THE _CONVERGENT CROSS MAPPING_ METHOD IN PALEOECOLOGY TIME SERIES

### Foreword

Our ultimate goal is to analyze the _structure_ of the networks of co-existence of plant species extracted from the time series of pollen abundance. To this aim, we need to _quantify_ the **causal relationships** between species.

In the paper of Ushio _et al._ [1] the authors use the **_Convergent Cross Mapping_** method (henceforth, CCM) to quantify the causal relationships between the time series of fish species' abundance (2 samples per month over 12 years = a total of 288 points).

### The CCM method

Ushio _et al._ [1] use the CCM method introduced by Sugihara _et al._ in [3]. The latter, in turn, introduce the method **not as a competitor** of the Granger causality (GC) [4] but, instead, as a method usable in those circumstances in which GC cannot be applied. What are the circumstances in which GC cannot be used?<br/>


* Non separability: namely that information about causative factors cannot be subtracted from the effects.
* Non linear systems with weak to modest interaction strengths.


The CCM estimates the causality, $\rho \in [0,1]$, between two observables $X$ and $Y$ by mapping their evolution into two manifolds (_i.e._, attractors) $M_X$ and $M_Y$. If the position of neighboring points in $M_X$ is also close in $M_Y$, then we can say that $X$ "cause" $Y$ (_i.e._, $\rho_{XY} \simeq 1$). It is worth noting that the relationship might not be symmetrical (_i.e._, $\rho_{XY} \neq \rho_{YX}$) because $X$ cause $Y$ but not the other way around. More in general, complex patterns of causality can exist. A complete set of possibilities can be seen in Fig. 4 of [3].

> [!WARNING]
> One caveat is that the CCM method is sensible to the length, $L$, of the time series used to estimate the causality.

### What is next?

In [1] the authors consider only one network (_i.e._, $L = 288$), but we need to consider several networks and our whole time series has only $140$ points. Hence, we have to split our time series into shorter ones. This poses a first question:

> [!IMPORTANT]
> <b>How long should be each subseries?</b>

which, in turn, calls for the second question:

> [!IMPORTANT]
> <b>Should we consider overlapped time series and, if so, how big should be the overlap?</b>

Leading to the last question:

> [!IMPORTANT]
> <b>How can we use Adam's code to generate the networks using CCM?</b>
 
> [!NOTE]
> In alternative to GC, one can use the so-called Transfer Entropy (TE) [5] which can be applied to nonlinear systems.


## BIBLIOGRAPHY

1. Ushio, M., Hsieh, C., Masuda, R., Deyle, E. R., Ye, H., Chang, C.-W., Sugihara, G., & Kondoh, M. (2018). _Fluctuating interaction network and time-varying stability of a natural fish community_. **Nature**, 554(7692), 360–363. (https://doi.org/10.1038/nature25504)


2. Park, J., Smith, C., Sugihara, G., Deyle, E., Saberski, E., & Ye, H. (2024) _rEDM: Empirical Dynamic Modeling ('EDM') package for R_. Available at: (https://cran.r-project.org/web/packages/rEDM/index.html) DOI:(https://doi.org/10.32614/CRAN.package.rEDM)


3. Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., & Munch, S. (2012). _Detecting Causality in Complex Ecosystems_. **Science**, 338(6106), 496–500. (https://doi.org/10.1126/science.1227079)


4. Granger, C. W. J., (1969). _Investigating causal relations by econometric models and cross-spectral methods_. **Econometrica** 37, 424.

    
5. Schreiber, T., (2000). _Measuring Information Transfer_. **Physical Review Letters**, 85, 461. (https://doi.org/10.1103/PhysRevLett.85.461).

