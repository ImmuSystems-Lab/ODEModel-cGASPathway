# ODEModel-cGASPathway-

This repository contains a computational model simulating the dynamics of the cGAS and subsequent JAK/STAT pathways. For more detail on the model, please refer to our associated [publication]( https://doi.org/10.1016/j.jtbi.2018.11.001).

In brief, this model simulates a response to a DNA transfection, which leads to the production of interferon, and activation of various feedback responses (e.g. IRF7). The model consists of 13 differential equations tracking intramolecular concentrations over time. 


$$
\begin{align}
    \frac{d\left[cGAS\right]}{dt}&=-k_{1f}\cdot cGAS\cdot DNA+k_{1r}\cdot cGASc\\
    \frac{d\left[DNA\right]}{dt}&=-k_{1f}\cdot cGAS\cdot DNA+k_{1r}\cdot cGASc-\frac{k_{cat2}\cdot TREX1\cdot DNA}{K_{m2}+DNA}\\
    \frac{d\left[cGAMP\right]}{dt}&=k_{3f}\cdot \left(cGAS_{tot}-cGAS\right)-k_{4f}\cdot cGAMP\cdot STING +k_{4r}\cdot \left(STING_{tot}-STING\right)-\tau_3\cdot cGAMP \\
    \frac{d\left[STING\right]}{dt}&=-k_{4f}\cdot cGAMP\cdot STING+k_{4r}\cdot \left(STING_{tot}-STING\right)\\
    \frac{d\left[IRF3\right]}{dt}&=-\frac{k_{cat5}\cdot IRF3\cdot \left(STING_{tot}-STING\right)}{K_{m5}+IRF3}+k_{5r}\cdot \left(IRF3_{tot}-IRF3\right)\\
    \frac{d\left[IFNβ m\right]}{dt}&=\frac{k_{cat6}\cdot \left(IRF3_{tot}-IRF3\right)}{K_{m6}+\left(IRF3_{tot}-IRF3\right)}+k_{6f}\cdot IRF7Pn-\tau_6\cdot IFNβ m\\
    \frac{d\left[IFNβ\right]}{dt}&=\frac{k_{cat7}\cdot IFNβ m}{K_{m7}+IFNβ m}-\tau_7\cdot IFNβ\\
    \frac{d\left[STATP2n\right]}{dt}&=\frac{k_{cat8}\cdot IFNβ}{{\rm Km}_8+IFNβ}\cdot \frac{1}{1+k_{8f}\cdot SOCSm}-\tau_8\cdot SOCSm\\
    \frac{d\left[SOCSm\right]}{dt}&=k_{9f}\cdot STATP2n-\tau_9\cdot SOCSm\\
    \frac{d\left[IRF7m\right]}{dt}&=k_{10f1}\cdot STATp2n+k_{10f2}\cdot IRF7Pn-\tau_{10}\cdot IRF7m\\
    \frac{d\left[TREX1m\right]}{dt}&=k_{11f}\cdot STATP2n-\tau_{11}\cdot TREX1m\\
    \frac{d\left[IRF7Pn\right]}{dt}&=k_{12f}\cdot IRF7m-\tau_{12}\cdot IRF7Pn\\
    \frac{d\left[TREX1\right]}{dt}&=k_{13f}\cdot TREX1m-\tau_{13}\cdot TREX1
\end{align}
$$

This model was trained on data from literature using MCMC to determine posterior probably distributions for all unknown parameters. The simulation file contains a parameter set that minimized model error. Use this file to simulate the model and generate plots of the dynamics.

![](DynamicsEx.png)
