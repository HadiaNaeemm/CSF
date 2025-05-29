# CSF
This repository contains all the codes used in the manuscript titled "[Hadia Naeem#, Chun-Wang Su#, Mei Ji, Nan Yao, Zi-Gang Huang*, and Celso Grebogi. Diverse Performance in Compressive Sensing-Based Reconstruction of an oscillatory dynamical system.Chaos Solitons & Fractals, 2025, 198:116508.]". The scripts and functions provided here were used for data processing, analysis, and figure generation.
### Code Overview
#### Section 2.1.2 (Figure 1)
- Phase Potrait of LSO (`LSO_Fig1`)
#### Section 3.1.1 (Figure 3)
- **SLSO_LO**, **SLSO_HO**, **SLSO_HO1**: Simulations for a single Landau-Stuart (LS) oscillator using:
  - Low-order basis functions (`SLSO_LO`)
  - Third-order basis functions (`SLSO_HO`)
  - Higher-order basis functions (`SLSO_HO1`)
- **CLSO_LO**, **CLSO_HO**, **CLSO_HO1**: Simulations for coupled LS oscillators using:
  - Low-order basis functions (`CLSO_LO`)
  - Third-order basis functions (`CLSO_HO`)
  - Higher-order basis functions (`CLSO_HO1`)

#### Section 3.2
- **HeatMap_SLSO** (Figure 4): Generates the heatmap of non-zero reconstruction error for a single oscillator across defined omega and alpha ranges.
- **LOX**, **LOY**, **LSO** (Figure 5): Compute reconstruction errors for the x and y time series of a single Lorenz oscillator and LS oscillator.
- **CLO**, **CLSO** (Figure 6): Compute reconstruction errors for the x₁ and x₂ time series of coupled Lorenz and LS oscillators.

#### Section 3.3.1 (Figure 8)
- **Osc1_**, **osc2**: Used to generate results for Figure 8.

#### Section 3.3.2 (Figure 10)
- **Osc1_**, **osc2**: Used to generate results for Figure 10.

#### Section 3.4 ( Figure 12)
- **C1_osc1**, **C1_osc2**, **C2_osc1**, **C2_osc2**: Used to produce results for Figure 12.

#### Figures 11, 13–17
- The corresponding `.mat` files for **Figure 11**, **Figure 13**, **Figure 14**, **Figure 15**, **Figure 16**, and **Figure 17** are used to reproduce the respective figures as presented in the manuscript.
