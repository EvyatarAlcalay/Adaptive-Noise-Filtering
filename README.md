# Adaptive Noise Filtering (LMS / NLMS / RLS)

A signal processing project focused on **adaptive noise filtering** using algorithms such as LMS, NLMS, and RLS.  
Includes MATLAB implementations, experiments on real audio files, performance plots, and NRdB evaluation.

## 📂 Repository Structure


## ⚙️ Requirements
- MATLAB R2020a or later  
- Built-in functions such as `audioread`, `audiowrite`, `toeplitz`  
- Example audio files placed in `data/raw` (airplane, cafe, city, vacuumcleaner)

## ▶️ How to Run
1. Open MATLAB in the repository root.  
2. Add source code to the path:
   ```matlab
   addpath(genpath('src/matlab'))
3. Run the experiment scripts:
   -adaptivepredict.m (Adaptive prediction function)

## 📊 Main Experiments

-Steepest Descent: Step-size effect and weight error convergence.
-LMS/NLMS: Convergence speed and misadjustment for different parameters.
-RLS: Effect of forgetting factor (λ) and initialization (δ).
-Audio Denoising: Instantaneous power plots, NRdB calculations, and side-by-side comparison of noisy vs denoised signals.

## 📈 Results
-Figures saved to results/figures/
-Denoised audio saved to results/audio/
   

