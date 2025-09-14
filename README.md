# Adaptive Noise Filtering (LMS / NLMS / RLS)

A signal processing project focused on **adaptive noise filtering** using algorithms such as LMS, NLMS, and RLS.  
Includes MATLAB implementations, experiments on real audio files, performance plots, and NRdB evaluation.

## 📂 Repository Structure
adaptive-noise-filtering/
├─ README.md
├─ docs/                  # Project PDF and additional documents
├─ data/
│  ├─ raw/                # Original audio files (e.g., airplane.wav, cafe.wav, city.wav, vacuumcleaner.wav)
│  └─ processed/          # Processed signals (created during experiments)
├─ src/
│  ├─ Adaptive_Filters.m      # Main script divided into sections (Q1–Q5)
│  ├─ adaptivepredict.m       # Prediction function for competition (Q6)


## ⚙️ Requirements
- MATLAB R2020a or later  
- Built-in functions such as `audioread`, `audiowrite`, `toeplitz`  
- Example audio files placed in `data/raw` (airplane, cafe, city, vacuumcleaner)  

## ▶️ How to Run
### 1) Running `Adaptive_Filters.m` by Sections:

- Open MATLAB and set the Current Folder to the repository root. Make sure the `.wav` files are accessible from this folder
- Open `Adaptive_Filters.m`
- The script is divided into Sections (marked with `%%`). Each section corresponds to a specific question. Run them one by one using **Run Section** (Ctrl+Enter)

#### Mapping Sections:

- `%% The code for section 4,5 in Q1` – Generates the AR(1) process, computes empirical statistics (mean/second moment), finds β, runs optimal filters for **L=1.5** and calculates NRdB.
- `%% 2A: calculation of the R matrix and its eignables values` – Builds the autocorrelation matrix 𝑅, computes eigenvalues, and analyzes λmax for step-size selection.
- `%%Q3` - LMS/NLMS experiments with different filter lengths (L) and step sizes (μ); produces plots and NRdB values.
- `%% Q4.1 **and** %% 4.2` - RLS experiments scanning δ and λ; plots weight error convergence and computes NRdB.
- `%% Question 5, %% 5(b) + (c) + (d), %% Tests for Q5` - Audio experiments using `airplane.wav`, `cafe.wav`, `city.wav`, `vacuumcleaner.wav`; instantaneous power plots with window size M = 10000; comparison to a trivial filter L = 1, w = [1]; parameter sweeps for λ, μ, and δ.
- `%% Helper Functions` – Utility functions for generating processes, computing R and p, etc. (not run directly).


### 2) Running `adaptivepredict.m`
- This function receives a vector of samples z_1,...,z_n (e.g. an audio segment) and z\hat_(n-1) returns the predicted next sample.
- It internally tests multiple parameter sets for RLS and LMS and selects the one that achieves the best NRdB.
   
   Example:

   ```
   Z = audioread('airplane.wav');   % or any other signal
   Z = Z(:,1);                      % if stereo, use one channel
   znext = adaptivepredict(Z);
   ```


   Inside the function:
   - **RLS** is tested with `L_values_RLS = [10 20 30], lambda_vals = [0.01, 0.5, 0.999], delta = 1`.
   - **LMS** is tested with `L_values_LMS = [5 10 15], mu_vals = [0.85 9 0.01]`.
   
   Note: some μ values (e.g., μ=9) are aggressive and may be unstable; they were included for experimentation.


## 📊 Main Experiments
- **Steepest Descent:** Step-size effect and weight error convergence.
- **LMS/NLMS:** Convergence speed and misadjustment for different parameters.
- **RLS:** Effect of forgetting factor (λ) and initialization (δ).
- **Audio Denoising:** Instantaneous power plots, NRdB calculations, and side-by-side comparison of noisy vs denoised signals.

## 📈 Results
- Figures saved to results/figures/
- Denoised audio saved to results/audio/
   

