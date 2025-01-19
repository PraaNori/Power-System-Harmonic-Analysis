# Power System Harmonic Analysis

Power System Harmonic Analysis is a Python-based web application built using Flask to analyze voltage signals and determine the presence of harmonics and interharmonics. The application supports visualization of data through various plots and enables users to upload CSV files for analysis.

---

## Features
- **CSV File Upload**: Upload voltage signal data in `.csv` format.
- **Harmonic and Interharmonic Analysis**: Determines whether the signal contains harmonics only or includes interharmonics.
- **Visualization**:
  - Downsampled signal plots (`xp` vs `yp`).
  - Rank Difference (RD) vs Rank Index (RDI) plot.
  - Frequency vs Amplitude plot.
  - Frequency spectrum visualization.
- **Algorithmic Analysis**:
  - **FFT Analysis** for harmonic detection.
  - **MEMO-ESPRIT Algorithm** for interharmonic analysis.

---

## Technologies Used
- **Backend**: Flask
- **Visualization**: Matplotlib
- **Data Processing**: NumPy, Pandas
- **Scientific Computing**: SciPy

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/PraaNori/Power-System-Harmonic-Analysis
   cd Power System Harmonic Analysis
   ```

2. Create a virtual environment (optional but recommended):
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```

3. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Run the Flask application:
   ```bash
   python app.py
   ```

5. Access the application in your browser at:
   [http://127.0.0.1:5000](http://127.0.0.1:5000)

---

## How to Use

1. Open the application in your browser.
2. Upload a `.csv` file containing a voltage signal dataset.
   - The file should have a single row of voltage values.
3. Specify the sampling frequency in the provided input field.
4. Submit the file for analysis.

### Outputs:
- If interharmonics are present:
  - **Plots**:
    - `xp` vs `yp` plot.
    - RD vs RDI plot.
    - Frequency vs Amplitude plot.
  - The message: "Interharmonics are present."
- If only harmonics are present:
  - **Plots**:
    - `xp` vs `yp` plot.
    - Frequency Spectrum plot.
  - The message: "Only harmonics are present."

---

## File Structure
- `app.py`: Main Flask application file.
- `templates/`: Contains the HTML templates for the web application.
  - `index.html`: Main page for uploading and displaying results.
- `static/images/`: Directory where plots are saved.
- `uploads/`: Directory for uploaded CSV files.

---

## Key Algorithms

### MEMO-ESPRIT Algorithm
- Constructs a Hankel matrix and computes the autocorrelation matrix.
- Performs eigenvalue decomposition to extract rank differences.
- Detects interharmonics by analyzing RD vs RDI.

### FFT Analysis
- Computes the Fast Fourier Transform of the signal to extract frequency components.
- Visualizes the frequency spectrum for harmonic analysis.

---

## Requirements
- Python 3.7+
- Libraries:
  - Flask
  - NumPy
  - Pandas
  - Matplotlib
  - SciPy

---

## Future Enhancements
- Add support for additional input formats (e.g., Excel).
- Enable real-time signal analysis.
- Provide a downloadable report of the analysis.
  
