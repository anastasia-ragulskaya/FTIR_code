# FTIR_Code

This GUI was written for the data analysis of the ["FTIR on Protein Secondary Structure" practical of the Institute of Applied Physics, Tübingen University](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/physik/institute/institut-fuer-angewandte-physik/studium/nanoscience-practical/). This practicum is a part of the Physics, Biomedical Technology, and Nanoscience studying programs.

## Disclaimer

This repository mainly serves as a public archive for the code used for the [FTIR practical](https://uni-tuebingen.de/securedl/sdl-eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpYXQiOjE2NzAyMTU0MTgsImV4cCI6MTY3MDMwNTQxNiwidXNlciI6MCwiZ3JvdXBzIjpbMCwtMV0sImZpbGUiOiJmaWxlYWRtaW5cL1VuaV9UdWViaW5nZW5cL0Zha3VsdGFldGVuXC9NYXRoZVBoeXNpa1wvSW5zdGl0dXRlXC9JQVBcL0ZvcnNjaHVuZ1wvWmhhbmdcL0ZUSVJfZW5nbGlzaF8yMDE4LnBkZiIsInBhZ2UiOjQzMzQ1fQ.TNpcm8kYLYcPkTrrntls6nanO_J0anspyI65YZYRqRs/FTIR_english_2018.pdf) and is not yet optimized for the use by other researchers. Future updates will improve the functionality and usability of the program.

## Main dependencies

To be able to run the code, a python 3.7 installation with the following dependencies is required:
* PySimpleGUI
* matplotlib
* numpy
* scipy.optimize

The github directory also contains "main.exe", which can be used to directly start the application on Windows.

## Data analysis
The raw data of the FTIR consists of two columns (wavenumber and absorption) and is saved as a .txt file. The data procedure is the following:
### Load the data
To load the data press the Browse button and choose the folder. Then the list of all .txt files in the folder will appear below. Once the file from the list is chosen, the corresponding spectrum appears on the right.
![Program window layout. The browse of the folder and the list of files in it are marked with yellow.](/step1.png)
### Extract Amide I band
The range [1600, 1700] cm<sup>-1</sup> of the spectrum corresponds to the Amide I band. The "Extract Amide I band" button display this part of the spectrum on the right. The obtained spectrum should be preprocessed before the following data analysis. A baseline correction can be done by pressing "Background Subtraction" button and is preformed by subtracting a linear background from the absorption spectrum between 1600 cm<sup>-1</sup> and 1700 cm<sup>-1</sup>. This approximately corresponds to a subtraction of the side-chains absorption. The corrected spectrum can be saved by pressing "Save data after baseline correction" button.
![](/step2.png)

### Fit

"Fit button" performs the fit of the absorption spectra using a multiple Gaussian fit (sum of five Gaussian functions). The starting values for central wavelengths are fixed and correspond to the band position of the  Amide I modes of different secondary structure elements for a protein dissolved in D<sub>2</sub>O. The heights, Ai, and half widths (FWHM), si, are free within the boundary conditions determined by the user (see slides in "Boundaries for fit parameters" section of the program). If the fit converges, the results appear in the table in the right side of the window.
![](/step3.png)

### Post analysis
The results obtained by this program can be used to determine the amount of secondary structures (β–sheets, not ordered, α–helix, and slopes) for the measured protein: 
* Integrate the obtained single Gaussian functions for a quantitative analysis of the secondary structure. 
* Determine the fraction of secondary structure components by comparing the obtained areas. This is done on the assumption that the extinction by the different secondary structures is approximately constant in the region of the main band. The Gaussian function which corresponds to the high frequency lateral bands of the β–sheets between 1682 cm<sup>-1</sup> and 1689 cm<sup>-1</sup> is neglected for the calculation of the secondary structure fractions. Otherwise the β–sheet contribution would be counted twice. Thus, only four out of the five fitted Gaussians are used to determine the specific fractions. 

# Authors
Dr. Anastasia Ragulskaya (Institut für Angewandte Physik, University of Tübingen)
