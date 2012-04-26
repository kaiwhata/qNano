Python code used for interacting with data generated by qNano.

-----
extractor.py - extracts .irbd files into .dat files [Current(nA) with time(s)]

execute with the following linux command: 
$ od -A d -t f8 -w32 Raw_data.irbd |./extractor.py > Raw_data.dat

od translates the binary into a form python can interface with, then pipes the output to the .dat file, 
which can later be read line-by-line

-----
locator.py - locates and returns a list of all 'blockade_1.csv' type files in a directory, 
copies and renames them with their parent directory's name and sticks them in a file called 'data' 
in the parent directory

------
xml_interpreter.py - extracts parameters from the .xml file accompanying each qNano recording.

options are:
  <SrAperture>A6006</SrAperture>
  <SrDilution>1</SrDilution>
  <SrElectrolyteID>SEB</SrElectrolyteID>
  <SrExperiment>Elf_charge_paper_review</SrExperiment>
  <SrFileName>c:\Izon Data\\Elf_charge_paper_review\B86_pH6_0.5V_30.1.12_b\B86_pH6_0.5V_30.1.12_b.irbd</SrFileName>
  <SrNotes />
  <SrPartNumber />
  <SrPressure>0</SrPressure>
  <SrRawConcentration>0</SrRawConcentration>
  <SrSampleName>B86_pH6_0.5V_30.1.12_b</SrSampleName>
  <SrSize>0</SrSize>
  <SrType>Sample</SrType>
  <SrZetaPotential>0</SrZetaPotential>
  <SrRecordTime>2012-01-30T15:32:12.5739991+13:00</SrRecordTime>
  <SrBandwidthFilter>0</SrBandwidthFilter>
  <SrBandwidthFilterOn>false</SrBandwidthFilterOn>

----
baseline_comparison.py - extracts baselines of all plots and overlays on the same graph. 
N.B. baselines vary wildly - but normalization (\frac{\Delta I}{I} of events seems to account for most effects.

----
bibliography_plot.py
Uses bibitex functionality to extract information including year and author from .bib library and plot grouping by shared authors.

-----
plotting_functions.py - performs a standard Cumulcative Count experimental plot including:
-cubic and parabolic fitting and inflection point printing to external file
-mean FWHM plotting and gaussian fit to extract mean and st-dev values
-baseline current plots overlaid with mean normalised blockades
-comparison plot between mean, media and mode FWHM fitting

N.B. this only works with blockade_1.csv files generated with the izon command line software, that have then been calibrated to convert a time into an applied pressure.
