# Modified_UNIFAC
The module can be used to predict activity coefficients using the Modified UNIFAC model. Uses Modified Raoult's law for binary VLE data generation and multicomponent VLE flash calculations 

A number of jupyter notebooks have been included that illustrate the use of these modules. "moduni.py" is the main module and us the heart of this system and the backend data required for it's working are "MODUNI.xlsx" and "Data.xlsx". Thus, these three are the essential files for all the other files. "flash_solver.py" is built on top of "moduni.py" and uses it to perform multicomponent VLE flash calculations using the Rachford-Rice algorithm. 

### What all can you do with these pieces of code 

1) Use "moduni.py" to find activity coefficients
2) Generate binary VLE data (isothermal and isobaric) using "moduni.py" for select compounds available in "Data.xlsx"
3) Perform multicomponent PT-flash calculation using "flash_solver.py" 

### Have a look at the jupyter notebooks to understand how to use the functions 

## Important points 
1) The code can be used to compute activity coefficients for any compound that you can construct using the subgroups and is not restricted to the select compounds in "Data.xlsx"
2) However, binary VLE data generation and flash calculations can be performed only if the compound is available in the list. This is because the program relies on vapor pressure data and the data for only a select compounds is stored in "Data.xlsx"
3) If you can add the vapor pressure data for your custom compounds to "Data.xlsx", VLE calculations are possible. However while doing so note that the code is written with a static reference and if you add a compound, you must update the range of excel cells within which the code searches for the compound name
4) The code makes use of Modified Raoult's law and gives good results only at lower and close to atmospheric pressures. High pressure VLE would require the introduction of fugacities which are currently not made available
5) The activity coefficients have been tested and are found to give good results for VLE calculations. The module has not been extensively tested on LLE data and is not recommended for LLE flash calculations. Also, the activity coefficients close to infinite dilution are not very accurate. Thus it is recommended to use this module only for low pressure VLE calculations as of now
6) The flash solver is not designed to handle non-condensables and supercritical components as of now
