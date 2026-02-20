# Modified_UNIFAC
The module can be used to predict activity coefficients using the Modified UNIFAC model. Uses Modified Raoult's law for binary VLE data generation and multicomponent VLE flash calculations 

A number of jupyter notebooks have been included that illustrate the use of these modules. "moduni.py" is the main module and us the heart of this system and the backend data required for it's working are "MODUNI.xlsx" and "Data.xlsx". Thus, these three are the essential files for all the other files. "flash_solver.py" is built on top of "moduni.py" and uses it to perform multicomponent VLE flash calculations using the Rachford-Rice algorithm. 

### What all can you do with these pieces of code 

1) Use "moduni.py" to find activity coefficients
2) Generate binary VLE data (isothermal and isobaric) using "moduni.py" for select compounds available in "Data.xlsx"
3) Perform multicomponent PT-flash calculation using "flash_solver.py" 
