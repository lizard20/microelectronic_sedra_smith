# microelectronic_sedra_smith
A library written in Python for Microelectronic circuits book, chapter 3, by Sedra-Smith. 7th. edition


## Getting started
Python version: 3.11.5

## Installing
1.- To use the library install the library `astropy`
```bash
pip install astropy 
```
2.- Download the repository
```bash
git clone https://github.com/lizard20/microelectronic_sedra_smith.git
```
## How to use it
- Put the `semiconductor.py` module in your working sub directory
- Example:

Exercise 3.1 

Calculate the intrinsic carrier density n<sub>i</sub> for silicon at T = 50 K and 350 K.
```python
'''
    exercise_3_1.py

    Author: Aldo Nunez Tovar
'''

import semiconductor as sc
from astropy import units as u 

T_1 = 50 * u.Kelvin
T_2 = 350 * u.Kelvin

n_i_1 = sc.n_i(T_1)
n_i_2 = sc.n_i(T_2)

print(f"Exercise 3.1:")
print(f"n_i({T_1:.0f}) = {n_i_1:0.2E}")
print(f"n_i({T_2:.0f}) = {n_i_2:0.2E}")
```
To execute
```bash
python exercise_3_1.py 
```
Outcome:
```
Exercise 3.1:
n_i(50 K) = 9.25E-39 1 / cm3
n_i(350 K) = 4.13E+11 1 / cm3
```

## Solved exercises and problems from chapter 3
We have included the files:

`chapter_3_exercises.py`

`chapter_3_problems.py`



