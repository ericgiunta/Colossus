This folder contains two types of script for generating person data and generating events.
The Gen_Base.py script produces a parquet file with specific distributions of biological sex, half-year intervals from 18 to 81, for 1 million individuals. Random uniform dose is assigned up to age 65.
Every other script applies a different type of risk model. Ranging from Exponential (EXP, LL0-4) up to more complicated linear and piecewise functions.
