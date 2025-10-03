#Given Data Generation

import numpy as np
import pandas as pd

num_row = int(100)
np.random.seed(37421016)

dose =  np.random.rand(num_row)
econ = np.random.randint(0, 4, size = num_row)
YOB = np.random.randint(0, 4, size = num_row)

df = pd.DataFrame({'dose':dose, 'econ':econ,'YOB':YOB})

df['econ0'] = (df.econ==0).astype(int)
df['econ1'] = (df.econ==1).astype(int)
df['econ2'] = (df.econ==2).astype(int)
df['econ3'] = (df.econ==3).astype(int)

df['YOB0'] = (df.YOB==0).astype(int)
df['YOB1'] = (df.YOB==1).astype(int)
df['YOB2'] = (df.YOB==2).astype(int)
df['YOB3'] = (df.YOB==3).astype(int)


df.to_csv("risk_table.csv")

#Given Value Calculation
import numpy as np
import pandas as pd

def df_refresh():
    df_T = pd.DataFrame({'ProductLinear':[0]*len(df.dose),'LogLinear':[0]*len(df.dose)})
    df_dT = pd.DataFrame({'0':[0]*len(df.dose), '1':[0]*len(df.dose), '2':[0]*len(df.dose), '3':[0]*len(df.dose), '4':[0]*len(df.dose), '5':[0]*len(df.dose), '6':[0]*len(df.dose)})
    df_ddT = pd.DataFrame({'00':[0]*len(df.dose), '01':[0]*len(df.dose), '02':[0]*len(df.dose), '03':[0]*len(df.dose), '04':[0]*len(df.dose), '05':[0]*len(df.dose), '06':[0]*len(df.dose),
                           '11':[0]*len(df.dose), '12':[0]*len(df.dose), '13':[0]*len(df.dose), '14':[0]*len(df.dose), '15':[0]*len(df.dose), '16':[0]*len(df.dose), '22':[0]*len(df.dose),
                           '23':[0]*len(df.dose), '24':[0]*len(df.dose), '25':[0]*len(df.dose), '26':[0]*len(df.dose), '33':[0]*len(df.dose), '34':[0]*len(df.dose), '35':[0]*len(df.dose),
                           '36':[0]*len(df.dose), '44':[0]*len(df.dose), '45':[0]*len(df.dose), '46':[0]*len(df.dose), '55':[0]*len(df.dose), '56':[0]*len(df.dose), '66':[0]*len(df.dose)})
    df_R = pd.DataFrame({'total':[0]*len(df.dose)})
    df_dR = pd.DataFrame({'0':[0]*len(df.dose), '1':[0]*len(df.dose), '2':[0]*len(df.dose), '3':[0]*len(df.dose), '4':[0]*len(df.dose), '5':[0]*len(df.dose), '6':[0]*len(df.dose)})
    df_ddR = pd.DataFrame({'00':[0]*len(df.dose), '01':[0]*len(df.dose), '02':[0]*len(df.dose), '03':[0]*len(df.dose), '04':[0]*len(df.dose), '05':[0]*len(df.dose), '06':[0]*len(df.dose),
                           '11':[0]*len(df.dose), '12':[0]*len(df.dose), '13':[0]*len(df.dose), '14':[0]*len(df.dose), '15':[0]*len(df.dose), '16':[0]*len(df.dose), '22':[0]*len(df.dose),
                           '23':[0]*len(df.dose), '24':[0]*len(df.dose), '25':[0]*len(df.dose), '26':[0]*len(df.dose), '33':[0]*len(df.dose), '34':[0]*len(df.dose), '35':[0]*len(df.dose),
                           '36':[0]*len(df.dose), '44':[0]*len(df.dose), '45':[0]*len(df.dose), '46':[0]*len(df.dose), '55':[0]*len(df.dose), '56':[0]*len(df.dose), '66':[0]*len(df.dose)})
    return df_T, df_dT, df_ddT, df_R, df_dR, df_ddR

def df_summary(df_T, df_dT, df_ddT, df_R, df_dR, df_ddR):
    accuracy = 4
    print("Subterm")
    print([round(sum(df_T[i]),accuracy)   for i in df_T.columns])
    print("Subterm first derivative")
    print([round(sum(df_dT[i]),accuracy)  for i in df_dT.columns])
    print("Subterm second derivative")
    print([round(sum(df_ddT[i]),accuracy) for i in df_ddT.columns])
    print("Risk")
    print([round(sum(df_R[i]),accuracy)   for i in df_R.columns])
    print("Risk first derivative")
    print([round(sum(df_dR[i]),accuracy)  for i in df_dR.columns])
    print("Risk second derivative")
    print([round(sum(df_ddR[i]),accuracy) for i in df_ddR.columns])
    return

df = pd.read_csv("risk_table.csv", index_col=0)

beta = [-0.1, 0.2, 0.3, 0.4, 0.2, 0.3, 0.4]
names = ['dose', 'econ1', 'econ2', 'econ3', 'YOB1', 'YOB2', 'YOB3']

df_T, df_dT, df_ddT, df_R, df_dR, df_ddR = df_refresh()

#Case 1, all loglinear
print("------------------------------------------------------------")
print("LogLinear Values")
for i in range(len(names)):
    a = beta[i]
    n = names[i]
    df_T.LogLinear = df_T.LogLinear + df[n]*a
df_T.LogLinear = np.exp(df_T.LogLinear)
#print(df_T)

for i in range(len(names)):
    n = names[i]
    df_dT[str(i)] = df_T.LogLinear * df[n]
#print(df_dT)

for i in df_ddT.columns:
    j = int(i[0])
    k = int(i[1])
    n0 = names[j]
    n1 = names[k]
    df_ddT[i] = df_dT[str(j)] * df[n1]
#print(df_ddT)

#
df_R.total = df_T.LogLinear
#print(df_R)

for i in range(len(names)):
    df_dR[str(i)] = df_dT[str(i)]
#print(df_dR)

for i in df_ddT.columns:
    df_ddR[i] = df_ddT[i]
#print(df_ddR)

df_summary(df_T, df_dT, df_ddT, df_R, df_dR, df_ddR)

df_T, df_dT, df_ddT, df_R, df_dR, df_ddR = df_refresh()

#Case 2, all product-linear
print("------------------------------------------------------------")
print("ProductLinear Values")
for i in range(len(names)):
    a = beta[i]
    n = names[i]
    df_T.ProductLinear = df_T.ProductLinear + df[n]*a
df_T.ProductLinear = df_T.ProductLinear + 1
#print(df_T)

for i in range(len(names)):
    n = names[i]
    df_dT[str(i)] = df[n]
#print(df_dT)

for i in df_ddT.columns:
    j = int(i[0])
    k = int(i[1])
    n0 = names[j]
    n1 = names[k]
    df_ddT[i] = 0
#print(df_ddT)

#
df_R.total = df_T.ProductLinear
#print(df_R)

for i in range(len(names)):
    df_dR[str(i)] = df_dT[str(i)]
#print(df_dR)

for i in df_ddT.columns:
    df_ddR[i] = df_ddT[i]
#print(df_ddR)

df_summary(df_T, df_dT, df_ddT, df_R, df_dR, df_ddR)

df_T, df_dT, df_ddT, df_R, df_dR, df_ddR = df_refresh()

tform = ['pl', 'll', 'll', 'll', 'll', 'll', 'll']
#Case 3, dose product-linear, rest loglinear
print("------------------------------------------------------------")
print("Combined Values")
for i in range(len(names)):
    a = beta[i]
    n = names[i]
    t = tform[i]
    if t=='ll':
        df_T.LogLinear = df_T.LogLinear + df[n]*a
    else:
        df_T.ProductLinear = df_T.ProductLinear + df[n]*a
df_T.ProductLinear = df_T.ProductLinear + 1
df_T.LogLinear = np.exp(df_T.LogLinear)
#print(df_T)

for i in range(len(names)):
    n = names[i]
    t = tform[i]
    if t=='ll':
        df_dT[str(i)] = df_T.LogLinear * df[n]
    else:
        df_dT[str(i)] = df[n]
#print(df_dT)

for i in df_ddT.columns:
    j = int(i[0])
    k = int(i[1])
    n0 = names[j]
    n1 = names[k]
    t0 = tform[j]
    t1 = tform[k]
    if (t0==t1) and (t0=='ll'):
        df_ddT[i] = df_dT[str(j)] * df[n1]
#print(df_ddT)

#
df_R.total = df_T.LogLinear * df_T.ProductLinear
#print(df_R)

for i in range(len(names)):
    t = tform[i]
    if t=='ll':
        df_dR[str(i)] = df_dT[str(i)] * df_T.ProductLinear
    else:
        df_dR[str(i)] = df_dT[str(i)] * df_T.LogLinear
#print(df_dR)

for i in df_ddT.columns:
    j = int(i[0])
    k = int(i[1])
    n0 = names[j]
    n1 = names[k]
    t0 = tform[j]
    t1 = tform[k]
    if (t0==t1):
        df_ddR[i] = df_ddT[i]
    else:
        df_ddR[i] = df_dT[str(j)] * df_dT[str(k)]
#print(df_ddR)

df_summary(df_T, df_dT, df_ddT, df_R, df_dR, df_ddR)


# Computing Risks

from logging import error
import numpy as np
import pandas as pd
import math

def make_risks(column_name, parameter_value, subterm, term_number, df, term_type):
  if not(len(column_name) == len(parameter_value) == len(term_type)):
    error("Number of columns and Parameters must be the same")
    return []  # Ensure the function exits if the check fails


  def compute_risk(row): # Should be updated with proper equation for R
    PL_term = 1
    LL_term = 1
    for i in range(0, len(parameter_value)):
      col = column_name[i]
      param = parameter_value[i]
      if term_type[i] == "PL":
        PL_term += (row[col] * param)
      elif term_type[i] == "LL":
        LL_term *= math.exp(row[col] * param)
    return PL_term * LL_term

  R_vals = df.apply(compute_risk, axis=1).tolist()
  return R_vals



#Note Subterm logic has not yet been setup, nor term number. Previously the pr

#Example Use Case:
# make_risks(column_name=["econ", "YOB0", "dose"], parameter_value=[0.1, 0.2, 0.4], subterm="LogLinear", term_number=4, df=df, term_type = ["LL", "PL", "PL"])
#print(len(R_vals))
#print(R_vals)