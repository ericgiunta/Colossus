import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random as rnd
import math as mth
import time as t

#Note: File Paths for data will need to be changed.

test_0_0 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_0_0.csv")
test_0_1 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_0_1.csv")
test_1_0 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_1_0.csv")
test_1_1 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_1_1.csv")
test_2_0 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_2_0.csv")
test_2_1 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_2_1.csv")
test_3_0 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_3_0.csv")
test_3_1 = pd.read_csv("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_3_1.csv")

def get_risk_ind(df, ntime):
  return (df.loc[(df['event'] == 1) & ((df['t1'] >= ntime) & (df['t0'] < ntime)) | (df['t0'] == ntime) & (df['t1'] == ntime)])["Unnamed: 0"].tolist()

def risk_group(df, ntime):
    return (df.loc[((df['t1'] >= ntime) & (df['t0'] < ntime)) | ((df['t0'] == ntime) & (df['t1'] == ntime))])["Unnamed: 0"].tolist()
def event_group(df, ntime):
    return (df.loc[(df['event'] == 1) & (df['t1'] == ntime)])["Unnamed: 0"].tolist()

def get_mult_risk_ind(df1, time_arr):
  res_group = []
  res_event = []
  for i in time_arr:
    temp = risk_group(df1, i)
    res_group += temp
    temp = event_group(df1, i)
    res_event += temp
  return res_group, res_event

df_end = test_0_0[test_0_0.event==1]
time_arr = np.unique(df_end.t1)

t_0_0 = t.time()
mult_val_0_0 = get_mult_risk_ind(test_0_0, time_arr)
t_0_0 = t.time() - t_0_0

df_end = test_0_1[test_0_1.event==1]
time_arr = np.unique(df_end.t1)

t_0_1 = t.time()
mult_val_0_1 = get_mult_risk_ind(test_0_1, time_arr)
t_0_1 = t.time() - t_0_1

df_end = test_1_0[test_1_0.event==1]
time_arr = np.unique(df_end.t1)

t_1_0 = t.time()
mult_val_1_0 = get_mult_risk_ind(test_1_0, time_arr)
t_1_0 = t.time() - t_1_0

df_end = test_1_1[test_1_1.event==1]
time_arr = np.unique(df_end.t1)

t_1_1 = t.time()
mult_val_1_1 = get_mult_risk_ind(test_1_1, time_arr)
t_1_1 = t.time() - t_1_1

df_end = test_2_0[test_2_0.event==1]
time_arr = np.unique(df_end.t1)

t_2_0 = t.time()
mult_val_2_0 = get_mult_risk_ind(test_2_0, time_arr)
t_2_0 = t.time() - t_2_0

df_end = test_2_1[test_2_1.event==1]
time_arr = np.unique(df_end.t1)

t_2_1 = t.time()
mult_val_2_1 = get_mult_risk_ind(test_2_1, time_arr)
t_2_1 = t.time() - t_2_1

df_end = test_3_0[test_3_0.event==1]
time_arr = np.unique(df_end.t1)

t_3_0 = t.time()
mult_val_3_0 = get_mult_risk_ind(test_3_0, time_arr)
t_3_0 = t.time() - t_3_0

df_end = test_3_1[test_3_1.event==1]
time_arr = np.unique(df_end.t1)

t_3_1 = t.time()
mult_val_3_1 = get_mult_risk_ind(test_3_1, time_arr)
t_3_1 = t.time() - t_3_1

#assert False

print("mult_val_3_0[0]: ", len(mult_val_3_0[0]), "\nmult_val_3_0[1]: ", len(mult_val_3_0[1]))
print("mult_val_3_1[0]: ", len(mult_val_3_1[0]), "\nmult_val_3_1[1]: ", len(mult_val_3_1[1]))


#Optional Data Visualization:

'''
dataset = ("test_0", "test_1", "test_2", "test_3")
times = {
    'a_b_0': (round(t_0_0, 2), round(t_1_0, 2), round(t_2_0, 2), round(t_3_0, 2)),
    'a_b_1': (round(t_0_1, 2), round(t_1_1, 2), round(t_2_1, 2), round(t_3_1, 2))
}

x = np.arange(len(dataset))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in times.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=2)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Execution Time (Seconds)')
ax.set_title('Lookup of both lists')
ax.set_xticks(x + (width / 2), dataset)
ax.legend(loc='upper left', ncols=2)
ax.set_ylim(0, max(t_3_0, t_3_1) * 1.2)

plt.savefig("example_test.png")
plt.show()

print("mult_val_3_0[0]: ", len(list(set(mult_val_3_0[0]))), "\nmult_val_3_0[1]: ", len(list(set(mult_val_3_0[1]))))
print("mult_val_3_1[0]: ", len(list(set(mult_val_3_1[0]))), "\nmult_val_3_1[1]: ", len(list(set(mult_val_3_1[1]))))

print("\n\nmult_val_3_0[0]: ", len(mult_val_3_0[0]), "\nmult_val_3_0[1]: ", len(mult_val_3_0[1]))
print("mult_val_3_1[0]: ", len(mult_val_3_1[0]), "\nmult_val_3_1[1]: ", len(mult_val_3_1[1]))

'''
