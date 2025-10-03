import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random as rnd
import math as mth
import time as t
import pyspark
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.window import Window
from pyspark.sql.types import *

# Initialize a new Spark session
spark = SparkSession.builder \
    .appName("RiskIndFinder") \
    .getOrCreate()

# Import csv to SparkSession.DataFrame
test_0_0_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_0_0.csv")
test_0_1_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_0_1.csv")
test_1_0_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_1_0.csv")
test_1_1_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_1_1.csv")
test_2_0_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_2_0.csv")
test_2_1_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_2_1.csv")
test_3_0_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_3_0.csv")
test_3_1_spark = spark.read.format("csv").option("header", "true").load("/content/drive/MyDrive/Colab Notebooks/Conversion_Testing_Data/test_3_1.csv")


# Find all of the distinct ending times in the dataframe
time_arr = test_0_0_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
# Create a broadcast variable for the time array, for synced data
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
#This needs to be done before each runthrough


# Function to filter DataFrame based on time array
def filter_between(df, time_arr_broadcast):
    time_arr = time_arr_broadcast.value
    return df.filter(
        (F.col('t0').isin(time_arr)) |
        (F.col('t1').isin(time_arr)) |
        ((F.col('t0') < F.array_max(F.array(*[F.lit(t) for t in time_arr]))) &
         (F.col('t1') >= F.array_min(F.array(*[F.lit(t) for t in time_arr]))))
    )

def filter_end(df, time_arr_broadcast):
  time_arr = time_arr_broadcast.value
  return df.filter(
      (F.col('t1').isin(time_arr))
  )





# Apply the filter function to each DataFrame
# Find all of the distinct ending times in the dataframe
time_arr = test_0_0_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
# Create a broadcast variable for the time array, for synced data
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
#This needs to be done before each runthrough

t_0_0 = t.time()
mult_val_0_0 = [filter_between(test_0_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_0_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_0_0 = t.time() - t_0_0



time_arr = test_0_1_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_0_1 = t.time()
mult_val_0_1 = [filter_between(test_0_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_0_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_0_1 = t.time() - t_0_1



time_arr = test_1_0_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_1_0 = t.time()
mult_val_1_0 = [filter_between(test_1_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_1_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_1_0 = t.time() - t_1_0



time_arr = test_1_1_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_1_1 = t.time()
mult_val_1_1 = [filter_between(test_1_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_1_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_1_1 = t.time() - t_1_1



time_arr = test_2_0_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_2_0 = t.time()
mult_val_2_0 = [filter_between(test_2_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_2_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_2_0 = t.time() - t_2_0



time_arr = test_2_1_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_2_1 = t.time()
mult_val_2_1 = [filter_between(test_2_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_2_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_2_1 = t.time() - t_2_1


time_arr = test_3_0_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_3_0 = t.time()
mult_val_3_0 = [filter_between(test_3_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_3_0_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_3_0 = t.time() - t_3_0
print("mult_val_3_0[0]: ", len(mult_val_3_0[0]), "\nmult_val_3_0[1]: ", len(mult_val_3_0[1]))



time_arr = test_3_1_spark.filter(F.col('event') == 1).select("t1").distinct().rdd.flatMap(lambda x: x).collect()
time_arr_broadcast = spark.sparkContext.broadcast(time_arr)
t_3_1 = t.time()
mult_val_3_1 = [filter_between(test_3_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect(),
                filter_end(test_3_1_spark, time_arr_broadcast).select("_c0").rdd.flatMap(lambda x: x).collect()]
t_3_1 = t.time() - t_3_1
print("mult_val_3_1[0]: ", len(mult_val_3_1[0]), "\nmult_val_3_1[1]: ", len(mult_val_3_1[1]))



#Plotting the results:

# Plotting the results
dataset = ("test_0", "test_1", "test_2", "test_3")
times = {
    'a_b_0': (round(t_0_0, 2), round(t_1_0, 2), round(t_2_0, 2), round(t_3_0, 2)),
    'a_b_1': (round(t_0_1, 2), round(t_1_1, 2), round(t_2_1, 2), round(t_3_1, 2))
}

x = np.arange(len(dataset))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

# Plot each set of measurements
for attribute, measurement in times.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=2)
    multiplier += 1

# Add labels, title, and custom x-axis tick labels
ax.set_ylabel('Execution Time (Seconds)')
ax.set_title('Lookup of both lists')
ax.set_xticks(x + (width / 2), dataset)
ax.legend(loc='upper left', ncols=2)
ax.set_ylim(0, max(t_3_0, t_3_1) * 1.2)

# Save and show the plot
plt.savefig("example_test.png")
plt.show()