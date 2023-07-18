import pandas as pd
import numpy as np
import pyspark.pandas as ps
from pyspark.sql import SparkSession, Window
from pyspark.sql import functions as F

from pyspark import SparkConf, SparkContext

from pyspark.sql.functions import udf, struct, when, rand, exp, pow, lit, monotonically_increasing_id
from pyspark.sql.types import IntegerType, NumericType

import os
os.environ["PYARROW_IGNORE_TIMEZONE"] = "1"

import time
#
if (False):
    spark = SparkSession.builder.config("spark.driver.memory", "20g").config("spark.driver.maxResultSize", 0).config("spark.executor.memory","20g").getOrCreate()
    spark.conf.set("spark.sql.execution.arrow.pyspark.enabled", True)
    sc=spark.sparkContext
else:
    conf = SparkConf().setAppName("AppBase")
    sc=SparkContext(conf=conf)
sc.setLogLevel("WARN")
#

n_group= int(1*10**6)

if (True):
    tstart = time.time()
    psdf0 = ps.DataFrame(
        {'UserID':[i+1 for i in range(n_group)],
         'age_entry': [18]*n_group,
         'age_exit': [18.5]*n_group,
         'sex': np.random.randint(0,2,n_group),
         'Urban': np.random.randint(0,2,n_group),
         'Work_Group': np.random.randint(0,5,n_group),
         'Smoking': np.random.randint(0,2,n_group)},
        index=[i+1 for i in range(n_group)]).to_spark()


    time0=32
    for i in range(7):
        psdf0 = psdf0.union(psdf0.withColumn("age_entry", psdf0.age_entry + time0*(0.5)**i).withColumn("age_exit", psdf0.age_exit + time0*(0.5)**i))
        
    psdf0 = psdf0.withColumn("annual_dose", when(psdf0.age_entry<65, 1).otherwise(0))
    psdf0 = psdf0.withColumn("annual_dose", rand(seed=3742) * psdf0.annual_dose)

    windowval = (Window.partitionBy('UserID').orderBy('age_entry')
                 .rowsBetween(Window.unboundedPreceding, 0))
    psdf0 = psdf0.withColumn('dose', F.sum('annual_dose').over(windowval))
    psdf0 = psdf0.withColumn("seq_id", monotonically_increasing_id())
    print("-------------------------------------------------------------------------3")
    print(time.time()-tstart)

    psdf0.select('UserID','seq_id','age_entry','age_exit','sex','Urban','Work_Group','Smoking','dose').write.parquet('5050.parquet','overwrite')
    print("-------------------------------------------------------------------------4")
    print(time.time()-tstart)

if (True):
    tstart = time.time()
    psdf0 = ps.DataFrame(
        {'UserID':[i+1 for i in range(n_group)],
         'age_entry': [18]*n_group,
         'age_exit': [18.5]*n_group,
         'sex': np.random.randint(0,4,n_group),
         'Urban': np.random.randint(0,2,n_group),
         'Work_Group': np.random.randint(0,5,n_group),
         'Smoking': np.random.randint(0,2,n_group)},
        index=[i+1 for i in range(n_group)]).to_spark()
    psdf0 = psdf0.withColumn("sex", when(psdf0.sex==0,0).otherwise(1))


    time0=32
    for i in range(7):
        psdf0 = psdf0.union(psdf0.withColumn("age_entry", psdf0.age_entry + time0*(0.5)**i).withColumn("age_exit", psdf0.age_exit + time0*(0.5)**i))
        
    psdf0 = psdf0.withColumn("annual_dose", when(psdf0.age_entry<65, 1).otherwise(0))
    psdf0 = psdf0.withColumn("annual_dose", rand(seed=3742) * psdf0.annual_dose)

    windowval = (Window.partitionBy('UserID').orderBy('age_entry')
                 .rowsBetween(Window.unboundedPreceding, 0))
    psdf0 = psdf0.withColumn('dose', F.sum('annual_dose').over(windowval))
    psdf0 = psdf0.withColumn("seq_id", monotonically_increasing_id())
    print("-------------------------------------------------------------------------3")
    print(time.time()-tstart)

    psdf0.select('UserID','seq_id','age_entry','age_exit','sex','Urban','Work_Group','Smoking','dose').write.parquet('2575.parquet','overwrite')
    print("-------------------------------------------------------------------------4")
    print(time.time()-tstart)

if (True):
    tstart = time.time()
    psdf0 = ps.DataFrame(
        {'UserID':[i+1 for i in range(n_group)],
         'age_entry': [18]*n_group,
         'age_exit': [18.5]*n_group,
         'sex': np.random.randint(0,4,n_group),
         'Urban': np.random.randint(0,2,n_group),
         'Work_Group': np.random.randint(0,5,n_group),
         'Smoking': np.random.randint(0,2,n_group)},
        index=[i+1 for i in range(n_group)]).to_spark()
    psdf0 = psdf0.withColumn("sex", when(psdf0.sex==0,1).otherwise(0))



    time0=32
    for i in range(7):
        psdf0 = psdf0.union(psdf0.withColumn("age_entry", psdf0.age_entry + time0*(0.5)**i).withColumn("age_exit", psdf0.age_exit + time0*(0.5)**i))
        
    psdf0 = psdf0.withColumn("annual_dose", when(psdf0.age_entry<65, 1).otherwise(0))
    psdf0 = psdf0.withColumn("annual_dose", rand(seed=3742) * psdf0.annual_dose)

    windowval = (Window.partitionBy('UserID').orderBy('age_entry')
                 .rowsBetween(Window.unboundedPreceding, 0))
    psdf0 = psdf0.withColumn('dose', F.sum('annual_dose').over(windowval))
    psdf0 = psdf0.withColumn("seq_id", monotonically_increasing_id())
    print("-------------------------------------------------------------------------3")
    print(time.time()-tstart)

    psdf0.select('UserID','seq_id','age_entry','age_exit','sex','Urban','Work_Group','Smoking','dose').write.parquet('7525.parquet','overwrite')
    print("-------------------------------------------------------------------------4")
    print(time.time()-tstart)

