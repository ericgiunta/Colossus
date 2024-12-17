import pandas as pd
import numpy as np
import pyspark.pandas as ps
from pyspark.sql import SparkSession, Window, SQLContext
from pyspark.sql import functions as F

from pyspark import SparkConf, SparkContext

from pyspark.sql.functions import udf, struct, when, rand, exp, pow, lit, monotonically_increasing_id, log
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
sqlContext = SQLContext(sc)
#
# spark-submit --master local[*] --driver-memory 20g --executor-memory 2g --verbose Gen_Base.py
n_group= int(1*10**4)

if (True):
    tstart = time.time()
    psdf0 = ps.DataFrame(
        {'UserID':[i+1 for i in range(n_group)],
         'age_entry': [18]*n_group,
         'age_exit': [19]*n_group,
         'sex': np.random.randint(0,2,n_group)},
        index=[i+1 for i in range(n_group)]).to_spark()


    time0=32
    for i in range(6):
        psdf0 = psdf0.union(psdf0.withColumn("age_entry", psdf0.age_entry + time0*(0.5)**i).withColumn("age_exit", psdf0.age_exit + time0*(0.5)**i))
        
    psdf0 = psdf0.withColumn("annual_dose", when(psdf0.age_entry<65, 1/30).otherwise(0))
    psdf0 = psdf0.withColumn("annual_dose", rand(seed=3742) * psdf0.annual_dose)

    windowval = (Window.partitionBy('UserID').orderBy('age_entry')
                 .rowsBetween(Window.unboundedPreceding, 0))
    psdf0 = psdf0.withColumn('dose', F.sum('annual_dose').over(windowval))
    psdf0 = psdf0.withColumn("seq_id", monotonically_increasing_id())
    print("-------------------------------------------------------------------------1")
    print(time.time()-tstart)

#    psdf0.select('UserID','seq_id','age_entry','age_exit','sex','dose').write.parquet('5050.parquet','overwrite')
    print("-------------------------------------------------------------------------2")
    psdf0 = psdf0.withColumn("sex1", when(psdf0.sex==1,1).otherwise(0))
    print("-------------------------------------------------------------------------3")
    print(time.time()-tstart)
    paras=[0.95, 0.1]
    seed=12132024
    psdf0=psdf0.withColumn("back_rate", exp(-0.03*psdf0.age_entry)*pow(psdf0.age_entry,2)*(1/200000))


    cov_cols = ['sex1']#, 'dose']

    psdf0=psdf0.withColumn("cov_rate", lit(0))
    psdf0=psdf0.withColumn("cov_rate", psdf0.cov_rate + psdf0[cov_cols[0]]*paras[1])
    psdf0=psdf0.withColumn("cov_rate", exp(psdf0.cov_rate))
    
    a0=paras[0]
    #
    psdf0=psdf0.withColumn("cov_rate", psdf0.cov_rate*(lit(1)+a0*psdf0.dose))
#    psdf0=psdf0.withColumn("cov_rate", when(psdf0.cov_rate < lit(1e-10),1e-10).otherwise(psdf0.cov_rate))

    print("-------------------------------------------------------------------------4")
    print(time.time()-tstart)

    psdf0=psdf0.withColumn("p_event", psdf0.cov_rate * psdf0.back_rate)
    psdf0=psdf0.withColumn("unirand", rand(seed))
    ##
    psdf0=psdf0.withColumn("exit_adjust", psdf0.age_entry-log(psdf0.unirand)/psdf0.p_event)
    psdf0=psdf0.withColumn("event", when(psdf0.age_exit>psdf0.exit_adjust,1).otherwise(0))
    psdf0=psdf0.withColumn("exit_adjust", when(psdf0.event>0,psdf0.exit_adjust).otherwise(psdf0.age_exit))
    ##
    psdf0=psdf0.withColumn("event_t",when(psdf0.event==1,psdf0.exit_adjust).otherwise(100))
    ##
    w = Window.partitionBy('UserID')
    psdf0 = psdf0.withColumn('minEt', F.min('event_t').over(w))
    psdf0=psdf0.withColumn("event", when(psdf0.minEt<psdf0.exit_adjust,2).otherwise(psdf0.event))
    psdf0=psdf0.withColumn("event", when(psdf0.age_entry>psdf0.exit_adjust,2).otherwise(psdf0.event))
    ##
    print("-------------------------------------------------------------------------5")
    print(time.time()-tstart)

#    psdf0.select('UserID','seq_id','age_entry','age_exit','sex','dose','exit_adjust','event').write.parquet("excess_test.parquet",'overwrite')
#    print("-------------------------------------------------------------------------6")
#    print(time.time()-tstart)
    psdf0.filter(psdf0.event < 2).select('UserID','seq_id','age_entry','exit_adjust','sex','dose','event','p_event').toPandas().to_csv("excess_test.csv", header=True)
    print("-------------------------------------------------------------------------7")
    print(time.time()-tstart)
