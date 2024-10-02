import pandas as pd
import numpy as np
import pyspark.pandas as ps
from pyspark.sql import SparkSession, Window, SQLContext
from pyspark.sql import functions as F

from pyspark import SparkConf, SparkContext

from pyspark.sql.functions import udf, struct, when, rand, exp, pow, lit
from pyspark.sql.types import IntegerType, NumericType

import os
os.environ["PYARROW_IGNORE_TIMEZONE"] = "1"

import time
#

conf = SparkConf().setAppName("AppGen")
sc=SparkContext(conf=conf)
sc.setLogLevel("WARN")
sqlContext = SQLContext(sc)
#


def Write_LL_Event(iname,oname,paras,seed=3742621):
    tstart = time.time()
    psdf0 = sqlContext.read.parquet(iname) #reads the file
    
    
    #Splits the factor covariates into their factor values
    psdf0 = psdf0.withColumn("sex1", when(psdf0.sex==1,1).otherwise(0))

    psdf0 = psdf0.withColumn("Urban1", when(psdf0.Urban==1,1).otherwise(0))

    psdf0 = psdf0.withColumn("Work_Group1", when(psdf0.Work_Group==1,1).otherwise(0))
    psdf0 = psdf0.withColumn("Work_Group2", when(psdf0.Work_Group==2,1).otherwise(0))
    psdf0 = psdf0.withColumn("Work_Group3", when(psdf0.Work_Group==3,1).otherwise(0))
    psdf0 = psdf0.withColumn("Work_Group4", when(psdf0.Work_Group==4,1).otherwise(0))

    psdf0 = psdf0.withColumn("Smoking1", when(psdf0.Smoking==1,1).otherwise(0))
    print("-------------------------------------------------------------------------1")
    print(time.time()-tstart)
    
    #we are applying a "cell-killing" curve
    psdf0=psdf0.withColumn("back_rate", (psdf0.age_exit-psdf0.age_entry)*exp(-0.035*psdf0.age_entry)*pow(psdf0.age_entry,2)*(1/24298))


    cov_cols = ['sex1', 'Urban1', 'Work_Group1', 'Work_Group2', 'Work_Group3', 'Work_Group4', 'Smoking1', 'dose']

    psdf0=psdf0.withColumn("cov_rate", lit(0))
    for i in range(len(cov_cols)): #combines to a weighted sum of covariates
        psdf0=psdf0.withColumn("cov_rate", psdf0.cov_rate + psdf0[cov_cols[i]]*paras[i])
    psdf0=psdf0.withColumn("cov_rate", exp(psdf0.cov_rate)) #takes exponential of sum

    print("-------------------------------------------------------------------------2")
    print(time.time()-tstart)

    psdf0=psdf0.withColumn("p_event", psdf0.cov_rate * psdf0.back_rate) #total event probability
    psdf0=psdf0.withColumn("unirand", rand(seed)) #reference uniform random value
    ##
    psdf0=psdf0.withColumn("event", when(psdf0.p_event>psdf0.unirand,1).otherwise(0)) #flag for rows with event
    psdf0=psdf0.withColumn("exit_adjust", when(psdf0.event>0,psdf0.age_entry+psdf0.unirand/psdf0.p_event).otherwise(psdf0.age_exit)) #shifts exit ages so higher probabilities give a shorter interval
    ##
    psdf0=psdf0.withColumn("event_t",when(psdf0.event==1,psdf0.exit_adjust).otherwise(100)) #gives event ages
    ##
    w = Window.partitionBy('UserID')
    psdf0 = psdf0.withColumn('minEt', F.min('event_t').over(w)) #first event age for every person, or a age above the upper limit
    psdf0=psdf0.withColumn("event", when(psdf0.minEt<psdf0.exit_adjust,2).otherwise(psdf0.event)) #for every individual, flags ages past event
    ##
    print("-------------------------------------------------------------------------3")
    print(time.time()-tstart)

    psdf0.select('seq_id','exit_adjust','event').write.parquet(oname,'overwrite') #minimum amount of columns to line up results and give new ages/events
    print("-------------------------------------------------------------------------4")
    print(time.time()-tstart)
    return



for i_f in ["5050","2575","7525"]: #the three base files used
    for i in [0,1,2]: #iterating over different levels
        if i==0:
            paras=[-0.69,0.1,0.69,0.1,-0.1,-0.69,0.6,-2/30.0]
        elif i==1:
            paras=[-0.1,0.1,0.69,0.1,-0.1,-0.69,0.6,0]
        elif i==2:
            paras=[0.0,0.1,0.69,0.1,-0.1,-0.69,0.6,2/30.0]
        iname="{}.parquet".format(i_f)
        oname="{}_LL_{}.parquet".format(i_f,i)
        print(oname)
        if os.path.exists(oname):
            print("exists")
        else:
            pass
            Write_LL_Event(iname,oname,paras)


