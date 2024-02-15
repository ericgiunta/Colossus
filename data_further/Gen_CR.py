                                ###################################################
                                ################ IMPORTS AND SETUP ################
                                ###################################################
import numpy as np
from typing import Optional
from pyspark.sql.types import *
from numpy.lib import math as M
from numpy import random as Random
from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession, Window, SQLContext
from pyspark.sql import functions as F
import os
import random
os.environ["PYARROW_IGNORE_TIMEZONE"] = "1"
from pyspark.sql.functions import udf, struct, when, exp, pow, lit, least, round as sqlround, sum as sqlsum
# round and sum are default functions, so the sqlsum and sqlround names are needed to differentiate
conf = SparkConf().setAppName("AppGen")
sc=SparkContext(conf=conf)
sc.setLogLevel("WARN")
sqlContext = SQLContext(sc)

                                ##################################################
                                ############# Event Func Definitions #############
                                ##################################################
def Time_Dependency(t: float, c: int, beta: float):
  return ((t**2) * M.exp((c + (beta/2) )* t))

class EventObj:
  def __init__(self, alpha=float,
               dose_Beta=float,
               smoke_Betas=[float, float, float],
               econ_Betas=[float, float, float, float, float, float, float],
               sex_Betas=[float, float], c=float):
    self.alpha = alpha
    self.dose_Beta = dose_Beta
    self.smoke_Beta = smoke_Betas
    self.econ_Beta = econ_Betas
    self.sex_Beta = sex_Betas
    self.c = c

def Event(t0, dose, smoke, econ, sex, event_params=EventObj):
  t1 = t0+1
  u = np.log(random.uniform(0, 1))
  cov = (event_params.c + (event_params.dose_Beta/2) ) # cov represents c'
  # The Following Array index calls get appropriate beta value.
  A = event_params.alpha * M.exp((event_params.dose_Beta * dose)+( event_params.smoke_Beta[smoke - 1])+(event_params.econ_Beta[econ - 1])+(event_params.sex_Beta[sex - 1]))
  F0 = ((cov*t0 * (cov*t0 - 2) + 2) * M.exp(t0 * cov)) / M.pow(cov, 3)
  F1 = ((cov*t1 * (cov*t1 - 2) + 2) * M.exp(t1 * cov)) / M.pow(cov, 3)
  x = (F1 - F0)-(-u / A)
  tl = t0
  th = t1
  Fl = ((cov*tl * (cov*tl - 2) + 2) * M.exp(tl * cov)) / M.pow(cov, 3)
  Fh = ((cov*th * (cov*th - 2) + 2) * M.exp(th * cov)) / M.pow(cov, 3)
  t1 = (tl+th)/2
  F1 = ((cov*t1 * (cov*t1 - 2) + 2) * M.exp(t1 * cov)) / M.pow(cov, 3)
  x = (F1 - F0)-(-u / A)
  count=0
  while (abs(x)> 0.01):
    count += 1
    Fl = ((cov*tl * (cov*tl - 2) + 2) * M.exp(tl * cov)) / M.pow(cov, 3)
    Fh = ((cov*th * (cov*th - 2) + 2) * M.exp(th * cov)) / M.pow(cov, 3)
    if (Fh-F0)<-u/A: #Upper estimate is too low
      tl = th #shifts the interval up
      th = tl + (t1-t0)
      t1 = (tl+th)/2
      F1 = ((cov*t1 * (cov*t1 - 2) + 2) * M.exp(t1 * cov)) / M.pow(cov, 3)
      x = (F1 - F0)-(-u / A)
    elif ((Fl-F0)<-u/A):
      #assert , "Lower estimate is also too high"
      #The lower bound is below and the upper bound is above, so we assume the zero lies between
      if (F1-F0)<-u/A: #The midpoint is also below
        tl=t1 #we move the lower bound up
        t1 = (tl+th)/2
        F1 = ((cov*t1 * (cov*t1 - 2) + 2) * M.exp(t1 * cov)) / M.pow(cov, 3)
        x = (F1 - F0)-(-u / A)
      else: #The midpoint is also above
        th=t1 #we move the upper bound down
        t1 = (tl+th)/2
        F1 = ((cov*t1 * (cov*t1 - 2) + 2) * M.exp(t1 * cov)) / M.pow(cov, 3)
        x = (F1 - F0)-(-u / A)
    else: return t1
  return t1

                                ##################################################
                                ############# DATASET GEN DEFINITION #############
                                ##################################################
def Test_Dataset_Values(set_size, e1_params=EventObj, e2_params=EventObj, e3_params=EventObj, output=None, decimal_precision=None, sort=None, ascending=None):
  """Creates a dataset of random individuals which are then censored via the hazard function definitions
  - set_size:   The number of individuals to represent
  - e1_params:   EventObj type representing the Beta values for Covariates in Event 1
  - e2_params:   EventObj type representing the Beta values for Covariates in Event 2
  - e3_params:   EventObj type representing the Beta values for Covariates in Event 3
  - output:   File path to write .parquet file to. "p" or "print" will print the dataset to console instead.
  - decimal_Precision:   The number of decimal places to round final dataset to. Leave empty for no rounding.
  - sort:   The column to order by
  - ascending:   True for ascending order, False for descending order
  """
  ##### Decimal Precision has a default of 12, thus returning the same values when Decimal_Precision is null
  if decimal_precision is None: decimal_precision = 12 # Set the decimal place precision for the output
  data = []

  # Loop through every individual which will need to be generated
  for i in range(0,  set_size):
    entrance_age = random.randint(18, 25)
    dose = random.random()
    age = entrance_age+1
    sex = random.randint(0, 1)
    smoker_status = random.randint(1, 3)
    econ_status = random.randint(1, 7)
    Included = 1
    event_1 = Event(dose=dose, smoke=smoker_status, econ=econ_status, sex=sex, t0=entrance_age,  \
                                            event_params=e1_params)
    event_2 = Event(dose=dose, smoke=smoker_status, econ=econ_status, sex=sex, t0=entrance_age,  \
                                            event_params=e2_params)
    event_3 = Event(dose=dose, smoke=smoker_status, econ=econ_status, sex=sex, t0=entrance_age,  \
                                            event_params=e3_params)
    min_event = min(event_3, event_2, event_1, 80.0)
    event = 0
    if min_event < age:
      Included = 0
      if (event_1<age):
        event=1
      elif event_2<age:
        event=2
      elif event_3<age:
        event=3
    data.append({"-ID": i,
                  "Included": Included,
                  "Smoker_Status": smoker_status,
                  "Entrance_Age": round(entrance_age, decimal_precision),
                  "Age": round(age, decimal_precision),
                  "Dosage": round(dose, decimal_precision),
                  "Sex": sex,
                  "Economic_Status": econ_status,
                  "Event1": round(event_1, decimal_precision), "Event2": round(event_2, decimal_precision), "Event3": round(event_3, decimal_precision),
                  "Earliest_Event_Age": round(min_event, decimal_precision),
                  "Event": event})
    #for each person loop through until their exit age occurs in the coming year
    while age < 80 and Included == 1:
      entrance_age = age
      age += 1
      dose += random.uniform(0, 1)
      event_1 = min(Event(dose=dose, smoke=smoker_status, econ=econ_status, sex=sex, t0=entrance_age,  \
                                              event_params=e1_params), 80.0)
      event_2 = min(Event(dose=dose, smoke=smoker_status, econ=econ_status, sex=sex, t0=entrance_age,  \
                                              event_params=e2_params), 80.0)
      event_3 = min(Event(dose=dose, smoke=smoker_status, econ=econ_status, sex=sex, t0=entrance_age,  \
                                              event_params=e3_params), 80.0)
      min_event = min(event_3, event_2, event_1, 80.0)
      if min_event < age:
        Included = 0
        if (event_1<age):
          event=1
        elif event_2<age:
          event=2
        elif event_3<age:
          event=3
      data.append({"-ID": i,
                    "Included": Included,
                    "Smoker_Status": smoker_status,
                    "Entrance_Age": round(entrance_age, decimal_precision),
                    "Age": round(age, decimal_precision),
                    "Dosage": round(dose, decimal_precision),
                    "Sex": sex,
                    "Economic_Status": econ_status,
                    "Event1": round(event_1, decimal_precision),
                    "Event2": round(event_2, decimal_precision),
                    "Event3": round(event_3, decimal_precision),
                    "Earliest_Event_Age": round(min_event, decimal_precision),
                    "Event": event})
  RDD = sc.parallelize(data)
  EventRDD = RDD.map(lambda x:
                        (x["-ID"],
                        x["Included"],
                        x["Smoker_Status"],
                        x["Entrance_Age"],
                        x["Age"],
                        x["Dosage"],
                        x["Sex"],
                        x["Economic_Status"],
                        x["Event1"],
                        x["Event2"],
                        x["Event3"],
                        x["Earliest_Event_Age"],
                        x["Event"]))
  df = EventRDD.toDF(["-ID", "Included", "Smoker_Status", "Entrance_Age", "Age", "Dosage", "Sex", "Economic_Status", "Event1", "Event2", "Event3", "Earliest_Event_Age","Event"])
  df = df.withColumn("Dosage", sqlround(df["Dosage"], decimal_precision))
  df = df.withColumn("Entrance_Age", sqlround(df["Entrance_Age"], decimal_precision))

  if sort != None:
    if ascending == None: ascending = False
    df = df.orderBy(sort, ascending=ascending)
  else:
    if ascending == None:
      ascending = False
      df = df.orderBy('-ID', ascending=True)
    if ascending == True and sort==None:
      df = df.orderBy('-ID', ascending=False)
  if output is None:
    print()
    return df
  elif output == "print" or output == 'p':
    df.show()
  else:
    df.write.parquet(output, 'overwrite')

                                ##########################################################
                                ################ RUNNING EVENT VARIATIONS ################
                                ##########################################################
                                
E_Late_High  = EventObj(alpha=(1/7500),  dose_Beta=(0.2/30),  smoke_Betas=([np.log(0.75),np.log(1),np.log(1.25)]), econ_Betas=(np.log(np.linspace(0.75,1.25,7))), sex_Betas=([0.0,np.log(0.8)]), c=-0.035)
E_Late_Low   = EventObj(alpha=(1/30000), dose_Beta=(0.15/30), smoke_Betas=([np.log(0.75),np.log(1),np.log(1.25)]), econ_Betas=(np.log(np.linspace(0.75,1.25,7))), sex_Betas=([0.0,np.log(0.8)]), c=-0.035)
E_Early_Low  = EventObj(alpha=(1/7500),  dose_Beta=(0.3/30),  smoke_Betas=([np.log(0.75),np.log(1),np.log(1.25)]), econ_Betas=(np.log(np.linspace(0.75,1.25,7))), sex_Betas=([0.0,np.log(0.8)]), c=-0.065)
E_Early_High = EventObj(alpha=(1/3000),  dose_Beta=(0.3/30),  smoke_Betas=([np.log(0.75),np.log(1),np.log(1.25)]), econ_Betas=(np.log(np.linspace(0.75,1.25,7))), sex_Betas=([0.0,np.log(0.8)]), c=-0.065)

E_Low  = EventObj(alpha=(1/40000), dose_Beta=(0.0), smoke_Betas=([0.0, 0.0, 0.0]), econ_Betas=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), sex_Betas=(0.0, 0.0), c=-0.035)
E_High = EventObj(alpha=(1/20000), dose_Beta=(0.0), smoke_Betas=([0.0, 0.0, 0.0]), econ_Betas=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), sex_Betas=(0.0, 0.0), c=-0.035)

E3 = E_High
k3 = "E3=E_High"
for k2, E2, in {"E2=E_LATE_HIGH": E_Late_High, "E2=E_LATE_LOW": E_Late_Low, "E2=E_EARLY_LOW": E_Early_Low, "E2=E_Early_High": E_Early_High}.items():
  for k1, E1, in {"E1=E_LATE_HIGH": E_Late_High, "E1=E_LATE_LOW": E_Late_Low, "E1=E_EARLY_LOW": E_Early_Low, "E1=E_Early_High": E_Early_High}.items():
    DF = Test_Dataset_Values(set_size=1000000, e1_params = E1, e2_params=E2, e3_params=E3, sort='-ID', ascending=False, decimal_precision=4, output="/homes/dtstutzman/MPS_Research/outputs/%s_%s_%s.parquet"%(k1, k2, k3))
    print("Completed: %s_%s_%s"%(k1, k2, k3))
