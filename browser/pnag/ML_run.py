# run ml model and pre. load from pickle file

import pickle
import sys

import numpy as np
import pandas as pd
import os
# import sklearn random forest regressor
from sklearn.ensemble import RandomForestRegressor


# load the model from disk
filename = './pnag/static/rf_optimized_model_mason.sav'
rf = pickle.load(open(filename, 'rb'))

# load the pre. data given in the python command
out_path = sys.argv[1]
#out_path = "./static/data/2024_11_18_16_12_24/b0185/outputs"

# load the data
data = pd.read_csv(out_path + "/saved_table_ml.csv")
# exclude MIC_pred column
#data = data.drop('MIC_pred', axis=1)

# run model on data
y_pred = rf.predict(data)

data['MIC_pred'] = np.exp2(y_pred)

# save the data
data.to_csv(out_path + "/saved_table_ml.csv", index=False)






