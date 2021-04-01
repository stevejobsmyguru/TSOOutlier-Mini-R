#!/usr/bin/env python3

import pandas as pd
from adtk.data import validate_series
from adtk.visualization import plot
from adtk.detector import LevelShiftAD

s_train = pd.read_csv("./Kohl_BB_Data.csv", index_col="time", parse_dates=True, squeeze=True)
s_train = validate_series(s_train)
#print(s_train)

plot(s_train)

level_shift_ad = LevelShiftAD(c=6.0, side='both', window=1) # This is almost matching to TSOutlier
anomalies_1 = level_shift_ad.fit_detect(s_train)


s_test_output = pd.concat([s_train,anomalies_1],axis=1)
print(s_test_output)
plot(s_train, anomaly=anomalies_1, anomaly_color='red');
