# Survival Analysis and Nomogram by R

This is a population based retrospective survival analysis on cancer patients from 1987 to 2018. The data involved patients' age at diagnosis, sex, date of diagnosis, status (alive or dead), date of death(if any), primary tumour site, tumour differentiation, and cause of death (Null, cancer-related, and non-cancer-related).

In this project, the impact of each factor on survival was analyzed by univariable and multivariable Cox proportional hazard test (survival package). Nomograms to predict one year, three year, five year and ten year survival were plotted using rms package. Internal validation and calibration were also performed.

machine learning methods, including random survival forest, and recursive partitioning analysis were also performed using ranger and rpart, respectively.
