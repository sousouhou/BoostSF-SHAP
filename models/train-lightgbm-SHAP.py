import os 
import shutil
import math

import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt
from scipy import stats

import lightgbm


Xtrainfile = "../data/pdbbind-2020refined-exclude-2016core-x41features.csv"
ytrainfile = "../data/pdbbind-2020refined-exclude-2016core-y.txt"
pdbnametrain = "../data/pdbbind-2020refined-exclude-2016core-pdbnamelist.txt"

generateModelFolder  = "lightgbm-model01"  # the folder to store result

flagConductTest = True  # whether using test data to evaluate the model

if flagConductTest == True :
    Xtestfile = "../data/pdbbind-2016core-x41features.csv"
    ytestfile = "../data/pdbbind-2016core-y.txt"
    pdbnametest = "../data/pdbbind-2016core-pdbnamelist.txt"   


if os.path.exists(generateModelFolder) == True : 
    shutil.rmtree(generateModelFolder)     # ensure folder does not exist.
os.mkdir(generateModelFolder)


Xtrain = pd.read_csv( Xtrainfile )
ytrain = np.loadtxt( ytrainfile )


# ------- initialize Regressor -------
LGBM_params = {      # parameters for lightgbm
    'n_estimators' : 1000 ,
    'max_depth' : 8 ,
    'learning_rate' : 0.02 ,
    'random_state' : 0
}


model = lightgbm.LGBMRegressor(** LGBM_params)
# fit model 
model.fit(Xtrain, ytrain)


# save model
model.booster_.save_model("%s/trained-model.txt"%generateModelFolder)


# ------- model information -------
strMI = """\
* This model is trained on the following data:
{ %s }
{ %s }

The number of samples is: %d.
The number of features is: %d.

* The training data are extracted from the following PDBs, which are listed in:
{ %s } \n\n
"""%(Xtrainfile, ytrainfile,  Xtrain.shape[0],  Xtrain.shape[1],  pdbnametrain)

pf = open("%s/model-information.txt"%generateModelFolder, "w")
pf.write(strMI)
pf.close()


# ------- SHAP figure -------
explainer = shap.TreeExplainer(model)
shap_values = explainer(Xtrain)


shap.plots.beeswarm(shap_values,  max_display=60, show=False) # not show
fig01 = plt.gcf()
ax01 = plt.gca()
ax01.set_title( "beeswarm figure for the training data\n\n Note: Feature C/C represents C(in receptor)/C(in ligand)." )
ax01.set_ylabel("Features", fontsize=14)
fig01.savefig("%s/SHAP-beeswarm.png"%generateModelFolder, dpi=600, bbox_inches='tight')
plt.clf()  # clear figure


shap.summary_plot(shap_values, Xtrain, max_display=60, plot_type="bar", show=False)
fig02 = plt.gcf()
ax02 = plt.gca()
ax02.set_title( "global importance of each feature for the training data\n\n Note: Feature C/C represents C(in receptor)/C(in ligand)." )
ax02.set_ylabel("Features", fontsize=14)
fig02.savefig("%s/SHAP-summary_plot.png"%generateModelFolder, dpi=600, bbox_inches='tight')



# ------- figure Pearson's correlation coefficient, SD  for training -------
ytrain_predicted = model.predict(Xtrain)

fig03, ax03 = plt.subplots()

ax03.scatter(ytrain, ytrain_predicted, s=8, linewidths=1, facecolors='none', edgecolors='C0')

ax03.set_xlabel("Experimental binding affinity")
ax03.set_ylabel("Predicted binding affinity")

Rp,pvalue = stats.pearsonr(ytrain, ytrain_predicted)

MSE = 0.0
N = ytrain.shape[0]
for i in range(0, N):
    MSE += (ytrain[i] - ytrain_predicted[i])**2
MSE = MSE*1.0/N
RMSE = math.sqrt(MSE)
SD = math.sqrt(N*1.0/(N-1))*RMSE

strtitle = "Correlation plot on the training data\n\nThe Pearson's correlation coefficient $R_p$ = %.3f\nThe standard deviation $SD$ = %.3f"%(Rp, SD) 
ax03.set_title(strtitle)

fig03.savefig("%s/correlation-plot-training.png"%generateModelFolder, dpi=300, bbox_inches='tight')


if flagConductTest == True :
    Xtest = pd.read_csv( Xtestfile )
    ytest = np.loadtxt( ytestfile )
    ytest_predicted = model.predict(Xtest)

    fig04, ax04 = plt.subplots()

    ax04.scatter(ytest, ytest_predicted, s=8, linewidths=1, facecolors='none', edgecolors='C0')

    ax04.set_xlabel("Experimental binding affinity")
    ax04.set_ylabel("Predicted binding affinity")

    Rp,pvalue = stats.pearsonr(ytest, ytest_predicted)

    MSE = 0.0
    N = ytest.shape[0]
    for i in range(0, N):
        MSE += (ytest[i] - ytest_predicted[i])**2
    MSE = MSE*1.0/N
    RMSE = math.sqrt(MSE)
    SD = math.sqrt(N*1.0/(N-1))*RMSE

    strtitle = "Correlation plot on the test data\n\nThe Pearson's correlation coefficient $R_p$ = %.3f\nThe standard deviation $SD$ = %.3f"%(Rp, SD) 
    ax04.set_title(strtitle)

    fig04.savefig("%s/correlation-plot-test.png"%generateModelFolder, dpi=300, bbox_inches='tight')

# ------- model information -------
    strMI = """--------------------\n
* This model is tested on the following data:
{ %s }
{ %s }

The number of samples is: %d.
The number of features is: %d.

* The test data are extracted from the following PDBs, which are listed in:
{ %s }
    """%(Xtestfile, ytestfile,  Xtest.shape[0],  Xtest.shape[1],  pdbnametest)

    pf = open("%s/model-information.txt"%generateModelFolder, "a+")
    pf.write(strMI)
    pf.close()
    