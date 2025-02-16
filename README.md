R code for the model estimation under "Contractual" and "Non-contractual" settings in the paper: 
"An Economic Model of Membership Subscription and Repeat Purchase Based on Continuous Information Tracking"

QUICK START:
run_demo() or run_demo_membership()
for a demonstration using CDNOW dataset of the estimation routine used by the authors for the model under "Contractual" and "Non-contractual"
settings, respectively.

CDNOW_data_augmentation() was used to synthesize the membership version of CDNOW data for this demonstration.

GENERAL GUIDE:
Use Estimate_Parameters_PDE(...) or Estimate_Parameters(...) to estimate the model's parameters given the training dataset.
Once the estimation parameters are obtained, use Make_Prediction(...) to evaluate the model on the test dataset.

The dataset consists of a list of (x,t1,t2,...,tx,T) for the "Non-contractual" setting, or (x,t1,t2,...,tx,Tm,T) for the "Contractual" setting, for each individual customers.
NOTE: In the current version, Estimate_Parameters(...) only supports the model under the "Non-contractual" setting.
