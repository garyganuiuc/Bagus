##########################################################
####################Bagus Simulation Code#################
##########################################################

Author: Lingrui Gan

##########################################################
####################Codes#########################
##########################################################
1.10fold-glasso-cv.R : 10-fold cross validation for Graphical Lasso
2.space-cv.R: 10-fold cross validation for SPACE method
3.FPandFN.R: Calculate the false positive/false negative/true positive/true negative rate for the selection
4.BIC_Bagus.R: Calculating the BIC criteria for Bagus
5.Bagus-algo.R: The main algorithm for Bagus method
6.tune-Bagus.R: Tuning for Bags method.
7.Bagus_simulation.R: The simulations for all the methods considered.


##########################################################
####################Instruction to run the code###########
##########################################################

1. Open Bagus_simulaiton.R
2. Preload other functions in the folder with source.
i.e Run all the source(###.R) line in the Bagus_simulation.R
2. Choose the destined simulation setting by identifying p and n and structure of C.
3. Run the simulation code section(50 iterations) to get the results.
Results will include methods: Glasso, SPACE, Bagus and CLIME.

##########################################################
####################Output################################
##########################################################
MCC*: Matthew’s Correlation Coefficient Criteria
Sens*: Sensitivity
Spec*: Specificity
Fnorm*: Difference between truth and the estimation in Frobenius Norm

Corresponding to each methods, the outputs are:
Glasso: MCC1,Sens1,Spec1,Fnorm1
Meinshausen and Buhlmann’s neighborhood selection method: MCC2,Sens2,Spec2,Fnorm2
SPACE: MCC3,Sens3,Spec3,Fnorm3
Bagus: MCC4,Sens4,Spec4,Fnorm4
CLIME: MCC5, Sens5,Spec5,Fnorm5
