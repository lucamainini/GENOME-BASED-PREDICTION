# GENOME-BASED PREDICTION OF BREAST CANCER CELL RESPONSE
> Given the genomic features of a cell line, can we predict if it is sensitive or resistant to a particular treatment? <br>
  Is it possible to identify a minimal set of features (genes) that predict the response?

<p align="center">
    <img src="./media/pcna.png" height="350" alt="pcna_nate"/>
    <img src="./media/fbl.png" height="350" alt="reg1a_luca"/>
</p>

<a name="description"/>

## Description
Here the results of a statistical analysis of breast cancer data are presented. This project was carried out during the `Applied Statistics` course at the `Politecnico di Milano`. The relevant datasets used for this project can be found through the [Cancer Cell Line Encyclopedia (CCLE)](https://depmap.org/portal/download/) database. The CCLE is a multi-institutional effort to develop a comprehensive database for varying types of cancer. **Our main objective was to infer which genomic aspects most influence the sensitivity to antitumoral drugs.**

<a name="motivation"/>

## Motivation and Methodology
The goal of this project is to understand the link between the genetic and physiological data of a patient and how efficiently cancer-treating drugs perform on their cells. 

<p align="center">
    <img src="./media/scope.png" height="350" alt="scope"/>
</p>

- First, patients were divided into groups based on the effectiveness of the treatments on their cells. 
- Then, we tried to explain these groups with their gene expression.

To simplify the project scope, the analysis is strictly limited to patients having breast cancer, as opposed to all cancers which were present in the main dataset. 

<a name="why_breast"/>

## Why Breast Cancer?
<p align="center">
    <img src="./media/why_breast.svg" height="350" alt="why breast?"/>
</p>

<a name="assumptions"/>

## Assumptions and terminology 
In the following sections, the terms "AUC" or "AUC score" are equivalent to "drug efficacy". AUC actually refers to the Area Under the dose-response Curve. These value should therefore not be understood as precise indicator of the effectiveness of treatments on patients, but as the result of in-vitro studies.

<a name="requirements"/>

## Requirements 
The project itself uses the following files from the CCLE database:
- *data_clinical_patient.txt*
- *data_clinical_sample.txt*
- *data_drug_treatment_auc.txt*
- *data_mrna_seq_rpkm.txt*

For seamless integration with the existing code, download and save these files in a directory denoted `Dataset` stored in the root. 

On the R side, loading the project in R Studio should prompt to load any missing packages you have.  

*Data preprocess:* treatments with insufficient data were removed, and finally missing values were filled.
You can find the obtained datasets in the `Dataset` folder.

<a name="results"/>

## Results

### Clustering cell lines by treatment response

The aim of this step is to divide the patients according to the effectiveness of the treatments. For example, we could try to obtain a group where all treatments are effective and one where treatments generally do not work. Remember we define the drug performance as in-vitro efficacy (measured by AUC scores). 

Using the implemented Shiny app, you can view the clusters obtained with a hierarchical approach or with a non-hierarchical one.

![APP](media/app_1.jpeg)
![APP](media/app_2.png)

To run this Shiny app, you can clone the git repository and then use `runApp()`:

```R
# First, clone the repository with git. 
# If you have cloned it into ~/clusters
setwd("~/clusters")
runApp()
```
Alternatively, you can use `runGitHub`.

```R
runUrl("https://github.com/lucamainini/GENOME-BASED-PREDICTION/archive/master.zip",
       subdir = "clusters/")
```

Based on the results obtained, we decided to use k-means for clustering. The groups obtained reflect the average effectiveness of the treatments. In particular, one group clearly contains patients with higher AUC.

### Influencial genes through Random Forest Classification
*Dimensionality reduction*
The main effort of this project was to predict the auc score on a given cell using the genetic and physiological information of the patient.  Before such analysis could be undertaken, however, it was necessary to reduce the dimensionality of the feature space since over 50000 gene expressions were included in the provided dataset. The first thing we did was to disregard those genes whose variability was too low between patients.

Then, using the previously obtained groups, we fitted a Random Forest Classification Model. The most influential gene was FBL which is known to be an independent marker of poor outcome in breast cancer[^fbl]. Some of the key genes are shown in next Figure.

<p align="center">
    <img src="./media/clustering.svg" height="550" alt="step1"/>
</p>

[^fbl]:Nguyen Van Long, F., Lardy-Cleaud, A., Car??ne, D. et al. Low level of Fibrillarin, a ribosome biogenesis factor, is a new independent marker of poor outcome in breast cancer. BMC Cancer 22, 526 (2022). https://doi.org/10.1186/s12885-022-09552-x

As can be seen in the graph below, the expression of these genes is particularly different in the various groups.

<p align="center">
    <img src="./media/expression_group.svg" height="350" alt="expression_of selected genes in groups"/>
</p>

### LDA Classification
Using a small number of selected genes and the previously determined groups, we trained an LDA classification algorithm. Using the app contained in `clusters` folder, you can interactively visualise the results. *We emphasise that LDA's assumptions are only partially verified.*

<p align="center">
    <img src="./media/LDA.png" height="450" alt="lda classification"/>
</p>

### Influencial genes through LASSO Regression
The second way in which we carried out features (genes) selection was via a LASSO regression of the average efficiency over these drugs on a cell using the genetic expression and extracted coefficients which are above a certain threshold.??

<p align="center">
    <img src="./media/lasso.svg" height="550" alt="influential genes by LASSO"/>
</p>

This gene pool was cross-referenced with the outcome of previous random forests, and the intersection of the chosen genes from each method was kept.  After this, exhaustive search was applied to the reduced gene pool to find the best selection of regressors for models with *p=1...n*-dimensional feature spaces. The identified genes were researched to see whether or not they have been previously identified as being influential in breast cancer expression.  In the end, the feature space was reduced to 8 significant genes (our regressors). 

<p align="center">
    <img src="./media/model_summary.png" height="450" alt="pcna_nate"/>
</p>

To account for the limited number of samples present in the cleaned dataset (removal of NA's and non-relevant cancer-treating drugs) k-fold cross validation was used. The (k=5) cross-validated R^2 values during the training phase for each fold were: 0.746, 0.770, 0.758, 0.745, and 0.786.

## Conclusion and Limits of our project
There is a high degree of diversity between and within breast tumours as well as amongst cancer-bearing individuals. Understanding the molecular features responsible for this heterogeneity is important to correctly diagnose and to obtain more efficient treatment strategies. To this end, our project statistically confirms results previously obtained from bio-chemical analyses.

In particular, it shows a relationship between the expression of a select subset of pivotal genes and the in-vitro treatment efficacy (measured with AUC scores). 
One should be cognizant that in order to confirm the validity of these results, the models here presented should be tested on a larger and more varied dataset. 
In addition, it is important to point out that the AUC score refers to an in-vitro analysis. The values predicted by our linear model should therefore be understood as a rough indicator for treatment efficacy. Having data available on actual patient therapies would allow for the real-life robustness of the model to be ascertained. Finally, since the variability among the AUC scores of the chosen treatments is relatively low, it would be prudent to repeat the analysis of the clusters on a larger, higher variance dataset. New relationships and understandings may emerge.

## Acknowledgement
*This project was realised together with `Jara--Mikolajczak Jimmy Aygalic, Montani Giulia, Nethercott Robert James Nate and Troilo Roberta`.<br>
In this git folder you will find the app (which I made myself) for displaying the obtained clusters and classification results. For any further information on our work, please contact me by e-mail. <br>
Last but not least, I would like to thank our supervisors `Dr. Pietro Pinoli` & `Prof. Francesca Ieva` for all their help and advice with this project.*
