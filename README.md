# GENOME-BASED PREDICTION OF BREAST CANCER CELL RESPONSE
> Given the genomic features of a cell line, can we predict if it is sensitive or resistant to a particular treatment? <br>
  Is it possible to identify a minimal set of features (genes) that predict the response?

<p align="center">
    <img src="./media/pcna.png" height="350" alt="pcna_nate"/>
    <img src="./media/reg1a.png" height="350" alt="reg1a_luca"/>
</p>

<a name="description"/>

## Description
Here the results of a statistical analysis of breast cancer data are presented. This project was carried out during the Applied Statistics course at the Politecnico di Milano. The relevant datasets used for this project can be found through the [Cancer Cell Line Encyclopedia (CCLE)](https://depmap.org/portal/download/) database. The CCLE is a multi-institutional effort to develop a comprehensive database for varying types of cancer. **Our main objective was to infer which genomic aspects most influence the sensitivity to antitumoral drugs.**

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

A collection of utilities and preprocessing scripts are included in the utilities folder of the project. 





<a name="results"/>

## Results

### Hierarchical Classification 
<p align="center">
    <img src="./media/hierarchical.png" height="350" alt="pcna_nate"/>
</p>

One vein of analysis was to consider the clusters arising from applying heirarchical clustering (euclidean distance, single linkage) to a reduced dataset containing only two cancers which in theory should be quite different in their manifestation.  We prove the merit of this approach using samples corresponding to breast and central nervous system cancers. 


