## MSINGB



### MSINGB：a novel computational method based on NGBoost for identifying microsatellite instability status from tumour mutation annotation data.

In this work, we comprehensively study the feature selection strategy, classification model, performance evaluation metrics and webserver/software usability,  proposing a learning framework termed MSINGB, which achieves a better prediction performance on both the 10-fold cross-validation test and independent test compared with state-of-the-art methods.

***



### Requirements

The following python packages are required when building the stacking model.

Python 3 >= 3.7

Mlxtend >= 0.17.3

 Sklearn >= 0.0 

NGBoost >= 0.3.6

intervaltree >= 3.1.0

pandas >= 1.1.2

***



### Running MSINGB
#### Step 1
download the file "simplerepeats. TXT" from the UCSC database and save it in the data folder. 
</br>***Download link: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz***
#### Step 2
open cmd in Windows or terminal in Linux, then cd to the MSINGB/codes folder which contains pre.py.
</br>`python msingb.py --maf [predicting data in maf format] -outdir [ predicting results in the folder]`</br>  </br>***Example:***
</br>`python msingb.py --maf ../data/example.maf --outdir ../out`</br>  

***



### Announcements

* In order to obtain the prediction results, please be sure to enter the file in MAF format. If you don't know the format of the MAF file, please download the example.maf file to check.
* The intermediate feature file and result file generated in the process of running return prediction are saved in the output folder required by the user.
* The code "getFeature. Py" and "getTag. Py" refer to the feature extraction in MSIpred (Wang and Liang, 2018), respectively.



### References
Wang, C. and Liang, C. MSIpred: a python package for tumor microsatellite instability classification from tumor mutation annotation data using a support vector machine. *Scientific Reports* 2018;8.
