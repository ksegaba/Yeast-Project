### 08/22/2023
The baseline RF models for SNPs, ORFs, CNVs must not be sklearn version 0.24.2,
unlike the FS models. At some point I remember upgrading the sklearn version
for some reason and ran the feature selection models with the upgraded version.

_The error I got:_
Traceback (most recent call last):
  File "SHAP_training.py", line 285, in <module>
    main()
  File "SHAP_training.py", line 159, in main
    my_model = joblib.load(args.model)
  File "/mnt/home/seguraab/.local/lib/python3.6/site-packages/joblib/numpy_pickle.py", line 587, in load
    obj = _unpickle(fobj, filename, mmap_mode)
  File "/mnt/home/seguraab/.local/lib/python3.6/site-packages/joblib/numpy_pickle.py", line 506, in _unpickle
    obj = unpickler.load()
  File "/opt/software/Python/3.6.4-foss-2018a/lib/python3.6/pickle.py", line 1050, in load
    dispatch[key[0]](self)
  File "/opt/software/Python/3.6.4-foss-2018a/lib/python3.6/pickle.py", line 1338, in load_global
    klass = self.find_class(module, name)
  File "/opt/software/Python/3.6.4-foss-2018a/lib/python3.6/pickle.py", line 1388, in find_class
    __import__(module, level=0)
ModuleNotFoundError: No module named 'sklearn.ensemble.forest'

_Solution:_
I created a conda environment called "shap" with scikit-learn v. 0.18.1 and
scipy v. 0.18.1 and numpy v. 1.19.2
pip install -U numpy-1.19.2 scikit-learn-0.18.1 scipy-0.18.1
Collecting scikit-learn==0.18.1
  Using cached scikit_learn-0.18.1-cp36-cp36m-manylinux1_x86_64.whl (11.8 MB)
Collecting scipy==0.18.1
  Using cached scipy-0.18.1-cp36-cp36m-manylinux1_x86_64.whl (42.5 MB)
Collecting numpy==1.11.3
  Downloading numpy-1.11.3-cp36-cp36m-manylinux1_x86_64.whl (15.7 MB)
     |████████████████████████████████| 15.7 MB 69.6 MB/s
Installing collected packages: scipy, scikit-learn, numpy
  Attempting uninstall: scipy
    Found existing installation: scipy 1.5.2
    Uninstalling scipy-1.5.2:
      Successfully uninstalled scipy-1.5.2
  Attempting uninstall: scikit-learn
    Found existing installation: scikit-learn 0.21.3
    Uninstalling scikit-learn-0.21.3:
      Successfully uninstalled scikit-learn-0.21.3
  Attempting uninstall: numpy
    Found existing installation: numpy 1.19.5
    Uninstalling numpy-1.19.5:
      Successfully uninstalled numpy-1.19.5
Successfully installed numpy-1.19.2 scikit-learn-0.18.1 scipy-0.18.1

Also installed pandas v. 0.18.1 to address this error:
  ImportError: this version of pandas is incompatible with numpy < 1.13.3
pip install -U pandas==0.18.1

### 08/02/2023
I ran YPDNACL15 trait with 1024 features through SHAP_training.py using the command:
scratch=/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
trait=YPDNACL15M
snp_feat=1024
python SHAP_training.py -df ${data}/geno.csv -sep , -df2 ${data}/pheno.csv \
    -y_name YPDNACL15M -test ${data}/Test.txt \
    -feat ${scratch}/feat*_rf_${trait}_top_${snp_feat} \
    -model ${scratch}/${trait}*_rf_${snp_feat}_models.pkl \
    -top 10 -save SNP_Figures/${trait}_snp_rf_${snp_feat}

_Package versions:_
datatable                                 1.0.0
matplotlib                                3.2.1
numpy                                     1.19.5
pandas                                    1.0.4
shap                                      0.40.0
scikit-learn                              0.21.3

_The error I got:_
Traceback (most recent call last):
  File "SHAP_training.py", line 285, in <module>
    main()
  File "SHAP_training.py", line 159, in main
    my_model = joblib.load(args.model)
  File "/mnt/home/seguraab/.local/lib/python3.6/site-packages/joblib/numpy_pickle.py", line 587, in load
    obj = _unpickle(fobj, filename, mmap_mode)
  File "/mnt/home/seguraab/.local/lib/python3.6/site-packages/joblib/numpy_pickle.py", line 506, in _unpickle
    obj = unpickler.load()
  File "/opt/software/Python/3.6.4-foss-2018a/lib/python3.6/pickle.py", line 1050, in load
    dispatch[key[0]](self)
  File "/opt/software/Python/3.6.4-foss-2018a/lib/python3.6/pickle.py", line 1338, in load_global
    klass = self.find_class(module, name)
  File "/opt/software/Python/3.6.4-foss-2018a/lib/python3.6/pickle.py", line 1388, in find_class
    __import__(module, level=0)
ModuleNotFoundError: No module named 'sklearn.ensemble._forest'

_Solution:_
My scikit-learn version is too old, so I installed scikit-learn 0.24.2 since that is the version of the RF model.