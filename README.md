# PlenoCalib
Methods for calibrating standard plenoptic cameras.

# Usage
1. Install  the toolbox using:  
```
    ins = setup.Setup;
    ins.install;
```  
It might be needed some manual download of the dependencies identified within the library.

2. Download Lytro datasets and white images from the following link:  
https://ieee-dataport.org/open-access/lytro-illum-calibration-dataset

3. Obtain white images database, decode calibration images and extract corner features using LFToolbox. The toolbox is added as a submodule.

# Citing
The appropriate citations for the calibration and modeling of viewpoint and microlens camera arrays are:  
```
@article{monteiro2019standard,
  title={Standard Plenoptic Cameras Mapping to Camera Arrays and Calibration based on DLT},
  author={Monteiro, Nuno Barroso and Barreto, Joao P and Gaspar, Jos{\'e} Ant{\'o}nio},
  journal={IEEE Transactions on Circuits and Systems for Video Technology},
  year={2019},
  publisher={IEEE}
}
```  
```
@inproceedings{monteiro2019generalized,
  title={Generalized Camera Array Model for Standard Plenoptic Cameras},
  author={Monteiro, Nuno Barroso and Gaspar, Jos{\'e} Ant{\'o}nio},
  booktitle={Iberian Robotics conference},
  pages={3--14},
  year={2019},
  organization={Springer}
}
```
