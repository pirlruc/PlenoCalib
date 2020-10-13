# PlenoCalib
Methods for calibrating standard plenoptic cameras.

# Usage
1. Install the toolbox using (it might be needed some manual download of the dependencies identified):  
```
    ins = setup.Setup;
    ins.install;
```  

2. Download Lytro datasets and white images from the following link:  
https://ieee-dataport.org/open-access/lytro-illum-calibration-dataset

3. Obtain white images database, decode calibration images and extract corner features using LFToolbox. The toolbox is added as a submodule.  

4. Calibrate the camera using the calibration procedure proposed. An example is provided [here](tcsvt_main.m "Calibration Example").

# Citing
The appropriate citations for the calibration and modeling of viewpoint and microlens camera arrays are:  
1. Standard Plenoptic Cameras Mapping to Camera Arrays and Calibration based on DLT, N. B. Monteiro, J. P. Barreto and J. Gaspar, IEEE Transactions on Circuits and Systems for Video Technology (TCSVT), November 2019. 
[[Journal](https://doi.org/10.1109/TCSVT.2019.2954305 "Journal URL")] 
[[PDF](http://vislab.isr.ist.utl.pt/wp-content/uploads/2020/01/nmonteiro-tcsvt2019.pdf "PDF File")]
```
@article{monteiro2019standard,
  title={Standard Plenoptic Cameras Mapping to Camera Arrays and Calibration based on DLT},
  author={Monteiro, Nuno Barroso and Barreto, Joao P and Gaspar, Jos{\'e} Ant{\'o}nio},
  journal={IEEE Transactions on Circuits and Systems for Video Technology},
  year={2019},
  publisher={IEEE}
}
```  
2. Generalized Camera Array Model for Standard Plenoptic Cameras, N. B. Monteiro and J. Gaspar, Iberian Robotics Conference, (ROBOT 2019). Porto, Portugal, 20-22 November 2019. 
[[Conference](https://doi.org/10.1007/978-3-030-36150-1_1 "Conference Proceedings")]
[[PDF](http://vislab.isr.ist.utl.pt/wp-content/uploads/2020/01/nmonteiro-robot2019.pdf "PDF File")]
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
