# DeepDynaForecast: Phylogenetic-informed graph deep learning for epidemic transmission dynamic prediction

During an outbreak or sustained epidemic, accurate prediction of patterns in transmission risk can reliably inform public health strategies. Projections indicating growth or decline of transmission for specific risk groups can significantly enhance the optimization of interventions, especially when resources are limited. To address this, we present DeepDynaForecast, a cutting-edge deep learning algorithm designed for forecasting pathogen transmission dynamics. Uniquely, DeepDynaForecast was trained on in-depth simulation data, classifying samples according to their dynamics (growth, static, or decline) with accuracy of 91.6%. We evaluated the model’s performance and application using simulated outbreak data and empirical, large-scale data from the HIV epidemic in Florida between 2012 and 2020. We conclude DeepDynaForecast represents a significant advancement in genomics-mediated pathogen transmission characterization and has the potential to catalyze new research directions within virology, molecular biology, and public health.


## Paper
This repository provides the official implementation of the model in the following paper:

**DeepDynaForecast: Phylogenetic-informed graph deep learning for epidemic transmission dynamic prediction**

Chaoyue Sun<sup>1</sup>, Ruogu Fang<sup>1,2,3</sup>, Marco Salemi<sup>4,5</sup>, Mattia Prosperi<sup>5,6\*</sup> and Brittany Rife
Magalis<sup>4,5\*</sup>

<sup>1</sup> Department of Electrical and Computer Engineering, Herbert Wertheim College of Engineering, University
of Florida, Gainesville, Florida, United States of America<br>
<sup>2</sup> J. Crayton Pruitt Family Department of Biomedical
Engineering, Herbert Wertheim College of Engineering, University of Florida, Gainesville, Florida, United
States of America<br>
<sup>3</sup> Center for Cognitive Aging and Memory, McKnight Brain Institute, University of Florida,
Gainesville, Florida, United States of America<br>
<sup>4</sup> Department of Pathology, Immunology, and Laboratory
Medicine, University of Florida, Gainesville, Florida, United States of America<br>
<sup>5</sup> Emerging Pathogens
Institute, University of Florida, Gainesville, Florida, United States of America<br>
<sup>6</sup> Department of Epidemiology,
University of Florida, Gainesville, Florida, United States of America<br>

## Dependencies
Please check the dependencies.txt.

## Datasets

- Agent-based simulation model [nosoi] to generate simulated trees and input datatables in R. Please see [sim] for the simulation codes.
- For Florida HIV epidemic data, please contact [FDOH] for availabiliy.

## Usage examples
 - Preparing dataset  
Both [original and preprocessed datasets] are avaliable. Please uzip them in the root folder. Jupyter notebook files in `aly` provide a detailed preprocessing step by step for edge feature. The provided preprocessed dataset is divided into train/validation/test sets.  
Before running the code, put the preprocessed dataset on a desired directory. By default, the data root is set as `data/preprocessed_data/split_rs123`.  
See: [dl/config.py]
```sh
data.add_argument("--ds_name", type=str, default="preprocessed_data", help="The name of dataset")
data.add_argument("--ds_dir", type=str, default="../data/", help="The base folder for data")
data.add_argument("--ds_split", type=str, default="split_rs123")
```

 - Run container  
The container is build by Docker and is avaliable at [here]. To run a shell within the container:
```sh
sudo docker run -it —gpus all -v ddf    
```

 - Train models  
Example commands 
```sh
# train a DeepDynaForecast model
python main.py --model ddf --model_num 0
# train a DeepDynaForecast model with specific setting
python main.py --model ddf --model_num 0 --batch_size 32 --init_lr 0.001 --min_lr 1e-6 --lr_decay_rate 0.1
```

## Results
To use the [well-trained neural network models], please download files and unzip it to `dl`.  
Please change the settings in [test/main_test.py] to run the test phase with different models.  
```sh
python main_test.py
```
Besides, scripts for generating figures in the main pages and supplementaries are avaliable at the `test` folder.


## Citation
If you use this code, please cite our papers:
```
@article{sun2024deepdynaforecast,
  title={DeepDynaForecast: Phylogenetic-informed graph deep learning for epidemic transmission dynamic prediction},
  author={Sun, Chaoyue and Fang, Ruogu and Salemi, Marco and Prosperi, Mattia and Rife Magalis, Brittany},
  journal={PLOS Computational Biology},
  volume={20},
  number={4},
  pages={e1011351},
  year={2024},
  publisher={Public Library of Science San Francisco, CA USA}
}
```

## Acknowledgement
The authors acknowledge University of Florida Research Computing for providing computational resources and support that have contributed to the research results reported in this publication.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [nosoi]: <https://github.com/slequime/nosoi>
   [here]: <https://www.dropbox.com/s/1mxzoiruldhacfx/DDF.tar?dl=0>
   [sim]: <https://github.com/lab-smile/DeepDynaForecast/tree/main/sim>
   [FDOH]: <Research@flhealth.gov>
   [original and preprocessed datasets]: <https://www.dropbox.com/s/ta7lsxtx04o3m61/cleaned_data.zip?dl=0>
   [dl/config.py]: <https://github.com/lab-smile/DeepDynaForecast/blob/main/dl/config.py>
   [well-trained neural network models]: <https://www.dropbox.com/s/aaozfa2wyhdkacg/saved_models.zip?dl=0>
   [test]: <https://github.com/lab-smile/DeepDynaForecast/tree/main/test>
   [test/main_test.py]: <https://github.com/lab-smile/DeepDynaForecast/blob/main/test/main_test.py>
