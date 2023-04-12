# DeepDynaForecast

This is a pytorch implementation of our research:  
DeepDynaForecast: Phylogenetic-informed graph deep learning for epidemic transmission dynamic prediction
If you find our work useful in your research or publication, please cite our work.

## Dependencies
Please check the dependencies.txt.

## Datasets

- Agent-based simulation model [nosoi] to generate simulated trees and input datatables in R.

## Usage examples

 - Run container  
The container is build by [singularity] and is avaliable at [here]. To run a shell within the container:
```sh
singularity shell --nv singularity_deepdynatree.sif    
```

 - Train models  
Example commands 
```sh
# train a PDGLSTM model
python main.py --model 'pdglstm' --model_num 0
# train a PDGLSTM model with specific setting
python main.py --model 'pdglstm' --model_num 0 --batch_size 32 --init_lr 0.001 --min_lr 1e-6 --lr_decay_rate 0.1
```

## Results
To use the [trained neural network models], please download files and unzip it to `dl`.  
Please change the settings in [test/main_test.py] to run the test phase with different models.  
```sh
python main_test.py
```
Besides, scripts for generating figures in the main pages and supplementaries are avaliable at [models/post_aly.ipynb] and [models/post_aly_all.ipynb] respectively.


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [nosoi]: <https://github.com/slequime/nosoi>
   [singularity]: <https://sylabs.io/guides/3.9/user-guide/>
   [here]: <https://genome.ufl.edu/download/ddt/singularity_deepdynatree.sif>
   [dl/config.py]: <https://github.com/salemilab/DeepDynaTree/blob/main/dl/config.py>
   [models/ml_models.ipynb]: <https://github.com/salemilab/DeepDynaTree/blob/main/models/ml_models.ipynb>
   [test/ml_test.ipynb]: <https://github.com/salemilab/DeepDynaTree/blob/main/test/ml_test.ipynb>
   [test/main_test.py]: <https://github.com/salemilab/DeepDynaTree/blob/main/test/main_test.py>
   [models/post_aly.ipynb]: <https://github.com/salemilab/DeepDynaTree/blob/main/models/post_aly.ipynb>
   [models/post_aly_all.ipynb]: <https://github.com/salemilab/DeepDynaTree/blob/main/models/post_aly_all.ipynb>
