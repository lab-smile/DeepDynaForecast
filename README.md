# DeepDynaForecast

This is a pytorch implementation of our research:  
DeepDynaForecast: Phylogenetic-informed graph deep learning for epidemic transmission dynamic prediction, which has been accepted by 
PLOS Computational Biology.
If you find our work useful in your research or publication, please cite our work.

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


## Coming soon...
A website for users to apply DeepDynaForecast on Phylogenetic Trees.

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
