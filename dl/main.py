#!/usr/bin/env
# -*- coding: utf-8 -*-

"""
main.py: 
"""
import sys 
sys.path.append("..")
import argparse
import logging
import importlib
from dl import logger, args
from dl.train import Trainer
from dl.info import *
import shutil
from dl.dataloader import Dataset,  collate_fn
from torch.utils.data import DataLoader
from dl.utils import *
import torch
import os.path as osp
import json
LOG_LEVELS = list(logging._levelToName.values())


def main(args):
    """Main train and evaluation function.
        Parameters
        ----------
        args: argparse.Namespace
            Arguments
    """
    set_seed(args.seed)
    logger.info(os.getcwd())
    logger.info('Process Id: {}'.format(os.getpid()))
    if args.num_gpus== 1:
        os.environ["CUDA_VISIBLE_DEVICES"]=f'{args.vis_cuda}'

    if args.seed is None:
        args.seed = random.randint(1, 10000)
    set_seed(args.seed)

    os.environ['CUDA_VISIBLE_DEVICES'] = f'{args.vis_cuda}'
    
    #rm_dump_folder(osp.join(SAVE_PATH, "exp", "{}_{}".format(args.model, "0")))
    # device = get_device(is_gpu=not args.no_cuda)
    exp_path = osp.join(SAVE_PATH, "exp", "{}_{}".format(args.model, args.model_num))
    os.makedirs(exp_path, exist_ok=True)
    logger.info("Root directory for saving and loading experiments: {}".format(exp_path))

    model_bk_f = osp.join(exp_path, "{}_{}.py".format(args.model.lower(), args.model_num))

    if args.mode == "eval" or args.restore:
        # tmp_model_path = osp.join(osp.dirname(script_dir), 'tmp')
        # os.makedirs(tmp_model_path, exist_ok=True)
        # shutil.copy(model_bk_f, osp.join(tmp_model_path))
        spec = importlib.util.spec_from_file_location("DLModel", model_bk_f)
        MODULE = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(MODULE)
        # MODULE = importlib.import_module('models.tmp.{}_{}'.format(args.model.lower(), args.model_num))
        train_config_path = osp.join(exp_path, "weight")
        args = load_config(train_config_path, 'train_config_0', args)

    else:
        spec = importlib.util.spec_from_file_location("DLModel", f"./models/{args.model.lower()}.py")
        MODULE = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(MODULE)
        # MODULE = importlib.import_module("models." + args.model.lower())
        shutil.copy(osp.join("models", "{}.py".format(args.model.lower())), model_bk_f)
        
    device = torch.device("cuda" if args.num_gpus > 0 else "cpu")
    model = MODULE.Net(args)
    #torch.cuda.set_device(args.local_rank)
    #torch.distributed.init_process_group(backend='nccl')
    dl_dict = dict()
    num_train_sample = 1
    for phase in ['train', 'valid', 'test']:

        ds = Dataset(args, phase=phase, device=device)

        if phase == 'train':
            shuffle = True
            num_train_sample = len(ds)
        else:
            shuffle = False

        dl_dict[phase] = DataLoader(ds, args.batch_size, num_workers=args.num_workers,
                                    shuffle=shuffle, collate_fn=collate_fn)

    trainable_count = count_parameters(model)
    logger.info(
        "The ratio is {} ({} / {})".format(trainable_count / num_train_sample, trainable_count, num_train_sample))

    trainer = Trainer(args, model, dl_dict, exp_path, device)
    if args.mode == "train":
        trainer.train()
        trainer.args.mode = "eval"
    #trainer.eval("valid")
    trainer.eval("test")
 
            


if __name__ == "__main__":
    script_dir = osp.realpath(__file__)
    os.chdir(osp.dirname(script_dir))
    print("Processing ID: {}".format(os.getpid()))
    main(args)
