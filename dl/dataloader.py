# -*- coding: utf-8 -*-
"""
@author: Suncy

dataloader for Graph NN Models
"""

import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from torchvision import transforms, utils
import os.path as osp
import dgl
import torch
import json
from dl import  logger
from dgl.data import DGLDataset
from collections import Counter
from torch.utils.data import Dataset

# dataloader for graph models
class Dataset(DGLDataset):
    def __init__(self, args, phase, device="cpu"):
        self.device = device
        ds_folder = osp.join(args.ds_dir, args.ds_name, args.ds_split)
        
        self.phase = phase
        if phase in ["train", "valid", "test"]:
            self.node_df = pd.read_csv(f"{ds_folder}/{phase}.csv", low_memory=False)
            self.edge_df = pd.read_csv(f"{ds_folder}/{phase}_edge.csv", low_memory=False)
        else:
            raise NotImplementedError
        self.tree_ids = self.node_df["sim"].unique()  # num of trees
        
        #self.node_df['dynamic_cat0'] = self.node_df['dynamic_cat'].map(lambda x: 1 if x == 0 else 0)
        #self.node_df['dynamic_cat1'] = self.node_df['dynamic_cat'].map(lambda x: 1 if x == 1 else 0)
        #self.node_df['dynamic_cat2'] = self.node_df['dynamic_cat'].map(lambda x: 1 if x == 2 else 0)
        #self.node_df['dynamic_cat3'] = 0

        #self.node_feat_cols = ['dynamic_cat0', 'dynamic_cat1', 'dynamic_cat2']
        self.edge_feat_cols = ['weight1_arsinh-norm','weight2_arsinh-norm']
        self.node_label_cols = 'dynamic_cat'
        
        cid_dict = {'Background':0, 'c1':1, 'c2':2, 'c3':3, 'c4':4, 'c5':5, 'c6':6}
        
        state_dict = {'A':0, 'B':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6}
        
        self.node_df['c_id'] = self.node_df['cluster_id'].map(cid_dict)
        
        self.node_df['s'] = self.node_df['state'].map(state_dict)
        
        self.node_feat_org = ['sim', 'c_id', 's', 'node']

        self.add_self_loop = args.add_self_loop
        self.bidirection = args.bidirection

    def process(self):        
        pass

    def __getitem__(self, index):

        tree_id = self.tree_ids[index]  # tree of index

        # dgl tree of index
        onetree_node_df = self.node_df[self.node_df['sim'] == tree_id]
        onetree_edge_df = self.edge_df[self.edge_df['sim'] == tree_id]
        src_ids = torch.tensor(onetree_edge_df['new_from'].values)
        dst_ids = torch.tensor(onetree_edge_df['new_to'].values)
        src_ids -= 1
        dst_ids -= 1
        g = dgl.graph((src_ids, dst_ids))  # create dgl
        sorted_onetree_node_df = onetree_node_df.sort_values(by='node')
        
        leaves = list(set(onetree_edge_df['new_to'].values)-set(onetree_edge_df['new_from'].values))
        leaves = [i-1 for i in leaves]
        
        # assign features and labels for background nodes
        #node_feat = np.zeros([len(sorted_onetree_node_df),20])
        #node_feat = np.random.rand(len(sorted_onetree_node_df),16)
        node_feat = np.zeros([len(onetree_edge_df)+1,16])+0.5
        node_label = np.zeros([len(sorted_onetree_node_df),])+3
        for leaf in leaves:
            #node_feat[leaf,:] = np.zeros([1,3])
            #if sorted_onetree_node_df['cluster_id'].values[leaf] == 'Background':
            node_label[leaf] = sorted_onetree_node_df[self.node_label_cols].values[leaf]
        num_nodes = node_feat.shape[0]
        num_feat = node_feat.shape[1]

        # assign features for nodes and edges, assign labels
        g.ndata["feat"] = torch.tensor(node_feat, dtype=torch.float32)
        g.ndata["org_feat"] = torch.tensor(sorted_onetree_node_df[self.node_feat_org].values, dtype=torch.float32)
        g.ndata["label"] = torch.tensor(node_label, dtype=torch.int64)
        g.edata["feat"] = torch.tensor(onetree_edge_df[self.edge_feat_cols].values, dtype=torch.float32)
        

        if self.add_self_loop:
            g = dgl.add_self_loop(g)  # TODO: Add self-loop with self-edge weight filled with zero
        if self.bidirection:
            g = dgl.add_reverse_edges(g, copy_ndata=True, copy_edata=True)

        g = g.to(self.device)
        return g
        
    def __len__(self):
        return len(self.tree_ids)   # number of trees



# create batch of trees(aggregate multiples trees to a single tree)
def collate_fn(batch_graphs):
    g = dgl.batch(batch_graphs)
    return g


def gen_label_weight(args):
    # Get the weights for the unbalanced sample based on the positive sample
    # weights inversely proportional to class frequencies in the training data
    '''
    ds_folder = osp.join(args.ds_dir, args.ds_name, args.ds_split)
    node_df = pd.read_csv(f'{ds_folder}/train.csv')

    node_label = node_df[args.node_label_cols].values
    label_counter = Counter(node_label)
    n_samples = len(node_label)
    n_classes = len(label_counter)

    label_weights = [n_samples / (n_classes * label_counter[i]) for i in range(n_classes)]
    '''
    #label_weights = [0.3410502395067024, 18.27036152196542, 76.06173713292357]
    label_weights = [0.3406808892621469, 51.519496230514754, 22.079112613356273]
    
    if args.loss_ignore_bg:
        label_weights.append(0)

    return label_weights





