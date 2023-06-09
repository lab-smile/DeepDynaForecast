# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 15:11:58 2021

@author: Suncy
"""

import numpy as np # linear algebra
import dgl
import torch
import torch.nn as nn
from torch.autograd import Variable
from dgl.nn.pytorch.conv import GraphConv
import torch.nn.functional as F
#from dl import feat_dict


######################################## Network Model ############################################

class Net(nn.Module):
    def __init__(self, args):
        super(Net, self).__init__()
        h_feat = 128
        n_iter = 20
        in_feats = 16
        num_classes = 4
        self.device = torch.device("cuda" if args.num_gpus > 0 else "cpu")
        
        self.conv1 = GraphConv(in_feats, h_feat)
        self.conv2 = GraphConv(h_feat, h_feat)
        self.conv3 = GraphConv(h_feat, h_feat)
        self.conv4 = GraphConv(h_feat, h_feat)
        self.conv5 = GraphConv(h_feat, h_feat)
        self.conv6 = GraphConv(h_feat, h_feat)
        self.conv7 = GraphConv(h_feat, h_feat)
        self.conv8 = GraphConv(h_feat, h_feat)
        self.conv9 = GraphConv(h_feat, h_feat)
        self.conv10 = GraphConv(h_feat, h_feat)
        self.conv11 = GraphConv(h_feat, h_feat)
        self.conv12 = GraphConv(h_feat, h_feat)
        self.conv13 = GraphConv(h_feat, h_feat)
        self.conv14 = GraphConv(h_feat, h_feat)
        self.conv15 = GraphConv(h_feat, h_feat)
        self.conv16 = GraphConv(h_feat, h_feat)
        self.conv17 = GraphConv(h_feat, h_feat)
        self.conv18 = GraphConv(h_feat, h_feat)
        self.conv19 = GraphConv(h_feat, h_feat)
        self.conv20 = GraphConv(h_feat, h_feat)
        self.m = nn.LeakyReLU()
        
        self.fc = nn.Linear(h_feat, num_classes)

    
    def forward(self, g):
        info = dict()
        # num of nodes & edges in batched g
        edge_feat = g.edata['feat'][:,0]
        node_feat = g.ndata['feat']
        
        num_nodes = g.ndata['feat'].shape[0]
        num_edges = g.edata['feat'].shape[0]
        
        h = self.conv1(g, node_feat, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv2(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv3(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv4(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv5(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv6(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv7(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv8(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv9(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv10(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv11(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv12(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv13(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv14(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv15(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv16(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv17(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv18(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv19(g, h, edge_weight=edge_feat)
        h = self.m(h)
        h = self.conv20(g, h, edge_weight=edge_feat)
        h = self.m(h)

        output = self.fc(h)

        
        return output, info
    
    def ce_loss(self, y_pred, y_true, weight=None):
        # print(y_pred.shape, y_true.shape, weight.shape)
        ce = F.cross_entropy(y_pred, y_true, weight=weight, size_average=None, reduce=None, reduction='mean')
        return {"loss": ce}
