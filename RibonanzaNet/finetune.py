import pandas as pd
import torch
import matplotlib.pyplot as plt
import numpy as np
import torch
import time
from torch.utils.data import Dataset, DataLoader
import random
import sys
from Network import *
import yaml

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

    
def load_dataset(datafile, length_limit = 850, length_min = 100):

    dataset = []
    file_ = open(datafile).readlines()

    for enum, line_ in enumerate(file_):

        if line_[0] == ">":

            if len(file_[enum+1][:-1]) > length_limit: continue
            if len(file_[enum+1][:-1]) < length_min: continue
            if len(file_[enum+1][:-1]) != len(file_[enum+2][:-2].split(" ")):
                print("missmatch")
                print(len(file_[enum+1][:-1]))
                print(file_[enum+2][:-1].split(" "))
                print(len(file_[enum+2][:-2].split(" ")))

                continue

            dataset.append([line_[:-1],file_[enum+1][:-1], file_[enum+2][:-2]])

    return dataset


class RNA_Dataset2(torch.utils.data.Dataset):
    def __init__(self,data):
        self.data=data
        self.tokens={nt:i for i,nt in enumerate('ACGT')}
        self.counter = 0

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence=[self.tokens[nt] for nt in self.data[idx][1]]
        sequence=np.array(sequence)
        sequence=torch.tensor(sequence)
        target = self.data[idx][2].split(" ")
        target = [float(t) for t in target  if t != ""]
        target=torch.tensor(target)

        mask = [1 if i=="A" or i =="C" else 0 for i in self.data[idx][1]]
        mask=torch.tensor(mask)

        if torch.isnan(target).sum() >0:
            target = torch.zeros(target.size())
            self.counter = self.counter + 1

        assert torch.isnan(target).sum() == 0, "NaN values found in target data!"

        return {'sequence':sequence, 'id':self.data[idx][0], 'target': target, 'mask':mask}


class RNA_Dataset(torch.utils.data.Dataset):
    def __init__(self,data):
        self.data=data
        self.tokens={nt:i for i,nt in enumerate('ACGU')}

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence=[self.tokens[nt] for nt in self.data.loc[idx,'sequence']]
        sequence=np.array(sequence)
        sequence=torch.tensor(sequence)

        return {'sequence':sequence}


class Config:
    def __init__(self, **entries):
        self.__dict__.update(entries)
        self.entries=entries

    def print(self):
        print(self.entries)

def load_config_from_yaml(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return Config(**config)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="finetuning RibonanzaNet")
    
    # Define arguments
    parser.add_argument("--config_path", type=str, required=True, help="Path to config file")
    parser.add_argument("--dataset_path", type=str, required=True, help="Path to dms dataset file")
    parser.add_argument("--input_model_path", type=str, required=True, help="Path to RibonanzaNet.pt model")
    parser.add_argument("--output_model_path", type=str, required=True, help="output path for final model")
    parser.add_argument("--test_dataset_path", type=str, required=False, help="output path for final model")
    parser.add_argument("--test_output_folder", type=str, required=False, help="output folder for test samples")
    # Parse arguments
    args = parser.parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    print("################################## start fine-tuning ####################################################")

    conf = load_config_from_yaml(args.config_path)
    model=RibonanzaNet(load_config_from_yaml(args.config_path)).to(device)
    dataset = load_dataset(args.dms_path)
    random.shuffle(dataset)
    train_dataset = RNA_Dataset2(dataset)

    model.load_state_dict(torch.load(args.input_model_path,map_location=device))
    model.train()

    optimizer = torch.optim.Adam(model.parameters(), lr= 0.00001)

    epochs = 1
    lossfull = 0
    sumsvals = 0
    epoch_loss = 0
    last_epoch_loss = np.inf


    for epoch in range(epochs):
        epoch_loss = 0

        for i in range(len(train_dataset)):
            example=train_dataset[i]
            sequence=example['sequence'].to(device).unsqueeze(0)
            id_ = example['id']
            target = example['target'].to(device)
            mask = example['mask'].to(device)

            assert torch.isnan(sequence).sum() == 0, "NaN values found in input data!"
            assert torch.isinf(sequence).sum() == 0, "Inf values found in input data!"
            assert torch.isnan(target).sum() == 0, "NaN values found in target data!"
            assert torch.isinf(target).sum() == 0, "Inf values found in target data!"
        
            res = model(sequence,torch.ones_like(sequence))
            loss = torch.sum(torch.nn.functional.mse_loss(res[..., 1:2].squeeze().squeeze(), target, reduction="none")*mask)/sum(mask)

            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
            print(f"loss: {loss.item()}")

        if epoch_loss < last_epoch_loss: torch.save(model.state_dict(), args.output_model_path)
        last_epoch_loss = epoch_loss

    print("################################## save test dataset ####################################################")
    if args.test_dataset_path != None: 
    
    
        if args.test_output_folder == None: 
            raise Exception("if you want to make use of a test set, please also set test_output_folder") 
    
        dataset2 = load_dataset(args.test_dataset_path)
        test_dataset=RNA_Dataset2(dataset2)
        test_preds=[]
        model.eval()
        for i in range(len(test_dataset)):
            example=test_dataset[i]
            sequence=example['sequence'].to(device).unsqueeze(0)
            id_ = example['id']
            with torch.no_grad():
                res = model(sequence,torch.ones_like(sequence).to(device))
            torch.save(res, args.test_output_folder + id_ + ".pkl")
