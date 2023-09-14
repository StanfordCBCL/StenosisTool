# File: train_nn.py
# File Created: Tuesday, 15th August 2023 1:30:50 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:42:45 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Trains a Neural network


from torch import nn
from torch import optim
import torch
import torch.utils.data as tdata
import pytorch_lightning as pl
import numpy as np
import os
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping
from pytorch_lightning.loggers import CSVLogger
from pathlib import Path


class BasicNN(nn.Module):
    """ Basic Neural Network """
    def __init__(self, input_neurons, output_neurons, hidden_layers, neurons_per_layer):
        super(BasicNN, self).__init__()
        self.input_layer = nn.Linear(input_neurons, neurons_per_layer )
        self.tanh = nn.Tanh()
        self.hidden = nn.Sequential()
        if hidden_layers < 1:
            raise ValueError('hidden layers must be > 0')
        else:
            for i in range(hidden_layers):
                self.hidden.append(nn.Linear(neurons_per_layer, neurons_per_layer))
                self.hidden.append(self.tanh)
        self.output_layer = nn.Linear(neurons_per_layer, output_neurons)
    
    def forward(self, x):
        x = self.input_layer(x)
        x = self.tanh(x)
        x = self.hidden(x)
        x = self.output_layer(x)
        return x

class LightningNN(pl.LightningModule):
    
    def __init__(self, model, lr, revert_map):
        super(LightningNN, self).__init__()
        self.model = model
        self.lr = lr
        self.revert_map = revert_map
    
    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.huber_loss(y_hat, y)
        self.log("train_loss", loss)
        return loss
    
    def validation_step(self, batch, batch_idx):
        """ Validation """
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.huber_loss(y_hat, y)
        self.log("val_loss", loss)
        return loss
    
    def test_step(self, batch, batch_idx):
        """ Testing """
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.huber_loss(y_hat, y)
        self.log("test_loss", loss)
        return loss
    
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        """ Prediction """
        x, y = batch
        y_hat = self.model(x)
        revert(y_hat, self.revert_map)
        revert(y, self.revert_map)
        
        return torch.stack((y, y_hat), dim = 1)
    
    def configure_optimizers(self):
        """ Optimizer """
        optimizer = optim.Adam(self.parameters(), lr=self.lr)# weight_decay=.00005)
        return {
        "optimizer": optimizer,
        "lr_scheduler": {
            "scheduler": optim.lr_scheduler.ReduceLROnPlateau(optimizer,mode = 'min', factor = .1, patience=5, min_lr=1e-8, verbose = True),
            "monitor": "val_loss",
            "frequency": 1
        },
    }
    
class Dataset0D(tdata.Dataset):
    """Dataset for input and outputs, performs a transformation"""
    
    def __init__(self, input_file, output_file, output_transformation = None, revert_map = None):
        self.input = np.load(input_file)
        self.output= np.load(output_file)
        self.revert_map = revert_map
        if output_transformation:
            self.output, self.revert_map = output_transformation(self.output, self.revert_map)
        
    def __len__(self):
        return len(self.input)
    
    def __getitem__(self, idx):
        return torch.from_numpy((self.input[idx])).float(), torch.from_numpy(self.output[idx]).float()

def normalization(output, revert_map = None):
    """ Performs a zscore normalization on targets.
    If revert_map is provided, then apply the revert_map rather than compute a new one (used for validation and test")"""
    
    # compute revert map
    if revert_map is None:
        revert_map = []
        for i in range(len(output[0])):
            std = output[:, i].std()
            mean = output[:, i].mean()
            output[:, i] = (output[:, i] - mean) / std if std != 0 else 0
            revert_map.append([mean, std])
        return output, revert_map
    
    # apply existing revert map
    else:
        for i in range(len(output[0])):
            std = revert_map[i][1]
            mean = revert_map[i][0]
            output[:, i] = (output[:, i] - mean) / std if std != 0 else 0
        return output, revert_map

# perform normalization revert
def revert(output, map_back):
    for i in range(len(output[0])):
        output[:, i] = (output[:, i] * map_back[i][1]) + map_back[i][0]
    return output


if __name__ == '__main__':
    
    #! Temp
    dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR')

    train_dataset = Dataset0D(dir / 'model_data' / 'train_data' / 'input.npy', dir / 'model_data' / 'train_data' / 'output.npy', normalization, revert_map=None)
    val_dataset = Dataset0D(dir / 'model_data' / 'val_data' / 'input.npy', dir / 'model_data' / 'val_data' / 'output.npy', normalization, revert_map=train_dataset.revert_map)
    test_dataset = Dataset0D(dir / 'model_data' / 'test_data' / 'input.npy', dir / 'model_data' / 'test_data' / 'output.npy', normalization, revert_map=train_dataset.revert_map)

    train_loader = tdata.DataLoader(train_dataset, batch_size = 128, shuffle = True,)
    val_loader = tdata.DataLoader(val_dataset, batch_size=128, shuffle = False, )
    test_loader = tdata.DataLoader(test_dataset, batch_size=128, shuffle = False)
    
    # retrieve first value of Dataset for sizes
    input_data, output_data = train_dataset[0]
    
    # construct model
    nnmodel = BasicNN(input_neurons=len(input_data), output_neurons=len(output_data), hidden_layers=3, neurons_per_layer=1000)
    litmodel = LightningNN(nnmodel, lr = 1e-3, revert_map = train_dataset.revert_map)
    
    all_results_folder = dir / 'training_results'
    if not os.path.exists(all_results_folder):
        os.mkdir(all_results_folder)
    
    cur_results_folder = all_results_folder / 'run1'
    if not os.path.exists(cur_results_folder):
        os.mkdir(cur_results_folder)
    
    # checkpointing
    # Init ModelCheckpoint callback, monitoring 'val_loss'
    checkpoint_callback = ModelCheckpoint(monitor="val_loss")
    early_stop = EarlyStopping(monitor="val_loss", mode="min",check_finite=True, patience=10)
    csv_logger = CSVLogger(cur_results_folder)
    
    # Trainer
    trainer = pl.Trainer( max_epochs=500, accelerator="auto", default_root_dir=cur_results_folder, callbacks=[checkpoint_callback, early_stop], logger = csv_logger, log_every_n_steps=5,)# fast_dev_run=True)
    trainer.fit(model=litmodel, train_dataloaders=train_loader, val_dataloaders=val_loader)
    
    # test and save test dataloader
    trainer.test(model=litmodel, dataloaders=test_loader, ckpt_path='best', verbose = True)

    # predict on the test loader and get normalized results
    rez = trainer.predict(model=litmodel, dataloaders=test_loader, ckpt_path="best", return_predictions=True)
    # retrieve x
    x = []
    for i in range(len(test_loader.dataset)):
        x.append(test_loader.dataset[i][0])
    x = torch.vstack(x)
    rez = torch.vstack(rez)
    
    # save results
    torch.save(x, dir / "training_results" /  "run1" / "lightning_logs" / f"version_{csv_logger.version}" / "predict_input.pt")
    torch.save(rez, dir / "training_results" /  "run1" / "lightning_logs" / f"version_{csv_logger.version}" / "predict_output.pt")