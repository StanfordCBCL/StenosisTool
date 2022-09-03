from torch import nn
from torch import optim
import torch
import torch.utils.data as tdata
import pytorch_lightning as pl
import numpy as np
import os
from pytorch_lightning.callbacks import ModelCheckpoint

from pathlib import Path


class BasicNN(nn.Module):
    
    def __init__(self, input_neurons, output_neurons, hidden_layers, neurons_per_layer):
        super(BasicNN, self).__init__()
        self.input_layer = nn.Linear(input_neurons, neurons_per_layer )
        
        self.input_relu = nn.ReLU()
        self.hidden = nn.Sequential()
        if hidden_layers < 1:
            raise ValueError('hidden layers must be > 0')
        else:
            for i in range(hidden_layers):
                self.hidden.append(nn.Linear(neurons_per_layer, neurons_per_layer))
                self.hidden.append(nn.ReLU())
        self.output_layer = nn.Linear(neurons_per_layer, output_neurons)
        self.output_relu = nn.ReLU()
    
    def forward(self, x):
        x = self.input_layer(x)
        x = self.input_relu(x)
        x = self.hidden(x)
        x = self.output_relu(self.output_layer(x))
        return x

class LightningNN(pl.LightningModule):
    
    def __init__(self, model, lr):
        super(LightningNN, self).__init__()
        self.model = model
        self.lr = lr
    
    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.mse_loss(y_hat, y)
        # Logging to TensorBoard by default
        self.log("train_loss", loss)
        return loss
    
    def validation_step(self, batch, batch_idx):
        x, y = batch
        print(y[0])
        y_hat = self.model(x)
        print(y_hat[0])
        print(y_hat.shape)
        loss = nn.functional.mse_loss(y_hat, y)
        self.log("val_loss", loss)
        return loss
    
    def test_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.mse_loss(y_hat, y)
        self.log("test_loss", loss)
        return loss
    
    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=self.lr)
        return optimizer
    
    
class Dataset0D(tdata.Dataset):
    
    def __init__(self, input_file, output_file):
        self.input = np.load(input_file)
        self.output= np.load(output_file)
        
    def __len__(self):
        return len(self.input)
    
    def __getitem__(self, idx):
        return torch.from_numpy(self.input[idx]).float(), torch.from_numpy(np.log(self.output[idx][1::2]+1)).float()



if __name__ == '__main__':
    
    #! Temp
    dir = Path('data/healthy/0080_0001/jc_solver_dir_0/artificial_stenosis/Manual_0')

    sim_dataset = Dataset0D(dir / 'training_data' / 'input.npy', dir / 'training_data' / 'output.npy')
    
    
    
    
    train_len = int(.8 * len(sim_dataset))
    val_len = int(.1 * len(sim_dataset))
    test_len = len(sim_dataset) - train_len - val_len
    train_dataset, val_dataset, test_dataset = tdata.random_split(sim_dataset, lengths=[train_len, val_len, test_len])
    
    
    train_loader = tdata.DataLoader(train_dataset, batch_size = 128, shuffle = True, num_workers=8)
    val_loader = tdata.DataLoader(val_dataset, batch_size=128, shuffle = False,num_workers=8 )
    test_loader = tdata.DataLoader(test_dataset, batch_size=128, shuffle = False, num_workers=8)
    
    # retrieve first value of Dataset for sizes
    input_data, output_data = sim_dataset[0]
    print(output_data)
    # Arbitrary
    nnmodel = BasicNN(input_neurons=len(input_data), output_neurons=len(output_data), hidden_layers=3, neurons_per_layer=1000)
    litmodel = LightningNN(nnmodel, lr = 1e-3)
    
    all_results_folder = dir / 'training_results'
    if not os.path.exists(all_results_folder):
        os.mkdir(all_results_folder)
    
    cur_results_folder = all_results_folder / 'run1'
    if not os.path.exists(cur_results_folder):
        os.mkdir(cur_results_folder)
    
    # checkpointing
    # Init ModelCheckpoint callback, monitoring 'val_loss'
    checkpoint_callback = ModelCheckpoint(monitor="val_loss")

    trainer = pl.Trainer( max_epochs=100, accelerator="gpu", default_root_dir=cur_results_folder, callbacks=[checkpoint_callback], log_every_n_steps=5)#, fast_dev_run=True)
    trainer.fit(model=litmodel, train_dataloaders=train_loader, val_dataloaders=val_loader)
    
    
    trainer.test(model=litmodel, dataloaders=test_loader, ckpt_path='best', verbose = True)