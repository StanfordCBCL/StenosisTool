
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
        x = self.output_layer(x)
        # x = self.output_relu(x)
        return x

class LightningNN(pl.LightningModule):
    
    def __init__(self, model, lr, revert_map):
        super(LightningNN, self).__init__()
        self.model = model
        self.lr = lr
        self.revert_map = revert_map
    
    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.huber_loss(y_hat, y)
        # Logging to TensorBoard by default
        self.log("train_loss", loss)
        return loss
    
    def validation_step(self, batch, batch_idx):
        x, y = batch
        #print(y[0])
        y_hat = self.model(x)
        #print(y_hat[0])
        #print(y_hat.shape)
        loss = nn.functional.huber_loss(y_hat, y)
        self.log("val_loss", loss)
        return loss
    
    def test_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = nn.functional.huber_loss(y_hat, y)
        self.log("test_loss", loss)
        return loss
    
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x, y = batch
        y_hat = self.model(x)
        revert(y_hat, self.revert_map)
        revert(y, self.revert_map)
        
        return torch.stack((y, y_hat), dim = 1)
    
    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=self.lr, weight_decay=.00005)
        return {
        "optimizer": optimizer,
        "lr_scheduler": {
            "scheduler": optim.lr_scheduler.ReduceLROnPlateau(optimizer,mode = 'min', factor = .1, patience=4, min_lr=1e-8, verbose = True),
            "monitor": "val_loss",
            "frequency": 1
            # If "monitor" references validation metrics, then "frequency" should be set to a
            # multiple of "trainer.check_val_every_n_epoch".
        },
    }
    
    
class Dataset0D(tdata.Dataset):
    
    def __init__(self, input_file, output_file, output_transformation = None):
        self.input = np.load(input_file)
        self.output= np.load(output_file)
        self.revert_map = None
        if output_transformation:
            self.output, self.revert_map = output_transformation(self.output)
        
    def __len__(self):
        return len(self.input)
    
    def __getitem__(self, idx):
        return torch.from_numpy((self.input[idx])).float(), torch.from_numpy(self.output[idx]).float()

# Normalization methods
def normalization(output):
    revert_map = []
    for i in range(len(output[0])):
        std = output[:, i].std()
        mean = output[:, i].mean()
        output[:, i] = (output[:, i] - mean) / std
        revert_map.append([mean, std])
    
    return output, revert_map

def revert(output, map_back):
    for i in range(len(output[0])):
        output[:, i] = (output[:, i] * map_back[i][1]) + map_back[i][0]
    return output


if __name__ == '__main__':
    
    #! Temp
    dir = Path('data/healthy/0080_0001/jc_solver_dir_0/artificial_stenosis/Manual_1')

    sim_dataset = Dataset0D(dir / 'training_data_old' / 'input.npy', dir / 'training_data_old' / 'output.npy', normalization)
    
    
    train_len = int(.8 * len(sim_dataset))
    val_len = int(.1 * len(sim_dataset))
    test_len = len(sim_dataset) - train_len - val_len
    
    train_dataset, val_dataset, test_dataset = tdata.random_split(sim_dataset, lengths=[train_len, val_len, test_len], generator= torch.Generator().manual_seed(42))
    
    
    train_loader = tdata.DataLoader(train_dataset, batch_size = 128, shuffle = True,)
    val_loader = tdata.DataLoader(val_dataset, batch_size=128, shuffle = False, )
    test_loader = tdata.DataLoader(test_dataset, batch_size=128, shuffle = False)
    # retrieve first value of Dataset for sizes
    input_data, output_data = sim_dataset[0]
    #print(output_data)
    # Arbitrary
    nnmodel = BasicNN(input_neurons=len(input_data), output_neurons=len(output_data), hidden_layers=3, neurons_per_layer=1000)
    litmodel = LightningNN(nnmodel, lr = 1e-3, revert_map = sim_dataset.revert_map)
    
    all_results_folder = dir / 'training_results'
    if not os.path.exists(all_results_folder):
        os.mkdir(all_results_folder)
    
    cur_results_folder = all_results_folder / 'run1'
    if not os.path.exists(cur_results_folder):
        os.mkdir(cur_results_folder)
    
    # checkpointing
    # Init ModelCheckpoint callback, monitoring 'val_loss'
    checkpoint_callback = ModelCheckpoint(monitor="val_loss")
    early_stop = EarlyStopping(monitor="val_loss", mode="min",check_finite=True, patience=10,  )

    csv_logger = CSVLogger(cur_results_folder)
    trainer = pl.Trainer( max_epochs=500, accelerator="gpu", default_root_dir=cur_results_folder, callbacks=[checkpoint_callback, early_stop], logger = csv_logger, log_every_n_steps=5)#, fast_dev_run=True)
    trainer.fit(model=litmodel, train_dataloaders=train_loader, val_dataloaders=val_loader)
    
    # test and save test dataloader
    trainer.test(model=litmodel, dataloaders=test_loader, ckpt_path='best', verbose = True)
    torch.save(test_loader, dir / "training_results" / "run1" / "lightning_logs" / f"version_{csv_logger.version}" / "test_dataloader.pt")
    
    # predict on the test loader and get normalized results
    rez = trainer.predict(model=litmodel, dataloaders=test_loader, ckpt_path="best", return_predictions=True)
    # retrieve x
    x = []
    for i in range(len(test_loader.dataset)):
        x.append(test_loader.dataset[i][0])
    x = torch.vstack(x)
    rez = torch.vstack(rez)
    
    torch.save(x, dir / "training_results" /  "run1" / "lightning_logs" / f"version_{csv_logger.version}" / "predict_input.pt")
    torch.save(rez, dir / "training_results" /  "run1" / "lightning_logs" / f"version_{csv_logger.version}" / "predict_output.pt")