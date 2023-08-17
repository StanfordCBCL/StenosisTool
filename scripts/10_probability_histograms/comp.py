# File: comp.py
# File Created: Thursday, 22nd September 2022 7:45:48 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 16th August 2023 11:42:49 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
#! Description: compares final output.


from torch import nn
from torch import optim
import torch
import torch.utils.data as tdata
import pytorch_lightning as pl
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from svinterface.utils.io import read_json
from svinterface.utils.misc import d2m

class BasicNN(nn.Module):
    
    def __init__(self, input_neurons, output_neurons, hidden_layers, neurons_per_layer):
        super(BasicNN, self).__init__()
        self.input_layer = nn.Linear(input_neurons, neurons_per_layer )
        self.relu = nn.Tanh()
        self.hidden = nn.Sequential()
        if hidden_layers < 1:
            raise ValueError('hidden layers must be > 0')
        else:
            for i in range(hidden_layers):
                self.hidden.append(nn.Linear(neurons_per_layer, neurons_per_layer))
                self.hidden.append(self.relu)
        self.output_layer = nn.Linear(neurons_per_layer, output_neurons)
        #self.output_relu = nn.ReLU()
    
    def forward(self, x):
        x = self.input_layer(x)
        x = self.relu(x)
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
        optimizer = optim.Adam(self.parameters(), lr=self.lr)# weight_decay=.00005)
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

# Normalization methods
def normalization(output, revert_map = None):
    if revert_map is None:
        revert_map = []
        for i in range(len(output[0])):
            std = output[:, i].std()
            mean = output[:, i].mean()
            output[:, i] = (output[:, i] - mean) / std if std != 0 else 0
            revert_map.append([mean, std])
        
        return output, revert_map
    
    else:
        for i in range(len(output[0])):
            std = revert_map[i][1]
            mean = revert_map[i][0]
            output[:, i] = (output[:, i] - mean) / std if std != 0 else 0
        return output, revert_map

def revert(output, map_back):
    for i in range(len(output[0])):
        output[:, i] = (output[:, i] * map_back[i][1]) + map_back[i][0]
    return output

class PredictDataset(tdata.Dataset):
    def __init__(self, inputs):
        self.inputs = inputs
        
    def __len__(self):
        ''' currently undefined '''
        return len(self.inputs)
    
    def __getitem__(self, idx):
        return torch.from_numpy((self.inputs[idx])).float()

class RepairDistribution():
    
    class RepairFixed():
        """A repair is fixed/completed and sampling always returns the fixed value"""
        def __init__(self, c):
            self.c = c
        
        def sample(self):
            return self.c
    
    class RepairRandom():
        """A repair is a random distribution"""
        class Category():
            def __init__(self, lower, upper):
                self.lower = lower
                self.upper = upper
            
            def sample(self):
                return np.random.uniform(self.lower, self.upper)
            
        class Fail(Category):
            def __init__(self, lower, upper):
                super().__init__(lower, upper)
    
        class Moderate(Category):
            def __init__(self, lower, upper):
                super().__init__(lower, upper)
        
        class Success(Category):
            def __init__(self, lower, upper):
                super().__init__(lower, upper)
                
        def __init__(self, category_probs):
            '''categories = (Fail prob, Moderate prob, Success prob)'''
            assert len(category_probs) == 3 # must be probabilities for (Fail, Moderate, Success)
            assert sum(category_probs) == 1 # must add to 1
            self.category_probs = category_probs
            
            # set up categories as 0-.2, .2-.8, .8-1.
            self.categories = (self.Fail(lower = 0,
                                         upper = .2),
                               self.Moderate(lower = .2,
                                             upper = .8),
                               self.Success(lower = .8,
                                            upper = 1)) # split 20, 60, 20
            
            self.frozen = False
        
        def sample(self):
            # if frozen, just return Fail.lower since it should be disabled
            if self.frozen:
                return self.categories[0].lower
            
            n = np.random.uniform()
            # check if in Fail
            if (n - self.category_probs[0] < 0):
                return self.categories[0].sample() # sample from Fail
            elif (n - self.category_probs[0] - self.category_probs[1] < 0):
                return self.categories[1].sample() # sample from Moderate
            else:
                return self.categories[2].sample() # sample from Success
    
    
    def __init__(self, num_repairs, category_probs = (.3, .4, .3)):
        
        self.category_probs = category_probs
        self.num_repairs = num_repairs
        self.repair_points = [self.RepairRandom(category_probs=self.category_probs) for i in range(num_repairs)]
    
    def __len__(self):
        return self.num_repairs
        
    def fixed(self, repair_idx, c):
        ''' change a particular vessel to fixed '''
        self.repair_points[repair_idx] = self.RepairFixed(c)
    
    def undo_fixed(self, repair_idx):
        ''' only for testing purposed... to undo a fixture.'''    
        self.repair_points[repair_idx] = self.RepairRandom(category_probs=self.category_probs)
        
    def freeze(self, to_freeze = []):
        ''' freeze particular vessels (not fixed, just not changing). If a vessel is already frozen or fixed, do nothing'''
        for repair_idx in to_freeze:
            if type(self.repair_points[repair_idx]) == self.RepairRandom:
                self.repair_points[repair_idx].frozen = True
    
    def unfreeze_all(self):
        ''' unfreeze all vessels that aren't fixed'''
        for repair_point in self.repair_points:
            if type(repair_point) == self.RepairRandom:
                repair_point.frozen = False

    def unfreeze(self, to_unfreeze = []):
        ''' unfreezes previously frozen vessels. If a vessel is already unfrozen, do nothing'''
        for repair_idx in to_unfreeze:
            if type(self.repair_points[repair_idx]) == self.RepairRandom:
                self.repair_points[repair_idx].frozen = False
        
    def sample(self):
        ''' generates a sample '''
        return np.array([repair_point.sample() for repair_point in self.repair_points])
    
    def create_dataset(self, num_samples):
        input_data = np.array([self.sample() for i in range(num_samples)])
        return PredictDataset(input_data)
    
    def get_histograms(self, results):
        return [(np.histogram(results[:, 2*i], bins = 'auto'), np.histogram(results[:, 2*i+1], bins = 'auto')) for i in range(len(self))]
    
    def save_histograms(self, results, path):
        np.save(path, self.get_histograms(results),allow_pickle=True)
            
    def animate_histograms(self, results, save):

        fig, ax = plt.subplots(1, 2, figsize = (16, 8))
        def animate(vessel):
            P_vessel = vessel * 2
            Q_vessel = vessel * 2 + 1 
            ax[0].clear()
            ax[1].clear()
            ax[0].set_title('Pressures (mmHg)')
            ax[0].set_xlabel('Pressures (mmHg)')
            ax[0].set_ylabel('count')
            ax[1].set_title('Flows (ml/s)')
            ax[1].set_xlabel('Flows (ml/s)')
            ax[1].set_ylabel('count')
            P_vessel = vessel * 2
            Q_vessel = vessel * 2 + 1
            ax[0].hist(d2m(results[:, P_vessel]), bins = 'auto')
            ax[1].hist(results[:, Q_vessel], bins = 'auto')
            fig.suptitle(f'Vessel {vessel}')
            
        anim = FuncAnimation(fig, animate, frames = range(len(results[0])// 2), interval = 200, repeat = True)
        anim.save(save, writer = 'pillow', fps = .5)
        
    def run_test(self, model, trainer, best_chkpt,  num_samples = 4096, batch_size = 4096):
        test_dataset = self.create_dataset(num_samples)
        test_loader = tdata.DataLoader(test_dataset, batch_size = batch_size)
        # predict on the test loader and get normalized results
        rez = trainer.predict(model=model, dataloaders=test_loader, ckpt_path=best_chkpt, return_predictions=True)
        # retrieve x
        x = []
        for i in range(len(test_loader.dataset)):
            x.append(test_loader.dataset[i][0])
        x = torch.vstack(x)
        rez = torch.vstack(rez)
        return x, rez

class LightningNNPredictor(LightningNN):
    
    def __init__(self, model, lr, revert_map):
        super().__init__(model, lr, revert_map)
        
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x = batch
        y_hat = self.model(x)
        revert(y_hat, self.revert_map)
        return y_hat

if __name__ == '__main__':
    
    #! Temp hardcode path
    dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR')

    # load training data to recompute normalization map (#TODO: SAVE THE NORM NEXT TIME.)
    train_dataset = Dataset0D(dir / 'model_data' / 'train_data' / 'input.npy', dir / 'model_data' / 'train_data' / 'output.npy', normalization, revert_map=None)

    # Construct model
    # retrieve first value of Dataset for sizes
    input_data, output_data = train_dataset[0]
    nnmodel = BasicNN(input_neurons=len(input_data), output_neurons=len(output_data), hidden_layers=3, neurons_per_layer=1000)
    litmodel = LightningNNPredictor(nnmodel, lr = 1e-3, revert_map = train_dataset.revert_map)
    
    
    # create probability distribution
    distribution = RepairDistribution(num_repairs=len(input_data), 
                                      category_probs=(.3,.4,.3))
    
    #! TO be changed
    best_checkpoint = dir / "training_results" / "run1" / "lightning_logs" / "version_0" / "checkpoints" / "epoch=498-step=127744.ckpt"
    trainer = pl.Trainer(logger = False, accelerator='auto')
    
    # create a dir to save data
    prob_dir = dir / 'probability_histograms'
    prob_dir.mkdir(exist_ok=True)
    
    
    # starting at random, no fixture (All possible distributions)
    all_dir = prob_dir / 'all_distributions'
    all_dir.mkdir(exist_ok=True)
    x, yhat = distribution.run_test(model=litmodel,
                          trainer=trainer,
                          best_chkpt=best_checkpoint,
                          num_samples=4096*24,
                          batch_size=4096)
    distribution.animate_histograms(yhat, save = str(all_dir / 'histograms.gif'))
    distribution.save_histograms(yhat, path=str(all_dir / 'histograms.npy') )
    
    
    # fix repair 0 (LPA_proximal) at 
    fixz_dir = prob_dir / 'fixeda1_conditional'
    fixz_dir.mkdir(exist_ok=True)
    distribution.fixed(repair_idx=0,
                       c=.8)
    x, yhat = distribution.run_test(model=litmodel,
                          trainer=trainer,
                          best_chkpt=best_checkpoint,
                          num_samples=4096*24,
                          batch_size=4096)
    distribution.animate_histograms(yhat, save = str(fixz_dir / 'histograms.gif'))
    distribution.save_histograms(yhat, path=str(fixz_dir / 'histograms.npy') )
    