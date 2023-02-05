# File: comp.py
# File Created: Thursday, 22nd September 2022 7:45:48 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 23rd January 2023 7:46:03 pm
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
from src.file_io import read_json
from src.misc import d2m
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
        #print(y[0])
        input_filter = (x< 1.0).any(dim = 1) == False
        x = x[input_filter]
        y = y[input_filter]
        
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
            output[:, i] = (output[:, i] - mean) / std
            revert_map.append([mean, std])
        
        return output, revert_map
    
    else:
        for i in range(len(output[0])):
            std = revert_map[i][1]
            mean = revert_map[i][0]
            output[:, i] = (output[:, i] - mean) / std
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

class StenosisDistribution():
    
    class StenosisFixed():
        def __init__(self, r):
            self.r = r
        
        def sample(self):
            return self.r
    
    class StenosisRandom():
        
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
                
        def __init__(self, category_probs, radius_range):
            '''categories = (Fail prob, Moderate prob, Success prob)'''
            assert len(category_probs) == 3 # must be probabilities for (Fail, Moderate, Success)
            assert sum(category_probs) == 1 # must add to 1
            self.category_probs = category_probs

            rnge = radius_range[1] - radius_range[0]
            
            
            self.categories = (self.Fail(lower = radius_range[0],
                                         upper = radius_range[0] + rnge * .2),
                               self.Moderate(lower = radius_range[0] + rnge * .2,
                                             upper = radius_range[0] + rnge * .8),
                               self.Success(lower = radius_range[0] + rnge * .8,
                                            upper = radius_range[0] + rnge)) # split 20, 60, 20
            
            self.frozen = False
        
        def sample(self):
            # if frozen, just return FAIL.lower since it should be 1
            if self.frozen:
                return self.categories[0].lower
            
            n = np.random.uniform()
            # check if in Fail
            if (n - self.category_probs[0] < 0):
                return self.categories[0].sample() # sample from FAIL
            elif (n - self.category_probs[0] - self.category_probs[1] < 0):
                return self.categories[1].sample() # sample from Moderate
            else:
                return self.categories[2].sample()
    
    
    def __init__(self, stenosis_parametrization):
        self.category_probs = (.3, .4, .3)
        self.stenosis_params = stenosis_parametrization
        self.stenosis_points = [self.StenosisRandom(category_probs=self.category_probs,
                                                  radius_range=sten_point) for sten_point in self.stenosis_params]
    
    def __len__(self):
        return len(self.stenosis_points)
        
    def fixed(self, sten_point, r_increase):
        ''' change a particular vessel to fixed '''
        self.stenosis_points[sten_point] = self.StenosisFixed(r_increase)
    
    def undo_fixed(self, sten_point):
        ''' only for testing purposed... to undo a fixture.'''    
        self.stenosis_points[sten_point] = self.StenosisRandom(category_probs=self.category_probs,
                                                                radius_range=self.stenosis_params[sten_point])
        
    def freeze(self, to_freeze = []):
        ''' freeze particular vessels (not fixed, just not changing). If a vessel is already frozen, do nothing'''
        for sten_idx in to_freeze:
            if type(self.stenosis_points[sten_idx]) == self.StenosisRandom:
                self.stenosis_points[sten_idx].frozen = True
    
    def unfreeze_all(self):
        ''' unfreeze all vessels that aren't fixed'''
        for sten_point in self.stenosis_points:
            if type(sten_point) == self.StenosisRandom:
                sten_point.frozen = False

    def unfreeze(self, to_unfreeze = []):
        ''' unfreezes previously frozen vessels. If a vessel is already unfrozen, do nothing'''
        for sten_idx in to_unfreeze:
            if type(self.stenosis_points[sten_idx]) == self.StenosisRandom:
                self.stenosis_points[sten_idx].frozen = False
        
    def sample(self):
        ''' generates a sample '''
        return np.array([sten_point.sample() for sten_point in self.stenosis_points])
    
    def create_dataset(self, num_samples):
        input_data = np.array([self.sample() for i in range(num_samples)])
        return PredictDataset(input_data)
    
    def get_histograms(self, results):
        return [(np.histogram(results[:, 2*i], bins = 'auto'), np.histogram(results[:, 2*i+1], bins = 'auto')) for i in range(len(self))]
            
    def plot_histograms(self, results, save):
        
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
        
def parametrize_stenosis(occlusion):
    ''' takes the artificial occlusion and computes the range of diameter changes for that particular vessel'''
    cur_area = 1 - occlusion
    area_increase = 1/cur_area
    radial_increase = np.sqrt(area_increase)
    # radial change goes from 1 (no repair) to (1/(1-occlusion))^(1/2) (complete repair)
    return np.dstack((np.ones_like(radial_increase), radial_increase)).squeeze()

class LightningNNPredictor(LightningNN):
    
    def __init__(self, model, lr, revert_map):
        super().__init__(model, lr, revert_map)
        
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x = batch
        y_hat = self.model(x)
        revert(y_hat, self.revert_map)
        return y_hat

if __name__ == '__main__':
    
    # path
    dir = Path('data/healthy/0080_0001/jc_solver_dir_0/artificial_stenosis/Manual_1')

    # load training data to recompute normalization map (#TODO: SAVE THE NORM NEXT TIME.)
    train_dataset = Dataset0D(dir / 'data' / 'training_data_small' / 'input.npy', dir / 'data' / 'training_data_small' / 'output.npy', normalization, revert_map=None)

    # Construct model
    # retrieve first value of Dataset for sizes
    input_data, output_data = train_dataset[0]
    nnmodel = BasicNN(input_neurons=len(input_data), output_neurons=len(output_data), hidden_layers=3, neurons_per_layer=1000)
    litmodel = LightningNNPredictor(nnmodel, lr = 1e-3, revert_map = train_dataset.revert_map)
    
    
    # create probability distribution
    stenosis_file = read_json(dir / 'stenosis_vessels.dat')
    stenosis_parametrizations = parametrize_stenosis(np.array(stenosis_file['occlusions']))
    distribution = StenosisDistribution(stenosis_parametrizations)
    
    
    best_checkpoint = dir / "training_results" / "run1" / "lightning_logs" / "version_11" / "checkpoints" / "epoch=147-step=37888.ckpt"
    trainer = pl.Trainer(logger = False, accelerator='gpu')
    
    '''
    # starting at random, no fixture
    x, yhat = distribution.run_test(model=litmodel,
                          trainer=trainer,
                          best_chkpt=best_checkpoint,
                          num_samples=4096*24,
                          batch_size=4096)
    distribution.plot_histograms(yhat, save = 'pptx/pptx_img/prob.random.gif')
    
    
    # fix 3, keep 1 random, freeze rest
    distribution.fixed(0, 1.5)
    distribution.fixed(1, 2.1)
    distribution.fixed(2, 1.9)
    # freeze all but 3
    distribution.freeze(to_freeze=list(range(4, len(distribution))))
    x, yhat = distribution.run_test(model=litmodel,
                          trainer=trainer,
                          best_chkpt=best_checkpoint,
                          num_samples=4096*24,
                          batch_size=4096)
    distribution.plot_histograms(yhat, save = 'pptx/pptx_img/prob.f3r1rf17.gif')
    
    '''
    # freeze all but 3
    distribution.freeze(to_freeze=list(range(1, len(distribution))))
    x, yhat = distribution.run_test(model=litmodel,
                          trainer=trainer,
                          best_chkpt=best_checkpoint,
                          num_samples=4096*24,
                          batch_size=4096)
    distribution.plot_histograms(yhat, save = 'pptx/pptx_img/prob.test.gif')