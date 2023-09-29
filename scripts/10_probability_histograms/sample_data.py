# File: comp.py
# File Created: Thursday, 22nd September 2022 7:45:48 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 25th September 2023 3:18:29 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
#! Description: compares final output.

from torch import nn
from torch import optim
import torch
import re
import tqdm
import torch.utils.data as tdata
import pytorch_lightning as pl
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import argparse

from svinterface.utils.io import read_json
from svinterface.utils.misc import d2m
from train_nn import *

class PredictDataset(tdata.Dataset):
    def __init__(self, inputs):
        self.inputs = inputs
        
    def __len__(self):
        ''' currently undefined '''
        return len(self.inputs)
    
    def __getitem__(self, idx):
        return torch.from_numpy((self.inputs[idx])).float()
    
def revert(yhat, revert_map):
    for i in range(len(yhat[0])):
        yhat[:, i] = (yhat[:, i] * revert_map[i][1]) + revert_map[i][0]
    return yhat

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
    
    def get_baseline(self, model, trainer, best_chkpt):
        ''' retrieves baseline sample (freezing all unfixed and unfrozen)'''
        to_freeze = []
        for idx, repair_point in enumerate(self.repair_points):
            if type(repair_point) == self.RepairRandom and repair_point.frozen == False:
                to_freeze.append(idx)
        self.freeze(to_freeze)
        # run base test
        x, yhat = self.run_test(model, trainer, best_chkpt,  num_samples = 1, batch_size = 1, p_var = 0, q_var = 0)
        self.unfreeze(to_unfreeze=to_freeze)
        return x, yhat
    
    def test_single(self, model, trainer, best_chkpt, c):
        data = PredictDataset(np.array([c, c, c, c]))
        test_loader = tdata.DataLoader(data, batch_size = 2, shuffle = False)
        # predict on the test loader and get normalized results
        rez = trainer.predict(model=model, dataloaders=test_loader, ckpt_path=best_chkpt, return_predictions=True)
        # retrieve x
        x = []
        for i in range(len(test_loader.dataset)):
            x.append(test_loader.dataset[i])
        x = torch.vstack(x)
        rez = torch.vstack(rez)
        return x, rez
    
    def get_histograms(self, yhat, points = 'all'):
        if points == 'all':
            points = range(0, len(yhat[0]),6)
        else:
            points = np.array(points) * 6
        return [[np.histogram(yhat[:, i + j], bins = 'auto', density=True) for j in range(3)]+[np.histogram(d2m(yhat[:, i + j]), bins = 'auto', density=True) for j in range(3,6)] for i in points]

    def to_mmHg(self, yhat):
        base_idx = np.arange(0, len(yhat[0]), 6)
        pressures_idx = np.concatenate([base_idx + 3, base_idx + 4, base_idx + 5])
        yhat[:, pressures_idx] = d2m(yhat[:, pressures_idx])
        return yhat
        
    
    def save_data(self, x, yhat, baseline, filepath):
        np.save(filepath, {'x': x,
                           'yhat': yhat,
                           'baseline':baseline},
                allow_pickle=True)
        
    def save_histograms(self, histograms, filepath):
        np.save(filepath, histograms, allow_pickle=True)
            
    def plot_histograms(self, histograms, baseline =None, path: Path= None):
        counter = 0
        for i in tqdm.tqdm(range(len(histograms)), desc='Generating Histograms'):
            fig = self.plot_single_histogram(histograms[i], baseline[counter:counter+6])
            fig.suptitle(f"Histograms at point {i +1}")
            fig.savefig(str(path / f'hist_point_{i+1}.png'))
            plt.close(fig)
            counter += 6

    def plot_single_histogram(self, h, baseline_yhat, p=True, q=True,):
        if p and q:
            fig, ax = plt.subplots(2, 3, figsize = (10, 5) )
            ax = ax.flatten()
            names = ['Diastolic Flow (mL/s)', 'Mean Flow (mL/s)', 'Systolic Flow (mL/s)', 'Diastolic Pressure (mmHg)','Mean Pressure (mmHg)','Systolic Pressure (mmHg)']
            for j in range(6):
                ax[j].bar(h[j][1][:-1], h[j][0], width=np.diff(h[j][1]),edgecolor="black", linewidth=.1, align="edge")
                ax[j].axvline(x=baseline_yhat[j], color = 'r', label = 'Baseline')
                ax[j].set_xlabel(names[j])
                ax[j].set_ylabel("Density")
            fig.legend(labels=['Baseline'], fontsize=plt.rcParams['font.size']-2)
            fig.tight_layout()
            return fig
        elif p:
            fig, ax = plt.subplots(1, 3, figsize = (10, 2.5) )
            ax = ax.flatten()
            names = ['Diastolic Pressure (mmHg)','Mean Pressure (mmHg)','Systolic Pressure (mmHg)']
            for j in range(3, 3+len(names)):
                ax[j-3].bar(h[j][1][:-1], h[j][0], width=np.diff(h[j][1]),edgecolor="black", linewidth=.1, align="edge")
                ax[j-3].axvline(x=baseline_yhat[j], color = 'r', label = 'Baseline')
                ax[j-3].set_xlabel(names[j-3])
                ax[j-3].set_ylabel("Density")
            fig.legend(labels=['Baseline'], fontsize=plt.rcParams['font.size']-2)
            fig.tight_layout()
            return fig
        elif q:
            fig, ax = plt.subplots(1, 3, figsize = (10, 2.5) )
            ax = ax.flatten()
            names = ['Diastolic Flow (mL/s)', 'Mean Flow (mL/s)', 'Systolic Flow (mL/s)']
            for j in range(len(names)):
                ax[j].bar(h[j][1][:-1], h[j][0], width=np.diff(h[j][1]),edgecolor="black",linewidth=.1, align="edge")
                ax[j].axvline(x=baseline_yhat[j], color = 'r', label = 'Baseline')
                ax[j].set_xlabel(names[j])
                ax[j].set_ylabel("Density")
            fig.legend(labels=['Baseline'], fontsize=plt.rcParams['font.size']-2)
            fig.tight_layout()
            return fig
    
    def add_uncertainty(self, yhat, p_var, q_var):
        base = np.array(range(0, len(yhat[0]), 6))
        pressures = np.concatenate([base + 3, base + 4, base + 5])
        flows = np.concatenate([base, base + 1, base + 2])
        yhat[:, pressures] += torch.from_numpy(np.random.normal(0, p_var, size=(len(yhat),1))).tile(len(pressures)).float()
        yhat[:, flows] +=  torch.from_numpy(np.random.normal(0, q_var, size=(len(yhat),1))).tile(len(flows)).float()
        return yhat
        
    def run_test(self, model, trainer, best_chkpt,  num_samples = 4096, batch_size = 4096, p_var = 1, q_var = 1):
        test_dataset = self.create_dataset(num_samples)
        test_loader = tdata.DataLoader(test_dataset, batch_size = batch_size, shuffle = False)
        # predict on the test loader and get normalized results
        rez = trainer.predict(model=model, dataloaders=test_loader, ckpt_path=best_chkpt, return_predictions=True)
        # retrieve x
        x = []
        for i in range(len(test_loader.dataset)):
            x.append(test_loader.dataset[i])
        x = torch.vstack(x)
        rez = torch.vstack(rez)
        rez = self.add_uncertainty(rez, p_var, q_var)
        return x, rez

class LightningNNPredictor(LightningNN):
    
    def __init__(self, model, lr, revert_map):
        super().__init__(model, lr, revert_map)
        
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x = batch
        y_hat = self.model(x)
        revert(y_hat, self.revert_map)
        return y_hat
    
def get_checkpoint(ckpt_dir: Path):
    for file in ckpt_dir.iterdir():
        return file
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Observe the results")
    parser.add_argument("-train_dir", default = 'data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR/training_results/run_32768', help = 'training data directory to use' )
    args = parser.parse_args()
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    dir = Path(args.train_dir)
    # load normalization map
    revert_map = torch.load(dir / 'lightning_logs' / 'version_0' / 'revert_map.pt', map_location=torch.device(device))
    
    # get checkpoint and model size
    best_ckpt = get_checkpoint(dir  / 'lightning_logs' / 'version_0' / 'checkpoints')
    ckpt = torch.load(str(best_ckpt), map_location=torch.device(device))
    input_size = ckpt['state_dict']['model.input_layer.weight'].size()[1]
    output_size = ckpt['state_dict']['model.output_layer.weight'].size()[0]

    # Construct model
    nnmodel = BasicNN(input_neurons=input_size, output_neurons=output_size, hidden_layers=3, neurons_per_layer=1000) 
    litmodel = LightningNNPredictor(nnmodel, lr = 1e-3, revert_map = revert_map)
    
    # trainer
    trainer = pl.Trainer(logger = False, accelerator='auto')
    
    # create a dir to save data
    prob_dir = dir / 'probability_histograms'
    prob_dir.mkdir(exist_ok=True)
    
#     distribution = RepairDistribution(num_repairs=input_size, 
#                                         category_probs=(.3,.4,.3))
    
#     x, yhat = distribution.test_single(model=litmodel,
#                             trainer=trainer,
#                             best_chkpt=str(best_ckpt),
#                             c=[0.7663, 0.2433, 0.9728]
# )
#     print(x)
#     print(yhat)
    # starting at random, no fixture (All possible distributions)
    
    all_dir = prob_dir / 'all_distributions'
    all_dir.mkdir(exist_ok=True)
    
    for std in (0, 1, 3):
        all_dir_std = all_dir / str(std)
        all_dir_std.mkdir(exist_ok=True)
        # create probability distribution
        distribution = RepairDistribution(num_repairs=input_size, 
                                        category_probs=(.3,.4,.3))
        
        x, yhat = distribution.run_test(model=litmodel,
                            trainer=trainer,
                            best_chkpt=str(best_ckpt),
                            num_samples=4096*24,
                            batch_size=4096,
                            p_var=1333.22 * std,
                            q_var=std)
        hist = distribution.get_histograms(yhat)
        
        baseline = distribution.get_baseline(model=litmodel,
                                            trainer=trainer,
                                            best_chkpt=str(best_ckpt))
        distribution.save_data(x, yhat, baseline, all_dir_std / 'data.npy')
        hist_dir = all_dir_std / 'histograms'
        hist_dir.mkdir(exist_ok=True)
        # only plot first 10
        points = 10
        distribution.plot_histograms(hist[:points],baseline[1][0][:6 * points], path = hist_dir)
        # distribution.save_histograms(hist, filepath=str(all_dir / 'histograms.npy') )
        
    
    
    # fix repair 0 (LPA_proximal) at 
    fixz_dir = prob_dir / 'fixeda3_conditional'
    fixz_dir.mkdir(exist_ok=True)
    
    for std in (0, 1, 3):
        fixz_dir_std = fixz_dir / str(std)
        fixz_dir_std.mkdir(exist_ok=True)
        
        distribution = RepairDistribution(num_repairs=input_size, 
                                        category_probs=(.3,.4,.3))
        distribution.fixed(repair_idx=2,
                        c=1)
        
        x, yhat = distribution.run_test(model=litmodel,
                            trainer=trainer,
                            best_chkpt=str(best_ckpt),
                            num_samples=4096*24,
                            batch_size=4096,
                            p_var=1333.22 * std,
                            q_var=std)
        hist = distribution.get_histograms(yhat)
        
        baseline = distribution.get_baseline(model=litmodel,
                                            trainer=trainer,
                                            best_chkpt=str(best_ckpt))
        distribution.save_data(x, yhat, baseline, fixz_dir_std / 'data.npy')
        hist_dir = fixz_dir_std / 'histograms'
        hist_dir.mkdir(exist_ok=True)
        # only plot first 10
        points = 10
        distribution.plot_histograms(hist[:points],baseline[1][0][:6 * points], path = hist_dir)
        # distribution.save_histograms(hist, filepath=str(fixz_dir / 'histograms.npy') )
    
