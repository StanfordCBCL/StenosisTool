
import matplotlib.pyplot as plt
import torch
from pathlib import Path
import torch.utils.data as tdata
import numpy as np

from svinterface.plotting.params import set_params


    
class Dataset0D(tdata.Dataset):
    """Dataset for input and outputs"""
    
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



if __name__ == '__main__':
    set_params()
    
    dir = Path('data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR')
    
    # load train and test dataset to get truths
    train_dataset = Dataset0D(dir / 'model_data' / 'train_data' / 'input.npy', dir / 'model_data' / 'train_data' / 'output.npy', normalization, revert_map=None)
    test_dataset = Dataset0D(dir / 'model_data' / 'test_data' / 'input.npy', dir / 'model_data' / 'test_data' / 'output.npy', normalization, revert_map=train_dataset.revert_map)

    # load 
    version = 'version_0'
    x = torch.load(dir / "training_results" /  "run1" / "lightning_logs" / version / "predict_input.pt")
    yhat = torch.load(dir / "training_results" /  "run1" / "lightning_logs" / version / "predict_output.pt")
    
    print(test_dataset.output.shape)
    print(yhat.shape())