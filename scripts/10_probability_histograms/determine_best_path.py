from sample_data import RepairDistribution, BasicNN, LightningNNPredictor, get_checkpoint

import pytorch_lightning as pl
import argparse
import torch
from pathlib import Path
import matplotlib.pyplot as plt

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
    
    
    distribution = RepairDistribution(num_repairs=3, 
                                        category_probs=(.3,.4,.3)) 
    
    distribution.fixed(2, 0)
    x, yhat = distribution.run_test(model=litmodel,
                        trainer=trainer,
                        best_chkpt=str(best_ckpt),
                        num_samples=4096*24,
                        batch_size=4096,
                        p_var=0,
                        q_var=0)
    hist = distribution.get_histograms(yhat)
    
    
    baseline = distribution.get_baseline(model=litmodel,
                                        trainer=trainer,
                                        best_chkpt=str(best_ckpt))
    baseline_yhat = distribution.to_mmHg(baseline[1])
    
    fig = distribution.plot_single_histogram(hist[0],baseline_yhat[0][:6], p = True, q = False)
    
    plt.show()
    