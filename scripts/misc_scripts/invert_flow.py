from svinterface.core.bc import Inflow
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Invert flow")
    parser.add_argument("-i", help = 'inflow waveform')
    parser.add_argument("-o", help = "outflow file of inverted inflow")
    
    args = parser.parse_args()
    
    i = Inflow.from_file(args.i, inverse = True, smooth = False)
    i.write_flow(args.o)
    print(f"Inverted inflow written to {args.o}")