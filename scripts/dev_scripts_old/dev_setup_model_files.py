
from src.parser import DevParser
from src.data_org import DataPath

if __name__ == '__main__':
    dev_parser = DevParser(desc='Sets up necessary files. And checks for missing information')
    args = dev_parser.parse_args()

    org = DataPath(args.root)
    