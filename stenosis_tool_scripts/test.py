import sys
print(sys.path)


from utils.data_org import DataPath


x = DataPath(root = '.')
print(x)