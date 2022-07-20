import sys
print(sys.path)


from src.data_org import DataPath


x = DataPath(root = '.')
print(x)