from svinterface.core.zerod.lpn import LPN



x = LPN.from_file("data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/LPN_DIR/AS1_SU0308.in")

tree = x.get_tree()

stack = []

def get_gens( node, gen = 0):
    
    stack.append(gen)
    for c in node.children:
        if type(c) == LPN.BranchNode:
            get_gens(c, gen=gen+1)
        else:
            get_gens(c, gen = gen)
            
get_gens(tree)
print(max(stack))