from baseManager import Manager


class svManager(Manager):
    
    def __init__(self, yaml_file = None):
        super().__init__(yaml_file)
    
    def register_model(self, model_name: str ):
        pass

    def register_simulation(self):
        pass
    
    def register_rom(self):
        pass

    def register_mesh(self):
        pass

    def register_centerlines(self):
        pass

    def register_flows(self):
        pass
    