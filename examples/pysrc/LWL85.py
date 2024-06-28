from sys import path
path.append("/mnt/ngs004/wayne/Project/VB/dNdS/kakscalculator2/examples/pysrc/")
from NG86 import Base

class LWL85(Base):
    def __init__(self):
        super().__init__()
        self.gamma = 0.0
        self.Pi = []
        self.Qi = []
        self.A = []
        self.B = []
        self.Si_temp = []
        self.Vi_temp = []
        self.K = []
    
    def run(self):
        pass
    
    def pre_process(self):
        pass
    
    def count_site_and_diff(self):
        pass
    
    def get_codon_class(self):
        pass
    
    def transition_transversion(self):
        pass
    


if __name__ == "__main__":
    lwl = LWL85()
       