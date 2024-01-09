# Tools specific to post-processing moist simulations

class MoistParams:
    """
    Class containing input parameters for moist convection,
    reading the input from the file `humid.in`
    """
    def __init__(self, folder):
        fname = folder+"/humid.in"
        with open(fname, "r") as f:
            flist = f.readlines()
            self.alpha = float(flist[4])
            self.beta = float(flist[7])
            self.gamma = float(flist[10])
            self.tau = float(flist[13])
            self.qfixN, self.qfixS = [int(x) for x in flist[20].split()]
            self.Sm = float(flist[23])