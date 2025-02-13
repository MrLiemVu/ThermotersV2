from sys import path as syspath
syspath.append("../")

import numpy as np
from utils.general_functions import multi_map, tensum, getDiNu
from utils.model_functions import getBricks

class ThermodynamicModel:
    def __init__(self, parameters):
        self.params = dict(parameters)

    def sequences2bricks(self, seqs, dinuCoordsAndValues=None):
        '''
        Converts a sequence of DNA to bricks for both forward and reverse complement strands.
        
        Parameters:
          seqs: a 2D numpy array of shape (n, L) where n is the number of sequences and L is the length of each sequence.
          dinuCoordsAndValues: a tuple of two lists. The first list contains the coordinates of 
                            dinucleotides in the sequence and the second list contains the values of the dinucleotides.
        
        Returns:
          bricks: a dictionary with keys "frw" and "rc" where the values are 2D numpy arrays
                  of shape (n, L) representing the bricks for the forward and reverse complement
                  strands respectively.
        '''
        # 1. Extract the dinucleotide coordinates and values
        if dinuCoordsAndValues is not None:
            dinuCoords, dinuValues = dinuCoordsAndValues
            
        # 2. Initialize the strands
        strands = ["frw"]
        if self.params["includeRC"]:
            strands += ["rc"]
        
        # 3. Create bricks for each strand
        bricks = {}
        for strand in strands:
            sq = seqs
            if strand=="rc":
                sq = 3 - sq[:, ::-1].copy(order="C")
            tmp = getBricks(
                self.params["matrices"],
                self.params["min.spacer"],
                self.params["sp.penalties"],
                sq,
                makeLengthConsistent=True).T
            tmp += -self.params["chem.pot"]

            if dinuCoordsAndValues is not None:
                global mp_getDiNu

                def mp_getDiNu(coord_):
                    return getDiNu(*coord_,
                                   n1=self.params["matrices"][0].shape[0],
                                   minSpacer=self.params["min.spacer"],
                                   n2=self.params["matrices"][1].shape[0],
                                   sequences=sq,
                                   nSpacer=len(self.params["sp.penalties"])).T

                tmpDn = np.array(multi_map(mp_getDiNu, dinuCoords, processes=14))
                #                 tmpDn = np.array([
                #                     getDiNu(*coord,
                #                         n1=mdl["matrices"][0].shape[0],
                #                         minSpacer=mdl["min.spacer"],
                #                         n2=mdl["matrices"][1].shape[0],
                #                         sequences=sq,
                #                         nSpacer=len(mdl["sp.penalties"])).T
                #                     for coord in dinuCoords])
                tmp += np.array(tensum(dinuValues, tmpDn))
            if strand=="rc":
                tmp = tmp[:, ::-1]
            bricks[strand] = tmp
        return bricks

    def bricks2pons(self, bricks):
        '''
        Converts bricks to Pons for a given set of bricks.
        
        Parameters:
            bricks: a dictionary with keys "frw" and "rc" where the values are 2D numpy arrays
                    of shape (n, L) representing the bricks for the forward and reverse complement
                    strands respectively.
        
        Returns:
            Pons: a 2D numpy array of shape (n, L) representing the Pons for the forward and reverse complement
                  strands respectively.
        '''
        # Check if the logClearanceRate parameter is present
        if "logClearanceRate" in self.params:
            R_ = np.exp(self.params["logClearanceRate"])
        else:
            R_ = 0
        
        # Extract the bricks for the forward strand
        bdni = bricks["frw"]
        thresholdPos = self.params["RBSthreshold"]
        
        # Check if the threshold position is less than or equal to 0
        if thresholdPos <= 0:
            thresholdPos = bdni.shape[1] + thresholdPos
        off = thresholdPos < bdni.shape[1]
        
        # Determine function to calculate the sum of the binding energies for the ON and OFF states
        bindMode_ = self.params["bindMode"]
        if bindMode_ == "add":
            bindsumF = lambda xi: np.sum(
                1. / (np.exp(xi) + R_),
                axis=tuple(range(1, xi.ndim))
            )
        elif bindMode_ == "max":
            bindsumF = lambda xi: np.exp(-np.min(xi, axis=tuple(range(1, xi.ndim))))
        
        # Calculate the sum of the binding energies for the ON state
        sumON_ = bindsumF(bdni[:, :thresholdPos])
        
        # Calculate the sum of the binding energies for the OFF state
        if off:
            sumOFF_ = bindsumF(bdni[:, thresholdPos:])
        else:
            sumOFF_ = 0.
        
        # Calculate the Pons
        if "rc" in bricks:
            bdni = bricks["rc"]
            rcOcclusion = self.params.get("rcOcclusion", np.arange(bdni.shape[1]))
            sumOFF_ += bindsumF(bdni[:, rcOcclusion])
        Pons_ = sumON_ / (1. + sumOFF_ + sumON_)
        
        return Pons_

    def __repr__(self):
        '''Returns string representation of model's parameters'''
        return self.params.__repr__()