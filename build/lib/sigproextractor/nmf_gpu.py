"""
Implementation of non-negative matrix factorization for GPU
"""

from datetime import datetime

from nimfa.methods.seeding import nndsvd
import numpy as np
import torch
import torch.nn
from torch import nn


class NMF:
    def __init__(self, V, rank, max_iterations=100000, tolerance=1e-8, test_conv=1000, gpu_id=0, seed=None,
                 init_method='random', floating_point_precision='double', min_iterations=2000):

        """
        Run non-negative matrix factorisation using GPU. Uses beta-divergence.

        Args:
          V: Matrix to be factorised
          rank: (int) number of latent dimensnions to use in factorisation
          max_iterations: (int) Maximum number of update iterations to use during fitting
          tolerance: tolerance to use in convergence tests. Lower numbers give longer times to convergence
          test_conv: (int) How often to test for convergnce
          gpu_id: (int) Which GPU device to use
          seed: random seed, if None (default) datetime is used
          init_method: how to initialise basis and coefficient matrices, options are:
            - random (will always be the same if seed != None)
            - NNDSVD
            - NNDSVDa (fill in the zero elements with the average),
            - NNDSVDar (fill in the zero elements with random values in the space [0:average/100]).
          floating_point_precision: (string or type). Can be `double`, `float` or any type/string which
              torch can interpret.
          min_iterations: the minimum number of iterations to execute before termination. Useful when using
              fp32 tensors as convergence can happen too early.
        """
        torch.cuda.set_device(gpu_id)

        if seed is None:
            seed = datetime.now().timestamp()

        if floating_point_precision == 'float':
            self._tensor_type = torch.FloatTensor
        elif floating_point_precision == 'double':
            self._tensor_type = torch.DoubleTensor
        else:
            self._tensor_type = floating_point_precision

        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)

        self.max_iterations = max_iterations
        self.min_iterations = min_iterations

        self._V = V.type(self._tensor_type).cuda()
        self._fix_neg = nn.Threshold(0., 1e-8)
        self._tolerance = tolerance
        self._prev_loss = None
        self._iter = 0
        self._test_conv = test_conv
        self._gpu_id = gpu_id
        self._rank = rank
        self._W, self._H = self._initialise_wh(init_method)

    def _initialise_wh(self, init_method):
        """
        Initialise baseis and coefficient matrices according to `init_method`
        """
        if init_method == 'random':
            W = torch.rand(self._V.shape[0], self._rank).type(self._tensor_type).cuda(self._gpu_id)
            H = torch.rand(self._rank, self._V.shape[1]).type(self._tensor_type).cuda(self._gpu_id)
            return W, H

        elif init_method == 'NNDSVD':
            nv = nndsvd.Nndsvd()
            vin = np.mat(self._V.cpu().numpy())
            W, H = nv.initialize(vin, self._rank, options={'flag': 0})

        elif init_method == 'NNDSVDa':
            nv = nndsvd.Nndsvd()
            vin = np.mat(self._V.cpu().numpy())
            W, H = nv.initialize(vin, self._rank, options={'flag': 1})

        elif init_method == 'NNDSVDar':
            nv = nndsvd.Nndsvd()
            vin = np.mat(self._V.cpu().numpy())
            W, H = nv.initialize(vin, self._rank, options={'flag': 2})

        W = torch.from_numpy(W).type(self._tensor_type).cuda(self._gpu_id)
        H = torch.from_numpy(H).type(self._tensor_type).cuda(self._gpu_id)
        return W, H

    @property
    def reconstruction(self):
        return self.W @ self.H

    @property
    def W(self):
        return self._W

    @property
    def H(self):
        return self._H

    @property
    def _kl_loss(self):
        return (self._V * (self._V / self.reconstruction).log()).sum() - self._V.sum() + self.reconstruction.sum()

    @property
    def _loss_converged(self):
        """
        Check if loss has converged
        """
        if not self._iter:
            self._loss_init = self._kl_loss
        elif ((self._prev_loss - self._kl_loss) / self._loss_init) < self._tolerance:
            return True
        self._prev_loss = self._kl_loss
        return False

    def fit(self, beta=1):
        """
        Fit the basis (W) and coefficient (H) matrices to the input matrix (V) using multiplicative updates and
            beta divergence
        Args:
          beta: value to use for generalised beta divergence. Default is 1 for KL divergence
            beta == 2 => Euclidean updates
            beta == 1 => Generalised Kullback-Leibler updates
            beta == 0 => Itakura-Saito updates
        """
        with torch.no_grad():
            def stop_iterations():
                stop = (self._iter % self._test_conv == 0) and self._loss_converged and (self._iter > self.min_iterations)
                if stop:
                    print("loss converged with {} iterations".format(self._iter))
                return stop
    
            if beta == 2:
                for self._iter in range(self.max_iterations):
                    self.H = self.H * (self.W.transpose(0, 1) @ self._V) / (self.W.transpose(0, 1) @ (self.W @ self.H))
                    self.W = self.W * (self._V @ self.H.transpose(0, 1)) / (self.W @ (self.H @ self.H.transpose(0, 1)))
                    if stop_iterations():
                        break
    
            # Optimisations for the (common) beta=1 (KL) case.
            elif beta == 1:
                ones = torch.ones(self._V.shape).type(self._tensor_type).cuda(self._gpu_id)
                for self._iter in range(self.max_iterations):
                    ht = self.H.t()
                    numerator = (self._V / (self.W @ self.H)) @ ht
                    denomenator = ones @ ht
                    self._W *= numerator / denomenator
    
                    wt = self.W.t()
                    numerator = wt @ (self._V / (self.W @ self.H))
                    denomenator = wt @ ones
                    self._H *= numerator / denomenator
                    if stop_iterations():
                        break
    
            else:
                for self._iter in range(self.max_iterations):
                    self.H = self.H * ((self.W.t() @ (((self.W @ self.H) ** (beta - 2)) * self._V)) /
                                       (self.W.transpose(0, 1) @ ((self.W @ self.H)**(beta-1))))
                    self.W = self.W * ((((self.W@self.H)**(beta-2) * self._V) @ self.H.transpose(0, 1)) /
                                       (((self.W @ self.H) ** (beta - 1)) @ self.H.transpose(0, 1)))
                    if stop_iterations():
                        break
