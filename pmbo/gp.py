#!/usr/bin/env python
# coding: utf-8

# In[1]:


import torch
import gpytorch
from gpytorch.kernels import ScaleKernel
from gpytorch.kernels import MaternKernel


# In[2]:


class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = ScaleKernel(
            MaternKernel(5/2, ard_num_dims=train_x.shape[1]))

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
def gp_model (train_x, train_y,iteration,error=1e-8):
    likelihood_optimization_criteria=error
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = ExactGPModel(train_x, train_y, likelihood)
    model.train()
    likelihood.train()
    # Use the adam optimizer
    optimizer_likelihood = torch.optim.Adam(model.parameters(), lr=0.1)  # Includes GaussianLikelihood parameters
    # "Loss" for GPs - the marginal log likelihood
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
    loss_prev = 0.1
    likelihood_optimization_iter_max=iteration
    for i in range(likelihood_optimization_iter_max):
        # Zero gradients from previous iteration
        optimizer_likelihood.zero_grad()
        # Output from model
        output =model(train_x)
        # Calc loss and backprop gradients
        loss = - mll(output, train_y)
        loss.backward()
        loss_residual = abs(loss.item() - loss_prev) / abs(loss_prev)
        loss_prev = loss.item()
        if (i +1)%1000==0:
            print('Iter %d/%d - Loss: %.3f  res: %.8f' % (
                i + 1, likelihood_optimization_iter_max,
                loss.item(),
                loss_residual
            ))
        if loss_residual < likelihood_optimization_criteria:
            print ("Step %d achieves convergence accuracy" % (i + 1))
            break
        optimizer_likelihood.step()
    print ("### complete ###")
    return model,likelihood
def gp_predict(model,likelihood,test_x):
    # Get into evaluation (predictive posterior) mode
    model.eval()
    likelihood.eval()
    # Test points are regularly spaced along [0,1]
    # Make predictions by feeding model through likelihood
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        observed_pred= likelihood(model(test_x))
    return observed_pred


# In[ ]:




