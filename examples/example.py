import pyenkf
import numpy as np
import matplotlib.pyplot as plt

# Example showing the use of pyenkf Ensemble Kalman square root filter
# the test function is linear in coefficients:
#       y = a + cos(t) + b*t + sin(t)
# the reader is supposed to be up-to-speed with the theory of the ensemble
# kalman square root filter

# define our test model
def test_model(t,a,b,c):
    return a*np.cos(t) +b*t+c*np.sin(t)

# true values of the coefficients
a_true=0.3
b_true=0.3
c_true=3

# independent variable
t=np.linspace(0,10,100)

# synthetic truth based on true values of the coefficients
truth = test_model(t,a_true,b_true,c_true)

# create observations with an uncertainty from the truth 
obs_uncertainty=0.1
obs=truth+np.random.normal(0,obs_uncertainty,np.shape(truth))

# create an outlier, which will be rejected based on a 3-sigma filter criterion
obs[3]+=4

# Uncorrelated observations, diagonal of covariance matrix
R=np.array([0.1]*len(obs))

# Create empty containers:
HPHR=np.zeros((len(obs),),float)
# List of integers. Will be returned by the filter.
# value of 0 means the observation is used in the fitting
# a value of 1 
rejected=np.zeros((len(obs),),int)

# Can we reject the observation?
may_reject=np.array([True]*len(obs))

# option to create a per-observation <number>-sigma criterion, 3 in the example below 
rejection_threshold=np.array([3]*len(obs))

# Define the prior guess 
a_prior=0.1
b_prior=0.1
c_prior=1.0

# prior uncertainty
a_prior_unc=0.5
b_prior_unc=0.5
c_prior_unc=0.5
x=np.array([a_prior,b_prior,c_prior])

# how many members to represent the ensemble?
ensemblesize=50

# deviations for each coefficient to be fitted 
X_prime=np.zeros((len(x),ensemblesize),float, order='F')

X_prime[0,1:]=np.random.normal(0,a_prior_unc,ensemblesize-1)
X_prime[1,1:]=np.random.normal(0,b_prior_unc,ensemblesize-1)
X_prime[2,1:]=np.random.normal(0,c_prior_unc,ensemblesize-1)

# Model results using prior mean
Hx=test_model(t,*x)

# Copy for plot, as Hx will be modified
prior_guess=Hx.copy()

# Model results for the deviations
HX_prime=np.zeros((len(obs),ensemblesize),float,order='F')
for member in range(ensemblesize):
    HX_prime[:,member]=test_model(t,*(x+X_prime[:,member]))-Hx

#
# Call the filter by parsing the data
# Note that above we defined 2 dimensional arrays as Fortran arrays which
# affects their memory layout. Not doing so results in errors
pyenkf.enkf_core.enksrf(obs, may_reject, rejection_threshold, Hx, HX_prime, X_prime, R, HPHR, x, rejected)

# create a figure that shows truth/observations/prior/posterior results
# and an annotation with the rejected sample. 
plt.plot(t,truth,label='truth %s %s %s' %(a_true,b_true,c_true),zorder=3,)
plt.plot(t,obs,marker='o',linestyle='', label='observations')

plt.annotate('outlier/rejected',xy=(t[3],obs[3]))
plt.plot(t,prior_guess,linestyle='-', label='prior guess %s %s %s' % (a_prior,b_prior,c_prior))
plt.plot(t,Hx,linestyle=':', label='posterior result a=%.2f b=%.2f c=%.2f' % (x[0],x[1],x[2]))

plt.xlabel('x')
plt.ylabel('y')

plt.legend()

plt.tight_layout()

plt.savefig('test.png')
