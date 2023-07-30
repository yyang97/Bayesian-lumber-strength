import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
@jax.jit



def sigmoid(x, a):
  # x = jnp.array(x)
  # a = jnp.array(a)
  return 0.5 * (jnp.tanh(x * a / 2) + 1)

def dmgmodel_py(y,alpha,l,c,s):
  y = jnp.array(y)
  return(y*sigmoid(c*y-l,s) + alpha*y*sigmoid(l-c*y,s))

def dmgmodel_root_py(y,alpha,l,c,s,ystar):
  return(dmgmodel_py(y,alpha,l,c,s) - ystar)



def dmginverse_py(ystar,alpha,l,c,s):
  ystar = jnp.array(ystar)
  bisec = jaxopt.Bisection(
    optimality_fun=dmgmodel_root_py,
    lower = 0,
    upper = 10000,
    check_bracket = False)
  return(bisec.run(alpha = alpha,l = l, c= c ,s = s,ystar = ystar).params)

def dmginvgrad_py(ystar,alpha,l,c,s):
  grad_func = jax.grad(dmginverse_py,0)
  return(jnp.abs(grad_func(ystar,alpha,l,c,s)))

def dmglik_py(ystar,alpha,l,c,s,mu,sigma):
  y =  dmginverse_py(ystar,alpha,l,c,s)
  return(jnp.log(jax.scipy.stats.norm.pdf(y,loc = mu,scale = sigma)) + 
jnp.log(dmginvgrad_py(ystar,alpha,l,c,s))
)


def dmglik_vmap(y_group,alpha,l,c,s,mu,sigma):
  y_group = jnp.array(y_group)
  lik = jax.vmap(lambda y_group: dmglik_py(ystar = y_group,
alpha = alpha,l = l, c= c,s =s,mu = mu, sigma=sigma))(y_group)
  return(jnp.sum(lik))

def dmglik_jax(y_group,alpha,l,c,s,mu,sigma):
  return(np.array(dmglik_vmap(y_group,alpha,l,c,s,mu,sigma)))

def negdmglik_jax(theta,y_obs_g1,y_obs_g2,y_obs_g3,l,s):
  mu = theta[0]
  sigma = theta[1]
  alpha = theta[2]
  c = theta[3]
  y_obs_g1 = jnp.array(y_obs_g1)
  y_obs_g3 = jnp.array(y_obs_g3)

  #lik1 = jnp.sum(jnp.log(jax.vmap(lambda y: jax.scipy.stats.norm.pdf(y,loc = mu, scale = sigma))(y_obs_g1)))
  lik1 = jnp.sum(jax.scipy.stats.norm.logpdf(y_obs_g1,loc = mu, scale = sigma))
  lik2 = y_obs_g2*jnp.log(
    jax.scipy.stats.norm.cdf(dmginverse_py(l,alpha,l,c,s), loc=mu, scale=sigma) - 
      jax.scipy.stats.norm.cdf(l, loc=mu, scale=sigma)
  )
  lik3 = dmglik_vmap(y_group = y_obs_g3,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)

  return(-lik1 - lik2-lik3)



def negdmglik_np(theta,y_obs_g1,y_obs_g2,y_obs_g3,l,s):
  return(np.array(negdmglik_jax(theta,y_obs_g1,y_obs_g2,y_obs_g3,l,s)))
  
  
#   return(jax.vmap(dmglik_py,
# in_axes = (0,None,None,None,None,None,None))(
#   y_group,alpha,l,c,s,mu,sigma
# ))
# def gktrans_py(z,theta):
#     """The transformation function to convert the std normal z ~ N(0,1)
#       to x, which follows the g-and-k distribution of the parameter theta = (A,B,g,k)
# 
#     Args:
#         z: The basic r.v. from the std normal N(0,1)
#         theta: the parameter of the gk dist theta = (A,B,g,k)
# 
#     Returns:
#        The r.v. from the gk dist of the parameter theta = (A,B,g,k)
# 
#     """
#     A = theta[0]
#     B = theta[1]
#     g = theta[2]
#     k = theta[3]
#     z = jnp.array(z)
#     expterm = jnp.exp(-g*z)
#     return (A + B*(1 + .8*(1 - expterm)/(1 + expterm))*(1 + z**2)**k*z)
# 
# def Gfun_jax(yobs_sumstat, z, z_quant, theta):
#     """The G matrix (moment condition) of the logel, based on the jax
# 
#     Args:
#         yobs_sumstat: the vector of size `4`, the summary statistics observations from the g-k dist:
#           [mean, 25% quantile, 50% quantile,75% quantile]
#         z: the matrix of size `N x m`, the basic r.v. from the std normal N(0,1)
#         z_quant: the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#         theta: the parameter of the gk dist theta = (A,B,g,k)
# 
#     Returns:
#        the G matrix of size `m x 4`
#        G[,0] the mean summary stat, is mean(y_obs) - colmean(simulated g-k r.v.)
#        G[,1] the 25% quantile summary stat, is ...
#        G[,2] the 50% quantile summary stat, is ...
#        G[,3] the 75% quantile summary stat, is ...
# 
#     """
#     # To convert it into the jnp.array so that it can be used by jax
#     yobs_sumstat = jnp.array(yobs_sumstat)
#     x = jnp.array(gktrans_py(z, theta))
#     # convert the quantiles of z to the quantiles of x
#     x_quant = jnp.array(gktrans_py(z_quant, theta))
#     # The mean summary statistics
#     mean_summary = jnp.array(jnp.mean(x,axis = 0) - yobs_sumstat[0])
#     # The 25% qunatile summary statistics
#     quant25_summary = jnp.array(x_quant[:,0] - yobs_sumstat[1])
#     # The 50% quantile summary statistics
#     quant50_summary = jnp.array(x_quant[:,1] - yobs_sumstat[2])
#     # THe 75% quantile summary statistics
#     quant75_summary = jnp.array(x_quant[:,2] - yobs_sumstat[3])
#   
#     # The inversed G matrix
#     G_inv = jnp.array(
#       [mean_summary,quant25_summary,quant50_summary,quant75_summary]
#     )
#     # The G matrix
#     G = jnp.transpose(G_inv)
#     return (G)
# 
# def Gfun(yobs_sumstat, z, z_quant, theta):
#     """
#     The wrapper function, to tranform the G matrix from python to R
#     """
#     return(np.array(Gfun_jax(yobs_sumstat, z, z_quant, theta)))
# 
# def Gfun_grad_jax(yobs_sumstat, z, z_quant, theta):
#     """
#     The gradient of G matrix (moment condition) w.r.t. theta, i.e., dG/dtheta, based on the jax
# 
#     Args:
#         yobs_sumstat: the vector of size `4`, the summary statistics observations from the g-k dist:
#           [mean, 25% quantile, 50% quantile,75% quantile]
#         z: the matrix of size `N x m`, the basic r.v. from the std normal N(0,1)
#         z_quant: the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#         theta: the parameter of the gk dist theta = (A,B,g,k)
# 
#     Returns:
#        the gradient matrix of size `m x 4 x 4`
# 
#     """
#     # To convert it into the jnp.array so that it can be used by jax
#     yobs_sumstat = jnp.array(yobs_sumstat)
#     z = jnp.array(z)
#     theta = jnp.array(theta)
#   
#     # The gradient function dG/dtheta. Set it for 3 so that it takes gradient w.r.t. the index of 3
#     G_grad_func = jax.jacobian(Gfun_jax,3)
#     G_grad = G_grad_func(yobs_sumstat, z, z_quant, theta)
#     return G_grad
# 
# 
# 
# def Gfun_grad(yobs_sumstat, z, z_quant, theta):
#     """
#     The wrapper function, to tranform the gradient from python to R
#     """
#     return np.array(Gfun_grad_jax(yobs_sumstat, z, z_quant, theta))
# 
