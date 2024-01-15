---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

<!-- #region -->
# Damage model 
The original model is 

\begin{equation}\label{eq: original model}
y^*  = y \cdot I\{c \cdot y > l\} + \alpha \cdot y \cdot I\{c \cdot y < l\}
\end{equation}



# Smooth function

Note that we have the indicator function in the model we might consider the smooth function:
\begin{equation}\label{eq: trans function}
S(x;s) = \frac{1}{1 + \exp(-s \cdot x)}
\end{equation},
where $s$ is the smoothing hyper-parameter to control the smoothness. 
<!-- #endregion -->

Then the model becomes 
$$
y^*  = y \cdot \frac{1}{1 + \exp(-s\cdot(cy-l))} + \alpha \cdot y \cdot \frac{1}{1 + \exp(-s\cdot(l-cy))}
$$

<!-- #region -->
We have $y \geq y^*$.

The lumber have three groups:

- Group 1: $y <l$, $y^* < l$, i.e., $y^* <y < l$. The lumber pieces are broken blow the proof loading.
- Group 2: $y >l$, $y^* < l$, i.e. $y^*<l<y$. The lumber pieces are broken during the proof loading process. This groups we only knows how many pieces. 
- Group 3: $y >l$, $y^* > l$, i.e. $l<y^*<y$. The lumber pieces survived in the proof-loading. And then we destruct them to test their strength.


# The PDF calculation 

Given $Y \sim N(\mu, \sigma^2)$, $Y^* = h(Y)$. Then the pdf of $Y^*$,
$$
f_{Y^*}(y^*) = f_{Y}(h^{-1}(y^*))|\frac{d}{dy^*}h^{-1}(y^*)|,
$$
where $f_Y()$ is the pdf of $Y$, i.e., normal. 

Following this, We need the numerical function of $h^{-1}(y^*)$, and its numerical gradient $\frac{d}{dy^*}h^{-1}(y^*)$. (The analytical form doesn't seem available.)
<!-- #endregion -->

```python
import jax
import jaxopt
import jax.numpy as jnp
import pyreadr
import projplot as pjp

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
```

```python
@jax.jit


def logit(x):
    return(jnp.log(x/(1-x)))


def expit(x):
    return 1/(1+jnp.exp(-x))




def sigmoid(x, s):
    # x = jnp.array(x)
    # a = jnp.array(a)
    return 0.5 * (jnp.tanh(x * s / 2) + 1)

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
    return(jax.scipy.stats.norm.logpdf(y,loc = mu,scale = sigma)+ 
           jnp.log(dmginvgrad_py(ystar,alpha,l,c,s))
          )


def dmglik_vmap(y_group,alpha,l,c,s,mu,sigma):
    y_group = jnp.array(y_group)
    lik = jax.vmap(lambda y_group: dmglik_py(ystar = y_group,
                                             alpha = alpha,l = l, c= c,s =s,mu = mu, sigma=sigma))(y_group)
    return(jnp.sum(lik))


```

```python
# # # read data 

# y_obs_g1 = pyreadr.read_r('y_obs_g1.Rds') # also works for RData
# y_obs_g1 = jnp.array(y_obs_g1[None])
# len_g1 = jnp.shape(y_obs_g1)[0]
# y_obs_g1 = jnp.reshape(y_obs_g1,newshape = (len_g1,))

# y_obs_g2 = pyreadr.read_r('y_obs_g2.Rds') # also works for RData
# y_obs_g2 = jnp.array(y_obs_g2[None])
# y_obs_g2 = y_obs_g2[0][0]

# y_obs_g3 = pyreadr.read_r('y_obs_g3.Rds') # also works for RData
# y_obs_g3 = jnp.array(y_obs_g3[None])
# len_g3 = jnp.shape(y_obs_g3)[0]

# y_obs_g3 = jnp.reshape(y_obs_g3,newshape = (len_g3,))
# y_obs_g3





```

```python
# the orignal sample size 
#N = 30000
N = 300
mu = 48
sigma = 19
l =  32
alpha = 0.25
c = 0.65
s = 10

N_marginal = 139




```

```python
# data generation 
key = jax.random.PRNGKey(0)
subkeys = jax.random.split(key, num=N)

y = sigma*jax.random.normal(subkeys[10], shape=(N, )) + mu

#y = y[y>0]

# g1
y_obs_g1  = y[y<l]


# g23_star
g23 = y[y>l]
g23_star = jax.vmap(lambda y: dmgmodel_py(y,alpha,l,c,s))(g23)


# g3
y_obs_g3 = g23_star[g23_star > l]

# g2
y_obs_g2 = N - len(y_obs_g1) - len(y_obs_g3)

#y_obs_g2 = jnp.shape(y)[0] - len(y_obs_g1) - len(y_obs_g3)

# g4, is all marginal data
y_obs_g4 = sigma*jax.random.normal(subkeys[0], shape=(N_marginal, )) + mu
```

```python
# negdmglik_jax(theta0)
#dmglik_vmap(y_group = y_obs_g3,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)

#dmglik_py(y_obs_g3[0][0],alpha,l,c,s,mu,sigma)

#dmglik_vmap(y_obs_g3,alpha,l,c,s,mu,sigma)

# s= 1
l = l
y_obs_g1 = y_obs_g1
y_obs_g2 = y_obs_g2
y_obs_g3 = y_obs_g3
y_obs_g4 = y_obs_g4
@jax.jit
def negdmglik_jax(theta):
    mu = theta[0]
    sigma = theta[1]
    alpha = theta[2]
    c = theta[3]
    lik1 = jnp.sum(jax.scipy.stats.norm.logpdf(y_obs_g1,loc = mu, scale = sigma))
    #lik1 = dmglik_vmap(y_group = y_obs_g1,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
    lik2 = y_obs_g2*jnp.log(
        jax.scipy.stats.norm.cdf(dmginverse_py(l,alpha,l,c,s), loc=mu, scale=sigma) - 
        jax.scipy.stats.norm.cdf(l, loc=mu, scale=sigma)
    )
    lik3 = dmglik_vmap(y_group = y_obs_g3,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
    lik4 = jnp.sum(jax.scipy.stats.norm.logpdf(y_obs_g4,loc = mu, scale = sigma))

    return(-lik1 - lik2-lik3-lik4)

theta0 = jnp.array([mu,sigma,alpha,c])



negdmglik_jax(theta0)
```

```python
theta0 = jnp.array([mu,sigma,alpha,c])


solver = jaxopt.ScipyMinimize(method = "Nelder-Mead",fun=negdmglik_jax,maxiter = 1000)
res = solver.run(theta0)

# solver = jaxopt.BFGS(fun=negdmglik_jax)
# res = solver.run(theta0)

res.params,res.state
theta0
```

```python
xsol = res.params
xrange = .5
theta_upper = xsol  + xrange *jnp.abs(xsol)
theta_lower = xsol - xrange * jnp.abs(xsol)
theta_rng_trans = jnp.array([theta_lower, theta_upper])

#res.params
theta_rng = jnp.transpose(theta_rng_trans)
theta_rng
# plot at the true value theta0
pjp.proj_plot(
    fun=negdmglik_jax,
    x_opt=xsol,
    x_lims=theta_rng, 
    n_pts=100, 
    opt_vlines=True
)
```

```python

```

# real data

```python
import pandas
import numpy as np

realdata = pandas.read_csv("bending-pl.csv", header = None)
realdata[0] = realdata[0]/1000*6.895

l = 4500/1000*6.895
# R100id = (realdata[1] == "R100")
# R100 = jnp.array(realdata[realdata[1] == "R100"][0])


R100 = jnp.array(realdata[realdata[1] == "R100"][0])
R20 = jnp.array(realdata[realdata[1] == "R20"][0])
R20R100 = jnp.array(realdata[realdata[1] == "R20R100"][0])

g2 = 97



s= 2
l = l
R100 = R100
R20 = R20
R20R100 = R20R100
g2 = 97
@jax.jit
def negdmglik_jax(theta):
    mu = theta[0]
    sigma = theta[1]
    alpha = theta[2]
    #alpha = jnp.exp(theta[2])
    c = theta[3]
    #s = theta[4]
    #s = jnp.exp(theta[4])
    #lik1 = jnp.sum(jnp.log(jax.vmap(lambda y: jax.scipy.stats.norm.pdf(y,loc = mu, scale = sigma))(y_obs_g1)))
    lik1 = jnp.sum(jax.scipy.stats.norm.logpdf(R20,loc = mu, scale = sigma))
    #lik1 = dmglik_vmap(y_group = R20,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
    lik2 = g2*jnp.log(
        jax.scipy.stats.norm.cdf(dmginverse_py(l,alpha,l,c,s), loc=mu, scale=sigma) - 
        jax.scipy.stats.norm.cdf(l, loc=mu, scale=sigma)
    )
    lik3 = dmglik_vmap(y_group = R20R100,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
    lik4 = jnp.sum(jax.scipy.stats.norm.logpdf(R100,loc = mu, scale = sigma))

    return(-lik1 - lik2-lik3- lik4)

#theta0 = jnp.array([mu,sigma,jnp.log(alpha),c,jnp.log(s)])

theta0 = jnp.array([mu,sigma,alpha,c])

negdmglik_jax(theta0)
```

```python
# solver = jaxopt.LBFGS(fun=negdmglik_jax)
# res = solver.run(theta0)



solver = jaxopt.ScipyMinimize(method = "Nelder-Mead",fun=negdmglik_jax, maxiter =1000)
res = solver.run(theta0)

res.params
theta0

# Array(1678.6377, dtype=float32)
negdmglik_jax(res.params)
```

```python
xsol = res.params
xrange = .5
theta_upper = xsol  + xrange *jnp.abs(xsol)
theta_lower = xsol - xrange * jnp.abs(xsol)
theta_rng_trans = jnp.array([theta_lower, theta_upper])

#res.params
theta_rng = jnp.transpose(theta_rng_trans)
theta_rng
# plot at the true value theta0
pjp.proj_plot(
    fun=negdmglik_jax,
    x_opt=xsol,
    x_lims=theta_rng, 
    n_pts=100, 
    x_names = ["mu","sigma","alpha","c"],
    opt_vlines=True
)
```

## treat s as an estimated parameter

```python
import pandas
import numpy as np

realdata = pandas.read_csv("bending-pl.csv", header = None)
realdata[0] = realdata[0]/1000*6.895

l = 4500/1000*6.895
# R100id = (realdata[1] == "R100")
# R100 = jnp.array(realdata[realdata[1] == "R100"][0])


R100 = jnp.array(realdata[realdata[1] == "R100"][0])
R20 = jnp.array(realdata[realdata[1] == "R20"][0])
R20R100 = jnp.array(realdata[realdata[1] == "R20R100"][0])

g2 = 97



#s= 1
l = l
R100 = R100
R20 = R20
R20R100 = R20R100
g2 = 97
@jax.jit
def negdmglik_jax_s(theta):
    mu = theta[0]
    sigma = theta[1]
   # alpha = theta[2]
    alpha = jnp.exp(theta[2])
    c = theta[3]
    #s = jnp.exp(theta[4])
    s = theta[4]
    #lik1 = jnp.sum(jnp.log(jax.vmap(lambda y: jax.scipy.stats.norm.pdf(y,loc = mu, scale = sigma))(y_obs_g1)))
    lik1 = jnp.sum(jax.scipy.stats.norm.logpdf(R20,loc = mu, scale = sigma))
    lik2 = g2*jnp.log(
        jax.scipy.stats.norm.cdf(dmginverse_py(l,alpha,l,c,s), loc=mu, scale=sigma) - 
        jax.scipy.stats.norm.cdf(l, loc=mu, scale=sigma)
    )
    lik3 = dmglik_vmap(y_group = R20R100,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
    lik4 = jnp.sum(jax.scipy.stats.norm.logpdf(R100,loc = mu, scale = sigma))

    return(-lik1 - lik2-lik3- lik4)

theta0 = jnp.array([mu,sigma,jnp.log(alpha),c,s])

#theta0 = jnp.array([mu,sigma,alpha,c,s])

negdmglik_jax_s(theta0)
```

```python
# solver = jaxopt.LBFGS(fun=negdmglik_jax)
# res = solver.run(theta0)



solver = jaxopt.ScipyMinimize(method = "Nelder-Mead",fun=negdmglik_jax_s, maxiter =1000)
res = solver.run(theta0)

res.params
theta0

# Array(1678.6377, dtype=float32)
negdmglik_jax_s(res.params)
```

```python
xsol = res.params
xrange = .5
theta_upper = xsol  + xrange *jnp.abs(xsol)
theta_lower = xsol - xrange * jnp.abs(xsol)
theta_rng_trans = jnp.array([theta_lower, theta_upper])

#res.params
theta_rng = jnp.transpose(theta_rng_trans)
theta_rng
# plot at the true value theta0
pjp.proj_plot(
    fun=negdmglik_jax_s,
    x_opt=xsol,
    x_lims=theta_rng, 
    n_pts=100, 
    x_names = ["mu","sigma","alpha","c","s"],
    opt_vlines=True
)
```

## damaged y and undamaged ystar

```python
import matplotlib.pyplot as plt
ystar = np.arange(0.01, 99,0.01)

#dmginverse_py
y1 = jax.vmap(lambda ystar:dmginverse_py(ystar,alpha = -50,l = l,c = 0.69,s = 0.5))(ystar)
y2 = jax.vmap(lambda ystar:dmginverse_py(ystar,alpha = 0.25,l = l,c = 0.69,s = 1))(ystar)

plt.plot(ystar, y1, label = "alpha = -50, s = 0.5", linestyle="-")
plt.plot(ystar, y2, label = "alpha = 0.25, s = 1", linestyle="-")
# plt.axvline(x = l/c, 
#            color = 'purple',
#            label = 'damage threshold',linestyle = ":")
plt.axvline(x = min(R20R100), 
           color = 'purple',
           label = 'min of R20R100',linestyle = ":")
plt.xlabel('ystar')
plt.ylabel('y')

plt.legend()
plt.show()
```

```python
y = np.arange(0.01, 99,0.01)

#dmginverse_py
ystar1 = jax.vmap(lambda y:dmgmodel_py(y,alpha = -50,l = l,c = 0.69,s = 0.5))(y)
ystar2 = jax.vmap(lambda y:dmgmodel_py(y,alpha = 0.25,l = l,c = 0.69,s = 1))(y)

plt.plot(y, ystar1, label = "alpha = -50, s = 0.5", linestyle="-")
plt.plot(y, ystar2, label = "alpha = 0.2, s = 1", linestyle="-")
# # plt.axvline(x = l/c, 
# #            color = 'purple',
# #            label = 'damage threshold',linestyle = ":")
# plt.axvline(x = l/.69, 
#            color = 'purple',
#            label = 'damage threshold',linestyle = ":")
plt.xlabel('y')
plt.ylabel('ystar')

plt.legend()
plt.show()




y = np.arange(0.01, 99,0.01)

#dmginverse_py
ystar1 = jax.vmap(lambda y:dmgmodel_py(y,alpha = -50,l = l,c = 0.69,s = 0.5))(y)
ystar2 = jax.vmap(lambda y:dmgmodel_py(y,alpha = 0.25,l = l,c = 0.69,s = 1))(y)

plt.plot(y, ystar1, label = "alpha = -50, s = 0.5", linestyle="-")
plt.plot(y, ystar2, label = "alpha = 0.2, s = 1", linestyle="-")
# # plt.axvline(x = l/c, 
# #            color = 'purple',
# #            label = 'damage threshold',linestyle = ":")
plt.axvline(x = l/.69, 
           color = 'purple',
           label = 'damage threshold',linestyle = ":")
plt.xlabel('y')
plt.ylabel('ystar')
plt.ylim(0, 100)

plt.legend()
plt.show()

```

```python

# dmginverse_py(ystar =  50,alpha = 0.2,l = l,c = 0.69,s = 1)
# dmginverse_py(ystar =  50,alpha = -50,l = l,c = 0.69,s = .5)

```

```python
# import pandas
# import numpy as np

# realdata = pandas.read_csv("bending-pl.csv", header = None)
# realdata[0] = realdata[0]/1000*6.895

# l = 4500/1000*6.895
# # R100id = (realdata[1] == "R100")
# # R100 = jnp.array(realdata[realdata[1] == "R100"][0])


# R100 = jnp.array(realdata[realdata[1] == "R100"][0])
# R20 = jnp.array(realdata[realdata[1] == "R20"][0])
# R20R100 = jnp.array(realdata[realdata[1] == "R20R100"][0])

# g2 = 97



# #s= 1
# l = l
# R100 = R100
# R20 = R20
# R20R100 = R20R100
# g2 = 97
# @jax.jit
# def negdmglik_jax_s(theta):
#     mu = theta[0]
#     sigma = theta[1]
#     alpha = theta[2]
#     alpha = jnp.exp(theta[2])
#     c = theta[3]
#     #s = jnp.exp(theta[4])
#     s = theta[4]
#     #lik1 = jnp.sum(jnp.log(jax.vmap(lambda y: jax.scipy.stats.norm.pdf(y,loc = mu, scale = sigma))(y_obs_g1)))
#     lik1 = dmglik_vmap(y_group = R20,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
#     lik2 = g2*jnp.log(
#         jax.scipy.stats.norm.cdf(dmginverse_py(l,alpha,l,c,s), loc=mu, scale=sigma) - 
#         jax.scipy.stats.norm.cdf(l, loc=mu, scale=sigma)
#     )
#     lik3 = dmglik_vmap(y_group = R20R100,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
#     lik4 = jnp.sum(jax.scipy.stats.norm.logpdf(R100,loc = mu, scale = sigma))

#     return(-lik1 - lik2-lik3- lik4)

# theta0 = jnp.array([mu,sigma,jnp.log(alpha),c,s])

# #theta0 = jnp.array([mu,sigma,alpha,c,s])

# negdmglik_jax_s(theta0)
```

```python
# # solver = jaxopt.LBFGS(fun=negdmglik_jax)
# # res = solver.run(theta0)



# solver = jaxopt.ScipyMinimize(method = "Nelder-Mead",fun=negdmglik_jax_s, maxiter =1000)
# res = solver.run(theta0)

# res.params
# theta0

# # Array(1678.6377, dtype=float32)
# negdmglik_jax_s(res.params)
```

```python
# xsol = res.params
# xrange = .5
# theta_upper = xsol  + xrange *jnp.abs(xsol)
# theta_lower = xsol - xrange * jnp.abs(xsol)
# theta_rng_trans = jnp.array([theta_lower, theta_upper])

# #res.params
# theta_rng = jnp.transpose(theta_rng_trans)
# theta_rng
# # plot at the true value theta0
# pjp.proj_plot(
#     fun=negdmglik_jax_s,
#     x_opt=xsol,
#     x_lims=theta_rng, 
#     n_pts=100, 
#     x_names = ["mu","sigma","alpha","c","s"],
#     opt_vlines=True
# )
```

# prediction check

```python
import pandas
import numpy as np

realdata = pandas.read_csv("bending-pl.csv", header = None)
realdata[0] = realdata[0]/1000*6.895

l = 4500/1000*6.895
# R100id = (realdata[1] == "R100")
# R100 = jnp.array(realdata[realdata[1] == "R100"][0])


R100 = jnp.array(realdata[realdata[1] == "R100"][0])
R20 = jnp.array(realdata[realdata[1] == "R20"][0])
R20R100 = jnp.array(realdata[realdata[1] == "R20R100"][0])

g2 = 97



#s= 1
l = l
R100 = R100
R20 = R20
R20R100 = R20R100
g2 = 97
@jax.jit
def negdmglik_jax_s(theta):
    mu = theta[0]
    sigma = theta[1]
    alpha = theta[2]
    #alpha = jnp.exp(theta[2])
    c = theta[3]
    #s = jnp.exp(theta[4])
    s = theta[4]
    #lik1 = jnp.sum(jnp.log(jax.vmap(lambda y: jax.scipy.stats.norm.pdf(y,loc = mu, scale = sigma))(y_obs_g1)))
    lik1 = jnp.sum(jax.scipy.stats.norm.logpdf(R20,loc = mu, scale = sigma))
    lik2 = g2*jnp.log(
        jax.scipy.stats.norm.cdf(dmginverse_py(l,alpha,l,c,s), loc=mu, scale=sigma) - 
        jax.scipy.stats.norm.cdf(l, loc=mu, scale=sigma)
    )
    lik3 = dmglik_vmap(y_group = R20R100,alpha = alpha,l = l, c = c,s = s, mu = mu, sigma = sigma)
    lik4 = jnp.sum(jax.scipy.stats.norm.logpdf(R100,loc = mu, scale = sigma))

    return(-lik1 - lik2-lik3- lik4)

theta0 = jnp.array([mu,sigma,alpha,c,s])

#theta0 = jnp.array([mu,sigma,alpha,c,s])

# solver = jaxopt.LBFGS(fun=negdmglik_jax)
# res = solver.run(theta0)



solver = jaxopt.ScipyMinimize(method = "Nelder-Mead",fun=negdmglik_jax_s, maxiter =1000)
res = solver.run(theta0)

res.params, res.state
theta0

# Array(1678.6377, dtype=float32)
negdmglik_jax_s(res.params)
```

```python
xsol = res.params
xrange = .5
theta_upper = xsol  + xrange *jnp.abs(xsol)
theta_lower = xsol - xrange * jnp.abs(xsol)
theta_rng_trans = jnp.array([theta_lower, theta_upper])

#res.params
theta_rng = jnp.transpose(theta_rng_trans)
theta_rng
# plot at the true value theta0
pjp.proj_plot(
    fun=negdmglik_jax_s,
    x_opt=xsol,
    x_lims=theta_rng, 
    n_pts=100, 
    x_names = ["mu","sigma","alpha","c","s"],
    opt_vlines=True
)
```

```python
import seaborn as sns

sns.distplot(R20R100, hist = False")
plt.legend(prop="real data")

```

```python
# data generation 
key = jax.random.PRNGKey(0)
key, *subkeys = jax.random.split(key, num=N)



# first 
y = res.params[1]*jax.random.normal(subkeys[0], shape=(N, )) + res.params[0]

y = y[y>0]

l = 4500/1000*6.895



# g23_star
g23 = y[y>l]
g23_star = jax.vmap(lambda y: dmgmodel_py(y,res.params[2],l,res.params[3],res.params[4]))(g23)


# g3
y_obs_g3_1 = g23_star[g23_star > l]


# second 




y = res.params[1]*jax.random.normal(subkeys[1], shape=(N, )) + res.params[0]

y = y[y>0]

l = 4500/1000*6.895



# g23_star
g23 = y[y>l]
g23_star = jax.vmap(lambda y: dmgmodel_py(y,res.params[2],l,res.params[3],res.params[4]))(g23)


# g3
y_obs_g3_2 = g23_star[g23_star > l]



# third 
y = res.params[1]*jax.random.normal(subkeys[2], shape=(N, )) + res.params[0]

y = y[y>0]

l = 4500/1000*6.895



# g23_star
g23 = y[y>l]
g23_star = jax.vmap(lambda y: dmgmodel_py(y,res.params[2],l,res.params[3],res.params[4]))(g23)


# g3
y_obs_g3_3 = g23_star[g23_star > l]


# fourth
y = res.params[1]*jax.random.normal(subkeys[3], shape=(N, )) + res.params[0]

y = y[y>0]

l = 4500/1000*6.895



# g23_star
g23 = y[y>l]
g23_star = jax.vmap(lambda y: dmgmodel_py(y,res.params[2],l,res.params[3],res.params[4]))(g23)


# g3
y_obs_g3_4 = g23_star[g23_star > l]
```

```python
import seaborn as sns

sns.distplot(R20R100, hist = False, label = "real data")

sns.distplot(y_obs_g3_1, hist = False, label =" simulated data 1",kde_kws={'linestyle':'--'})
sns.distplot(y_obs_g3_2, hist = False, label =" simulated data 2",kde_kws={'linestyle':'--'})
sns.distplot(y_obs_g3_3, hist = False, label =" simulated data 3",kde_kws={'linestyle':'--'})

plt.legend()

```
