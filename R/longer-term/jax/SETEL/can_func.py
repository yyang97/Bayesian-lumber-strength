# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

N = 300
R = 60

# +
# use self-defined incomplete gamma function
import jax
import jax.scipy as jsp
import jax.numpy as jnp
from jax import lax, jit
import jaxopt
import relax as rx
import relax.logel
import relax.newton
import relax.logstar
import ott
import logetel


@jax.jit


def gammainc_series(a, x):
    """
    Incomplete Gamma function by series expansion.
    
    Args:
        a: Shape parameter (scalar).
        x: Quantile parameter (scalar).
        n_terms: Number of terms in the expansion.
        
    Returns:
        An estimate of `jax.scipy.gammainc(a=a, x=x)`.
    """
    n_terms = 100
    C = jnp.cumprod(x / (a + jnp.arange(n_terms) + 1.))
    S = 1. + jnp.sum(C)
    f = jnp.power(x, a) * jnp.exp(-x) / jsp.special.gamma(a + 1.)
    return f * S


def gammainc_confrac(a, x):
    """
    Incomplete Gamma function by continued fraction.

    Args:
        a: Shape parameter (scalar).
        x: Quantile parameter (scalar).
        n_terms: Number of terms in the expansion.
        
    Returns:
        An estimate of `jax.scipy.gammainc(a=a, x=x)`.
    """
    n_terms = 100

    def f(carry, xs):
        A0, B0, A1, B1 = carry
        an, bn = xs
        A2 = bn * A1 + an * A0
        B2 = bn * B1 + an * B0
        # avoid large / large fractions
        const = jnp.sqrt((A1*A1 + A2*A2 + B1*B1 + B2*B2) * .25)
        A1 = A1/const
        A2 = A2/const
        B1 = B1/const
        B2 = B2/const
        return (A1, B1, A2, B2), None
    n = jnp.arange(1, n_terms) + 1.
    an = (n - 1.) * (a - n)
    bn = 2*n - a + x
    init = (1., x, x + 1., x * (2. - a + x))
    carry, _ = jax.lax.scan(f=f, init=init, xs=(an, bn))
    _, _, An, Bn = carry
    S = An/Bn
    f = jnp.power(x, a) * jnp.exp(-x) / jsp.special.gamma(a)
    return 1. - f * S


def gammainc(a, x):
    """
    Incomplete Gamma function using series or continued fraction representation.

    Args:
        a: Shape parameter (scalar).
        x: Quantile parameter (scalar).
        n_terms: Number of terms in the expansion.
        
    Returns:
        An estimate of `jax.scipy.gammainc(a=a, x=x)`.
    """
    return jax.lax.cond(
        jnp.logical_or(jnp.less_equal(x, 1.), jnp.less(x, a)),
        lambda a, x: gammainc_series(a=a, x=x),
        lambda a, x: gammainc_confrac(a=a, x=x),
        # operands
        a,
        x
    )




def gamma(z):
    return jnp.exp(jax.scipy.special.gammaln(z))

def sigmoid(x, a):
    return 0.5 * (jnp.tanh(x * a / 2) + 1)


def logit(x):
    return(jnp.log(x/(1-x)))


def expit(x):
    return 1/(1+jnp.exp(-x))



def Teqn(Ts, pars, k):
    """
    The T equation to solve `Ts`.

    Args:
        Ts: The evaluated time.
        pars: The parameter vector, [(a, b, c, n, sigma_0)].
        k: The ramp-loading rate `k`.

    Returns:
        A number, the output of the T equation, which is solved to get the `Ts`.

    """
    A = pars[0]
    b = pars[1]
    C = pars[2]
    n = pars[3]
    s0 = pars[4]
    mu = 1

    As = A * k * Ts
    Cs = C * k * Ts

    intfac = 1 / mu * jnp.power(Cs, n) * Ts / \
    (n + 1) * jnp.power(1 - s0, n + 1)
    bnfr = (b + 1) / (n + 1)
    result = jnp.power(As, b) / jnp.power(Cs, n * bnfr) 
    result = result * jnp.power(mu / Ts * (n + 1), (b - n) / (n + 1)) 
    result = result *gammainc(bnfr, intfac) * gamma(bnfr) - jnp.exp(-intfac)
    return result


def const_Tc(Ts, T0, pars, k, uprbd):
    """
    To compute the constants for the constant time-to-failure Tc.

    Args:
        Ts: The evaluated time.
        T0: The time `T0` after which the load is held as constant.
        pars: The parameter vector, [a, b, c, n, s0].
        k: The ramp-loading ratio.

    Returns:
        The time-to-failure observation under the stage of constant load.

    """
    A = pars[0]
    b = pars[1]
    C = pars[2]
    n = pars[3]
    s0 = pars[4]
    mu = 1

    As = A * k * Ts
    Cs = C * k * Ts

    indicator = 0.0000001
    extraload = jnp.where(T0 / Ts - s0 >= 0, T0 / Ts - s0, indicator)

    # extraload = jnp.where(T0/Ts-s0>=0,T0/Ts-s0,0)
    intfac = jnp.where(extraload == indicator, 1, 1 / mu *
                       jnp.power(Cs, n) * Ts / (n + 1) * jnp.power(extraload, n + 1))

    bnfr = (b + 1) / (n + 1)
    res = jnp.zeros(4)

    # res[0] good
    res = res.at[0].set(1 / mu * jnp.power(As * extraload, b))
    res = res.at[1].set(1 / mu * jnp.power(Cs * extraload, n))
    result = jnp.exp(-res[1] * T0 + intfac) /jnp.power(As,b)
    result = result /jnp.power(Cs,n * bnfr)*jnp.power(mu / Ts * (n + 1),(b - n) / (n + 1))
    result = result * gammainc(bnfr, intfac) * gamma(bnfr) 
    res = res.at[2].set(result)
#     res = res.at[2].set(jnp.exp(-res[1] * T0 + intfac) * jnp.power(As,
#                                                                    b) / jnp.power(Cs,
#                                                                                   n * bnfr) * jnp.power(mu / Ts * (n + 1),
#                                                                                                         (b - n) / (n + 1)) * jsp.special.gammainc(bnfr,
#                                                                                                                                                   intfac) * gamma(bnfr))
    res = res.at[3].set(jnp.exp(-res[1] * T0))
    Tc = -1 / res[1] * jnp.log((res[0] / res[1] *
                                res[3] + res[2]) / (1 + res[0] / res[1]))

    return jnp.where(Tc > 0, Tc, uprbd)



def Teqnroot(pars, k):
    """
    The root-finding for Teqn.
    """
    bisec = jaxopt.Bisection(
        optimality_fun=Teqn,
        lower=0.00001,
        upper=0.1,
        check_bracket=False)
    return bisec.run(pars=pars, k=k).params


def canmodel_trans(eps, theta, sigma_b_s, k, tau_c, uprbd):
    """
    The function transforms the standard normal noise to time-to-failure observations from the Canadian model.

    Args:
        eps: The noise ~ N(0,1), with len = 5.
        theta: The parameter.
        sigma_b_s: the sigma parameter sigma_b, sigma_C,sigma_n,sigma_s
        k: The rate in the ramp load.
        tau_c: The level of constant load.
        N: The sample size.
        upbrd: All censored data set as upbrd.

    Returns:
        One single time-to-failure observation, a number.
    """
    # generate the noise
    # theta = theta.at[5].set(jnp.exp(theta[5]))
    # generate the result object
    # eps_obs = jnp.array(jax.random.normal(key, shape=[N,5]))
    # simualate the five random effects from the parameters
    A = jnp.exp(theta[0] + jnp.exp(theta[5]) * eps[0])
    b = jnp.exp(theta[1] + sigma_b_s[0] * eps[1])
    C = jnp.exp(theta[2] + sigma_b_s[1] * eps[2])
    n = jnp.exp(theta[3] + sigma_b_s[2] * eps[3])
    # use transformations so that s0 is bounded at [0,1]
    eta = jnp.exp(theta[4] + sigma_b_s[3] * eps[4])

    s0 = eta / (1 + eta)

    pars = jnp.array([A, b, C, n, s0])
    # calculate the T0, after which the load is held as a constant
    T0 = tau_c / k
    # calculate the two boundary values of Teqs, preparing to solve Teqs
#     Teqnupper = Teqn(Ts = 0.1, pars = pars, k = k)
#     Teqnlower = Teqn(Ts = 0.00001, pars = pars, k = k)
    z = Teqnroot(pars, k)
    res = jnp.select(
        [(z < 0.00001) | (z > 0.1), (T0 < 0) | (z < T0), T0 / z < pars[4]],
        [uprbd, z, uprbd], default=const_Tc(z, T0, pars, k, uprbd))
    # return res
#     res = jnp.select(
#     [res == uprbd , res!= uprbd],
#         [uprbd*(1+jnp.mean(jax.vmap(expit)(theta))),res]
#     )
    return res

# change uprbd to 

def sumstat(y_obs, quant, t_c):
    """
    The function that calculates summary statatics.

    Args:
        y_obs: The observed dataset.
        quant: The quantile values.
        t_c: The censored time.

    Returns:
    The log of quantiles of uncensored observations and the censoring proportion.

    """
    y_obs_cond = jnp.where(y_obs < t_c, jnp.log(y_obs), jnp.nan)
    sumstat_quant = jax.numpy.nanquantile(y_obs_cond, quant)
    #y_obs_nan = jnp.where(y_obs > t_c, 1, 0)
    prop_cen = jnp.mean(jnp.where(y_obs > t_c, 1, 0))
    #prop_cen = jnp.mean(jax.vmap(lambda y_obs: sigmoid(x = (y_obs - t_c), a = 1))(y_obs))
    #jnp.nanmean(y_obs_cond) = jnp.mean(jnp.where(y_obs < t_c, jnp.log(y_obs), 0))/(1-prop_cen)
    return jnp.append(sumstat_quant, 
                      jnp.array([jnp.mean(jnp.where(y_obs < t_c, jnp.log(y_obs), 0))/(1-prop_cen),
#                                  jnp.nanstd(y_obs_cond),
                                 prop_cen]))
#     return jnp.append(sumstat_quant, 
#                       jnp.array([jnp.nanmean(y_obs_cond),jnp.nanstd(y_obs_cond),
#                                  prop_cen]))

# def sumstat(y_obs,quant, t_c):
#     """
#     Continous summary statistics. 
#     For quantile values, use soft-soft from ott.tool;
#     For censoring proportion, using \sum(I(y_obs>t_c))/N \propto \sum(sigmoid(y_obs-t_c))/N
    
#     """
#     prop_cen = jnp.mean(jax.vmap(lambda y_obs: sigmoid(x = (y_obs - t_c), a = 1))(y_obs))
    
# #     sumstat_quant = jax.vmap(lambda quant: ott.tools.soft_sort.quantile(jnp.log(jnp.where(y_obs < t_c, y_obs, t_c)), 
# #                                                                         level=quant))(quant*(1-prop_cen))
#     sumstat_quant = ott.tools.soft_sort.quantile(jnp.log(jnp.where(y_obs < t_c, y_obs, t_c)), 
#                                                                         quant*(1-prop_cen))
#     return jnp.append(sumstat_quant,
#                       prop_cen)


def canmodel_gen_obs(theta, sigma_b_s, k, tau_c, uprbd, key):
    """
    The function generates the time-to-failure observations from the Canadian model.

    Args:
        theta: The parameter.
        sigma_b_s: the sigma parameter sigma_b, sigma_C,sigma_n,sigma_s.
        N: the sample size.
        k: The rate in the ramp load.
        tau_c: The level of constant load.
        N: The sample size.
        upbrd: All censored data set as upbrd.
        key: The seed to control the result.

    Returns:
        The time-to-failure observations, len = `N`.
    """
    eps_obs = jnp.array(jax.random.normal(key, shape=[N, 5]))
    return jax.vmap(
        canmodel_trans,
        in_axes=(
            0,
            None,
            None,
            None,
            None,
            None),
        out_axes=0)(
            eps_obs,
            theta,
            sigma_b_s,
            k,
            tau_c,
        uprbd)





def Gmat_adjust(G):
    """
    It transform the adjusted G matrix, which gives the adjusted log empirical likelihood.
    
    Args:
        G: the original G matrix
    
    Returns:
        A new matrix, which gives the adjusted log empirical likelihood.
    """
    R = jnp.shape(G)[0]
    return jnp.row_stack([G, -jnp.sum(G, axis=0) / R * jnp.maximum(1, jnp.log(R) / 2)])

def canmodel_gen(eps_obs, theta, sigma_b_s, k, tau_c, uprbd):
    """
    The function generates the time-to-failure observations from the Canadian model.

    Args:
        theta: The parameter.
        sigma_b_s: the sigma parameter sigma_b, sigma_C,sigma_n,sigma_s.
        N: the sample size.
        k: The rate in the ramp load.
        tau_c: The level of constant load.
        N: The sample size.
        upbrd: All censored data set as upbrd.
        key: The seed to control the result.

    Returns:
        The time-to-failure observations, len = `N`.
    """
    return jax.vmap(
        canmodel_trans,
        in_axes=(
            0,
            None,
            None,
            None,
            None,
            None),
        out_axes=0)(
            eps_obs,
            theta,
            sigma_b_s,
            k,
            tau_c,
        uprbd)


def eps_gen(key):
    return (jax.random.normal(key, shape=(R,N, 5)))


def canmodel_G(theta, sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj, key):
    """
    The G matrx from the Canadian model.

    Args:
        theta: The parameter.
        sigma_b_s: The sigma parameter sigma_b, sigma_C,sigma_n,sigma_s.
        N: The sample size.
        R: The number of simulated dataset.
        quant: The quantile values of the summary statistics, with len = q.
        y_obs_sum: The observed summary statistics.
        t_c: the censored time point, after which all data are censored.
        k: The rate in the ramp load.
        tau_c: The level of constant load.
        N: The sample size.
        upbrd: All censored data set as upbrd.
        key: The seed to control the result.

    Returns:
        The G matrix `R x (q+1)`.
    """
    eps_obs = eps_gen( key)
    t_sim = jax.vmap(
        canmodel_gen,
        in_axes=(
            0,
            None,
            None,
            None,
            None,
            None),
        out_axes=0)(
            eps_obs,
            theta,
            sigma_b_s,
            k,
            tau_c,
        uprbd)
    t_sumstat = jax.vmap(
        sumstat,
        in_axes=(
            0,
            None,
            None),
        out_axes=0)(
            t_sim,
            quant,
        t_c)
    G = jax.vmap(
        jnp.subtract,
        in_axes=(
            0,
            None),
        out_axes=0)(
            t_sumstat,
        y_obs_sum)
    if supp_adj == False:
        return(G)
    else:
        return(Gmat_adjust(G))
#     return (jnp.select([supp_adj ==False, supp_adj == True],
#                       [G,Gmat_adjust(G)]))


def canmodel_logel(theta, sigma_b_s, quant, y_obs_sum,
                   k, t_c, tau_c, uprbd,
                   supp_adj, n_steps, key):
    """
    The logEL from the Canadian model.

    Args:
        theta: The parameter.
        sigma_b_s: The sigma parameter sigma_b, sigma_C,sigma_n,sigma_s.
        N: The sample size.
        R: The number of simulated dataset.
        quant: The quantile values of the summary statistics.
        y_obs_sum: The observed summary statistics.
        t_c: the censored time point, after which all data are censored.
        k: The rate in the ramp load.
        tau_c: The level of constant load.
        N: The sample size.
        upbrd: All censored data set as upbrd.
        key: The seed to control the result.
        supp_adj: If True, calculated the adjusted empirical likelihood
    Returns:
        A number, the log empirical likelihood of the Canadian model.
    """
    G = canmodel_G(
        theta=theta,
        sigma_b_s=sigma_b_s,
        quant=quant,
        y_obs_sum=y_obs_sum,
        k=k,
        t_c=t_c,
        tau_c=tau_c,
        uprbd=uprbd,
        supp_adj = supp_adj,
        key=key)
    logEL,omega,_,_,_ = rx.logel._logel(G = G, n_steps = 100)
    return(jnp.select(
        [jnp.abs(jnp.sum(omega)-1) < 1e-3, jnp.abs(jnp.sum(omega)-1) > 1e-3],
        [logEL/G.shape[0], -jnp.inf],default = jnp.inf))





def canmmodel_logetel(theta, sigma_b_s, quant, y_obs_sum,
                   k, t_c, tau_c, uprbd,
                   supp_adj, n_steps, key):
    G = canmodel_G(
        theta=theta,
        sigma_b_s=sigma_b_s,
        quant=quant,
        y_obs_sum=y_obs_sum,
        k=k,
        t_c=t_c,
        tau_c=tau_c,
        uprbd=uprbd,
        supp_adj = supp_adj,
        key=key)
        # _logetel is calcualing 1/n*\sum{\log(\omega)}
    return(logetel._logETEL(G))





def can_G_row(params,eps, sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj):
    t_sim_row = canmodel_gen(eps, params, sigma_b_s, k, tau_c, uprbd)
    t_sumstat_row = sumstat(t_sim_row,quant,t_c)
    G_row = t_sumstat_row - y_obs_sum
    return(G_row)

# def can_G_row_jac(eps,params,sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj):
#     return(jax.jacobian(can_G_row,0)(params,eps, sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj))

can_G_row_jac = jax.jacobian(can_G_row,0)

def can_G_jac(params,sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj,key):
    """
    The calculates dG/dtheta, which should gives the same result compared to 
    jax.jacobian(SV_G,0)(params0,z_init,quant,x_obs_sum,supp_adj,key)
    
    Example:
    dGdt = SV_G_jac(params0,z_init,quant,x_obs_sum,supp_adj,key)
    # It is already checked that
    # dGdt equals to jax.jacobian(SV_G,0)(params0,z_init,quant,x_obs_sum,supp_adj,key)
    """
    eps = eps_gen(key)
    dGdt = jax.vmap(can_G_row_jac,
                    in_axes = (
                        None,
                        0,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None),
                    out_axes = 0)(
        params,eps,sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj
    )
    return(dGdt)


# #dGdt = SV_G_jac(res.params,z_init,quant,x_obs_sum,supp_adj,key)
# # It is already checked that
# # dGdt equals to jax.jacobian(SV_G,0)(params0,z_init,quant,x_obs_sum,supp_adj,key)#



def solve_V(V):
    """
    Calculate the inverse of V using choleksy decomposition, where V is positive definite.
    Args:
        V: The matrix.
    Returns:
        The inverse matrix of V.
    """
    return jax.scipy.linalg.cho_solve(jax.scipy.linalg.cho_factor(V),
                                        jnp.identity(jnp.shape(V)[0]))

def var_sandwich(G, dGdt, omega):
    """
    Calulate the standard error and the confidence interval of the parameter.
    Args:
        G: The g function matrix, g(x,theta), at theta0.
        dGdt: The jacobian of G w.r.t., theta, at theta0 
        omega: The implied probability in EL, at theta0
    Returns:
        The asymptotic Variance evaluated at params0.
    """
    ggt = jax.vmap(jnp.outer,in_axes=(
        0,
        0))(G,G)
    middle_inv = jnp.average(ggt,axis = 0, weights= omega)
    middle = jax.scipy.linalg.cho_solve(jax.scipy.linalg.cho_factor(middle_inv),
                                        jnp.identity(jnp.shape(middle_inv)[0]))
    #middle_part = 1/jnp.average(jnp.sum(jnp.square(G),axis = 1),weights = omega)
    # the gradient dG/dtheta

    outside = jnp.average(dGdt, axis = 0, weights = omega)
    # now the inverse of the matrix
    V_inv = jnp.matmul(jnp.matmul(jnp.transpose(outside),middle),
                       outside
                      )
    #V_inv = middle_part*jnp.matmul(jnp.transpose(outside_part),outside_part)
    #chol = jax.scipy.linalg.cho_factor(V_inv)
    return(solve_V(V_inv))



def V_corr(V,i,j):
    """
    Calculate the correlation between theta[i] and theta[j], given the covariance matrix of theta
    Args:
        V: The covariance matrix of theta.
        i: The index of theta[i].
        j: The index of theta[j].
    Returns: 
        The correlation between theta[i] and theta[j].
    """
    return V[i,j]/jnp.sqrt(V[i,i]*V[j,j])



def solve_beta(theta, sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj, key):
    """
    Given theta and key, find the beta [tau,eta,keppa,theta], of which the dim is 1, N_g, N_g, N_theta, respectively
    Args:
        theta: The parameter for ETEL.
        key: The key.
    Returns:
     beta [tau,eta,keppa,theta], of which the dim is 1, N_g, N_g, N_theta
     g1: tau_i for each x_i
    """
    G = canmodel_G(theta, sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj, key)
    eta_init = jnp.zeros(G.shape[1])
    eta,_ = logetel.newton_solve_ETEL(G,eta_init)

    dGdt = can_G_jac(theta,sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj,key)

    g1 = jnp.exp(jnp.dot(G, eta))
    tau = jnp.mean(g1)
    ggt = jax.vmap(jnp.outer,in_axes=(
    0,
    0))(G,G)
    keppa_inv_factor = jnp.average(ggt,axis = 0, weights= g1/tau)*jnp.sum(g1/tau)/G.shape[0]
    keppa = - jnp.matmul(solve_V(keppa_inv_factor),jnp.mean(G,axis= 0))
    beta = jnp.concatenate((theta,eta,keppa))
    beta = jnp.append(beta,tau)
    #beta = jnp.concatenate((jnp.array([tau]),keppa,eta, theta))
    return beta

# def etel_var_mr_row(beta,eps,N_g,N_theta):
#     """
#     The asymptotic varince of ETEL under misspecified model.
#     Args:
#         beta: [tau,keppa,eta,theta], of which the dim is 1, N_g, N_g, N_theta, respectively
#         G: The moment functions.
#         eta: The langrange multiplier
#     """
# #     tau = beta[0]
# #     keppa = beta[1:(N_g+1)]
# #     eta = beta[(N_g +1): (2*N_g+1)]
# #     theta = beta[-N_theta:]
#     theta = beta[:N_theta]
#     eta = beta[(N_theta):(N_theta + N_g)]
#     keppa = beta[(N_theta + N_g): (N_theta + 2*N_g)]
#     tau = beta[-1:]
    
#     g = can_G_row(theta,eps, sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj)
#     dgdt = can_G_row_jac(eps,theta,sigma_b_s, quant, y_obs_sum, t_c, k, tau_c, uprbd,supp_adj)
#     g1 = jnp.exp(jnp.dot(eta,g))
#     phi_tau = g1 - tau
#     phi_keppa = g1*g
#     phi_eta = (tau - g1) * g + g1 *jnp.outer(g,g) @ keppa
#     dgdt_eta = dgdt.T @ eta
#     phi_theta = g1 * dgdt.T @ keppa + jnp.outer(g1 * dgdt_eta,g)@keppa
#     phi_theta = phi_theta  + (tau - g1) * dgdt_eta
#     #phi = jnp.append(phi_tau,jnp.concatenate((phi_eta,phi_keppa,phi_theta)))
#     phi = jnp.concatenate((phi_theta,phi_eta,phi_keppa))
#     phi = jnp.append(phi,phi_tau)
#     #return(phi_keppa)
#     #return(phi_tau)
#     return(phi)



# # def etel_var_mr_row_jac(eps,beta,N_g,N_theta):
# #     return(jax.jacobian(etel_var_mr_row,0)(beta,eps,N_g,N_theta))
# etel_var_mr_row_jac = jax.jacobian(etel_var_mr_row, 0)

# def etel_var_mr(beta,key):
#     """
#     The asymptotic varince of ETEL under misspecified model.
#     Args:
#         beta: [tau,eta,keppa,theta], of which the dim is 1, N_g, N_g, N_theta, respectively
#         G: The moment functions.
#         eta: The langrange multiplier
#     """
#     eps = eps_gen(key)
#     theta = beta[:N_theta]
   
#     eta = beta[(N_theta):(N_theta + N_g)]

#     phi_x = jax.vmap(etel_var_mr_row,
#                      in_axes = (
#                          None,
#                          0,
#                          None,None),
#                      out_axes = 0)(beta,eps,N_g,N_theta)
#     dphidb = jax.vmap(etel_var_mr_row_jac,
#                     in_axes = (
#                          None,0,
#                          None,None),
#                      out_axes = 0)(beta,eps,N_g,N_theta)
#     dphidb_mean = jnp.mean(dphidb,axis = 0)
#     phixphix_t =  jax.vmap(jnp.outer,in_axes=(
#         0,
#         0))(phi_x,phi_x)
#     Phi = jnp.mean(phixphix_t,axis = 0)
#     #return(dphidb)
#     return(dphidb_mean,Phi,phi_x,dphidb)




