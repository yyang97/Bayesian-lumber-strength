import relax.newton
import relax as rx
from jax import lax, jit
import jax.scipy as jsp
import jax.numpy as jnp
import jax


@jax.jit

def newton_step_ETEL(G,eta):
    """
    Single step of the Newton-Raphson solver.

    Args:
        G: Matrix of moment conditions.
        eta: Initial of the Langrange multiplier.
        norm_weights: Vector of weights normalized to sum to one.

    Returns:
        Final value of the Lagrange multiplier.
    """
    Gl =  jnp.dot(G, eta)


    def calc_Q12(gl, g):
        # w is weights
        # Q1 takes the gradient of lambda, so it times g.
        # Q2 takes the second gradient of lambda, so it times g*g.
        Q1 = jnp.exp(gl) * g
        Q2 = jnp.exp(gl) * jnp.outer(g, g)
        return Q1, Q2

    Q1, Q2 = jax.vmap(calc_Q12)(Gl, G)
    Q1 = jnp.sum(Q1, axis=0)
    Q2 = jnp.sum(Q2, axis=0)
    return eta - relax.newton.chol_solve(Q2, Q1)


def newton_solve_ETEL(G,eta):
    """
    Newton-Raphson solver.
    eta is the langrange multiplier
    TODO:
    - Full documentation.
    - Early exit when tolerance is reached.
    - Converge check.
    """

    def step(eta, _):
        eta = newton_step_ETEL(G, eta)
        return eta, eta
    n_steps = 100
    eta, eta_iter = jax.lax.scan(step, eta, jnp.arange(n_steps))
    return eta,eta_iter


def get_omega(G, eta):
    """
    Calculate EL probabilities.

    TODO: Full documentation.
    """
    g1 = jnp.exp(jnp.dot(G, eta))
    omega = g1/jnp.sum(g1)
    return omega


def logETEL(G):
    """
    The log ETEL value of G. armin_{\theta} n^{-1}\sum -log(nw)
    Returns: The objective value,n^{-1}\sum -log(nw), note that this is value that you wanna minimize.
    """
    eta_init = jnp.zeros(G.shape[1])
    eta,eta_iter = newton_solve_ETEL(G,eta_init)
    omega = get_omega(G, eta)
    n = G.shape[0]
    ll = 1/n*jnp.sum(-jnp.log(n*omega))
    # convergence check 
    # check the L2-norm of eta_{n_step} - eta_{n_step-1}
    # the closer error is to 0, the better the result is 
    error = jnp.linalg.norm(eta_iter[-2] - eta_iter[-1])
    return ll,omega,eta,error

def _logETEL(G):
    ll,_,_,error = logETEL(G)
    tol = 1e-1
    return(jnp.select([error < tol,error >= tol],[ll,jnp.inf],default = jnp.inf))


