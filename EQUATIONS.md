# PRISMS PhaseField: Snow Crystal Growth

*Mathematical derivation by Professor Pal, based on the model presented by Demange et al. Transcribed from the supplied three-page derivation, preserving its notation and calculation sequence.*

> This document records the original derivation. The PRISMS-PF 4.0.0 application uses configurable flux conventions, anisotropy regularization, and natural boundary conditions; the derivation below is reproduced without changing those equations.

## Kinetics

The first PDE is given by:

```math
A(\mathbf{n})^2\partial_t\phi
=
-f'(\phi)+\lambda B(\mathbf{n})g'(\phi)u
+\frac{1}{2}\nabla_\Gamma\cdot
\left(
|\nabla\phi|^2\frac{\partial[A(\mathbf{n})^2]}{\partial\nabla\phi}
+A(\mathbf{n})^2\nabla_\Gamma\phi
\right)
```

Defining $\nabla_\Gamma=\Gamma\nabla$, we get

```math
A(\mathbf{n})^2\frac{\partial\phi}{\partial t}
=
\left(-f'(\phi)+\lambda B(\mathbf{n})g'(\phi)u\right)
+\nabla\cdot\frac{\Gamma^T}{2}
\left(
|\nabla\phi|^2\frac{\partial[A(\mathbf{n})^2]}{\partial\nabla\phi}
+A(\mathbf{n})^2\Gamma\nabla\phi
\right)
```

The second PDE is given by:

```math
\partial_t u
=
\tilde{D}\nabla_\Gamma\cdot(q(\phi)\nabla_\Gamma)u
-\frac{L_{\mathrm{sat}}}{2}B(\mathbf{n})\partial_t\phi
```

or

```math
\frac{\partial u}{\partial t}
=
\tilde{D}\nabla\cdot
\left(\Gamma^T\Gamma q(\phi)\nabla u\right)
-\frac{L_{\mathrm{sat}}}{2}B(\mathbf{n})\frac{\partial\phi}{\partial t}
```

Thus, the PDEs are:

```math
A(\mathbf{n})^2\frac{\partial\phi}{\partial t}
=
\left(-f'(\phi)+\lambda B(\mathbf{n})g'(\phi)u\right)
+\nabla\cdot\frac{\Gamma^T}{2}
\left(
|\nabla\phi|^2\frac{\partial[A(\mathbf{n})^2]}{\partial\nabla\phi}
+A(\mathbf{n})^2\Gamma\nabla\phi
\right)
```

```math
\frac{\partial u}{\partial t}
=
\tilde{D}\nabla\cdot
\left(\Gamma^T\Gamma q(\phi)\nabla u\right)
-\frac{L_{\mathrm{sat}}}{2}B(\mathbf{n})\frac{\partial\phi}{\partial t}
```

In abbreviated form, they can be rewritten as

```math
A(\nabla\phi)^2\frac{\partial\phi}{\partial t}
=f_1(\phi,\nabla\phi,u)+\nabla\cdot F_1(\nabla\phi)
```

```math
\frac{\partial u}{\partial t}
=\nabla\cdot F_2(\phi,\nabla u)
+f_2(\nabla\phi)\frac{\partial\phi}{\partial t}
```

where:

```math
f_1(\phi,\nabla\phi,u)
=-f'(\phi)+\lambda B(\mathbf{n})g'(\phi)u
```

```math
F_1(\nabla\phi)
=\frac{\Gamma^T}{2}
\left(
|\nabla\phi|^2\frac{\partial[A(\mathbf{n})^2]}{\partial\nabla\phi}
+A(\mathbf{n})^2\Gamma\nabla\phi
\right)
```

```math
f_2(\nabla\phi)=-\frac{L_{\mathrm{sat}}}{2}B(\mathbf{n})
```

```math
F_2(\phi,\nabla u)=\tilde{D}\Gamma^T\Gamma q(\phi)\nabla u
```

## Time discretization

To do time discretization, we rewrite our equations in the following form.

```math
\frac{\partial\phi}{\partial t}
=\frac{1}{A(\nabla\phi)^2}
\left(\xi_1(\phi,\nabla\phi,u)\right)
```

```math
\frac{\partial u}{\partial t}
=\nabla\cdot F_2(\phi,\nabla u)
+\frac{f_2(\nabla\phi)}{A(\nabla\phi)^2}
\left(\xi_1(\phi,\nabla\phi,u)\right)
```

```math
\xi_1(\phi,\nabla\phi,u)
=f_1(\phi,\nabla\phi,u)+\nabla\cdot F_1(\nabla\phi)
```

Now we can write the stepping schemes:

```math
\xi_1^n=f_1^n+\nabla\cdot F_1^n
```

```math
\phi^{n+1}=\phi^n+\frac{\Delta t}{A_n^2}\xi_1^n
```

```math
u^{n+1}=u^n+\Delta t(\nabla\cdot F_2^n)
+\Delta t\frac{f_2^n}{A_n^2}\xi_1^n
```

## Weak formulation

To get the weak form, we multiply the above stepping schemes with $w$ and integrate over the volume $\Omega$.

```math
\int_\Omega w\xi_1^n\,dV
=\int_\Omega wf_1^n\,dV
+\int_\Omega w\nabla\cdot F_1^n\,dV
```

```math
\int_\Omega w\phi^{n+1}\,dV
=\int_\Omega w
\left(\phi^n+\frac{\Delta t}{A_n^2}\xi_1^n\right)dV
```

```math
\int_\Omega wu^{n+1}\,dV
=\int_\Omega w
\left(u^n+\Delta t\frac{f_2^n}{A_n^2}\xi_1^n\right)dV
+\Delta t\int_\Omega w\nabla\cdot F_2^n\,dV
```

Applying the divergence theorem, and ensuring $w=0$ at the boundary (since $u,\phi$ are known on $\partial\Omega$):

```math
\int_\Omega w\xi_1^n\,dV
=\int_\Omega wf_1^n\,dV
-\int_\Omega(\nabla w)\cdot F_1^n\,dV
```

```math
\int_\Omega w\phi^{n+1}\,dV
=\int_\Omega w
\left(\phi^n+\frac{\Delta t}{A_n^2}\xi_1^n\right)dV
```

```math
\int_\Omega wu^{n+1}\,dV
=\int_\Omega w
\left(u^n+\Delta t\frac{f_2^n}{A_n^2}\xi_1^n\right)dV
-\Delta t\int_\Omega(\nabla w)\cdot F_2^n\,dV
```

This gives us the following residuals:

```math
r_\xi=f_1^n=-f'(\phi^n)+\lambda B_ng'(\phi^n)u^n
```

```math
r_{\xi x}=-F_1^n=-\frac{\Gamma^T}{2}
\left(
|\nabla\phi^n|^2\frac{\partial[A_n^2]}{\partial\nabla\phi^n}
+A_n^2\Gamma\nabla\phi^n
\right)
```

```math
r_\phi=\phi^n+\frac{\Delta t}{A_n^2}\xi_1^n
```

```math
r_u=u^n+\Delta t\frac{f_2^n}{A_n^2}\xi_1^n
=u^n-\Delta t\frac{\xi_1^n}{A_n^2}\frac{L_{\mathrm{sat}}}{2}B_n
```

```math
r_{ux}=(-\Delta t)F_2^n
=(-\Delta t)\tilde{D}\Gamma^T\Gamma q(\phi^n)\nabla u^n
```

The term involving the partial derivative with respect to $\nabla\phi$ is expressed by:

```math
\frac{\partial[A_n^2]}{\partial\nabla\phi^n}
=
\frac{\partial[A_n^2]}{\partial\theta}
\frac{\partial\theta}{\partial\nabla\phi^n}
+
\frac{\partial[A_n^2]}{\partial\psi}
\frac{\partial\psi}{\partial\nabla\phi^n}
```

```math
\frac{\partial\theta}{\partial\nabla\phi}
=
\left(
\frac{\partial\theta}{\partial(\frac{\partial\phi}{\partial x})},
\frac{\partial\theta}{\partial(\frac{\partial\phi}{\partial y})},
\frac{\partial\theta}{\partial(\frac{\partial\phi}{\partial z})}
\right)
=
\left(
\frac{\partial\theta}{\partial\phi_x},
\frac{\partial\theta}{\partial\phi_y},
\frac{\partial\theta}{\partial\phi_z}
\right)
```

```math
\tan\theta=\frac{\phi_y}{\phi_x},
\qquad
\tan\psi=\frac{\sqrt{\phi_x^2+\phi_y^2}}{-\phi_z}
```

Therefore:

```math
\frac{\partial\theta}{\partial\nabla\phi}
=\frac{1}{\phi_x}
\left(-\cos^2\theta\tan\theta,\ \cos^2\theta,\ 0\right)
```

```math
\frac{\partial\psi}{\partial\nabla\phi}
=\frac{1}{\phi_z^2\tan\psi}
\left(
\phi_x\cos^2\psi,\ \phi_y\cos^2\psi,\
-\phi_z\tan^2\psi\cos^2\psi
\right)
```
