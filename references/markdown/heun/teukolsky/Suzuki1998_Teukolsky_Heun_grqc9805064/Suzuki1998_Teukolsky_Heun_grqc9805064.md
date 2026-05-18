## Perturbations of Kerr-de Sitter Black Hole and Heun's Equations

Hisao Suzuki ∗

Department of Physics, Hokkaido University Sapporo, Hokkaido 060 Japan Eiichi Takasugi †

and

Hiroshi Umetsu ‡ Department of Physics, Osaka University Toyonaka, Osaka 560 Japan

#### Abstract

It is well known that the perturbation equations of massless fields for the Kerr-de Sitter geometry can be written in the form of separable equations. The equations have five definite singularities so that the analysis has been expected to be difficult. We show that these equations can be transformed to Heun's equations, for which we are able to use known technique for the analysis of the solutions. We reproduce results known previously for the Kerr geometry and de Sitter geometry in the confluent limits of the Heun's functions. Our analysis applies can be extended to Kerr-Newman-de Sitter geometry for massless fields with spin 0 and 12 .

<sup>∗</sup> e-mail address: hsuzuki@particle.sci.hokudai.ac.jp

<sup>†</sup> e-mail address: takasugi@phys.wani.osaka-u.ac.jp

e-mail address: umetsu@phys.wani.osaka-u.ac.jp

# 1 Introduction

One of the most non-trivial aspects of the perturbation equations for Kerr geometry is the separability of the radial and the angular parts. It was Carter[1] who first found that the scalar wave function is separable in the Kerr-Newman-de Sitter geometries. Later, this observation has been extended for spin 1/2, electromagnetic fields, gravitational perturbations and gravitino for the Kerr geometries and even for the Kerr-de Sitter class of geometries. The equations are called Teukolsky equations[2]. Except for the electromagnetic and gravitational perturbations, the separability persists even for the Kerr-Newman-de Sitter solutions. An important application of this fact was the proof of the stability of the Kerr black hole[3].

Though the Teukolsky equation for Kerr geometries is separable, both angular and radial equations has two regular and one irregular singularities so that the solutions cannot be written by a single form of any special functions, but they are expressed as series of special functions whose coefficients satisfy the three term recurrence relations[4-7]. The solution of the angular equation is expressed in the form of series of Jacobi functions[4]. The solution of the radial equation is rather complicated because we need the solution which is valid in the entire region extending from the outer horizon to infinity. The solution is written in the form of series of confluent hypergeometric functions[5] which is convergent around infinity and the hypergeometric functions[6-7] which is convergent around the outer horizon. By matching these solutions in the region where both solutions are convergent, we can obtain the solution which is valid from the outer horizon to infinity[6-8]. The great benifit of this kind of solutions is that the coefficients of the series are obtained as the Post Minkowskian expansion[6-7]. This technique has been successfully applied for the Post Newtonian expansion of gravitational waves from a particle in circular orbits around a rotating black hole[9-10].

Strictly speaking, the technique of constructing solutions as series of special functions is not new, but extremely old. The similar expansion has been considered in the case of Heun's equation[11-12]. It may not be difficult to see the correspondense between the Teukolsky equations for the Kerr black hole and the Heun's equation. That is, the Teukolsky equation with two regular and one irregular singularities will be given as a confluent case of the Heun's equation which has four definite singularities as we can see from Eq.(2.2) in Ref.6 and page 27 of Ref.12.

It is known that the separability holds true for the Kerr-de Sitter geometries and even

for the Kerr-Newman-de Sitter geometries except for electormagnetic and gravitational fields. Therefore, it may be meaningful to consider whether the Teukolsky equations for these geometries can be transformed into the Heun's equations, although these Teukolsky equations have five regular singularities. A usuful suggestion for this problem can be seen in the radial part of Teukolsky equation for the de Sitter geometries[13]. We showed that although this equation has four regular singularities, the one of the singularities can be eliminated by a suitable change of variables and by the redefinition of the radial function, so that the solutions of the equation are given by the hypergeometric function. From this experience, we expect that the similar non-trivial transformations may exist even for the Kerr(-Newman)-de Sitter geometries.

Our main aim of this short letter is to show that the Teukolsky equations for Kerr-Newman-de Sitter geometries can be transformed into the Heun's equations. Here, we stress again that our discussions for the Kerr-de Sitter geometries appliy to massless particles, but those for the Kerr-Newman-de Sitter geometries do apply to massless particles except for electromagnetic and gravitational fields. Then, we can use very old technique[12] for the analysis of the solutions. The previous solutions for the Kerr geometries[6-7] and de Sitter geometries[13] can be obtained by confluent limits of the Heun's functions.

In section 2, we provide a short review of the Teukolsky equations for the Kerr-Newman-de Sitter geometires. In section 3, we show that both angular and radial equations can be transformed into the Heun's equations. In section 4, we analize the solutions of the angular and radial equations by using the results of the Heun's equation and show that the previously obtained solutions for the Kerr geometries and de Sitter geometries can be obtained by the confluent limits of the soluitions for Kerr(-Newman)-de Sitter geometries. The last section is devoted to some discussions and further possible applications.

# 2 Teukolsky equation for the Kerr(-Newman)-de Sitter geometry

We consider the Teukolsky equations for the Kerr-Newman-de Sitter geometries. Since the electromagnetic field and the gravitational field couple through electric charge of black hole in this geometry and the equation is not separable so far, we deal with massless fields with any spins except for these fields. In the limit of null electric charge, Q → 0, the Kerr-de Sitter geometries are realized where the Teukolsky equation is valid for all massless fields. In the Boyer-Lindquist coordinates the Kerr-Newman-de Sitter metric has the form,

$$ds^{2} = -\rho^{2} \left( \frac{dr^{2}}{\Delta_{r}} + \frac{d\theta^{2}}{\Delta_{\theta}} \right) - \frac{\Delta_{\theta} \sin^{2} \theta}{(1+\alpha)^{2} \rho^{2}} [adt - (r^{2} + a^{2})d\varphi]^{2}$$

$$+ \frac{\Delta_{r}}{(1+\alpha)^{2} \rho^{2}} (dt - a\sin^{2} \theta d\varphi)^{2}, \qquad (2.1)$$

where

$$\Delta_{r} = (r^{2} + a^{2}) \left( 1 - \frac{\alpha}{a^{2}} r^{2} \right) - 2Mr + Q^{2} = -\frac{\alpha}{a^{2}} (r - r_{+})(r - r_{-})(r - r'_{+})(r - r'_{-}),$$

$$\Delta_{\theta} = 1 + \alpha \cos^{2} \theta, \qquad \alpha = \frac{\Lambda a^{2}}{3},$$

$$\bar{\rho} = r + ia \cos \theta, \qquad \rho^{2} = \bar{\rho} \bar{\rho}^{*}.$$
(2.2)

Here  $\Lambda$  is the cosmological constant, M is the mass of the black hole, aM its angular momentum and Q its charge. The electromagnetic field due to the charge of the black hole is given by

$$A_{\mu}dx^{\mu} = -\frac{Qr}{(1+\alpha)^{2}\rho^{2}}(dt - a\sin^{2}\theta d\varphi). \tag{2.3}$$

In particular we adopt the following vectors as the null tetrad,

$$l^{\mu} = \left(\frac{(1+\alpha)(r^{2}+a^{2})}{\Delta_{r}}, 1, 0, \frac{a(1+\alpha)}{\Delta_{r}}\right),$$

$$n^{\mu} = \frac{1}{2\rho^{2}}\left((1+\alpha)(r^{2}+a^{2}), -\Delta_{r}, 0, a(1+\alpha)\right),$$

$$m^{\mu} = \frac{1}{\bar{\rho}\sqrt{2\Delta_{\theta}}}\left(ia(1+\alpha)\sin\theta, 0, \Delta_{\theta}, \frac{i(1+\alpha)}{\sin\theta}\right), \quad \bar{m}^{\mu} = m^{*\mu}.$$
(2.4)

Assuming that the time and azimuthal dependence of the fields has the form  $e^{-i(\omega t - m\varphi)}$ , the tetrad components of the derivative and the electromagnetic field are

$$l^{\mu}\partial_{\mu} = \mathcal{D}_{0}, \qquad n^{\mu}\partial_{\mu} = -\frac{\Delta_{r}}{2\rho^{2}}\mathcal{D}_{0}^{\dagger},$$

$$m^{\mu}\partial_{\mu} = \frac{\sqrt{\Delta_{\theta}}}{\sqrt{2}\bar{\rho}}\mathcal{L}_{0}^{\dagger}, \qquad \bar{m}^{\mu}\partial_{\mu} = \frac{\sqrt{\Delta_{\theta}}}{\sqrt{2}\bar{\rho}^{*}}\mathcal{L}_{0},$$

$$l^{\mu}A_{\mu} = -\frac{Qr}{\Delta_{r}}, \qquad n^{\mu}A_{\mu} = -\frac{Qr}{2\rho^{2}},$$

$$m^{\mu}A_{\mu} = \bar{m}^{\mu}A_{\mu} = 0,$$

$$(2.5)$$

where

$$\mathcal{D}_n = \partial_r - \frac{i(1+\alpha)K}{\Delta_r} + n\frac{\partial_r \Delta_r}{\Delta_r},$$

$$\mathcal{D}_{n}^{\dagger} = \partial_{r} + \frac{i(1+\alpha)K}{\Delta_{r}} + n\frac{\partial_{r}\Delta_{r}}{\Delta_{r}},$$

$$\mathcal{L}_{n} = \partial_{\theta} + \frac{(1+\alpha)H}{\Delta_{\theta}} + n\frac{\partial_{\theta}(\sqrt{\Delta_{\theta}}\sin\theta)}{\sqrt{\Delta_{\theta}}\sin\theta},$$

$$\mathcal{L}_{n}^{\dagger} = \partial_{\theta} - \frac{(1+\alpha)H}{\Delta_{\theta}} + n\frac{\partial_{\theta}(\sqrt{\Delta_{\theta}}\sin\theta)}{\sqrt{\Delta_{\theta}}\sin\theta},$$
(2.6)

<span id="page-4-0"></span>with 
$$K = \omega(r^2 + a^2) - am$$
 and  $H = -a\omega \sin \theta + \frac{m}{\sin \theta}$ .

Using the Newman-Penrose formalism it is known that perturbation equations in the Kerr-de Sitter geometry are separable for massless spin  $0, \frac{1}{2}, 1, \frac{3}{2}$  and 2 fields. Similarly in the Kerr-Newman-de Sitter space those for spin  $0, \frac{1}{2}$  fields are also separable. The separated equations for fields with spin s and charge e are given by

$$\left[\sqrt{\Delta_{\theta}} \mathcal{L}_{1-s}^{\dagger} \sqrt{\Delta_{\theta}} \mathcal{L}_{s} - 2(1+\alpha)(2s-1)a\omega \cos \theta -2\alpha(s-1)(2s-1)\cos^{2}\theta + \lambda\right] S_{s}(\theta) = 0, \quad (2.7)$$

$$\left[\Delta_{r} \mathcal{D}_{1} \mathcal{D}_{s}^{\dagger} + 2(1+\alpha)(2s-1)i\omega r - \frac{2\alpha}{a^{2}}(s-1)(2s-1) + \frac{-2(1+\alpha)eQKr + iseQr\partial_{r}\Delta_{r} + e^{2}Q^{2}r^{2}}{\Delta_{r}} - 2iseQ - \lambda\right] R_{s}(r) = 0. \quad (2.8)$$

# 3 Transformation of Teukolsky equation to Heun's equation

In this section we show that the Teukolsky equations (2.7) and (2.8) can be transformed to the Heun's equations by factoring out a single regular singularity.

## 3.1 Angular Teukolsky equation

From eq.(2.7), the angular Teukolsky equation become

$$\left\{ \frac{d}{dx}(1+\alpha x^{2})(1-x^{2})\frac{d}{dx} + \lambda - s(1-\alpha) + \frac{(1+\alpha)^{2}}{\alpha}\xi^{2} - 2\alpha x^{2} + \frac{1+\alpha}{1+\alpha x^{2}} \left[ 2s(\alpha m - (1+\alpha)\xi)x - \frac{(1+\alpha)^{2}}{\alpha}\xi^{2} - 2m(1+\alpha)\xi + s^{2} \right] - \frac{(1+\alpha)^{2}m^{2}}{(1+\alpha x^{2})(1-x^{2})} - \frac{(1+\alpha)(s^{2}+2smx)}{1-x^{2}} \right\} S(x) = 0,$$
(3.1)

where  $x = \cos \theta$  and  $\xi = a\omega$ . This equation has five regular singularities at  $\pm 1$ ,  $\pm \frac{i}{\sqrt{\alpha}}$  and  $\infty$ . We also note that the angular equation has no dependence on M and Q. By choosing

<span id="page-5-0"></span>the variable z

$$z = \frac{1 - \frac{i}{\sqrt{\alpha}}}{2} \frac{x + 1}{x - \frac{i}{\sqrt{\alpha}}},$$

then eq.(3.1) takes the following form,

$$\left\{ \frac{d^{2}}{dz^{2}} + \left[ \frac{1}{z} + \frac{1}{z-1} + \frac{1}{z-z_{s}} - \frac{2}{z-z_{\infty}} \right] \frac{d}{dz} - \left( \frac{m-s}{2} \right)^{2} \frac{1}{z^{2}} - \left( \frac{m+s}{2} \right)^{2} \frac{1}{(z-1)^{2}} + \left( \frac{1+\alpha}{2\sqrt{\alpha}} \xi - \frac{\sqrt{\alpha}m + is}{2} \right)^{2} \frac{1}{(z-z_{s})^{2}} + \frac{2}{(z-z_{\infty})^{2}} + \left[ -\frac{m^{2}}{2} \left( 1 + \frac{4\alpha}{(1+i\sqrt{\alpha})^{2}} \right) + \frac{s^{2}}{2} \left( \frac{1-i\sqrt{\alpha}}{1+i\sqrt{\alpha}} \right)^{2} + \frac{2ims\sqrt{\alpha}}{1+i\sqrt{\alpha}} + \frac{\lambda - s(1-\alpha) - 2\alpha + 2(1+\alpha)(m+s)\xi}{(1+i\sqrt{\alpha})^{2}} \right] \frac{1}{z} + \left[ \frac{m^{2}}{2} \left( 1 + \frac{4\alpha}{(1-i\sqrt{\alpha})^{2}} \right) - \frac{s^{2}}{2} \left( \frac{1+i\sqrt{\alpha}}{1-i\sqrt{\alpha}} \right)^{2} - \frac{2ims\sqrt{\alpha}}{1-i\sqrt{\alpha}} - \frac{\lambda - s(1-\alpha) - 2\alpha + 2(1+\alpha)(m-s)\xi}{(1-i\sqrt{\alpha})^{2}} \right] \frac{1}{z-1} + \left[ -m^{2} \frac{2i\alpha\sqrt{\alpha}}{(1+\alpha)^{2}} - s^{2} \frac{i\sqrt{\alpha}(1-\alpha)}{(1+\alpha)^{2}} - ms\frac{\alpha}{1+\alpha} + \frac{i\sqrt{\alpha}(\lambda - s(1-\alpha) + 2)}{(1+\alpha)^{2}} + \frac{(2i\sqrt{\alpha}m + (\alpha-1)s)\xi}{1+\alpha} \right] \frac{4}{z-z_{s}} - \frac{8i\sqrt{\alpha}}{1+\alpha} \frac{1}{z-z_{\infty}} \right\} S(z) = 0, \tag{3.2}$$

where  $z_s = -\frac{i(1+i\sqrt{\alpha})^2}{4\sqrt{\alpha}}$  and  $z_{\infty} = -\frac{i(1+i\sqrt{\alpha})}{2\sqrt{\alpha}}$ . It turns out that the singularity at  $z = z_{\infty}$  which corresponds to  $x = \infty$  can be factored out by the transformation,

$$S(z) = z^{A_1}(z-1)^{A_2}(z-z_s)^{A_3}(z-z_\infty)f(z), \tag{3.3}$$

where  $A_1 = \frac{|m-s|}{2}$ ,  $A_2 = \frac{|m+s|}{2}$  and  $A_3 = \pm \frac{i}{2} \left( \frac{1+\alpha}{\sqrt{\alpha}} \xi - \sqrt{\alpha}m - is \right)$ . Now f(z) satisfies the equation

$$\left\{ \frac{d^2}{dz^2} + \left[ \frac{2A_1 + 1}{z} + \frac{2A_2 + 1}{z - 1} + \frac{2A_3 + 1}{z - z_s} \right] \frac{d}{dz} + \frac{\rho_+ \rho_- z + u}{z(z - 1)(z - z_s)} \right\} f(z) = 0,$$
(3.4)

where

$$\rho_{\pm} = A_1 + A_2 + A_3 \pm A_3^* + 1, \tag{3.5}$$

$$u = \frac{-i}{4\sqrt{\alpha}} \left\{ \lambda - s(1-\alpha) - 2i\sqrt{\alpha} + 2(1+\alpha)(m+s)\xi - (1+i\sqrt{\alpha})^2 (2A_1A_2 + A_1 + A_2) - 4i\sqrt{\alpha}(2A_1A_3 + A_1 + A_3) - \frac{m^2}{2} \left[ 4\alpha + (1+i\sqrt{\alpha})^2 \right] + \frac{s^2}{2} (1-i\sqrt{\alpha})^2 + 2ims\sqrt{\alpha}(1+i\sqrt{\alpha}) \right\}. \tag{3.6}$$

<span id="page-6-0"></span>Equation (3.4) is called the Heun's equation which has four regular singularities. The f(z) is determined by requiring non-singular behaviors at z = 0 and 1. We can take either one of signs of  $A_3$  to find the solution S(z).

#### 3.2 Radial Teukolsky equation

From eq.(2.8), the radial Teukolsky equation is explicitly written by

$$\left\{ \Delta_r^{-s} \frac{d}{dr} \Delta_r^{s+1} \frac{d}{dr} + \frac{1}{\Delta_r} \left[ (1+\alpha)^2 \left( K - \frac{eQr}{1+\alpha} \right)^2 - is(1+\alpha) \left( K - \frac{eQr}{1+\alpha} \right) \frac{d\Delta_r}{dr} \right] + 4is(1+\alpha)\omega r - \frac{2\alpha}{a^2} (s+1)(2s+1)r^2 + 2s(1-\alpha) - 2iseQ - \lambda \right\} R = 0.$$
(3.7)

This equation has five regular singularities at  $r_{\pm}$ ,  $r'_{\pm}$  and  $\infty$  which are assigned such that  $r_{\pm} \to M \pm \sqrt{M^2 - a^2 - Q^2} \equiv r_{\pm}^0$  and  $r'_{\pm} \to \pm \frac{a}{\sqrt{\alpha}}$  in the limit  $\alpha \to 0$  ( $\Lambda \to 0$ ). By using the new variable

$$z = \left(\frac{r_{+} - r'_{-}}{r_{+} - r_{-}}\right) \left(\frac{r - r_{-}}{r - r'_{-}}\right),$$

eq (3.7) becomes an equation which has regular singularities at  $0, 1, z_r, z_\infty$  and  $\infty$ ,

$$z_r = \left(\frac{r_+ - r'_-}{r_+ - r_-}\right) \left(\frac{r'_+ - r_-}{r'_+ - r'_-}\right) , \qquad z_\infty = \frac{r_+ - r'_-}{r_+ - r_-}.$$

Again we can factor out the singularity at  $z=z_{\infty}$  by the transformation as

$$R(z) = z^{B_1}(z-1)^{B_2}(z-z_r)^{B_3}(z-z_\infty)^{2s+1}g(z),$$
(3.8)

where

$$B_{1} = \frac{1}{2} \left\{ -s \pm i \left[ \frac{2(1+\alpha)a^{2} \left( \omega(r_{-}^{2}+a^{2}) - am - \frac{eQr_{-}}{1+\alpha} \right)}{\alpha(r'_{+} - r_{-})(r'_{-} - r_{-})(r_{+} - r_{-})} - is \right] \right\},$$

$$B_{2} = \frac{1}{2} \left\{ -s \pm i \left[ \frac{2(1+\alpha)a^{2} \left( \omega(r_{+}^{2}+a^{2}) - am - \frac{eQr_{+}}{1+\alpha} \right)}{\alpha(r'_{+} - r_{+})(r'_{-} - r_{+})(r_{-} - r_{+})} - is \right] \right\},$$

$$B_{3} = \frac{1}{2} \left\{ -s \pm i \left[ \frac{2(1+\alpha)a^{2} \left( \omega(r'_{+}^{2}+a^{2}) - am - \frac{eQr'_{+}}{1+\alpha} \right)}{\alpha(r_{-} - r'_{+})(r'_{-} - r'_{+})(r_{+} - r'_{+})} - is \right] \right\}.$$

Then, g(z) satisfies the Heun's equation as

$$\left\{ \frac{d^2}{dz^2} + \left[ \frac{2B_1 + s + 1}{z} + \frac{2B_2 + s + 1}{z - 1} + \frac{2B_3 + s + 1}{z - z_r} \right] \frac{d}{dz} + \frac{\sigma_+ \sigma_- z + v}{z(z - 1)(z - z_r)} \right\} g(z) = 0.$$
(3.9)

where

$$\sigma_{\pm} = B_{1} + B_{2} + B_{3} + 2s + 1 + \frac{1}{2} \left\{ -s \pm i \left[ \frac{2(1+\alpha)a^{2} \left( \omega(r'_{-}^{2} + a^{2}) - am - \frac{eQr'_{-}}{1+\alpha} \right)}{\alpha(r_{+} - r'_{-})(r_{-} - r'_{-})(r'_{+} - r'_{-})} - is \right] \right\},$$

$$v = \frac{2a^{4}(1+\alpha)^{2}}{\alpha^{2}D} \frac{(r_{+} - r'_{+})^{2}(r_{+} - r'_{-})^{2}(r_{-} - r'_{-})(r'_{+} - r'_{-})}{r_{+} - r_{-}} \left\{ -\omega^{2}r_{-}^{3}(r_{+}r_{-} - 2r_{+}r'_{+} + r_{-}r'_{+}) + 2a\omega(a\omega - m)r_{-}(r_{+}r'_{+} - r_{-}^{2}) - a^{2}(a\omega - m)^{2}(2r_{-} - r_{+} - r'_{+}) + \frac{eQ}{1+\alpha} \left[ \omega r_{-}^{2}(r_{+}r_{-} + r_{-}^{2} - 3r_{+}r'_{+} + r_{-}r'_{+}) - a(a\omega - m)(r_{+}r_{-} - 3r_{-}^{2} + r_{+}r'_{+} + r_{-}r'_{+}) \right] + \left\{ \frac{eQ}{1+\alpha} \right\}^{2} r_{-}(-r_{-}^{2} + r_{+}r'_{+}) + \frac{eQ}{1+\alpha} \left[ \frac{\omega(r_{-}r'_{-} + a^{2}) - am - \frac{eQ}{1+\alpha} \frac{r_{-}+r'_{-}}{2}}{2} \right] - am - \frac{eQ}{1+\alpha} \frac{r_{-}+r'_{-}}{2} + r_{-}r'_{-}(r'_{+} - r_{-})(r'_{+} - r'_{-})(r'_{-} - r'_{-}) + (s+1)(2s+1) \left[ \frac{2r'_{-}^{2}}{(r_{+} - r_{-})(r'_{+} - r'_{-})} - z_{\infty} \right] - 2B_{1}(z_{r}B_{2} + B_{3}) - (s+1) \left[ (1+z_{r})B_{1} + z_{r}B_{2} + B_{3} \right] - \frac{a^{2}}{\alpha(r_{+} - r_{-})(r'_{+} - r'_{-})} \left[ -\lambda - 2iseQ + 2s(1-\alpha) \right].$$
(3.10)

Here  $\mathcal{D}$  is the discriminant of  $\Delta_r = 0$ ,

$$\mathcal{D} = (r_{+} - r_{-})^{2} (r_{+} - r'_{+})^{2} (r_{+} - r'_{-})^{2} (r_{-} - r'_{+})^{2} (r_{-} - r'_{-})^{2} (r'_{+} - r'_{-})^{2}$$

$$= \frac{16a^{10}}{\alpha^{5}} \left\{ (1 - \alpha)^{3} \left[ M^{2} - (1 - \alpha)(a^{2} + Q^{2}) \right] + \frac{\alpha}{a^{2}} \left[ -27M^{4} + 36(1 - \alpha)M^{2}(a^{2} + Q^{2}) - 8(1 - \alpha)^{2}(a^{2} + Q^{2})^{2} \right] - \frac{16\alpha^{2}}{a^{4}} (a^{2} + Q^{2})^{3} \right\}.$$

The sign ambiguity in  $B_2$  or  $B_3$  are related to the boundary condition at the horizon or at the de Sitter horizon, respectively. We can take either one of signs of  $B_1$ .

## 4 Solution of Teukolsky equation

In the former section, we have shown that both the angular and radial Teukolsky equations can be transformed into Heun's equations by factoring out the singularity at  $z_{\infty}$ . Thus, both equations can be solved in a similar way by taking into account the boundary conditions. First we briefly explain the formal prescription solving the Heun's equation by a series expansion of the hypergeometric functions and then consider each case including the boundary conditions.

#### <span id="page-8-0"></span>4.1 Solution of Heun's equation

The standard form of the Heun's equation with regular singularities at 0, 1,  $a_H$  and  $\infty$  is given by [12]

$$\left\{ \frac{d^2}{dz^2} + \left[ \frac{\gamma}{z} + \frac{\delta}{z - 1} + \frac{\epsilon}{z - a_H} \right] \frac{d}{dz} + \frac{\alpha \beta z - q}{z(z - 1)(z - a_H)} \right\} y_{\nu}(z) = 0, \tag{4.1}$$

where parameters should satisfy the condition,  $\alpha+\beta+1=\gamma+\delta+\epsilon$  and we assume  $|a_H|>1$ . Solutions which are analytic in some domain including two singularities are called Heun's functions. We particularly explain a formal construction of a Heun's function relative to the points 0, 1 by means of a series of hypergeometric functions.

We first construct a solution valid in the neighborhood of 0. The solution is analytic in the interior of the ellipse with foci at 0, 1 if the condition for augmented convergence presented bellow is satisfied. To the end we write  $y_{\nu}(z)$  in the following form,

$$y_{\nu}(z) = \sum_{n=-\infty}^{+\infty} c_n^{\nu} u_{\nu+n}(z),$$
 (4.2)

$$u_{\nu}(z) = F(-\nu, \nu + \omega; \gamma; z) , \qquad (4.3)$$

where  $\omega = \gamma + \delta - 1 = \alpha + \beta - \epsilon$ .  $\omega$  is usually used for the definition of Heun's equation and should not be confused with the frequency which we use in this paper frequently.

By using the following recurrence relations[12],

$$zu_{\nu}(z) = -\frac{(\nu+\omega)(\nu+\gamma)}{(2\nu+\omega)(2\nu+\omega+1)}u_{\nu+1}(z) + \frac{2\nu(\nu+\omega)+\gamma(\omega-1)}{(2\nu+\omega+1)(2\nu+\omega-1)}u_{\nu}(z) - \frac{\nu(\nu+\delta-1)}{(2\nu+\omega)(2\nu+\omega-1)}u_{\nu-1}(z), \tag{4.4}$$

$$z(z-1)\frac{d}{dz}u_{\nu}(z) = -\frac{\nu(\nu+\omega)(\nu+\gamma)}{(2\nu+\omega)(2\nu+\omega+1)}u_{\nu+1}(z) + \frac{\nu(\nu+\omega)(\gamma-\delta)}{(2\nu+\omega+1)(2\nu+\omega-1)}u_{\nu}(z) - \frac{\nu(\nu+\omega)(\nu+\delta-1)}{(2\nu+\omega)(2\nu+\omega-1)}u_{\nu-1}(z),$$
(4.5)

it is shown that  $y_{\nu}(z)$  becomes a solution if coefficients  $c_n$ 's satisfy the following three term recurrence relations,

$$\alpha_n^{\nu} c_{n+1}^{\nu} + \beta_n^{\nu} c_n^{\nu} + \gamma_n^{\nu} c_{n-1}^{\nu} = 0, \tag{4.6}$$

where

$$\alpha_n^{\nu} = -\frac{(\nu + n + 1)(\nu + n + \omega - \alpha + 1)(\nu + n + \omega - \beta + 1)(\nu + n + \delta)}{(2\nu + 2n + \omega + 2)(2\nu + 2n + \omega + 1)}, \quad (4.7)$$

$$\beta_n^{\nu} = \frac{J_n^{\nu}}{(2\nu + 2n + \omega + 1)(2\nu + 2n + \omega - 1)} - a_H(\nu + n)(\nu + n + \omega) - q, \quad (4.8)$$

$$\gamma_n^{\nu} = -\frac{(\nu + n + \alpha - 1)(\nu + n + \beta - 1)(\nu + n + \gamma - 1)(\nu + n + \omega - 1)}{(2\nu + 2n + \omega - 2)(2\nu + 2n + \omega - 1)}, \quad (4.9)$$

<span id="page-9-0"></span>and

$$J_n^{\nu} = \epsilon(\nu + n)(\nu + n + \omega)(\gamma - \delta) + [(\nu + n)(\nu + n + \omega) + \alpha\beta][2(\nu + n)(\nu + n + \omega) + \gamma(\omega - 1)]. \quad (4.10)$$

By defining the continued fractions

$$R_n(\nu) = \frac{c_n^{\nu}}{c_{n-1}^{\nu}}, \qquad L_n(\nu) = \frac{c_n^{\nu}}{c_{n+1}^{\nu}},$$
 (4.11)

they satisfy

$$R_n(\nu) = -\frac{\gamma_n^{\nu}}{\beta_n^{\nu} + \alpha_n^{\nu} R_{n+1}(\nu)} , \qquad L_n(\nu) = -\frac{\alpha_n^{\nu}}{\beta_n^{\nu} + \gamma_n^{\nu} L_{n-1}(\nu)} . \tag{4.12}$$

Now we can evaluate the coefficients  $c_n^{\nu}$  by using either series  $R_n(\nu)$  or  $L_n(\nu)$  with appropriate initial data. The convergence of the series requires a following transcendental equation,

$$R_n(\nu)L_{n-1}(\nu) = 1. (4.13)$$

This equation determines either the eigenvalue  $\lambda$  for the angular equation, where the characteristic exponent  $\nu$  (the renormalized angular momentum) can be taken zero, or  $\nu$  for the radial cases, respectively.

It will be worthwhile to note that if  $\nu$  is a solution of eq.(4.13), then  $-\nu - \omega$  is also the solution. This can be shown as followings: First, we notice the algebraic relations,  $\alpha_{-n}^{-\nu-\omega} = \gamma_n^{\nu}$ ,  $\beta_{-n}^{-\nu-\omega} = \beta_n^{\nu}$  and  $\gamma_{-n}^{-\nu-\omega} = \alpha_n^{\nu}$ . From these relations, we find that  $c_{-n}^{-\nu-\omega}$ 's satisfy the same recurrence relations as these for  $c_n^{\nu}$ 's. If we choose  $c_0^{\nu} = c_0^{-\nu-\omega}$ , we have  $c_n^{\nu} = c_{-n}^{-\nu-\omega}$  so that

$$R_n(-\nu - \omega)L_{n-1}(-\nu - \omega) = R_n(\nu)L_{n-1}(\nu) = 1.$$
(4.14)

From eq.(4.3), we see that  $y_{-\nu-\omega}(z)$  agrees with  $y_{\nu}(z)$  by choosing  $c_0^{\nu}=c_0^{-\nu-\omega}$ .

## 4.2 Solution of the angular equation

By comparing Heun's equation (4.1) with eq.(3.4), parameters to define Heun's equation,  $\alpha, \beta, \gamma, \delta, \epsilon, \omega, q$  and  $a_H$  are given by

$$\alpha = \rho_{+} , \qquad \beta = \rho_{-} , \qquad \gamma = 2A_{1} + 1 ,$$

$$\delta = 2A_{2} + 1 , \quad \epsilon = 2A_{3} + 1 , \quad \omega = 2(A_{1} + A_{2}) + 1 , \qquad (4.15)$$

$$q = -u , \qquad a_{H} = z_{s} ,$$

We require that f(z) is a regular function at z = 0, 1, which is satisfied by taking  $\nu = 0$ . Then the solution is given in the form of the series of the Jacobi polynomials extending from n = 0 to  $\infty$ , similarly to the Kerr geometry case,

$$f(z) = \sum_{n=0}^{\infty} a_n u_n(z), \tag{4.16}$$

$$u_n(z) = F(-n, n+\omega; \gamma; z) = (-)^n \frac{\Gamma(2n+\omega)n!}{\Gamma(n+\gamma)} P_n^{(\omega-\gamma,\gamma-1)}(2z-1). \tag{4.17}$$

The coefficients  $a_n$ 's are determined by solving the equation for  $L_n(\nu)$  in eq.(4.6) with  $L_{-1}(\nu) = 0$ , where  $a_0$  gives the overall normalization. The eigenvalue  $\lambda$  can be obtained by transcendental equation (4.13) with the condition  $R_{\infty}^{\nu} = L_{-1}^{\nu} = 0$ . We give  $\lambda$  up to  $O(\xi^3)$  and the coefficient of each power of  $\xi$  up to  $O(\alpha^2)$  as

$$\lambda = l(l+1) - s^{2} + s$$

$$+\alpha \left[ -l(l+1) + s^{2} - s + 2m^{2} + \frac{2m^{2}s^{2}}{l(l+1)} - l^{2}H(l) + (l+1)^{2}H(l+1) \right]$$

$$+ \left\{ -2m\left(1 + \frac{s^{2}}{l(l+1)}\right) + \left(1 + \frac{s^{2}}{l(l+1)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) + \left(1 + \frac{s^{2}}{l(l+2)}\right) +$$

where

$$H(l) = \frac{2(l^2 - m^2)(l^2 - s^2)^2}{(2l - 1)l^3(2l + 1)},$$
(4.19)

and l is the angular momentum which takes an integer or half integer number satisfying  $l \geq \max(|m|, |s|)$ . From the result it seems that  $\lambda - s(1 - \alpha)$  is an even function of s. Indeed this holds true because eq.(4.13) which determines  $\lambda$  contains  $\alpha_n, \beta_n$  and  $\gamma_n$  in the forms of  $\alpha_{n-1}\gamma_n$  and  $\beta_n + \frac{i}{4\sqrt{\alpha}}(\lambda - s(1 - \alpha))$ , which are invariant under  $s \to -s$ , respectively.

In the rest of this section, we consider two limits, the Kerr-Newman limit and the Reissner-Nordström-de Sitter limit. It may be worthwhile to note that the angular Teukolsky equation (3.1) does not contain mass and charge parameters of the black hole so that the above limits are only meaningful ones.

The Kerr-Newman limit:  $\Lambda = (3/a^2)\alpha \rightarrow 0$ 

In the limit, parameters behave as

$$\rho_{+} = A_{1} + A_{2} \pm s + 1 + O(\sqrt{\alpha}), \quad \rho_{-} = \pm \frac{i}{\sqrt{\alpha}} \xi + O(1), \quad A_{3} = \pm \frac{i}{2\sqrt{\alpha}} \xi + O(1),$$

$$u = -\frac{i}{4\sqrt{\alpha}} \left[ \lambda - 2A_{1}A_{2} - A_{1} - A_{2} + 2(m + s \mp (2A_{1} + 1)) \xi - \frac{m^{2} - s^{2}}{2} - s \right], \quad (4.20)$$

$$z_{s} = -\frac{i}{4\sqrt{\alpha}} + O(1), \qquad z_{\infty} = -\frac{i}{2\sqrt{\alpha}} + O(1),$$

then parameters in the recurrence relation (4.7) -(4.9) become by omitting the superscript  $\nu$  as

$$\alpha_{n} = \pm \frac{i}{\sqrt{\alpha}} \xi \frac{(n+1)(n+A_{1}+A_{2}\mp s+1)(n+2A_{2}+1)}{2(2n+2A_{1}+2A_{2}+3)(n+A_{1}+A_{2}+1)}, \tag{4.21}$$

$$\beta_{n} = \frac{i}{\sqrt{\alpha}} \left\{ \pm \xi \frac{j_{n}}{2(n+A_{1}+A_{2}+1)(n+A_{1}+A_{2})} + \frac{n(n+2A_{1}+2A_{2}+1)}{4} - \frac{1}{4} \left[ \lambda - 2A_{1}A_{2} - A_{1} - A_{2} + 2(m+s\mp (2A_{1}+1))\xi - \frac{m^{2}-s^{2}}{2} - s \right] \right\},$$

$$\gamma_{n} = \mp \frac{i}{\sqrt{\alpha}} \xi \frac{(n+A_{1}+A_{2}\pm s)(n+2A_{1})(n+2A_{1}+2A_{2})}{2(2n+2A_{1}+2A_{2}-1)(n+A_{1}+A_{2})}, \tag{4.23}$$

where

$$j_n = n(n+2A_1+2A_2+1)(A_1-A_2) + (A_1+A_2\pm s+1)[n(n+2A_1+2A_2+1)+(2A_1+1)(A_1+A_2)]$$

Thus, the recurrence relations for  $a_n$ 's are obtained from those for  $c_n$ 's in eq.(4.6) by dropping the infinite factor  $\frac{i}{\sqrt{\alpha}}$ . This recurrence relation agrees with those in Ref.4 for the Kerr geometries. In this limit, the Heun's equation (3.4) reduces to a confluent Heun's equation which has two regular singularities at z = 0, 1 and an irregular singularity at  $\infty$ . The angular wave function S(z) is given by using f(z) as

$$S(z) \longrightarrow \left(-\frac{i}{4\sqrt{\alpha}}\right)^{\frac{i}{2\sqrt{\alpha}}\xi} z^{A_1} (z-1)^{A_2} e^{\mp 2\xi z} f(z) . \tag{4.24}$$

The Reissner-Nordström-de Sitter limit:  $a \rightarrow 0$ 

<span id="page-12-0"></span>In this limit, we have to take  $\alpha \to 0$  simultaneously, in order to keep  $\Lambda = 3\alpha/a^2$  fixed. Then, parameters become the following forms;

$$\rho_{+} = A_{1} + A_{2} \pm s + 1 , \quad \rho_{-} = A_{1} + A_{2} \pm \sqrt{\frac{3}{\Lambda}}\omega + 1 , \quad A_{3} = \pm \frac{i}{2}\sqrt{\frac{3}{\Lambda}}\omega , 
u = \frac{-i}{4a}\sqrt{\frac{3}{\Lambda}}\left(\lambda - 2A_{1}A_{2} - A_{1} - A_{2} - \frac{m^{2} - s^{2}}{2} - s\right) , \qquad (4.25)$$

$$z_{s} = \frac{-i}{4a}\sqrt{\frac{3}{\Lambda}} , \qquad z_{\infty} = \frac{-i}{2a}\sqrt{\frac{3}{\Lambda}} .$$

Then,  $\alpha_n^{\nu}$  and  $\gamma_n^{\nu}$  are O(1), while  $\beta_n^{\nu}$  are  $O(1/a) \to \infty$ ,

$$\beta_n^{\nu} = \frac{i}{4a} \sqrt{\frac{3}{\Lambda}} \left[ n(n+2A_1+2A_2+1) - \lambda + 2A_1A_2 + A_1 + A_2 + \frac{m^2 - s^2}{2} + s \right]. \tag{4.26}$$

In order that coefficients  $a_n$  have definite limits,  $\beta_n^{\nu}$  must be 0 for a specific value of n. The eigenvalue  $\lambda$  is determined by this condition as

$$\lambda = (l - s + 1)(l + s), \tag{4.27}$$

where  $l = n + A_1 + A_2 \ge \max(|m|, |s|)$ . Therefore the angular wave function has the following form,

$$S(z) \rightarrow \frac{i}{2a} \sqrt{\frac{3}{\Lambda}} \left( \frac{i}{4a} \sqrt{\frac{3}{\Lambda}} \right)^{\pm \frac{i}{2} \sqrt{\frac{3}{\Lambda}} \omega} z^{A_1} (z-1)^{A_2}$$

$$(4.28)$$

$$\times F(-l + A_1 + A_2, l + A_1 + A_2 + 1; 2A_1 + 1; z) . \tag{4.29}$$

In this limit, the Heun's equation (3.4) reduces to a hypergeometric equation. These results coincide with that of [14].

In the former subsection, we comment that if  $\nu$  is a solution of eq.(4.13), then  $-\nu - \omega$  is a solution. For the angular solution, we chose  $\nu = 0$  and we know that  $\omega$ =integer so that this transformation of the exponent does not lead any new solution.

#### 4.3 Solution of the radial equation

By comparing eq.(4.1) with eq.(3.9), we find parameters in Heun's equation,  $\alpha, \beta, \gamma, \delta, \epsilon, q$  and  $a_H$  is given by

$$\alpha = \sigma_{+},$$
  $\beta = \sigma_{-},$   $\gamma = 2B_{1} + s + 1, \quad \delta = 2B_{2} + s + 1, \quad \epsilon = 2B_{3} + s + 1,$   $q = -v,$   $a_{H} = z_{r}.$  (4.30)

The solution convergent in the ellipse with foci at z = 0, 1 which correspond to  $r = r_-, r_+$  respectively is given by

$$g_{\nu}(z) = \sum_{n=-\infty}^{+\infty} b_n^{\nu} u_{\nu+n}(z),$$
 (4.31)

$$u_{\nu}(z) = F(-\nu, \nu + 2(B_1 + B_2 + s) + 1; 2B_1 + s + 1; z).$$
 (4.32)

where  $\omega = \gamma + \delta - 1 = 2(B_1 + B_2 + s) + 1$ . The expansion coefficients  $b_n^{\nu}$ 's are determined by the recurrence relation (4.6) in which  $c_n$ 's are replaced with  $b_n^{\nu}$ 's. The radial solution is expressed by a series where n runs from  $-\infty$  to  $+\infty$ . This is because the eigenvalue  $\lambda$ is fixed from the angular solution and thus the characteristic exponent (the renormalized momentum)  $\nu$  is determined by the transcendental equation (4.13) which guarantees the convergence of series. Remind that if  $g_{\nu}(z)$  is a solution then  $g_{-\nu-\omega}(z)$  is a solution too. the Kerr-Newman limit:  $\Lambda = (3/a^2)\alpha \to 0$ 

In this limit,  $r_{\pm}$  and  $r'_{\pm}$  are

$$r_{\pm} \longrightarrow r_{\pm}^{0} \left[ 1 + \frac{r_{\pm}^{0}(r_{\pm}^{0}{}^{2} + a^{2})}{2a^{2}(r_{\pm}^{0} - M)} \alpha \right] + O(\alpha^{2}),$$
 (4.33)

$$r'_{\pm} \longrightarrow \pm \frac{a}{\sqrt{\alpha}} \left[ 1 \pm \frac{M}{a} \sqrt{\alpha} + \frac{Q^2 - 3M^2}{2a^2} \alpha \right] + O(\alpha),$$
 (4.34)

then parameters are

$$\begin{split} B_1 &= \frac{1}{2} \left[ -s \mp i (-\tilde{\epsilon} + \tilde{\tau} + is) \right] + O(\sqrt{\alpha}) \equiv B_1^0 + O(\sqrt{\alpha}), \\ B_2 &= \frac{1}{2} \left[ -s \pm i (\tilde{\epsilon} + \tilde{\tau} - is) \right] + O(\sqrt{\alpha}) \equiv B_2^0 + O(\sqrt{\alpha}), \\ B_3 &= \mp \frac{ia\omega}{2\sqrt{\alpha}} + \frac{1}{2} \left[ -s \pm i (-\tilde{\epsilon} - is) \right] + O(\sqrt{\alpha}), \\ \sigma_+ &= B_1^0 + B_2^0 + s + 1 \pm i [-\tilde{\epsilon} - is] + O(\alpha), \\ \sigma_- &= \mp \frac{ia\omega}{\sqrt{\alpha}} + B_1^0 + B_2^0 + s + 1 + O(\sqrt{\alpha}), \\ v &= -\frac{a}{2\sqrt{\alpha}(r_+^0 - r_-^0)} \left\{ -\lambda + 2s - 2iseQ + 2B_1^0 B_2^0 + (s + 1) \left( B_1^0 + B_2^0 \right) \right. \\ &\left. - i\epsilon\tilde{\kappa} \left[ s \pm (2B_1^0 + s + 1) \right] + 2\epsilon^2 - \epsilon^2\tilde{\kappa} + eQ\epsilon(\tilde{\kappa} - 3) - \epsilon\tilde{q} - \frac{1}{2} \left( \tilde{\epsilon}^2 + \tilde{\tau}^2 - 2is\tilde{\epsilon} \right) \right\} \\ &\left. + O(\sqrt{\alpha}), \\ z_r &= \frac{a}{2\sqrt{\alpha}(r_+^0 - r_-^0)} + O(1), \qquad z_\infty = \frac{a}{\sqrt{\alpha}(r_+^0 - r_-^0)} + O(1), \end{split} \tag{4.36}$$

where  $\tilde{\epsilon} = 2M\omega - eQ = \epsilon - eQ$ ,  $\epsilon = 2M\omega$ ,  $\tilde{\kappa} = \sqrt{1 - \frac{a^2 + Q^2}{M^2}}$ ,  $\tilde{q} = \frac{am + Q^2\omega}{M}$  and  $\tilde{\tau} = \frac{\tilde{\epsilon} - \tilde{q}}{\tilde{\kappa}}$ . The parameter  $\omega$  in the above expressions is used for the frequency. The singularity at

 $\infty$  becomes irregular because of  $z_r \to \infty$ . Thus the radial equation (3.9) becomes to a confluent Heun's equation as the angular one does. In the Kerr limit ( $\Lambda \to \infty$  and Q = 0), coefficients of the recurrence relation (4.6) coincide with those in Ref.6, if we divide them by  $-\frac{a}{2\sqrt{\alpha}(r_+^0 - r_-^0)}$ . the de Sitter limit: M, Q and  $a \to 0$ 

We rearrange the locations of the singularities so that  $z = 0, 1, z'_r$ , and  $\infty$  correspond to  $r = r'_+, r_-, r_+$  and  $r'_-$  respectively;

$$z = \frac{(r_{-} - r'_{-})}{(r_{-} - r'_{+})} \frac{(r - r'_{+})}{(r - r'_{-})},$$
(4.37)

thus  $z'_r = \frac{(r_- - r'_-)}{(r_- - r'_+)} \frac{(r_+ - r'_+)}{(r_+ - r'_-)}$ . The effect of this rearrangement is only to replace  $r_-, r_+, r'_+$  and  $r'_-$  with  $r'_+, r_-, r_+$  and  $r'_-$ , respectively, in the radial equation. In this limit we find

$$r_{\pm} \to \pm ia, \qquad r'_{\pm} \to \pm \sqrt{\frac{3}{\Lambda}},$$
 (4.38)

and

$$B_{1} = \frac{1}{2} \left[ -s \mp i \left( \sqrt{\frac{3}{\Lambda}} \omega + is \right) \right],$$

$$B_{2} = \frac{1}{2} \left[ -s \pm (m+s) \right], \qquad B_{3} = \frac{1}{2} \left[ -s \pm (-m+s) \right],$$

$$\sigma_{\pm} = B_{1} + B_{2} + B_{3} + 2s + 1 + \frac{1}{2} \left[ -s \pm i \left( \sqrt{\frac{3}{\Lambda}} \omega - is \right) \right],$$

$$v = -\lambda + 2s - (s+1)(2s+1) - 2B_{1}(B_{2} + B_{3}) - (s+1)(2B_{1} + B_{2} + B_{3}) + is\sqrt{\frac{3}{\Lambda}} \omega,$$

$$z'_{r} = 1, \qquad z_{\infty} = -1.$$

$$(4.39)$$

The resulting radial equation is the following with the three singularities,

$$\left\{ \frac{d^2}{dz^2} + \left[ \frac{2B_1 + s + 1}{z} + \frac{2(B_1 + B_2 + s + 1)}{z - 1} \right] \frac{d}{dz} + \frac{\sigma_+ \sigma_- + v}{z(z - 1)^2} \right\} g(z) = 0.$$
(4.40)

By using  $\lambda = (l - s + 1)(l + s)$  (4.27) and setting  $g(z) = (z - 1)^{-B_2 - B_3 - s - \frac{1}{2} + (l + \frac{1}{2})} g_{\text{deS}}(z)$ , this equation reduces to a hypergeometric equation thus we obtain

$$R(z) = z^{B_1}(z-1)^{l-s}(z+1)^{2s+1}g_{deS}(z)$$

$$= z^{B_1}(z-1)^{l-s}(z+1)^{2s+1}F(l+1-s,l+1+i\sqrt{\frac{3}{\Lambda}}\omega;1-s+i\sqrt{\frac{3}{\Lambda}}\omega;z).(4.41)$$

This radial function coincides with that in Ref.13 except for  $\omega \to -\omega$  because of a difference between our choice of the null tetrad and theirs.

# 5 Conclusions and discussions

In this note, we have shown that equations for perturbations of the Kerr-de Sitter black hole for massless particles with any spins are reduced to the Heun's equation and thus we can obtain solutions in the form of series of hypergeometric functions. This method has been extended to the cases of Kerr-Newman-de Sitter geometries except for photons and gravitons for which the technique to separate angular and radial parts are not succeeded so far. We have explicitly constructed solutions for angular and radial equations and examined some limiting cases. We have shown that these limits reproduce the former known results.

The radial solutions we obtained in the text are valid inside a ellipse with foci at z = 0, 1 and are not valid at de Sitter horizon. In order to apply for some physical situations, we have to construct solutions which are valid in the entire regions of r and satisfy some specific boundary condition. For this, we should obtain solutions valid around the de Sitter horizon and match solutions with different convergence regions in the region where both solutions are convergent. This procedure was made for the Kerr black hole case in Ref.6,7. Also, it may be interesting to see the Teukolsky-Starobinsky identities analytically. We will answer these questions in the future publication.

# References

- [1] B.Carter, Comm. Math. Phys. 10, 280 (1968).
- [2] S.A.Teukolsky, Astorophys. J. 185, 635 (1973).
- [3] B.F.Whiting, J. Math. Phys. 30,1301 (1989).
- [4] E.D.Fackerell and R.G.Grossman, J. Math. Phys. 18,1849 (1977).
- [5] E.W.Leaver, J. Math. Phys. 27,1238 (1986).
- [6] S.Mano, H.Suzuki and E.Takasugi, Prog. Theor. Phys. 96,1079 (1996).
- [7] S.Mano and E.Takasugi, Prog. Theor. Phys. 97,213 (1997).
- [8] V.S.Otchik, in Quantum Systems; Newtends and Methods, edited by A.I.Barnt et.al.(World Scientific,1995), p.154.
- [9] T.Tanaka, H.Tagoshi and M.Sasaki, Prog. Theor. Phys. 96, 1087 (1996).
- [10] H.Tagoshi, S.Mano and E.Takasugi, Prog. Theor. Phys. 98,829 (1997).
- [11] A.Erd´elyi, Quart.J.of Math. 15, 62 (1944).
- [12] Heun's Differential Equations edited by A.Ronveaux (Oxford Science Publications).
- [13] H.Suzuki and E.Takasugi,Mod. Phys. Lett., A11, 431 (1996).
- [14] J.N.Goldberg, A.J.Macfarlane, E.T.Newman, F.Rohrlich and E.C.G,Sudarshan, J. Math. Phys. 8,2155 (1967).