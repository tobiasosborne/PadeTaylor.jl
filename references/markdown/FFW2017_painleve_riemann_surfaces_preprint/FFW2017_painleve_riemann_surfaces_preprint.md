## Methods for the computation of the multivalued Painlev´e transcendents on their Riemann surfaces

Marco Fasondini[1](#page-0-0) , Bengt Fornberg[2](#page-0-1) and J.A.C. Weideman[3](#page-0-2)

Abstract. We extend the numerical pole field solver (B. Fornberg and J.A.C. Weideman, J. Comput. Phys. 230:5957–5973, 2011) to enable the computation of the multivalued Painlev´e transcendents, which are the solutions to the third, fifth and sixth Painlev´e equations, on multiple sheets of their Riemann surfaces. We display, for the first time we believe, solutions to these equations on multiple Riemann sheets. We also provide numerical evidence for the existence of solutions to the sixth Painlev´e equation that have pole-free sectors, known as tronqu´ee solutions.

Key words. Painlev´e transcendents, PIII equation, P<sup>V</sup> equation, PVI equation, Pad´e approximation, Riemann surface

# 1 Introduction

The Painlev´e equations are six second order nonlinear ODEs, denoted by PI–PVI, that were singled out in the early 1900s because their solutions, known as the Painlev´e transcendents, have a special singularity structure. A singularity of an ODE solution is movable if its location depends on the initial conditions (ICs). The Painlev´e equations possess the Painlev´e property, which means that their solutions are free from movable branch point singularities. However, some Painlev´e transcendents may have essential singularities or branch points at certain fixed locations as well as movable poles, the latter by far the most common type of singularity. Furthermore, the Painlev´e transcendents, as their name implies, generally cannot be expressed in terms of other known functions such as the elementary functions or the classical special functions. The Painlev´e equations have not only attracted attention because of their analytical properties but also for their appearance in many applications, particularly in mathematical physics. For example, PIII, P<sup>V</sup> and PVI feature in, respectively, the quantum sine-Gordon model [\[22\]](#page-19-0), interaction models of fermions [\[7\]](#page-18-0), and the Ising model [\[13\]](#page-19-1). The history, properties and applications of the Painlev´e transcendents are discussed in more detail in [\[3,](#page-18-1) [6,](#page-18-2) [12\]](#page-19-2).

The pole field solver (PFS) presented in [\[9\]](#page-19-3) was the first numerical method that enabled the efficient and accurate computation of the pole fields of the Painlev´e transcendents on extended regions in the complex plane. The PFS was used to survey the solution space of the Painlev´e equations with single-valued solutions, viz. P<sup>I</sup> [\[9\]](#page-19-3), PII [\[10,](#page-19-4)[11\]](#page-19-5) and PIV [\[24](#page-19-6)[–26\]](#page-20-0). In these surveys new solution features and pole field patterns were discovered in solutions that had not been studied before. Our aims are to extend the PFS to the computation of the generally multivalued solutions of PIII, P<sup>V</sup> and PVI and to eventually survey, if not the entire solution spaces, then at least certain classes of PIII, P<sup>V</sup> and PVI solutions.

For a discussion of computational methods for the Painlev´e equations that preceded the PFS, we refer to [\[9\]](#page-19-3). After the introduction of the PFS, another method for computing Painlev´e transcendents was presented by Abramov and Yukhno [\[1\]](#page-18-3). Their method avoids singularities by making certain changes of variables in the neighborhoods of poles. While [\[1\]](#page-18-3) and the earlier

<span id="page-0-0"></span><sup>1</sup> University of the Free State, Department of Mathematics and Applied Mathematics, Bloemfontein 9300, South Africa (fasondinim@ufs.ac.za).

<span id="page-0-1"></span><sup>2</sup>University of Colorado, Department of Applied Mathematics, 526 UCB, Boulder, CO 80309, USA (fornberg@colorado.edu).

<span id="page-0-2"></span><sup>3</sup>University of Stellenbosch, Department of Mathematical Sciences, Private Box X1, Matieland 7602, South Africa (weideman@sun.ac.za).

methods, viz. the pole vaulting method [4] and the method in [30], should in theory be capable of computing Painlevé solutions over extended regions of the complex plane, no such results have been presented for either the multivalued or single-valued Painlevé transcendents, most likely due to a combination of cost and complexity.

The following is an outline of the paper. First we show how the PFS can be extended to accommodate the singularity structure of  $P_{\rm III}$  and  $P_{\rm V}$  solutions. The efficiency and accuracy of our numerical methods are then tested experimentally. This is followed by illustrations of different approaches to the computation of  $P_{\rm VI}$  solutions. We then use these methods to illustrate more examples of  $P_{\rm III}$ ,  $P_{\rm V}$  and  $P_{\rm VI}$  solutions on multiple Riemann sheets after which we make some concluding remarks.

### 2 Computing the multivalued Painlevé transcendents

The  $P_{III}$ ,  $P_{V}$  and  $P_{VI}$  equations are as follows:

$$\begin{aligned} & P_{\text{III}}: & \frac{d^{2}u}{dz^{2}} = \frac{1}{u} \left(\frac{du}{dz}\right)^{2} - \frac{1}{z} \frac{du}{dz} + \frac{\alpha u^{2} + \beta}{z} + \gamma u^{3} + \frac{\delta}{u}, \\ & P_{\text{V}}: & \frac{d^{2}u}{dz^{2}} = \left(\frac{1}{2u} + \frac{1}{u - 1}\right) \left(\frac{du}{dz}\right)^{2} - \frac{1}{z} \frac{du}{dz} + \frac{(u - 1)^{2}}{z^{2}} \left(\alpha u + \frac{\beta}{u}\right) + \frac{\gamma u}{z} + \frac{\delta u(u + 1)}{u - 1}, \end{aligned}$$

and

$$P_{VI}: \frac{d^{2}u}{dz^{2}} = \frac{1}{2} \left( \frac{1}{u} + \frac{1}{u-1} + \frac{1}{u-z} \right) \left( \frac{du}{dz} \right)^{2} - \left( \frac{1}{z} + \frac{1}{z-1} + \frac{1}{u-z} \right) \frac{du}{dz} + \frac{u(u-1)(u-z)}{z^{2}(z-1)^{2}} \left( \alpha + \frac{\beta z}{u^{2}} + \frac{\gamma(z-1)}{(u-1)^{2}} + \frac{\delta z(z-1)}{(u-z)^{2}} \right),$$

where  $\alpha, \beta, \gamma$  and  $\delta$  are arbitrary constants. Our computational approaches to these equations are determined by the possible locations of the branch points of their solutions. We know from the Painlevé property that the solutions to these equations are free from movable branch points. This implies that branch points can only occur at the fixed singularities of the equations. We see that in the finite plane  $P_{III}$  and  $P_{V}$  both have a fixed singularity at z=0 and  $P_{VI}$  has fixed singularities at z=0 and z=1. Hence, our computational methods for  $P_{III}$  and  $P_{V}$  are the same while  $P_{VI}$  requires a different approach.

## 2.1 Computing $P_{III}$ and $P_V$ solutions

For  $P_{III}$  and  $P_V$  one can use an exponential transformation to map the fixed singularity at z=0 out of the finite complex plane and then obtain equations whose solutions are meromorphic and thus single-valued. Specifically, setting  $z=e^{\zeta/2}$  and  $u(z)=e^{-\zeta/2}w(\zeta)$  in  $P_{III}$  results in a modified third Painlevé equation,

$$\widetilde{P}_{\text{III}}: \frac{d^2w}{d\zeta^2} = \frac{1}{w} \left(\frac{dw}{d\zeta}\right)^2 + \frac{1}{4} \left(\alpha w^2 + \gamma w^3 + \beta e^{\zeta} + \frac{\delta e^{2\zeta}}{w}\right),$$

whose solutions are meromorphic [15]. Likewise, setting  $z = e^{\zeta}$  and  $u(z) = w(\zeta)$  in  $P_V$  yields a modified fifth Painlevé equation,

$$\widetilde{P}_{V}: \qquad \frac{d^{2}w}{d\zeta^{2}} = \left(\frac{1}{2w} + \frac{1}{w-1}\right) \left(\frac{dw}{d\zeta}\right)^{2} + (w-1)^{2} \left(\alpha w + \frac{\beta}{w}\right) + \gamma e^{\zeta}w + \frac{\delta e^{2\zeta}w(w+1)}{w-1},$$

whose solutions are meromorphic [14]. One may thus obtain multivalued  $P_{III}$  and  $P_V$  solutions u in the z-plane by computing the single-valued solutions w of  $\widetilde{P}_{III}$  and  $\widetilde{P}_V$  in the  $\zeta$ -plane. We now give a brief overview of the PFS, see [9] for details, before describing how we extend it to allow for the computation of  $\widetilde{P}_{III}$  and  $\widetilde{P}_V$  solutions in the  $\zeta$ -plane.

#### <span id="page-2-1"></span>2.1.1 The PFS

The PFS computes pole fields using a two-stage approach:

Stage 1: A node set is generated on the region where the solution is to be computed; for the single-valued Painlevé transcendents, a uniform two-dimensional coarse grid on a rectangle in the complex plane was used. Starting from the point where the ICs w and w' are supplied, the Taylor coefficients of the solution are generated recursively. Then the coefficients of the type  $(\nu, \nu)$  Padé approximant are computed, where  $\nu = n/2$  and n is even, i.e., the coefficients  $a_0, \ldots, a_{\nu}$  and  $b_1, \ldots, b_{\nu}$  are calculated such that

<span id="page-2-0"></span>
$$w(\zeta + h) = \frac{a_0 + a_1 h + \dots + a_{\nu} h^{\nu}}{1 + b_1 h + \dots + b_{\nu} h^{\nu}} + \mathcal{O}\left(h^{n+1}\right). \tag{1}$$

A target node is randomly chosen and the Padé approximant is evaluated in five directions at a constant distance |h| from the current point: one pointing straight at the target node as well as 15° and 30° on either side of this direction. The Padé step is taken in the direction in which the modulus of the solution is smallest among the five directions; this is to prevent the loss of numerical stability that may occur close to a pole. Padé steps are taken in this fashion until the target node is within a distance of |h|; all the while w, w' and the Padé coefficients are stored along the path. A different target node is chosen randomly and, starting from the closest point where the Padé coefficients are available, Padé steps are taken along the minimum modulus directions until the |h|-neighborhood of that target node is reached. This is repeated until paths have been run to the neighborhoods of all the nodes of the Stage 1 node set.

Stage 2: The solution is computed at all points on a fine grid. This is accomplished as follows. Let  $\zeta_i$ ,  $1 \leq i \leq N$  denote the points where Padé coefficients are available from Stage 1. Padé steps are taken from each  $\zeta_i$  to the points on the fine grid to which it is the closest point among the  $\zeta_i$ . As mentioned in [9], these Padé steps can be taken very rapidly by vectorizing the evaluation of the Padé approximants with software such as MATLAB.

It should be noted that the conversion at each step from a recursively generated Taylor expansion to Padé rational form (1), proposed in [30], is crucial for both stages. Without it, step sizes would be limited to a small fraction of each local Taylor expansion's radius of convergence, vastly adding to the computational cost.

#### 2.1.2 PFS enhancements

As mentioned above, a uniform Stage 1 node set and constant Padé step lengths sufficed for the computation of the single-valued solutions of  $P_{II}$ ,  $P_{II}$  and  $P_{IV}$ . However, as we shall illustrate below, the solutions of  $\widetilde{P}_{III}$  and  $\widetilde{P}_{V}$  (as well as solutions of  $P_{III}$  and  $P_{V}$ ) generally have highly non-uniform pole densities. Therefore two enhancements of the PFS are required: a variable density Stage 1 node set and a variable step size Padé method.

The non-uniform Stage 1 node set should reflect at least the general trend in the pole densities of  $\widetilde{P}_{III}$  and  $\widetilde{P}_{V}$  solutions. We have not found any results concerning the pole densities of  $P_{III}$  and  $P_{V}$  solutions (which could be translated to those of  $\widetilde{P}_{III}$  and  $\widetilde{P}_{V}$  solutions) but

our experiments indicate that they are not only highly non-uniform functions of  $\zeta$  but also of the parameters and ICs. This makes it difficult to choose a Stage 1 node set that conforms to the pole densities of all  $\widetilde{P}_{III}$  and  $\widetilde{P}_{V}$  solutions. Nevertheless, it is to be expected from the exponential transformations used to arrive at  $\widetilde{P}_{III}$  and  $\widetilde{P}_{V}$  that the pole density will increase rapidly on the region Re  $\zeta \gg 0$ . For simplicity we have therefore chosen Stage 1 node sets with a node separation function  $R(\zeta)$  that decreases linearly with Re  $\zeta$ . The node separation function specifies the distance from a node at  $\zeta$  to its neighboring nodes. For example, the first column of Figure 1 shows a Stage 1 node set in a rectangular domain on the  $\zeta$ -plane with a node separation function given by  $R(\zeta) = (8-\text{Re }\zeta)/20$ . Our node sets are generated using the node placement algorithm introduced in [8].

We consider two variable step size methods, both of which are applicable only to the first stage of the PFS: one in which the step sizes vary in proportion to the node separation function of the Stage 1 grid, which we refer to as the prescribed step size method, and an adaptive step size method. For the latter method we need an estimate of the error incurred by a Padé step. Denote the numerator polynomial in (1) by p(h), thus  $p(h) = a_0 + a_1 h + \cdots + a_{\nu} h^{\nu}$ , and the denominator polynomial by q(h) and assume u has the formal power series expansion  $w(\zeta + h) = \sum_{k=0}^{\infty} c_k h^k$  (it is explained in [9] how the ODE can be used to generate the coefficients  $c_k$ ). Then it follows from (1) that

$$q(h)w(\zeta+h) - p(h) = \sum_{k=n+1}^{\infty} \epsilon_k h^k, \qquad \epsilon_k = c_k + \sum_{r=1}^{\nu} b_r c_{k-r}.$$

We estimate the relative error of the Padé step by

$$\left| \frac{w(\zeta + h) - p(h)/q(h)}{p(h)/q(h)} \right| = \left| \frac{1}{p(h)} \sum_{n=1}^{\infty} \epsilon_r h^r \right| \approx \left| \frac{\epsilon_{n+1} h^{n+1}}{p(h)} \right| := T(h).$$

Suppose T(h) is greater than some specified tolerance Tol, then we must rescale the step size, h := qh, and find the scaling factor q such that  $T(qh) \leq Tol$ . We have that

$$Tol \ge T(qh) = q^{n+1} \frac{p(h)}{p(qh)} T(h) \approx q^{n+1} T(h),$$

and since we have made two crude approximations we choose q conservatively, i.e.,

<span id="page-3-0"></span>
$$q = \left(\frac{k \cdot Tol}{T(h)}\right)^{1/(n+1)},\tag{2}$$

where k is a small positive constant.

We use the variable step size method only in Stage 1 of the PFS, as follows. After a Padé step is taken in the minimum modulus direction in the manner described above, the error is estimated using T(h). If T(h) > Tol, then h is replaced by qh, the solution is again computed in the five directions, the minimum modulus direction is found again and T(h) is computed. This is repeated until  $T(h) \le Tol$ . If  $T(h) \le Tol$ , then w, w', the Padé coefficients and the scaled step length |qh| are stored at the point. The initial step length is always the scaled step length stored at the current point; if the current point is the initial point, a user-specified step length is used. Variable step size Padé steps are taken in this manner until the target node is within a distance of |qh|. Paths are run in this fashion to reach close to all the Stage 1 nodes after which Stage 2 is implemented as described above. As an example, the second column of Figure 1 shows the Padé steps taken by the adaptive step size method in the first stage of the PFS. At each of the 2701 points in the second column, w, w' and the Padé coefficients

are available and each of the 3041 nodes in the first column is within a distance of |qh| to a point in the second column. To compute the solution shown in the third column, we used a  $571\times140$  uniform grid with spacing 1/15 on the  $\zeta$ -plane domain as our Stage 2 fine grid. Thus, the Padé coefficients at the 2701 points in the second column were used to compute Padé approximations to the solution w at the  $571\times140=79940$  points of the fine grid.

If we use the prescribed step size method, the length of each Padé step in Stage 1 is  $|h(\zeta)| = cR(\zeta)$ , where c is a positive constant and  $R(\zeta)$  is the node separation function. Paths are run to within a distance  $|h(\zeta)|$  of each Stage 1 node. Otherwise, the implementation of the prescribed step size and adaptive step size methods are the same.

![](_page_4_Figure_2.jpeg)

<span id="page-4-0"></span>Figure 1: A  $P_{III}$  solution on three sheets of its Riemann surface, as computed by our enhanced PFS method. Left to right: The Stage 1 node set with node separation function  $R(\zeta) = (8-\text{Re }\zeta)/20$ ; the Padé steps taken by the adaptive step size method in Stage 1; a  $P_{III}$  solution in the  $\zeta$ -plane and on the corresponding sheets of its Riemann surface in the z-plane  $(z=e^{\zeta/2})$ . This solution has parameter values  $\alpha=-1/2=\beta,\ \gamma=1=-\delta$  and ICs (u(1),u'(1))=(1/4,1) in the z-plane.

The third column of Figure 1 is a plot of the modulus of a  $P_{III}$  solution, i.e., a plot of  $|e^{-\zeta/2}w(\zeta)| = |u|$ , where w is the solution computed on the Stage 2 grid mentioned above. We map the  $\tilde{P}_{III}$  solution w on the strip  $-2\pi + 4\pi s$  < Im  $\zeta \le 2\pi + 4\pi s$ ,  $s \in \mathbb{Z}$  to the  $P_{III}$  solution u on the s-th sheet of the Riemann surface in the z-plane according to  $u(z) = u(e^{\zeta/2}) = e^{-\zeta/2}w(\zeta)$ . The strips with s = -1, 0, 1 are indicated by dashed horizontal lines in Figure 1. The fourth

column is another plot of |u|, but mapped to the s=-1,0,1 sheets of the Riemann surface. According to Table 1 all the poles of the solution in the figure are first order with residue  $+\sqrt{\gamma}=+1$  or  $-\sqrt{\gamma}=-1$  in the z-plane, indicated by red and yellow circles, respectively. Similarly, all the zeros are simple and u' either has the value  $\sqrt{-\delta}=+1$  (purple squares) or  $-\sqrt{-\delta}=-1$  (light blue squares) at each zero in the z-plane. Note how the lengths of the Padé steps taken in the second column by the adaptive step size method conform to the pole density of the solution in the third column. The modulus of the solution in the figure has an up-down symmetry in the  $\zeta$ -plane since the solution in the upper and lower half planes are conjugate; this is a consequence of the real valued parameters and real ICs on the real  $\zeta$ -axis. The distinctive spirals of poles whirling around the branch point z=0 in the fourth column is a pole field pattern that has not been observed before: the pole fields of the single-valued Painlevé transcendents shown in [9–11,25,26] have very different characteristics.

<span id="page-5-0"></span>Table 1: The labels indicating the poles and zeros of the  $P_{III}$  solution in Figure 1. For every solution displayed in this paper we give a table similar to this one to describe the types of poles and zeros of the solution. The coefficients  $c_k$  in this table and the similar Tables 3–6 are derived as follows. In a neighborhood of a pole or zero at  $z_0$  we have  $u \approx c_k(z-z_0)^k$  with k < 0 or k > 0, respectively. Making this substitution in the relevant equation and taking the limit  $z \to z_0$  readily yields the order k of the pole or zero as well as the leading order coefficient  $c_k$ . We assume that  $z_0$  is not a fixed singularity of the equation. The coefficients  $c_k$  in the transformed plane, e.g., the  $\zeta$ -plane in Figure 1, can be derived similarly or by applying the appropriate transformation to the poles and zeros in the z-plane.

|                                     |                | Poles                                     | Zeros |                                        |  |  |
|-------------------------------------|----------------|-------------------------------------------|-------|----------------------------------------|--|--|
| $P_{\rm III}, \gamma \delta \neq 0$ | z-plane        | $c_{-1} = +1/\sqrt{\gamma}$               | •     | $c_1 = +\sqrt{-\delta}$                |  |  |
|                                     |                | $c_{-1} = -1/\sqrt{\gamma}$               | 0     | $c_1 = -\sqrt{-\delta}$                |  |  |
|                                     | $\zeta$ -plane | $c_{-1} = +2e^{-\zeta_0/2}/\sqrt{\gamma}$ | •     | $c_1 = +e^{\zeta_0/2}\sqrt{-\delta}/2$ |  |  |
|                                     |                | $c_{-1} = -2e^{-\zeta_0/2}/\sqrt{\gamma}$ | 0     | $c_1 = -e^{\zeta_0/2}\sqrt{-\delta}/2$ |  |  |

One could compute  $P_{III}$  and  $P_V$  solutions in the z-plane instead of the  $\zeta$ -plane by making the Padé steps run around the possible branch point at z=0 in clockwise and counterclockwise directions (and thus onto the  $s \leq 0$  and  $s \geq 0$  sheets, respectively). We shall illustrate this approach for the computation of  $P_{VI}$  solutions in section 2.2. We have found that the pole densities of  $P_{III}$  and  $P_V$  solutions can also be highly non-uniform in the z-plane, as Figure 1 illustrates for  $P_{III}$  (solutions of  $P_V$  in the z-plane will be shown in section 3). Thus, variable density Stage 1 node sets and a variable step size Padé method are also required in the z-plane. We have not found that there is any advantage, in terms of accuracy or speed, to computing  $P_{III}$  and  $P_V$  solutions in the z-plane as opposed to the  $\zeta$ -plane. If anything, the implementation of our method in the z-plane is more complicated because of the need to impose certain directions on the integration paths.

#### <span id="page-5-1"></span>2.1.3 Experiments

In practice we choose specific parameters in our numerical method (e.g.,  $R(\zeta)$ , the value of c for the prescribed step sizes  $|h(\zeta)| = cR(\zeta)$ , the value of k in (2), the order n of the Padé approximations, see (1)) based on experimentation for which error estimates are essential.

Solutions of <sup>P</sup>eIII and <sup>P</sup>e<sup>V</sup> with real parameter values and real ICs on the real axis satisfy w(ζ) = w(ζ). If we do not make use of this symmetry in our numerical method but, instead, compute the Pad´e steps in the upper and lower halves of the ζ-plane independently, as is the case in the second column of Figure [1,](#page-4-0) then we can estimate the relative numerical error by calculating E(ζ) = |w(ζ) − w(ζ)|/|w(ζ)|. By plotting E(ζ) and recording the execution time for different choices of the parameters we can experimentally determine parameter choices that give satisfactory results.

In Table [2](#page-6-1) we compare the performance of the adaptive and prescribed step size methods for two different situations: when high accuracy is required (Experiment 1), e.g., for the verification of theoretical asymptotic formulae and known closed-form solutions, and when efficiency is essential (Experiment 2), e.g., in a survey of a large number of solutions. As expected, the adaptive step size method is more accurate but slower than the prescribed step size method since it incorporates error control. The estimated error is orders of magnitude greater than T ol for the adaptive step size method since error control is only applied in Stage 1 but not in Stage 2. In addition, T ol is only a bound on the local error which accumulates with the number of Pad´e steps on the domain, which also accounts for the fact that the error on sheets +1 and −1 is larger than on the 0th sheet. The efficiency of both methods can be increased by using the up-down symmetry in the ζ-plane (for real-valued ICs and parameter values) and by parallelizing the implementation of Stage 1 and Stage 2. Stage 1 can be parallelized by partitioning the domain and the parallelization of Stage 2 is trivial since the evaluations of the Pad´e approximations from the available Pad´e coefficients are independent.

<span id="page-6-1"></span>Table 2: Statistics for the computation of the solution shown in Figure [1.](#page-4-0) In Experiment 1 and Experiment 2 both methods used a Stage 1 node set with a node separation function given by R1(ζ) = (8 − Re ζ)/40 and R2(ζ) = (8 − Re ζ)/10, respectively, and both used 30th order Pad´e steps. The adaptive step size method used k = 10<sup>−</sup><sup>3</sup> (see [\(2\)](#page-3-0)) in both experiments and T ol = 10<sup>−</sup><sup>14</sup> in Experiment 1 and T ol = 10<sup>−</sup><sup>11</sup> in Experiment 2. The step sizes of the prescribed step size method were |h(ζ)| = 1.85R1(ζ) in Experiment 1 and |h(ζ)| = 0.6R2(ζ) in Experiment 2. The third column refers to the number of Pad´e steps taken in Stage 1; the fine grid used for Stage 2, a 571×140 uniform grid with spacing 1/15 on the ζ-plane domain, were the same in both experiments. The execution times include Stage 1 and Stage 2 and were recorded on a machine with a clock speed of 3.6 GHz using only one core. The first and second relative errors in each row are for the solution on the regions −2π < Im ζ ≤ 2π (corresponding to the 0th sheet) and |Im ζ − 4π| ≤ 2π (−1th and +1th sheets), respectively. These errors were estimated using the symmetry-based method discussed above.

|              |                       | Number of steps | Time (seconds) | Relative error |
|--------------|-----------------------|-----------------|----------------|----------------|
| Experiment 1 | Adaptive step sizes   | 4259            | 5.43           | 4e-10, 4e-8    |
|              | Prescribed step sizes | 4270            | 4.52           | 8e-6, 4e-3     |
| Experiment 2 | Adaptive step sizes   | 1324            | 2.27           | 1e-6, 8e-2     |
|              | Prescribed step sizes | 1314            | 1.48           | 2e-6, 4e-2     |

### <span id="page-6-0"></span>2.2 Computing PVI solutions

An exponential transformation mapped the fixed singular point of PIII and P<sup>V</sup> at z = 0 out of the finite plane since the exponential function is entire and never assumes the value zero. This allowed us to avoid the possible branch point of  $P_{\rm III}$  and  $P_{\rm V}$  solutions at z=0. In general we cannot avoid the branch points of  $P_{\rm VI}$  solutions because the  $P_{\rm VI}$  equation has fixed singular points at z=0 and z=1 and, by Picard's theorem, no non-constant entire function can avoid two different values, as would be required to map these points out of the finite plane. Thus it is generally necessary to steer the integration paths of the PFS around branch points of  $P_{\rm VI}$  solutions. However, it is possible to use transformations to avoid branch points on restricted parts of the Riemann surfaces of  $P_{\rm VI}$  solutions. We illustrate both of these approaches, starting with the latter.

#### 2.2.1 Avoiding branch points

The transformation  $z=e^{\zeta}$  maps the fixed singularity of  $P_{VI}$  at z=0 out of the finite  $\zeta$ -plane and it maps the fixed singularity at z=1 to  $\zeta=2i\pi k,\ k\in\mathbb{Z}$ . These points are the fixed singularities of the equation

<span id="page-7-0"></span>
$$\frac{d^2w}{d\zeta^2} = \frac{1}{2} \left( \frac{1}{w} + \frac{1}{w-1} + \frac{1}{w-e^{\zeta}} \right) \left( \frac{dw}{d\zeta} \right)^2 - \left( \frac{e^{\zeta}}{e^{\zeta} - 1} + \frac{e^{\zeta}}{w - e^{\zeta}} \right) \frac{dw}{d\zeta} + \frac{w(w-1)(w-e^{\zeta})}{(e^{\zeta} - 1)^2} \left( \alpha + \frac{\beta e^{\zeta}}{w^2} + \frac{\gamma(e^{\zeta} - 1)}{(w-1)^2} + \frac{\delta e^{\zeta}(e^{\zeta} - 1)}{(w-e^{\zeta})^2} \right), \tag{3}$$

obtained by setting  $u(z) = w(\zeta)$ ,  $z = e^{\zeta}$  in  $P_{VI}$ . Another exponential transformation,  $\zeta = e^{\eta}$ , maps the fixed singularity at  $\zeta = 0$  out of the finite  $\eta$ -plane and it maps the remaining fixed singularities to

<span id="page-7-1"></span>
$$\eta = \log|2\pi k| + i\arg(2\pi i k), \qquad |k| \ge 1,\tag{4}$$

which are the fixed singularities of the equation obtained by setting  $w(\zeta) = v(\eta)$ ,  $\zeta = e^{\eta}$  in (3), i.e.,

<span id="page-7-2"></span>
$$\frac{d^{2}v}{d\eta^{2}} = \frac{1}{2} \left( \frac{1}{v} + \frac{1}{v-1} + \frac{1}{v-e^{e^{\eta}}} \right) \left( \frac{dv}{d\eta} \right)^{2} - \left( \frac{e^{\eta}e^{e^{\eta}}}{e^{e^{\eta}} - 1} + \frac{e^{\eta}e^{e^{\eta}}}{v - e^{e^{\eta}}} - 1 \right) \frac{dv}{d\eta} + \frac{v(v-1)(v-e^{e^{\eta}})e^{2\eta}}{(e^{e^{\eta}} - 1)^{2}} \left( \alpha + \frac{\beta e^{e^{\eta}}}{v^{2}} + \frac{\gamma(e^{e^{\eta}} - 1)}{(v-1)^{2}} + \frac{\delta e^{e^{\eta}}(e^{e^{\eta}} - 1)}{(v-e^{e^{\eta}})^{2}} \right).$$
(5)

We conclude from (4) that the region Re  $\eta < \log 2\pi$  is branch point-free. The first column of Figure 2 shows a  $P_{VI}$  solution computed within this region. This solution was computed by applying the method discussed above, i.e., the PFS method with a non-uniform node set and variable step sizes in Stage 1, to (5). The solution within the region  $-3\pi \leq \text{Im } \eta \leq -\pi$ , not shown in the figure, is the conjugate of the displayed solution contained in  $\pi \leq \text{Im } \eta \leq 3\pi$ ; as before, this follows from the real-valued parameters and real ICs on the real  $\eta$ -axis. The  $\eta$ -plane region in the figure corresponds to the following two rectangular sheets of the Riemann surface that winds around the branch point at  $\zeta = 0$  in the second column of the figure:

$$\{\zeta\in\mathbb{C}\ :\ \log(1/100)\leq\mathrm{Re}\ \zeta\leq\log(10), -\pi\leq\mathrm{Im}\ \zeta\leq\pi,\ -\pi\leq\arg\zeta\leq3\pi\}.$$

The asymptotic behaviors of the solution close to its poles and zeros are given in Table 3.

#### 2.2.2 Circumambulating branch points

Unlike the branch point at  $\zeta = 0$ , the branch point at  $\zeta = 2\pi i$  in the second column is not mapped out of the finite  $\eta$ -plane. Hence, this branch point cannot be avoided by computing in the  $\eta$ -plane. The approach we follow to compute the solution in the neighborhood of  $\zeta = 2\pi i$ 

<span id="page-8-0"></span>Table 3: The poles and zeros of the  $P_{VI}$  solutions in Figures 2, 3 and 5.

|                              |                | Poles                                                                                                                       |   | Zeros                                                                    |  |  |  |  |
|------------------------------|----------------|-----------------------------------------------------------------------------------------------------------------------------|---|--------------------------------------------------------------------------|--|--|--|--|
| $P_{VI}, \alpha\beta \neq 0$ | $\eta$ -plane  | $c_{-1} = +e^{-\eta_0} (e^{e^{\eta_0}} - 1) / \sqrt{2\alpha}$ $c_{-1} = -e^{-\eta_0} (e^{e^{\eta_0}} - 1) / \sqrt{2\alpha}$ | • | $c_1 = +e^{\eta_0} e^{e^{\eta_0}} \sqrt{-2\beta}/(e^{e^{\eta_0}} - 1)$   |  |  |  |  |
|                              |                | $c_{-1} = -e^{-\eta_0} (e^{e^{\eta_0}} - 1) / \sqrt{2\alpha}$                                                               | 0 | $c_1 = -e^{\eta_0} e^{e^{\eta_0}} \sqrt{-2\beta} / (e^{e^{\eta_0}} - 1)$ |  |  |  |  |
|                              | $\zeta$ -plane | $c_{-1} = +(e^{\zeta_0} - 1)/\sqrt{2\alpha}$ $c_{-1} = -(e^{\zeta_0} - 1)/\sqrt{2\alpha}$                                   | • | $c_1 = +e^{\zeta_0} \sqrt{-2\beta}/(e^{\zeta_0} - 1)$                    |  |  |  |  |
|                              |                | $c_{-1} = -(e^{\zeta_0} - 1)/\sqrt{2\alpha}$                                                                                | 0 | $c_1 = -e^{\zeta_0} \sqrt{-2\beta}/(e^{\zeta_0} - 1)$                    |  |  |  |  |
|                              | z-plane        | $c_{-1} = +z_0(z_0 - 1)/\sqrt{2\alpha}$                                                                                     | • | $c_1 = +\sqrt{-2\beta}/(z_0 - 1)$                                        |  |  |  |  |
|                              |                | $\begin{vmatrix} c_{-1} = +z_0(z_0 - 1)/\sqrt{2\alpha} \\ c_{-1} = -z_0(z_0 - 1)/\sqrt{2\alpha} \end{vmatrix}$              | 0 | $c_1 = -\sqrt{-2\beta}/(z_0 - 1)$                                        |  |  |  |  |

is to apply our enhanced PFS method to (3). The only modification of our method that is required is to steer the Padé paths in the appropriate directions around the branch points in the  $\zeta$ -plane, as shown in the third column of the figure. In the bottom-right frame the Stage 1 Padé paths run in clockwise and counterclockwise directions around the branch point at  $\zeta = 0$ . Note that, as required, none of the paths overstep the branch cuts, indicated by solid lines, on a given sheet. To move onto the sheet in the top-right frame the paths move through the branch cut by running only in a counterclockwise direction around  $\zeta = 0$ , as indicated by the arrows. The paths then run in both directions around the branch point at  $\zeta = 2\pi i$ . The dots on which the paths are superimposed are the points of the Stage 1 node set. We chose a node set that becomes increasingly dense close to the branch points since our numerical experiments have shown that there can be high pole densities in the neighborhoods of the branch points, which is evidently the case for the solution in the top-center frame of the figure. Note from the bottom-right frame that the adaptive step size method chose a few, relatively large steps on the pole-free sheet while the step sizes in the top-right frame reflect the increasing pole densities close to the branch points.

Figure 3 again depicts the  $P_{VI}$  solution in the  $\zeta$ -plane but also on the corresponding sheets in the z-plane. The second and third columns show phase portraits [28] of the solution, i.e., plots of  $w(\zeta)/|w(\zeta)| \in [-\pi,\pi]$  and  $u(z)/|u(z)| \in [-\pi,\pi]$ , respectively. The phase of the solution is indicated according to the color wheel at the top of the figure and so, for example, positive real solution values are indicated by red. The branch cuts in the second column are unmistakable and clearly indicate the manner in which the two sheets in the  $\zeta$ -plane are connected through the branch cut or, equivalently, the directions in which the PFS integration paths must have run. The traversal of the paths through the branch cut in the bottom-center frame corresponds to a traversal through the branch cut  $z \in (0,1)$  on the 0th sheet in the z-plane. The movement of the paths across the dotted line in the top-center frame corresponds to a movement across the branch cut on the negative real axis on the 1st sheet in the z-plane. Hence, the paths in the  $\zeta$ -plane correspond to counterclockwise revolutions around the branch points at z = 0 and z = 1 and thus the sheets in the z-plane are parametrized by

sheet 
$$k = \{z \in \mathbb{C} : 1/100 \le |z| \le 10, \ z = re^{i\theta_k}, \ z = 1 + \rho e^{i\phi_k}, \ \rho \ge 1/100\},$$

where

$$-\pi < \theta_0 \le \pi, \qquad -\pi < \phi_0 \le \pi, 
-\pi < \theta_1 \le \pi, \qquad \pi < \phi_1 \le 3\pi, 
\pi < \theta_2 \le 3\pi, \qquad \pi < \phi_2 \le 3\pi.$$

The most striking feature of the solution in Figures 2 and 3 is its large pole-free sector. Tronquée solutions are characterized by a pole-free sector on which it satisfies a divergent asymptotic expansion near infinity. Tronquée solutions of  $P_{I}$ — $P_{V}$  have been studied, see [2, 5, 17–21], but we are not aware of any results concerning tronquée solutions of  $P_{VI}$ . Figures 2 and 3 provide numerical evidence for the existence of tronquée  $P_{VI}$  solutions. The ICs of the solution, given in the caption to Figure 2, were obtained by substituting the formal expansion  $u = \sum_{n=0}^{\infty} a_n z^{-n}$  into  $P_{VI}$ , generating the coefficients  $a_n$  recursively and evaluating the expansion by optimal truncation far out on the positive real axis.

![](_page_9_Figure_1.jpeg)

<span id="page-9-0"></span>Figure 2: Two approaches to the computation of a  $P_{VI}$  solution: computing on a region in the  $\eta$ -plane that is branch point-free (first column) and computing in the  $\zeta$ -plane by making the PFS paths run around the branch points (second and third columns). This  $P_{VI}$  solution has parameters  $(\alpha, \beta, \gamma, \delta) = (4, -4, 8, -8)$  and ICs u(10) = 0.429534600325223 and u'(10) = -1.61713114374804e-3 in the z-plane. The symmetry-based error estimates for the  $\eta$ -plane method are 4e-7 and 3e-6 for the solution contained in the regions  $-\pi < \text{Im } \eta \le \pi$  and  $\pi < \text{Im } \eta \le 3\pi$ , respectively; those for the  $\zeta$ -plane method are 6e-7, for the solution shown in the bottom-centre frame, and 6e-4 and 5e-4 for the solution within the strips  $-\pi < \text{Im } \zeta \le \pi$  and  $\pi < \text{Im } \zeta \le 3\pi$ , respectively, in the top-center frame.

![](_page_10_Figure_0.jpeg)

<span id="page-10-0"></span>Figure 3: The phase of the PVI solution shows the structure of the Riemann surfaces in the ζ and z planes (recall that z = e ζ ). The phase of the solution is depicted according to the color wheel, taken from <http://dlmf.nist.gov/help/vrml/aboutcolor>. The pole-free 0th sheet in the z-plane only has a branch cut on z ∈ (0, 1) whereas the other sheets have branch cuts on z ∈ (0, 1) and the negative real z-axis.

# <span id="page-11-0"></span>3 More examples of $P_{III}$ , $P_{V}$ and $P_{VI}$ solutions

In this section we use the methods discussed above to compute more examples of the multivalued Painlevé transcendents on multiple sheets of their Riemann surfaces. In addition, we derive a condition number which shows that in pole-free regions the solution can be sensitive to perturbations. In these regions we enhance the method by solving the equation by boundary-value techniques, a remedy originally suggested in [9].

### 3.1 A tronquée P<sub>V</sub> solution

Figure 4 shows a  $P_V$  solution computed in the same manner as the  $P_{III}$  solution in Figure 1. The solution on sheets -1 and -2 (not shown) are the conjugates of the solution on sheets +1 and +2, respectively. According to Table 4, all the zeros of the solution have double multiplicity. They can be identified on the phase portraits in the third column as points around which each color is assumed twice in the order indicated by the color wheel above Figure 3 (red $\rightarrow$ yellow $\rightarrow$ green etc. for a counterclockwise traversal around the point).

It is shown in [2] that for a given set of parameters  $\alpha, \beta, \gamma, \delta$  with  $\delta \neq 0$ , there is a unique tronquée  $P_V$  solution with  $u(z) \sim -1$ ,  $z \to \infty$  for  $-\pi < \arg z < \pi$ . The solution in Figure 4 is such a solution and its ICs were obtained by substituting  $u = \sum_{n=0}^{\infty} a_n z^{-n}$  into  $P_V$  and evaluating the truncated expansion far out on the positive real axis.

|                                   |                | Poles                                                         | Zeros |                          |  |
|-----------------------------------|----------------|---------------------------------------------------------------|-------|--------------------------|--|
|                                   | z-plane        | $c_{-1} = +z_0/\sqrt{2\alpha}$ $c_{-1} = -z_0/\sqrt{2\alpha}$ | •     | $c_2 - k \in \mathbb{C}$ |  |
| $P_{V}, \alpha \neq 0, \beta = 0$ |                | $c_{-1} = -z_0/\sqrt{2\alpha}$                                | 0     | $C_2 - n \in \mathbb{C}$ |  |
| $[1, \alpha \neq 0, \beta = 0]$   | $\zeta$ -plane | $c_{-1} = +1/\sqrt{2\alpha}$                                  | •     | $c_2 = e^{2\zeta_0} k$   |  |
|                                   |                | $c_{-1} = +1/\sqrt{2\alpha}$ $c_{-1} = -1/\sqrt{2\alpha}$     | 0     |                          |  |

<span id="page-11-1"></span>Table 4: The poles and zeros of the P<sub>V</sub> solution in Figure 4.

### 3.2 A tronquée P<sub>III</sub> solution

The first and second columns of Figure 5 show a tronquée  $P_{III}$  solution that is pole-free on the region  $-3\pi/4 < \arg z < 9\pi/4$  (column 1), which in the  $\zeta$ -plane corresponds to the region  $-3\pi/2 < \operatorname{Im} \zeta < 9\pi/2$  (column 2). In the z-plane the asymptotic behavior of the solution on the pole-free region is  $u \sim \sqrt[3]{z}, z \to \infty$ , which in the  $\zeta$ -plane is equivalent to  $u \sim (e^{\zeta/2})^{1/3} = e^{\zeta/6}$ , Re  $\zeta \to +\infty$ . This solution is an example of the following set of tronquée  $P_{III}$  solutions, the existence of which is proved in [21]: for  $\gamma = 0$ ,  $\alpha = 1 = -\delta$ ,  $\beta$  arbitrary and any of the branches of  $z^{1/3}$ , there is a unique tronquée  $P_{III}$  solution with the behavior  $u \sim z^{1/3}$ ,  $z \to \infty$  on a certain region with angular width  $3\pi$ .

As discussed in [9], the accurate computation of the solution of a Painlevé equation on a pole-free region requires that it be solved as a boundary value problem (BVP) on that region. For the tronquée solutions in Figures 2 and 4, we achieved satisfactory accuracy without the use of a BVP solver. However, the tronquée solution in Figure 5 has a much larger pole-free sector and thus a BVP solver greatly improves the accuracy of the solution. Indeed, if our enhanced PFS method is used without a BVP solver to compute the tronquée  $P_{\rm III}$  solution,

then the error is on the order of 10<sup>−</sup><sup>1</sup> . If we use a BVP solver, then we achieve the much smaller errors reported in the caption of Figure [5.](#page-13-0) We now describe how we used a BVP solver to compute the tronqu´ee PIII solution, how the error of the computed solution was estimated and why the PFS method is unstable on the smooth region.

According to a result in [\[21\]](#page-19-12), in the z-plane the tronqu´ee solution in Figure [5](#page-13-0) satisfies

$$u \sim \sqrt[3]{z} \left[ 1 + \sum_{n=1}^{\infty} a_n (\sqrt[3]{z})^{-2n} \right], \qquad z \to \infty, \qquad -\frac{3\pi}{4} < \arg z < \frac{9\pi}{4}.$$

Substituting this expansion into the PIII equation, generating the coefficients a<sup>n</sup> recursively and evaluating the optimally truncated expansion (and its derivative) at two points close to the boundary of the pole-free sector at arg z = −3π/4 and arg z = 9π/4, we find that

![](_page_12_Figure_4.jpeg)

<span id="page-12-0"></span>Figure 4: A tronqu´ee P<sup>V</sup> solution in the ζ and z planes with parameters (α, β, γ, δ) = (1, 0, 1/4, −1/2) and ICs u(30) = −1.05294551349665 and u 0 (30) = 2.47019460566845e-3 in the z-plane. The third column is a phase portrait of the solution shown in the second column. The symmetry-based error estimates for the solution on sheets 0–2 are 3e-10, 7e-7 and 1e-6, respectively.

![](_page_13_Figure_0.jpeg)

<span id="page-13-0"></span>Figure 5: A tronquée  $P_{III}$  solution (left) and a plot of  $\log_{10}$  (min $\{\kappa_r, \kappa_a\}$ ), where  $\kappa_r$  is defined in (7) and  $\kappa_a = \kappa_r |\widetilde{w}''|$ ;  $\kappa_r$  and  $\kappa_a$  are interpreted as relative and absolute condition numbers of the  $\widetilde{P}_{III}$  equation, indicating its sensitivity to numerical errors. This solution has parameters  $(\alpha, \beta, \gamma, \delta) = (1, -1/20, 0, -1)$  and its ICs in the z-plane are given in (6). The error estimates for the solution on each of the strips, indicated by the dashed horizontal lines, are, from bottom to top: 4e-6, 3e-7, 2e-8, 3e-9, 4e-8.

<span id="page-14-1"></span>
$$z_{1} = 30e^{9\pi i/4 - i\pi/12}, z_{2} = 30e^{-3\pi i/4 + i\pi/12}, u(z_{1}) = -2.000735432319 + 2.376177147900i, u(z_{2}) = 2.384379236170 - 1.993845650158i, u'(z_{1}) = -5.939523100e-3 + 3.402038641e-2i, u'(z_{2}) = 6.050817704e-3 + 3.398020750e-2i, (6)$$

for the parameter values given in the caption of Figure 5. We translate these ICs at  $z_1$  and  $z_2$  to ICs for the  $\widetilde{P}_{III}$  equation at the corresponding points in the  $\zeta$ -plane: the two points indicated by crosses in the second column of Figure 5. Our enhanced PFS method is launched from these points to compute the solution everywhere on the rectangular domain except on the region inscribed on the pole-free sector in the second column. The solution values computed on the two curved boundaries of the inscribed region are used as boundary conditions for the BVP solver that is used to compute the solution on the pole-free sector. As in [9], we use the DMSUITE package [29] to implement a Chebyshev spectral collocation method as BVP solver.

We estimate the error on the pole-free region (computed with the BVP solver) and on the rest of the domain (computed with the enhanced PFS method) using different methods. The error of the PFS-computed solution in Figure 5 cannot be estimated using the symmetry-based method discussed in section 2.1.3 since the solution does not have the up-down symmetry in the  $\zeta$ -plane. Instead, we use the following method. Recall from section 2.1.1 that the PFS method selects the target nodes in Stage 1 in a random order. Hence, if we compute the same solution twice, different paths will be run in Stage 1, resulting in solutions that should differ by approximately the numerical error. We therefore use the relative difference between two independently calculated solutions as an error estimate. The error of the solution computed by the BVP solver can be estimated using the method discussed in [9]: by measuring the difference between the derivative values of the PFS-computed solution and the BVP solvercomputed solution at the boundaries. We increase the number of collocation points of the BVP solver until the derivative values match to the desired accuracy. We also use the difference between BVP solver-computed solutions with different numbers of collocation points as an error estimate. For the solution in Figure 5, we increased the number of collocation points until error estimates for the BVP solver-computed solution were smaller than the error estimate for the PFS-computed solution.

The instability of the computation of  $\widetilde{P}_{III}$  solutions as initial-value problems on pole-free regions can be demonstrated as follows<sup>4</sup>. Let the exact solution to the  $\widetilde{P}_{III}$  equation be w and let  $\widetilde{w}$  be the approximate (numerical) solution. Let  $w(\zeta) \approx \widetilde{w}(\zeta) + \epsilon$ , where  $\epsilon$  is constant. Making these substitutions in the  $\widetilde{P}_{III}$  equation, we find that

<span id="page-14-0"></span>
$$\left| \frac{w'' - \widetilde{w}''}{\widetilde{w}''} \right| \approx \left| \frac{1}{\widetilde{w}''} \left[ -\frac{(\widetilde{w}')^2}{\widetilde{w}} + \frac{1}{4} \left( 2\alpha \widetilde{w}^2 + 3\gamma \widetilde{w}^3 - \frac{\delta e^{2\zeta}}{\widetilde{w}} \right) \right] \right| \left| \frac{\epsilon}{\widetilde{w}} \right| + \mathcal{O}\left( \epsilon^2 \right)$$

$$:= \kappa_r \left| \frac{\epsilon}{\widetilde{w}} \right| + \mathcal{O}\left( \epsilon^2 \right), \qquad \kappa_r = \left| \frac{1}{\widetilde{w}''} \left[ -\frac{(\widetilde{w}')^2}{\widetilde{w}} + \frac{1}{4} \left( 2\alpha \widetilde{w}^2 + 3\gamma \widetilde{w}^3 - \frac{\delta e^{2\zeta}}{\widetilde{w}} \right) \right] \right|. (7)$$

Thus, we interpret  $\kappa_r$  as a 'relative condition number' of the  $\widetilde{P}_{III}$  equation. It gives the approximate factor with which the relative error of the solution  $(|\epsilon/\widetilde{w}|)$  is amplified to give the relative error in the evaluation of the right-hand side of the  $\widetilde{P}_{III}$  equation  $(|w'' - \widetilde{w}''|/|\widetilde{w}|)$ . Since the  $\widetilde{P}_{III}$  equation is used to generate the Taylor coefficients that are converted to Padé coefficients in the PFS method, we expect the error of the PFS-computed solution to grow rapidly on regions where  $\kappa_r$  is large. That is, we expect the PFS method to be unstable

<span id="page-14-2"></span><sup>&</sup>lt;sup>4</sup>Although we only consider this instability for the  $\widetilde{P}_{III}$  equation, it is also present in the computation of  $P_{III}$  solutions in the z-plane. Recall that we obtain  $P_{III}$  solutions u by computing  $\widetilde{P}_{III}$  solutions w (with  $w = e^{\zeta/2}u$ ).

where  $\kappa_r$  is large. On the pole-free sector of the solution in the second column of Figure 5,  $\widetilde{w} \approx w = e^{\zeta/2} u \sim e^{\zeta/2} e^{\zeta/6} = e^{2\zeta/3}$ , Re  $\zeta \to +\infty$ , in which case  $\kappa_r$  simplifies to

$$\kappa_r \sim \frac{27}{16} e^{2\text{Re }\zeta/3}, \quad \text{Re }\zeta \to +\infty,$$

for the parameter values of the solution in Figure 5. This shows that computation on the smooth region is exponentially unstable. For the  $P_{\rm III}$  solution in Figure 5, Re  $\zeta \leq 2\log 30$  and thus the maximum value of  $\kappa$  is approximately  $27/16(30)^{4/3} \approx 157$ . By contrast, the equation is well-conditioned in the neighborhood of a pole or zero. For the solution in the figure, the behavior in the neighborhood of a pole at  $\zeta_0$  is (see Table 5)  $\widetilde{w} \approx w = e^{\zeta/2}u \sim e^{\zeta_0/2}c_{-2}(\zeta - \zeta_0)^{-2}$ . Similarly, in the neighborhood of a zero,  $\widetilde{w} \approx w = e^{\zeta/2}u \sim e^{\zeta_0/2}c_1(\zeta - \zeta_0)$ . For the parameter values of the solution in the figure,  $\kappa_r$  and the absolute condition number,  $\kappa_a = \kappa_r |\widetilde{w}''|$ , simplify to

$$\kappa_r \sim \frac{e^{2\operatorname{Re}\zeta_0}}{1536} |\zeta - \zeta_0|^6 \quad \text{and} \quad \kappa_a \sim \frac{e^{2\operatorname{Re}\zeta_0}}{32} |\zeta - \zeta_0|^2, \quad \zeta \to \zeta_0,$$

in the neighborhood of a pole and

$$\kappa_r \sim 10e^{2\operatorname{Re}\zeta_0}|\zeta - \zeta_0|^2$$
 and  $\kappa_a \sim \frac{e^{2\operatorname{Re}\zeta_0}}{8}|\zeta - \zeta_0|^2$ ,  $\zeta \to \zeta_0$ ,

in the neighborhood of a zero. The third column of Figure 5 shows a plot of the condition number for the tronquée  $P_{\rm III}$  solution. The plot clearly shows the exponential increase of the condition number on the pole-free region with Re  $\zeta$  as well as the comparatively well-conditioned nature of the pole fields. However, we also observe isolated points inside and between the pole fields with large condition numbers. These are the points at which the second derivative is small.

Poles Zeros  $c_{1} = -\sqrt{-\delta/2}$   $c_{2} = 2/\alpha$   $c_{1} = +\sqrt{-\delta/2}$   $c_{1} = -\sqrt{-\delta/2}$   $c_{1} = -\sqrt{-\delta/2}$   $c_{1} = -\sqrt{-\delta/2}$   $c_{2} = 8e^{-\zeta_{0}/2}/\alpha$   $c_{1} = -e^{\zeta_{0}/2}\sqrt{-\delta/2}$   $c_{1} = -e^{\zeta_{0}/2}\sqrt{-\delta/2}$ 

<span id="page-15-0"></span>Table 5: The poles and zeros of the solution in Figure 5.

### 3.3 Generic $P_V$ and $P_{VI}$ solutions

We close with  $P_V$  and  $P_{VI}$  solutions in Figures 6 and 7, respectively. These solutions were computed by applying the enhanced PFS method to the  $\widetilde{P}_V$  equation and the transformed  $P_{VI}$  equation (5). It follows from Table 3 that for a  $P_{VI}$  solution with a pole at  $z_0$ , where  $|z_0| > 1$ , the residue of the pole in the z-plane is larger than the corresponding residues in the  $\eta$  and  $\zeta$  planes by a factor of at least  $|z_0|$ . We therefore found it necessary to plot  $\log_{10}|u|$  in the z-plane in Figure 7 (column 3), instead of |u|, which is what we plot in the  $\eta$  and  $\zeta$  planes (columns 1 and 2, respectively).

The generic solutions in Figures 6 and 7 and the generic  $P_{III}$  solution in the third column of Figure 1 share a common feature: they have poles and zeros along oblique lines in the transformed planes. These sloping lines of poles and zeros are mapped to spirals in z-plane (for the  $P_{III}$  and  $P_{V}$  solutions) or  $\zeta$ -plane (for the  $P_{VI}$  solution). As we remarked above, spirals of poles and zeros were not observed in the single-valued Painlevé transcendents in [9–11,25,26].

<span id="page-16-0"></span>Table 6: The poles and zeros of the solutions in Figures 6 and 7.

|                         |                              |                                             | Poles                          |   | Zeros                       |  |
|-------------------------|------------------------------|---------------------------------------------|--------------------------------|---|-----------------------------|--|
| Figure 6 P <sub>V</sub> | $P_V, \alpha\beta \neq 0$    | z-plane                                     | $c_{-1} = +z_0/\sqrt{2\alpha}$ | • | $c_1 = +\sqrt{-2\beta}/z_0$ |  |
|                         |                              |                                             | $c_{-1} = -z_0/\sqrt{2\alpha}$ | 0 | $c_1 = -\sqrt{-2\beta}/z_0$ |  |
|                         |                              | (-plane                                     | $c_{-1} = +1/\sqrt{2\alpha}$   | • | $c_1 = +\sqrt{-2\beta}$     |  |
|                         |                              |                                             | $c_{-1} = -1/\sqrt{2\alpha}$   | 0 | $c_1 = -\sqrt{-2\beta}$     |  |
| Figure 7                | $P_{VI}, \alpha\beta \neq 0$ | $\eta, \zeta \text{ and } z \text{ planes}$ | See Table 3                    |   |                             |  |

![](_page_16_Figure_2.jpeg)

<span id="page-16-1"></span>Figure 6: A generic  $P_V$  solution in the  $\zeta$  and z planes  $(z=e^{\zeta})$ . The solution has parameter values  $(\alpha,\beta,\gamma,\delta)=(1,-1,1,-1/2)$  and ICs  $z_0=1,\ u(z_0)=2$  and  $u'(z_0)=-1$ . The error estimates for the solution on sheets 0–2 are 3e-9, 7e-6 and 2e-5, respectively.

![](_page_17_Figure_0.jpeg)

<span id="page-17-0"></span>Figure 7: A generic PVI solution in the η, ζ and z planes (ζ = e η , z = e ζ ). The solution has parameters (α, β, γ, δ) = (1, −1, 3/4, −3/2) and ICs z<sup>0</sup> = 2, u(z0) = 3/2 and u 0 (2) = −1. The error estimates for the solution on the strips indicated in the η-plane are, from bottom to top: 5e-8, 4e-4 and 8e-4.

# 4 Conclusions

To our knowledge we have presented the first numerical method for computing multivalued PIII, P<sup>V</sup> and PVI solutions on multiple sheets of their Riemann surfaces. In the process we have displayed, for the first time as far as we are aware, pole field patterns of generic and tronqu´ee PIII and P<sup>V</sup> solutions. For both these equations tronqu´ee solutions had been studied only theoretically before. We also displayed a generic PVI solution as well as what appears to be a tronqu´ee PVI solution, which to the best of our knowledge has not been proposed in the literature as yet.

We extended the capabilities of the PFS method for the single-valued Painlev´e transcendents. In particular, our enhanced PFS method can compute highly non-uniform pole fields accurately and efficiently and it can move onto the desired sheets of the Riemann surfaces of the multivalued Painlev´e transcendents by following appropriate paths around the branch points. In forthcoming studies we intend to use our methods to systematically explore multivalued Painlev´e transcendents. Our enhanced PFS method has already proven to be a valuable tool for the exploration of a family of tronqu´ee PIII solutions. We shall present these results elsewhere.

The methods we have presented in this paper are applicable to any ODE that possesses the Painlev´e property. Our methods could be used to explore ODEs that do not possess the Painlev´e property if they are used in conjunction with methods for singularity detection, see [\[4,](#page-18-4) [16,](#page-19-13) [27\]](#page-20-4), so that movable branch points can be identified. However, our method of Pad´e approximation at each step would need to be modified if the solution has essential singularities. Furthermore, our methods are not applicable beyond natural boundaries, which are closed curves in the complex plane beyond which the solution cannot be analytically continued. The Chazy equation is a well-known example of an ODE whose solutions have a movable natural boundary.

# References

- <span id="page-18-3"></span>[1] A. A. Abramov and L. F. Yukhno. A method for calculating the Painlev´e transcendents. Appl. Numer. Math., 93:262–269, 2015.
- <span id="page-18-5"></span>[2] F.V. Andreev and A.V. Kitaev. Exponentially small corrections to divergent asymptotic expansions of solutions of the fifth Painlev´e equation. Math. Res. Lett., 4(5):741–759, 1997.
- <span id="page-18-1"></span>[3] P. A. Clarkson. Painlev´e equations–nonlinear special functions. In F. Marcell´an and W. van Assche, editors, Orthogonal Polynomials and Special Functions: Computation and Application, volume 1883 of Lecture Notes in Math., pages 331–411. Springer, Berlin, 2006.
- <span id="page-18-4"></span>[4] G. F. Corliss. Integrating ODEs in the complex plane—pole vaulting. Math. Comp., 35(152):1181–1189, 1980.
- <span id="page-18-6"></span>[5] D. Dai and L. Zhang. On tronqu´ee solutions of the first Painlev´e hierarchy. J. Math. Anal. Appl., 368(2):393–399, 2010.
- <span id="page-18-2"></span>[6] NIST Digital Library of Mathematical Functions. http://dlmf.nist.gov/, Release 1.0.10 of 2015-08-07. Online companion to [\[23\]](#page-19-14).
- <span id="page-18-0"></span>[7] F.H.L. Essler, H. Frahm, A. R. Its, and V. E. Korepin. Painlev´e transcendent describes quantum correlation function of the XXZ antiferromagnet away from the free-fermion point. J. Phys. A Math. Gen., 29(17):5619–5626, 1996.

- <span id="page-19-9"></span>[8] B. Fornberg and N. Flyer. Fast generation of 2-D node distributions for mesh-free PDE discretizations. Comput. Math. Appl., 69(7):531–544, 2015.
- <span id="page-19-3"></span>[9] B. Fornberg and J.A.C. Weideman. A numerical methodology for the Painlev´e equations. J. Comput. Phys., 230:5957–5973, 2011.
- <span id="page-19-4"></span>[10] B. Fornberg and J.A.C. Weideman. A computational exploration of the second Painlev´e equation. Found. Comput. Math., 14(5):985–1016, 2014.
- <span id="page-19-5"></span>[11] B. Fornberg and J.A.C. Weideman. A computational overview of the solution space of the imaginary Painlev´e II equation. Physica D, 309:108–118, 2015.
- <span id="page-19-2"></span>[12] V.I. Gromak, I. Laine, and S. Shimomura. Painlev´e Differential Equations in the Complex Plane. Walter de Gruyter, Berlin, 2002.
- <span id="page-19-1"></span>[13] S. Hassani and J.M. Maillard. Scaling functions in the square Ising model. J. Phys. A Math. Theor., 48(11):115205, 2015.
- <span id="page-19-8"></span>[14] A. Hinkkanen and I. Laine. Solutions of a modified fifth Painlev´e equation are meromorphic. Report. Univ. Jyv¨askyl¨a, 83:133–146, 2001.
- <span id="page-19-7"></span>[15] A. Hinkkanen and I. Laine. Solutions of a modified third Painlev´e equation are meromorphic. J. Anal. Math., 85:323–337, 2001.
- <span id="page-19-13"></span>[16] C. Hunter and B. Guerrieri. Deducing the properties of singularities of functions from their Taylor series coefficients. SIAM J. on Appl. Math., 39(2):248–263, 1980.
- <span id="page-19-11"></span>[17] N. Joshi. Tritronqu´ee solutions of perturbed first Painlev´e equations. Theoret. and Math. Phys., 137(2):1515–1519, 2003.
- [18] N. Joshi and A.V. Kitaev. On Boutroux's tritronqu´ee solutions of the first Painlev´e equation. Stud. Appl. Math., 107(3):253–291, 2001.
- [19] N. Joshi and M. Mazzocco. Existence and uniqueness of tri-tronqu´ee solutions of the second Painlev´e hierarchy. Nonlinearity, 16(2):427, 2002.
- [20] N. Joshi and T. Morrison. Existence and uniqueness of tronqu´ee solutions of the fourthorder Jimbo–Miwa second Painlev´e equation. Proc. Amer. Math. Soc., 137(6):2005–2014, 2009.
- <span id="page-19-12"></span>[21] Y. Lin, D. Dai, and P. Tibboel. Existence and uniqueness of tronqu´ee solutions of the third and fourth Painlev´e equations. Nonlinearity, 27(2):171–186, 2014.
- <span id="page-19-0"></span>[22] S. L. Lukyanov. Critical values of the Yang–Yang functional in the quantum sine-Gordon model. Nucl. Phys. B, 853(2):475–507, 2011.
- <span id="page-19-14"></span>[23] F. W. J. Olver, D. W. Lozier, R. F. Boisvert, and C. W. Clark, editors. NIST Handbook of Mathematical Functions. Cambridge University Press, New York, NY, 2010. Print companion to [\[6\]](#page-18-2).
- <span id="page-19-6"></span>[24] J. A. Reeger. A Computational Study of the Fourth Painlev´e Equation and a Discussion of Adams Predictor-Corrector Methods. PhD thesis, University of Colorado, 2013.
- <span id="page-19-10"></span>[25] J. A. Reeger and B. Fornberg. Painlev´e IV with both parameters zero: A numerical study. Stud. Appl. Math., 130(2):108–133, 2013.

- <span id="page-20-0"></span>[26] J. A. Reeger and B. Fornberg. Painlev´e IV: A numerical study of the fundamental domain and beyond. Physica D, 280–281:1–13, 2014.
- <span id="page-20-4"></span>[27] Y. Tourigny and M. Grinfeld. Deciphering singularities by discrete methods. Math. Comput., 62(205):155–169, 1994.
- <span id="page-20-2"></span>[28] E. Wegert. Visual Complex Functions: An Introduction with Phase Portraits. Birkh¨auser/Springer Basel AG, Basel, 2012.
- <span id="page-20-3"></span>[29] J. A. C. Weideman and S. C. Reddy. A MATLAB differentiation matrix suite. ACM TOMS, 26(4):465–519, 2000.
- <span id="page-20-1"></span>[30] I. M. Willers. A new integration algorithm for ordinary differential equations based on continued fraction approximations. Comm. ACM, 17:504–508, 1974.