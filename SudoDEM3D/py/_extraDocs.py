# encoding: utf-8
# 2010 © Václav Šmilauer <eudoxos@arcig.cz>
#
# This module is imported at startup. It is meant to update
# docstrings of wrapper classes, which are not practical to document
# in the c++ source itself, due to the necessity of writing
# \n for newlines and having everything as "string".
#
# PLEASE:
#
# 1. provide at least brief description of the class
#    in the c++ code (for those who read it) and 
#
# 2. Add something like
#
#    "Full documentation of this class is in py/_extraDocs.py."
#
#    to the c++ documentation.

import wrapper

# Update docstring of your class/function like this:
#
#	wrapper.YourClass.__doc__="""
#		This class is documented from _extraDocs.py. Yay!
#
#		.. note::
#			The c++ documentation will be overwritten by this string.
#	"""


wrapper.Peri3dController.__doc__=r'''
Class for controlling independently all 6 components of "engineering" :yref:`stress<Peri3dController.stress>` and :yref:`strain<Peri3dController.strain>` of periodic :yref:`Cell`. :yref:`goal<Peri3dController.goal>` are the goal values, while :yref:`stressMask<Peri3dController.stressMask>` determines which components prescribe stress and which prescribe strain.

If the strain is prescribed, appropriate strain rate is directly applied. If the stress is prescribed, the strain predictor is used: from stress values in two previous steps the value of strain rate is prescribed so as the value of stress in the next step is as close as possible to the ideal one. Current algorithm is extremly simple and probably will be changed in future, but is roboust enough and mostly works fine.

Stress error (difference between actual and ideal stress)  is evaluated in current and previous steps ($\mathrm{d}\sigma_i,\mathrm{d}\sigma_{i-1}$). Linear extrapolation is used to estimate error in the next step

.. math:: \mathrm{d}\sigma_{i+1}=2\mathrm{d}\sigma_i - \mathrm{d}\sigma_{i-1}

According to this error, the strain rate is modified by :yref:`mod<Peri3dController.mod>` parameter

.. math:: \mathrm{d}\sigma_{i+1}\left\{\begin{array}{c} >0 \rightarrow \dot{\varepsilon}_{i+1} = \dot{\varepsilon}_i - \max(\mathrm{abs}(\dot{\boldsymbol{\varepsilon}}_i))\cdot\mathrm{mod} \\ <0 \rightarrow \dot{\varepsilon}_{i+1} = \dot{\varepsilon}_i + \max(\mathrm{abs}(\dot{\boldsymbol{\varepsilon}}_i))\cdot\mathrm{mod} \end{array}\right.

According to this fact, the prescribed stress will (almost) never have exact prescribed value, but the difference would be very small (and decreasing for increasing :yref:`nSteps<Peri3dController.nSteps>`. This approach works good if one of the dominant strain rates is prescribed. If all stresses are prescribed or if all goal strains is prescribed as zero, a good estimation is needed for the first step, therefore the compliance matrix is estimated (from user defined estimations of macroscopic material parameters :yref:`youngEstimation<Peri3dController.youngEstimation>` and :yref:`poissonEstimation<Peri3dController.poissonEstimation>`) and respective strain rates is computed form prescribed stress rates and compliance matrix (the estimation of compliance matrix could be computed autamatically avoiding user inputs of this kind).

The simulation on rotated periodic cell is also supported. Firstly, the `polar decomposition <http://en.wikipedia.org/wiki/Polar_decomposition#Matrix_polar_decomposition>`_ is performed on cell's transformation matrix :yref:`trsf<Cell.trsf>` $\mathcal{T}=\mat{U}\mat{P}$, where $\mat{U}$ is orthogonal (unitary) matrix representing rotation and $\mat{P}$ is a positive semi-definite Hermitian matrix representing strain. A logarithm of $\mat{P}$ should be used to obtain realistic values at higher strain values (not implemented yet). A prescribed strain increment in global coordinates $\mathrm{d}t\cdot\dot{\boldsymbol{\varepsilon}}$ is properly rotated to cell's local coordinates and added to $\mat{P}$

.. math:: \mat{P}_{i+1}=\mat{P}+\mat{U}^{\mathsf{T}}\mathrm{d}t\cdot\dot{\boldsymbol{\varepsilon}}\mat{U}

The new value of :yref:`trsf<Cell.trsf>` is computed at $\mat{T}_{i+1}=\mat{UP}_{i+1}$. From current and next :yref:`trsf<Cell.trsf>` the cell's velocity gradient :yref:`velGrad<Cell.velGrad>` is computed (according to its definition) as

.. math:: \mat{V} = (\mat{T}_{i+1}\mat{T}^{-1}-\mat{I})/\mathrm{d}t

Current implementation allow user to define independent loading "path" for each prescribed component. i.e. define the prescribed value as a function of time (or :yref:`progress<Peri3dController.progress>` or steps). See :yref:`Paths<Peri3dController.xxPath>`.

Examples :ysrc:`examples/test/peri3dController_example1.py` and :ysrc:`examples/test/peri3dController_triaxialCompression.py` explain usage and inputs of Peri3dController, :ysrc:`examples/test/peri3dController_shear.py` is an example of using shear components and also simulation on rotated cell.
'''


wrapper.Ig2_Sphere_Sphere_L3Geom.__doc__=r'''Functor for computing incrementally configuration of 2 :yref:`Spheres<Sphere>` stored in :yref:`L3Geom`; the configuration is positioned in global space by local origin $\vec{c}$ (contact point) and rotation matrix $\mat{T}$ (orthonormal transformation matrix), and its degrees of freedom are local displacement $\vec{u}$ (in one normal and two shear directions); with :yref:`Ig2_Sphere_Sphere_L6Geom` and :yref:`L6Geom`, there is additionally $\vec{\phi}$. The first row of $\mat{T}$, i.e. local $x$-axis, is the contact normal noted $\vec{n}$ for brevity. Additionally, quasi-constant values of $\vec{u}_0$ (and $\vec{\phi}_0$) are stored as shifted origins of $\vec{u}$ (and $\vec{\phi}$); therefore, current value of displacement is always $\curr{\vec{u}}-\vec{u}_0$.

Suppose two spheres with radii $r_i$, positions $\vec{x}_i$, velocities $\vec{v}_i$, angular velocities $\vec{\omega}_i$.

When there is not yet contact, it will be created if $u_N=|\curr{\vec{x}}_2-\curr{\vec{x}}_1|-|f_d|(r_1+r2)<0$, where $f_d$ is :yref:`distFactor<Ig2_Sphere_Sphere_L3Geom.distFactor>` (sometimes also called 
\`\`interaction radius''). If $f_d>0$, then $\vec{u}_{0x}$ will be initalized to $u_N$, otherwise to 0. In another words, contact will be created if spheres enlarged by $|f_d|$ touch, and the \`\`equilibrium distance'' (where $\vec{u}_x-\vec{u}-{0x}$ is zero) will be set to the current distance if $f_d$ is positive, and to the geometrically-touching distance if negative.

Local axes (rows of $\mat{T}$) are initially defined as follows:

* local $x$-axis is $\vec{n}=\vec{x}_l=\normalized{\vec{x}_2-\vec{x}_1}$;
* local $y$-axis positioned arbitrarily, but in a deterministic manner: aligned with the $xz$ plane (if $\vec{n}_y<\vec{n}_z$) or $xy$ plane (otherwise);
* local $z$-axis $\vec{z}_l=\vec{x}_l\times\vec{y}_l$.

If there has already been contact between the two spheres, it is updated to keep track of rigid motion of the contact (one that does not change mutual configuration of spheres) and mutual configuration changes. Rigid motion transforms local coordinate system and can be decomposed in rigid translation (affecting $\vec{c}$), and rigid rotation (affecting $\mat{T}$), which can be split in rotation $\vec{o}_r$ perpendicular to the normal and rotation $\vec{o}_t$ (\`\`twist'') parallel with the normal:

.. math:: \pprev{\vec{o}_r}=\prev{\vec{n}}\times\curr{\vec{n}}.

Since velocities are known at previous midstep ($t-\Dt/2$), we consider mid-step normal

.. math:: \pprev{\vec{n}}=\frac{\prev{\vec{n}}+\curr{\vec{n}}}{2}.

For the sake of numerical stability, $\pprev{\vec{n}}$ is re-normalized after being computed, unless prohibited by :yref:`approxMask<Ig2_Sphere_Sphere_L3Geom.approxMask>`. If :yref:`approxMask<Ig2_Sphere_Sphere_L3Geom.approxMask>` has the appropriate bit set, the mid-normal is not compute, and we simply use $\pprev{\vec{n}}\approx\prev{\vec{n}}$.

Rigid rotation parallel with the normal is

.. math:: \pprev{\vec{o}_t}=\pprev{\vec{n}}\left(\pprev{\vec{n}}\cdot\frac{\pprev{\vec{\omega}}_1+\pprev{\vec{\omega}}_2}{2}\right)\Dt.

*Branch vectors* $\vec{b}_1$, $\vec{b}_2$ (connecting $\curr{\vec{x}}_1$, $\curr{\vec{x}}_2$ with $\curr{\vec{c}}$ are computed depending on :yref:`noRatch<Ig2_Sphere_Sphere_L3Geom.noRatch>` (see :yref:`here<Ig2_Sphere_Sphere_ScGeom.avoidGranularRatcheting>`).

.. math::
	:nowrap:

	\begin{align*}
		\vec{b}_1&=\begin{cases} r_1 \curr{\vec{n}} & \mbox{with \texttt{noRatch}} \\ \curr{\vec{c}}-\curr{\vec{x}}_1 & \mbox{otherwise} \end{cases} \\
		\vec{b}_2&=\begin{cases} -r_2\curr{\vec{n}} & \mbox{with \texttt{noRatch}} \\ \curr{\vec{c}}-\curr{\vec{x}}_2 & \mbox{otherwise} \end{cases} \\
	\end{align*}

Relative velocity at $\curr{\vec{c}}$ can be computed as 

.. math:: \pprev{\vec{v}_r}=(\pprev{\vec{\tilde{v}}_2}+\vec{\omega}_2\times\vec{b}_2)-(\vec{v}_1+\vec{\omega}_1\times\vec{b}_1)

where $\vec{\tilde{v}}_2$ is $\vec{v}_2$ without mean-field velocity gradient in periodic boundary conditions (see :yref:`Cell.homoDeform`). In the numerial implementation, the normal part of incident velocity is removed (since it is computed directly) with $\pprev{\vec{v}_{r2}}=\pprev{\vec{v}_r}-(\pprev{\vec{n}}\cdot\pprev{\vec{v}_r})\pprev{\vec{n}}$.

Any vector $\vec{a}$ expressed in global coordinates transforms during one timestep as

.. math:: \curr{\vec{a}}=\prev{\vec{a}}+\pprev{\vec{v}_r}\Dt-\prev{\vec{a}}\times\pprev{\vec{o}_r}-\prev{\vec{a}}\times{\pprev{\vec{t}_r}}

where the increments have the meaning of relative shear, rigid rotation normal to $\vec{n}$ and rigid rotation parallel with $\vec{n}$. Local coordinate system orientation, rotation matrix $\mat{T}$, is updated by rows, i.e.

.. math:: \curr{\mat{T}}=\begin{pmatrix} \curr{\vec{n}_x} & \curr{\vec{n}_y} & \curr{\vec{n}_z} \\ \multicolumn{3}{c}{\prev{\mat{T}_{1,\bullet}}-\prev{\mat{T}_{1,\bullet}}\times\pprev{\vec{o}_r}-\prev{\mat{T}_{1,\bullet}}\times\pprev{\vec{o}_t}} \\ \multicolumn{3}{c}{\prev{\mat{T}_{2,\bullet}}-\prev{\mat{T}_{2,\bullet}}\times\pprev{\vec{o}_r}-\prev{\mat{T}_{,\bullet}}\times\pprev{\vec{o}_t}} \end{pmatrix}

This matrix is re-normalized (unless prevented by :yref:`approxMask<Ig2_Sphere_Sphere_L3Geom.approxMask>`) and mid-step transformation is computed using quaternion spherical interpolation as

.. math:: \pprev{\mat{T}}=\mathrm{Slerp}\,\left(\prev{\mat{T}};\curr{\mat{T}};t=1/2\right).

Depending on :yref:`approxMask<Ig2_Sphere_Sphere_L3Geom.approxMask>`, this computation can be avoided by approximating $\pprev{\mat{T}}=\prev{\mat{T}}$.

Finally, current displacement is evaluated as 

.. math:: \curr{\vec{u}}=\prev{u}+\pprev{\mat{T}}\pprev{\vec{v}_r}\Dt.

For the normal component, non-incremental evaluation is preferred, giving

.. math:: \curr{\vec{u}_x}=|\curr{\vec{x}_2}-\curr{\vec{x}_1}|-(r_1+r_2)

If this functor is called for :yref:`L6Geom`, local rotation is updated as 

.. math:: \curr{\vec{\phi}}=\prev{\vec{\phi}}+\pprev{\mat{T}}\Dt(\vec{\omega}_2-\vec{\omega}_1)

.. note: TODO: ``distFactor`` is not yet implemented as described above; some formulas mix values at different times, should be checked carefully.

'''

wrapper.LawTester.__doc__='''Prescribe and apply deformations of an interaction in terms of local mutual displacements and rotations. The loading path is specified either using :yref:`path<LawTester.path>` (as sequence of 6-vectors containing generalized displacements $u_x$, $u_y$, $u_z$, $\phi_x$, $\phi_y$, $\phi_z$) or :yref:`disPath<LawTester.disPath>` ($u_x$, $u_y$, $u_z$) and :yref:`rotPath<LawTester.rotPath>` ($\phi_x$, $\phi_y$, $\phi_z$). Time function with time values (step numbers) corresponding to points on loading path is given by :yref:`pathSteps<LawTester.pathSteps>`. Loading values are linearly interpolated between given loading path points, and starting zero-value (the initial configuration) is assumed for both :yref:`path<LawTester.path>` and :yref:`pathSteps<LawTester.pathSteps>`. :yref:`hooks<LawTester.hooks>` can specify python code to run when respective point on the path is reached; when the path is finished, :yref:`doneHook<LawTester.doneHook>` will be run.


LawTester should be placed between :yref:`InteractionLoop` and :yref:`NewtonIntegrator` in the simulation loop, since it controls motion via setting linear/angular velocities on particles; those velocities are integrated by :yref:`NewtonIntegrator` to yield an actual position change, which in turn causes :yref:`IGeom` to be updated (and :yref:`contact law<LawFunctor>` applied) when :yref:`InteractionLoop` is executed. Constitutive law generating forces on particles will not affect prescribed particle motion, since both particles have all :yref:`DoFs blocked<State.blockedDOFs>` when first used with LawTester.

LawTester uses, as much as possible, :yref:`IGeom` to provide useful data (such as local coordinate system), but is able to compute those independently if absent in the respective :yref:`IGeom`:


=================== ===== ==================================
:yref:`IGeom`       #DoFs LawTester support level
=================== ===== ==================================
:yref:`L3Geom`      3     full
:yref:`L6Geom`      6     full
:yref:`ScGeom`      3     emulate local coordinate system
:yref:`ScGeom6D`    6     emulate local coordinate system
=================== ===== ==================================

Depending on :yref:`IGeom`, 3 ($u_x$, $u_y$, $u_z$) or 6 ($u_x$, $u_y$, $u_z$, $\phi_x$, $\phi_y$, $\phi_z$) degrees of freedom (DoFs) are controlled with LawTester, by prescribing linear and angular velocities of both particles in contact. All DoFs controlled with LawTester are orthogonal (fully decoupled) and are controlled independently.

When 3 DoFs are controlled, :yref:`rotWeight<LawTester.rotWeight>` controls whether local shear is applied by moving particle on arc around the other one, or by rotating without changing position; although such rotation induces mutual rotation on the interaction, it is ignored with :yref:`IGeom` with only 3 DoFs. When 6 DoFs are controlled, only arc-displacement is applied for shear, since otherwise mutual rotation would occur. 

:yref:`idWeight<LawTester.idWeight>` distributes prescribed motion between both particles (resulting local deformation is the same if ``id1`` is moved towards ``id2`` or ``id2`` towards ``id1``). This is true only for $u_x$, $u_y$, $u_z$, $\phi_x$ however ; bending rotations $\phi_y$, $\phi_z$ are nevertheless always distributed regardless of ``idWeight`` to both spheres in inverse proportion to their radii, so that there is no shear induced.

LawTester knows current contact deformation from 2 sources: from its own internal data (which are used for prescribing the displacement at every step), which can be accessed in :yref:`uTest<LawTester.uTest>`, and from :yref:`IGeom` itself (depending on which data it provides), which is stored in :yref:`uGeom<LawTester.uGeom>`. These two values should be identical (disregarding numerical percision), and it is a way to test whether :yref:`IGeom` and related functors compute what they are supposed to compute.

LawTester-operated interactions can be rendered with :yref:`GlExtra_LawTester` renderer.

See :ysrc:`scripts/test/law-test.py` for an example.

'''
