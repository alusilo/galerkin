import json
import logging
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as la
from scipy.interpolate import griddata

from .nodalDG import NodalDG2D, Vandermonde2D, rk4a, rk4b

logger = logging.getLogger(__name__)

# np.set_printoptions(precision=4,suppress=True)

# standard element
SV1 = np.array([[-1.0], [-1.0]])
SV2 = np.array([[1.0], [-1.0]])
SV3 = np.array([[-1.0], [1.0]])


def Grad2D(info, u):
    ur = np.dot(u, info.Dr.T)
    us = np.dot(u, info.Ds.T)
    ux = info.rx * ur + info.sx * us
    uy = info.ry * ur + info.sy * us

    return ux, uy


class WaveDrive2D(NodalDG2D):
    """docstring for WaveDrive2D"""

    def __init__(self, **kwargs):
        super(WaveDrive2D, self).__init__(**kwargs)
        self.src_file = kwargs.pop("src_file", False)
        self.src_freq = kwargs["src_freq"]
        self.src_delay = kwargs["src_delay"]

        PROJECT_NAME = kwargs["project"]
        # Get project root (two levels up from src/galerkin/)
        DG_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        LOCAL_PROJECT_DIR = os.path.join("resources", "output", PROJECT_NAME)
        self.PROJECT_DIR = os.path.join(DG_ROOT, LOCAL_PROJECT_DIR)
        if not os.path.exists(self.PROJECT_DIR):
            logger.info("Starting project -> %s...", PROJECT_NAME)
            os.makedirs(self.PROJECT_DIR)
        else:
            logger.warning("Directory %s already exist.", LOCAL_PROJECT_DIR)
            logger.warning("This process will overwrite all project files.")

        # resolution
        self.pixel_size = kwargs["pixel_size"]
        self.tpf = kwargs["tpf"]

        # abc parameters
        self.pml_coef = kwargs["pml_coef"]
        # read source momentum
        mxx, myy, mxy = kwargs.get("stress", (0, 0, 0))
        self.stress = (mxx, myy, mxy)
        self.displacement = kwargs.get("displacement", (0, 0))

        # data arrays and initial conditions
        self.x = np.array([])
        self.y = np.array([])
        self.density = np.array([])
        self.mapSource = np.array([], dtype=int)
        self.srcArea = np.array([])
        for obj in self.pinfo:
            self.x = np.append(self.x, obj.x)
            self.y = np.append(self.y, obj.y)
            self.density = np.append(self.density, obj.rho)
            self.mapSource = np.append(self.mapSource, obj.mapSource)
            self.srcArea = np.append(self.srcArea, obj.srcArea)

        self.dist_loc = []
        sx, sy = self.src_position
        point2src = (self.x[self.mapSource] - sx) ** 2 + (self.y[self.mapSource] - sy) ** 2
        Np_ref = 6
        Np = int((self.source_order + 1) * (self.source_order + 2) / 2)
        logger.info("Creating spatial support on mesh...")
        if self.src_cells == 1:
            # source location element
            ielem = self.surrounding[0]
            # element area
            area = self.area_glob[ielem]
            # obtain barycentric coordinates
            l1, l2, l3 = self.barycentric(ielem, (sx, sy))
            # transform global coord to local coord
            sr, ss = l2 * SV1 + l3 * SV2 + l1 * SV3
            Vout = Vandermonde2D(self.source_order, sr, ss)
            # inverse of Vandermonde matrix
            invV = self.pinfo[self.source_order - 1].invV
            # interpolation polynomials
            imatrix = np.dot(Vout, invV)[0]
            # mass matrix
            mmatrix = self.pinfo[self.source_order - 1].massM
            # jacobian
            jacobian = self.pinfo[self.source_order - 1].J[0]
            JM = jacobian * mmatrix
            self.dist_loc = np.dot(imatrix, la.inv(JM).T)
            logger.info(
                "Exitation will be projected inside element %s of %s[m^2] of area.",
                ielem,
                area,
            )
            logger.info("Spatial support spreaded on %s elements.", self.src_cells)
        else:
            self.betatot = np.sum(np.exp(-(point2src) / (self.src_smooth**2)))
            self.dist_loc = (
                (Np / Np_ref)
                * np.exp(-(point2src) / (self.src_smooth**2))
                / (self.betatot * self.srcArea)
            )
        logger.info(
            "Spatial support covering %s[m^2] of total area.",
            np.sum(self.srcArea),
        )
        logger.info("Spatial support spreaded on %s elements.", self.src_cells)

        np.save(
            os.path.join(self.PROJECT_DIR, "initial_source_n{}.npy".format(self.Nmax)),
            np.concatenate(
                ([self.x[self.mapSource]], [self.y[self.mapSource]], [self.dist_loc]), axis=0
            ),
        )
        logger.info("Source exitation projection was stored in initial_source_n%s.npy.", self.Nmax)

        logger.info("Initializing fields...")
        self.Vx = np.zeros(len(self.x))
        self.Vy = np.zeros(len(self.x))
        self.Sxx = np.zeros(len(self.x))
        self.Syy = np.zeros(len(self.x))
        self.Sxy = np.zeros(len(self.x))

        # Set simulation time
        self.FinalTime = kwargs["duration"]
        self.src_dt = kwargs.get("src_dt", 0.0)
        # Solve the problem
        self.WavePNonCon2D()

    def interpolationMatrix(self, k, sx, sy):
        # obtain barycentric coordinates
        l1, l2, l3 = self.barycentric(k, (sx, sy))
        # transform global coord to local coord
        sr, ss = l2 * SV1 + l3 * SV2 + l1 * SV3
        Vout = Vandermonde2D(self.source_order, sr, ss)
        # inverse of Vandermonde matrix
        invV = self.pinfo[self.source_order - 1].invV
        # return interpolation matrix
        return np.dot(Vout, invV)[0]

    def updateSource(self, n):
        # Excite stresses
        self.Sxx[self.mapSource] += self.stress[0] * self.wavelet[n] * self.dist_loc * self.dt
        self.Syy[self.mapSource] += self.stress[1] * self.wavelet[n] * self.dist_loc * self.dt
        self.Sxy[self.mapSource] += self.stress[2] * self.wavelet[n] * self.dist_loc * self.dt
        # Excite displacements
        self.Vx[self.mapSource] += (
            self.displacement[0]
            * self.wavelet[n]
            * self.dist_loc
            * self.dt
            / (self.density[self.mapSource])
        )
        self.Vy[self.mapSource] += (
            self.displacement[1]
            * self.wavelet[n]
            * self.dist_loc
            * self.dt
            / (self.density[self.mapSource])
        )

    def abc(self, x, y):
        res = np.ones(x.shape)
        ids = np.nonzero(x <= self.xmin + self.pml_layer[0])
        res[ids] *= np.exp(-((self.pml_coef * (self.xmin + self.pml_layer[0] - x[ids])) ** 2))

        ids = np.nonzero(x >= self.xmax - self.pml_layer[1])
        res[ids] *= np.exp(-((self.pml_coef * (self.xmax - self.pml_layer[1] - x[ids])) ** 2))

        ids = np.nonzero(y <= self.ymin + self.pml_layer[2])
        res[ids] *= np.exp(-((self.pml_coef * (self.ymin + self.pml_layer[2] - y[ids])) ** 2))

        ids = np.nonzero(y >= self.ymax - self.pml_layer[3])
        res[ids] *= np.exp(-((self.pml_coef * (self.ymax - self.pml_layer[3] - y[ids])) ** 2))

        return res

    def WavePNonCon2D(self):
        logger.info("Initializing Runge-Kutta parameters...")
        # Runge-Kutta residual storage
        resVx = np.zeros(len(self.Vx))
        resVy = np.zeros(len(self.Vx))
        resSxx = np.zeros(len(self.Vx))
        resSyy = np.zeros(len(self.Vx))
        resSxy = np.zeros(len(self.Vx))
        # stability parameter
        logger.info("Calculating stability parameter...")
        # Kaser and Dumbster Heuristic Criterion
        self.dt = np.min((1 / (2 * self.Nmax + 1)) * 2 * self.inscribed_r / self.vp)
        self.dt = 0.000932

        if self.dt > self.FinalTime or self.dt < 0:
            logger.error(
                "Simulation time must be greater than sampling time, and this must be positive."
            )
            exit(1)

        # total simulation time steps
        time_steps = int(np.ceil(self.FinalTime / self.dt) + 1)
        # simulation time array
        time = np.array(
            [
                i * self.dt if i * self.dt < self.FinalTime else self.FinalTime
                for i in range(time_steps)
            ]
        )

        # Total simulation snaphots
        total_snaps = int(np.ceil(self.FinalTime / self.tpf))
        # frames between snapshots
        fps = int(np.ceil(time_steps / total_snaps))
        json.dump(
            {
                "duration": self.FinalTime,
                "limits": (self.xmin, self.xmax, self.ymin, self.ymax),
                "pml_layer": self.pml_layer,
                "source": {
                    "order": self.source_order,
                    "elements": self.src_cells,
                    "position": self.src_position,
                },
                "stress": self.stress,
                "displacement": self.stress,
                "gather": self.gather,
                "mesh": self.mesh_file,
                "dt": self.dt,
                "total_snaps": total_snaps,
                "frames_sampling": self.tpf,
                "time_steps": time_steps,
                "src_smooth": self.src_smooth,
            },
            open(os.path.join(self.PROJECT_DIR, "model.param"), "w"),
        )

        # number of pixels in each snapshot (cap to avoid huge memory use)
        MAX_GRID_PIXELS = 2000
        h_pixels = int(np.ceil((self.xmax - self.xmin) / self.pixel_size))
        v_pixels = int(np.ceil((self.ymax - self.ymin) / self.pixel_size))
        if h_pixels > MAX_GRID_PIXELS or v_pixels > MAX_GRID_PIXELS:
            # coarsen so grid is at most MAX_GRID_PIXELS per dimension
            scale = max(h_pixels / MAX_GRID_PIXELS, v_pixels / MAX_GRID_PIXELS)
            h_pixels = int(np.ceil(h_pixels / scale))
            v_pixels = int(np.ceil(v_pixels / scale))
            logger.info(
                "Wave-field output grid capped to %sx%s (increase pixel_size for large domains to avoid this).",
                h_pixels,
                v_pixels,
            )

        # interval to print out simulation progress
        print_step = int(np.ceil(time_steps / 100))

        # read source information
        if not self.src_file:
            self.wavelet = self.createWavelet(time)
        else:
            self.wavelet = self.ReadSource(time)

        # array to storage the information of perturbation in several points
        tracesx = np.zeros((len(self.gather), time_steps))
        tracesz = np.zeros((len(self.gather), time_steps))
        # movie array
        moviex = []
        moviez = []
        logger.info("###########################################################")
        logger.info("################   Simulation parameters   ################")
        logger.info("###########################################################")
        logger.info("Simulation time -> %s[sec]", self.FinalTime)
        logger.info("dt -> %s[sec]", self.dt)
        logger.info("Total time steps = %s", time_steps)
        logger.info("Number of simulation snapshots = %s", total_snaps)
        logger.info(
            "Computational domain -> [(%s,%s), (%s,%s)]",
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
        )
        logger.info("Minimum order of approximation -> %s", self.Nmin)
        logger.info("Maximum order of approximation -> %s", self.Nmax)
        logger.info("Maximum P-wave velocity -> %s[m/s]", np.max(self.vp))
        logger.info("Minimum P-wave velocity -> %s[m/s]", np.min(self.vp))
        logger.info("Maximum S-wave velocity -> %s[m/s]", np.max(self.vs))
        logger.info("Minimum S-wave velocity -> %s[m/s]", np.min(self.vs))
        logger.info("Maximum density value -> %s[Kg/m^3]", np.max(self.rho))
        logger.info("Minimum density value -> %s[Kg/m^3]", np.min(self.rho))
        logger.info("Maximum source frequency -> %s[Hz]", self.src_freq)
        logger.info("Source position -> %s", self.src_position)
        logger.info("Number of cells source is spreaded -> %s", self.src_cells)
        logger.info("Gather(s) position -> %s", self.gather)
        logger.info("###########################################################")
        logger.info("Calculating fields...")

        t0 = 0.0
        # outer time step loop
        for n, t in enumerate(time):
            if t != t0:
                self.dt = t - t0
            # update perturbation
            self.updateSource(n)
            # Runge-Kutta method
            for INTRK in range(5):
                # compute right hand side of TM-mode Maxwell's equations
                rhsSxx, rhsSyy, rhsSxy, rhsVx, rhsVy = self.WavePNonConRHS2D()
                # initiate and increment Runge-Kutta residuals
                resSxx = rk4a[INTRK] * resSxx + self.dt * rhsSxx
                resSyy = rk4a[INTRK] * resSyy + self.dt * rhsSyy
                resSxy = rk4a[INTRK] * resSxy + self.dt * rhsSxy
                resVx = rk4a[INTRK] * resVx + self.dt * rhsVx
                resVy = rk4a[INTRK] * resVy + self.dt * rhsVy

                # update fields
                self.Sxx = self.Sxx + rk4b[INTRK] * resSxx
                self.Syy = self.Syy + rk4b[INTRK] * resSyy
                self.Sxy = self.Sxy + rk4b[INTRK] * resSxy
                self.Vx = self.Vx + rk4b[INTRK] * resVx
                self.Vy = self.Vy + rk4b[INTRK] * resVy

            # ABC conditions
            self.Sxx[self.vmapPML] *= self.abc(self.x[self.vmapPML], self.y[self.vmapPML])
            self.Syy[self.vmapPML] *= self.abc(self.x[self.vmapPML], self.y[self.vmapPML])
            self.Sxy[self.vmapPML] *= self.abc(self.x[self.vmapPML], self.y[self.vmapPML])
            self.Vx[self.vmapPML] *= self.abc(self.x[self.vmapPML], self.y[self.vmapPML])
            self.Vy[self.vmapPML] *= self.abc(self.x[self.vmapPML], self.y[self.vmapPML])

            # print progress
            if n % print_step == 0 or n == time_steps - 1:
                sys.stdout.write("\r")
                sys.stdout.write(
                    "[%-30s] %d%%"
                    % ("=" * int(30 * n / (time_steps - 1)), 100 * (n / (time_steps - 1)))
                )
                sys.stdout.flush()

            # store simulation snapshot
            if n % fps == 0 or n == time_steps - 1:
                # performe a linear interpolation funtion of the whole computational domain
                points = np.array([self.x, self.y]).T
                X, Y = np.mgrid[
                    self.xmin : self.xmax : h_pixels * 1j, self.ymin : self.ymax : v_pixels * 1j
                ]
                Z0 = griddata(points, self.Vx, (X, Y), method="cubic", fill_value=0)
                Z1 = griddata(points, self.Vy, (X, Y), method="cubic", fill_value=0)
                moviex.append(Z0)
                moviez.append(Z1)

            # store gather
            for i, gth in enumerate(self.gather):
                points = np.array([self.x[self.vmapG[i]], self.y[self.vmapG[i]]]).T
                tracesx[i, n] = griddata(
                    points, self.Vx[self.vmapG[i]], (gth[0], gth[1]), method="cubic", fill_value=0
                )
                tracesz[i, n] = griddata(
                    points, self.Vy[self.vmapG[i]], (gth[0], gth[1]), method="cubic", fill_value=0
                )
            t0 = t

        sys.stdout.write("\n")
        sys.stdout.flush()
        moviex = np.array(moviex)
        moviez = np.array(moviez)
        logger.info("Writing wave field and receivers response...")
        np.save(os.path.join(self.PROJECT_DIR, "movieVx.npy"), moviex)
        np.save(os.path.join(self.PROJECT_DIR, "movieVy.npy"), moviez)
        np.save(os.path.join(self.PROJECT_DIR, "tracesVx.npy"), np.concatenate(([time], tracesx)))
        np.save(os.path.join(self.PROJECT_DIR, "tracesVy.npy"), np.concatenate(([time], tracesz)))
        logger.info(
            "Wave Field Files: -> (%s, %s)",
            os.path.join(self.PROJECT_DIR, "movieVx.npy"),
            os.path.join(self.PROJECT_DIR, "movieVy.npy"),
        )
        logger.info(
            "Seismogram Files -> (%s, %s)",
            os.path.join(self.PROJECT_DIR, "tracesVx.npy"),
            os.path.join(self.PROJECT_DIR, "tracesVy.npy"),
        )

        return

    def ReadSource(self, time):
        # reading source values from file
        if os.path.exists(self.src_file):
            with open(self.src_file) as f:
                data = f.readlines()
                source = np.array([np.array(line).astype(float) for line in data])
                ts = np.array([i * self.src_dt for i in range(len(source))])
            logger.info("Source file %s loadded successfully.", self.src_file)
        else:
            logger.error("Especified file '%s' does not exists.", self.src_file)
            exit(1)

        wout = np.zeros(len(time))
        # time iteration
        logger.info("Interpolating source data...")
        for n, t in enumerate(time):
            # check time sub domain associated to points in file
            for p in range(1, np.size(source, 0)):
                if t >= ts[p - 1] and t <= ts[p]:
                    t0 = ts[p - 1]
                    t1 = ts[p]
                    w0 = source[p - 1]
                    w1 = source[p]
                    break
                else:
                    t0 = 0.0
                    t1 = 0.0
            # Linear interpolation of the source values
            if t1 - t0 == 0.0:
                wout[n] = 0.0
            else:
                wout[n] = w0 + (w1 - w0) * (t - t0) / (t1 - t0)
        logger.info("Plotting source data...")
        plt.grid()
        plt.xlabel(r"$t$ [$s$]")
        plt.ylabel("Amplitud")
        plt.plot(time, wout)
        # plt.show()

        return wout

    def createWavelet(self, time):
        logger.info("Creating source data...")
        dispv = 0
        half = 0.5
        two = 2.0
        tau = time - self.src_delay
        twopi = two * np.pi
        a1 = -2.0
        a2 = -half * half * (twopi * self.src_freq) ** 2
        temp = a2 * tau * tau

        if dispv == 0:
            """
            -------------------------------------------
            --- Choix source, sortie en deplacement ---
            -------------------------------------------
            """
            # *** Gaussienne, sortie en deplacement
            wout = np.exp(temp)
            # *** derivee premiere de Gaussienne, sortie en deplacement
            # wout = a1*( tau )*np.exp(temp)
            # *** Ricker (derivee seconde de gaussienne), sortie en deplacement
            # wout = a1*( half + temp )*np.exp(temp)
        else:
            a1 = two * a1
            """
			---------------------------------------
			--- Choix source, sortie en vitesse ---
			---------------------------------------
			"""
            # *** Gaussienne, sortie en vitesse
            wout = a1 * (tau) * np.exp(temp)
            # *** derivee premiere de Gaussienne, sortie en vitesse
            # wout = a1*( half + temp )*np.exp(temp)
            # *** Ricker (derivee seconde de gaussienne), sortie en vitesse
            # wout = a1*a2*tau*np.exp(temp)*(3.0 + 2.0*temp)
        # wout = 4.513*np.exp(-(self.src_freq**2)*(time-2*self.src_delay)**2)
        wfile = os.path.join(self.PROJECT_DIR, "wavelet.src")
        with open(wfile, "w") as f:
            for value in wout:
                f.write(str(value))
                f.write("\n")
            logger.info("Source data saved in %s", "wavelet.src")
        logger.info("Plotting source data...")
        plt.figure(figsize=(10, 4))
        plt.grid()
        plt.xlabel(r"$t$ [$s$]")
        plt.ylabel("Amplitud")
        plt.plot(time, wout)
        # plt.show()
        return wout

    def WavePNonConRHS2D(self):
        # Initialize storage for right hand side residuals
        rhsVx = np.zeros(len(self.Vx))
        rhsVy = np.zeros(len(self.Vx))
        rhsSxx = np.zeros(len(self.Vx))
        rhsSyy = np.zeros(len(self.Vx))
        rhsSxy = np.zeros(len(self.Vx))

        # For each possible polynomial order
        for N in range(len(self.pinfo)):
            # Extract information for this polynomial order
            pinf = self.pinfo[N]
            K = pinf.K
            # Check to see if any elements of this order exist
            if K > 0:
                # Find location of N'th order nodes
                ids = pinf.ids
                Fmask = pinf.Fmask
                # Extract N'th order nodes
                SxxN = self.Sxx[ids].reshape(ids.shape)
                SyyN = self.Syy[ids].reshape(ids.shape)
                SxyN = self.Sxy[ids].reshape(ids.shape)
                VxN = self.Vx[ids].reshape(ids.shape)
                VyN = self.Vy[ids].reshape(ids.shape)
                xmuN = pinf.xmu.reshape(ids.shape)
                xlambdaN = pinf.xlambda.reshape(ids.shape)
                rhoN = pinf.rho.reshape(ids.shape)

                # Extract '-' traces of N'th order nodal data
                SxxM = SxxN[:, Fmask.flat]
                SyyM = SyyN[:, Fmask.flat]
                SxyM = SxyN[:, Fmask.flat]
                VxM = VxN[:, Fmask.flat]
                VyM = VyN[:, Fmask.flat]
                xmuM = xmuN[:, Fmask.flat]
                xlambdaM = xlambdaN[:, Fmask.flat]
                rhoM = rhoN[:, Fmask.flat]

                # Storage for '+' traces
                SxxP = np.zeros(SxxM.shape)
                SyyP = np.zeros(SyyM.shape)
                SxyP = np.zeros(SxyM.shape)
                VxP = np.zeros(VxM.shape)
                VyP = np.zeros(VyM.shape)
                xmuP = np.zeros(xmuM.shape)
                xlambdaP = np.zeros(xlambdaM.shape)
                rhoP = np.zeros(xlambdaM.shape)

                # For each possible order
                for N2 in range(len(self.pinfo)):
                    # Check to see if any neighbor nodes of this order were located
                    if len(pinf.fmapM[N2]) > 0:
                        # L2 project N2'th order neighbor data onto N'th order trace space
                        interp = pinf.interpP[N2]
                        fmapM = pinf.fmapM[N2]
                        vmapP = pinf.vmapP[N2]
                        SxxP.flat[fmapM] = np.dot(self.Sxx[vmapP].reshape(vmapP.shape), interp.T)
                        SyyP.flat[fmapM] = np.dot(self.Syy[vmapP].reshape(vmapP.shape), interp.T)
                        SxyP.flat[fmapM] = np.dot(self.Sxy[vmapP].reshape(vmapP.shape), interp.T)
                        VxP.flat[fmapM] = np.dot(self.Vx[vmapP].reshape(vmapP.shape), interp.T)
                        VyP.flat[fmapM] = np.dot(self.Vy[vmapP].reshape(vmapP.shape), interp.T)

                # Compute jumps of trace data at faces
                dSxx = SxxM - SxxP
                dSyy = SyyM - SyyP
                dSxy = SxyM - SxyP
                dVx = VxM - VxP
                dVy = VyM - VyP
                dxmu = xmuM - xmuP
                dxlambda = xlambdaM - xlambdaP
                drho = rhoM - rhoP

                # Apply PEC boundary condition at wall boundary faces
                dSxx.flat[pinf.mapW] = 2 * SxxM.flat[pinf.mapW]
                dSyy.flat[pinf.mapW] = 2 * SyyM.flat[pinf.mapW]
                dSxy.flat[pinf.mapW] = 2 * SxyM.flat[pinf.mapW]
                dVx.flat[pinf.mapW] = 0
                dVy.flat[pinf.mapW] = 0

                # Impose free surface boundary conditions in the boundary with normal vectors [nx ny]^T
                # dSxx.flat[pinf.mapF] = 2*SxxM.flat[pinf.mapF]
                dSyy.flat[pinf.mapF] = 2 * SyyM.flat[pinf.mapF]
                dSxy.flat[pinf.mapF] = 2 * SxyM.flat[pinf.mapF]
                dSxx.flat[pinf.mapF] = 0.0
                # displacement continuity
                dVx.flat[pinf.mapF] = 0.0
                dVy.flat[pinf.mapF] = 0.0

                # evaluate upwind fluxes
                fluxSxx = -(pinf.ny * dxlambda * dVy + pinf.nx * (dxlambda + 2 * dxmu) * dVx)
                fluxSyy = -(pinf.nx * dxlambda * dVx + pinf.ny * (dxlambda + 2 * dxmu) * dVy)
                fluxSxy = -dxmu * (pinf.ny * dVx + pinf.nx * dVy)
                fluxVx = -(1 / drho) * (pinf.nx * dSxx + pinf.ny * dSxy)
                fluxVy = -(1 / drho) * (pinf.ny * dSyy + pinf.nx * dSxy)

                # local derivatives of fields
                Sxx_x, Sxx_y = Grad2D(pinf, SxxN)
                Syy_x, Syy_y = Grad2D(pinf, SyyN)
                Sxy_x, Sxy_y = Grad2D(pinf, SxyN)
                Vx_x, Vx_y = Grad2D(pinf, VxN)
                Vy_x, Vy_y = Grad2D(pinf, VyN)

                # compute right hand sides of the PDE's
                rhsSxx[ids] = (
                    (xlambdaN + 2 * xmuN) * Vx_x
                    + xlambdaN * Vy_y
                    + np.dot(pinf.Fscale * fluxSxx, pinf.LIFT) / 2
                )
                rhsSyy[ids] = (
                    xlambdaN * Vx_x
                    + (xlambdaN + 2 * xmuN) * Vy_y
                    + np.dot(pinf.Fscale * fluxSyy, pinf.LIFT) / 2
                )
                rhsSxy[ids] = xmuN * (Vx_y + Vy_x) + np.dot(pinf.Fscale * fluxSxy, pinf.LIFT) / 2
                rhsVx[ids] = (1 / rhoN) * (Sxx_x + Sxy_y) + np.dot(
                    pinf.Fscale * fluxVx, pinf.LIFT
                ) / 2
                rhsVy[ids] = (1 / rhoN) * (Sxy_x + Syy_y) + np.dot(
                    pinf.Fscale * fluxVy, pinf.LIFT
                ) / 2

        return rhsSxx, rhsSyy, rhsSxy, rhsVx, rhsVy
