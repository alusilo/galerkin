import logging
import os
import operator

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

logger = logging.getLogger(__name__)


class Element(object):
    """docstring for Element"""

    def __init__(self):
        super(Element, self).__init__()
        self.elem = np.array([])
        self.vx = np.array([])
        self.vy = np.array([])
        self.rho = np.array([])
        self.mapW = np.array([])
        self.vmapW = np.array([])
        self.mapF = np.array([])
        self.vmapF = np.array([])
        self.vmapPML = np.array([])
        self.xmu = np.array([])
        self.xlambda = np.array([])
        self.srcArea = np.array([])
        self.mapSource = np.array([], dtype=int)
        self.rx = np.array([])
        self.sx = np.array([])
        self.ry = np.array([])
        self.sy = np.array([])
        self.nx = np.array([])
        self.ny = np.array([])
        self.sJ = np.array([])
        self.J = np.array([])
        self.x = np.array([])
        self.y = np.array([])
        self.Fscale = np.array([])
        self.K = 0
        self.ks = np.array([])
        self.tri = np.array([])
        self.ids = np.array([])
        self.area = np.array([])

    def calc_area(self):
        len1 = np.sqrt((self.vx[0] - self.vx[1]) ** 2 + (self.vy[0] - self.vy[1]) ** 2)
        len2 = np.sqrt((self.vx[1] - self.vx[2]) ** 2 + (self.vy[1] - self.vy[2]) ** 2)
        len3 = np.sqrt((self.vx[2] - self.vx[0]) ** 2 + (self.vy[2] - self.vy[0]) ** 2)
        sper = (len1 + len2 + len3) / 2.0
        self.area = np.sqrt(sper * (sper - len1) * (sper - len2) * (sper - len3))


class MeshReader(object):
    """docstring for MeshReader"""

    def __init__(self, **kwargs):
        super(MeshReader, self).__init__()
        self.plot = kwargs.pop("plot", True)
        self.mesh_file = kwargs.pop("mesh_file", None)
        self.src_smooth = kwargs["src_smooth"]
        # find source position
        self.src_position = kwargs["src_position"]
        # find gather position
        self.gather = kwargs["gather"]

        if self.mesh_file is None:
            logger.error(
                "To use this code you need a mesh file and define mesh_file in the input parameters!"
            )
            logger.info(
                "Ensure the mesh file is indentified with its correct extension (i.e. .msh, .neu)."
            )
            exit(-1)

        self.mesh_xt = os.path.splitext(self.mesh_file)[1][1:]
        self.filename = os.path.split(self.mesh_file)[1]

        if self.mesh_xt == "neu":
            logger.info("mesh file indentified -> Gambit Neutral File Format.")
            self.Nv, self.VX, self.VY, self.K, self.e2v = self.__neu()
        elif self.mesh_xt == "msh":
            logger.info("mesh file indentified -> GMSH File Format.")
            self.Nv, self.VX, self.VY, self.K, self.e2v = self.__msh()
        elif self.mesh_xt == "ele":
            self.ele_file = self.mesh_file
            self.node_file = os.path.join(
                os.path.split(self.mesh_file)[0], os.path.split(self.mesh_file)[1][:-3] + "node"
            )
            self.filename = self.node_file
            logger.info("mesh file indentified -> Triangle File Format.")
            self.Nv, self.VX, self.VY, self.K, self.e2v = self.__tri()
        else:
            logger.error("can't find this mesh format.")
            exit(-1)

        self.xmin = np.min(self.VX)
        self.xmax = np.max(self.VX)
        self.ymin = np.min(self.VY)
        self.ymax = np.max(self.VY)

        self.maps = self.FindSource()

        self.mapg = self.FindGather()

        # Build source rounding elements
        self.surrounding = self.getElementsNearSource()

        self.src_cells = len(self.surrounding)
        self.area_glob = self.get_area()
        self.inscribed_r = self.inscribed_radius()
        if self.plot:
            self.__plot()

    def get_values(self):
        return self.K, self.Nv, self.VX, self.VY, self.e2v

    def get_area(self):
        len1 = np.sqrt(
            (self.VX[self.e2v[0]] - self.VX[self.e2v[1]]) ** 2
            + (self.VY[self.e2v[0]] - self.VY[self.e2v[1]]) ** 2
        )
        len2 = np.sqrt(
            (self.VX[self.e2v[1]] - self.VX[self.e2v[2]]) ** 2
            + (self.VY[self.e2v[1]] - self.VY[self.e2v[2]]) ** 2
        )
        len3 = np.sqrt(
            (self.VX[self.e2v[2]] - self.VX[self.e2v[0]]) ** 2
            + (self.VY[self.e2v[2]] - self.VY[self.e2v[0]]) ** 2
        )
        sper = (len1 + len2 + len3) / 2.0

        return np.sqrt(sper * (sper - len1) * (sper - len2) * (sper - len3))

    def inscribed_radius(self):
        len1 = np.sqrt(
            (self.VX[self.e2v[0]] - self.VX[self.e2v[1]]) ** 2
            + (self.VY[self.e2v[0]] - self.VY[self.e2v[1]]) ** 2
        )
        len2 = np.sqrt(
            (self.VX[self.e2v[1]] - self.VX[self.e2v[2]]) ** 2
            + (self.VY[self.e2v[1]] - self.VY[self.e2v[2]]) ** 2
        )
        len3 = np.sqrt(
            (self.VX[self.e2v[2]] - self.VX[self.e2v[0]]) ** 2
            + (self.VY[self.e2v[2]] - self.VY[self.e2v[0]]) ** 2
        )

        sper = (len1 + len2 + len3) / 2.0
        area = np.sqrt(sper * (sper - len1) * (sper - len2) * (sper - len3))

        return area / sper

    def __plot(self):
        logger.info("Plotting mesh...")
        fig, ax = plt.subplots()
        mesh = tri.Triangulation(self.VX / 1000, self.VY / 1000)
        plt.triplot(mesh, lw=0.5, color="blue")
        sx, sy = self.src_position
        plt.scatter(sx / 1000, sy / 1000, color="red", marker="*", s=80)
        plt.ylabel(r"$z$ [$km$]")
        plt.xlabel(r"$x$ [$km$]")
        for gth in self.gather:
            gx, gy = gth
            plt.scatter(gx / 1000, gy / 1000, color="blue", marker="v", s=50)

        for k in self.surrounding:
            cx = np.sum(self.VX[self.e2v[:, k]]) / 3
            cy = np.sum(self.VY[self.e2v[:, k]]) / 3
            ax.text(
                cx / 1000,
                cy / 1000,
                ("{}").format(k),
                fontsize=9,
                bbox={"facecolor": "red", "alpha": 0.5, "pad": 10},
            )

        plt.gca().invert_yaxis()
        ax.set_aspect("equal")
        # plt.show()

    def getElementsNearSource(self):
        sx, sy = self.src_position
        logger.info("Indentifying elements on spatial support...")
        distance = np.array(
            list(
                map(
                    lambda a: np.sqrt(
                        (np.sum(self.VX[a]) / 3 - sx) ** 2 + (np.sum(self.VY[a]) / 3 - sy) ** 2
                    ),
                    self.e2v.T,
                )
            )
        )
        ids = np.nonzero(distance <= 2 * self.src_smooth)[0]
        if len(ids) == 0:
            ids = np.array(self.maps)

        return ids

    def barycentric(self, k, p):
        # get vertices of the elements
        v1 = np.array([self.VX[self.e2v[0, k]], self.VY[self.e2v[0, k]]])
        v2 = np.array([self.VX[self.e2v[1, k]], self.VY[self.e2v[1, k]]])
        v3 = np.array([self.VX[self.e2v[2, k]], self.VY[self.e2v[2, k]]])
        a = v1 - v2
        b = v3 - v2
        s = p - v2
        l1 = np.cross(s, a) / np.cross(b, a)
        l3 = np.cross(s, b) / np.cross(a, b)
        l2 = 1.0 - l1 - l3
        return l1 + 0, l2 + 0, l3 + 0

    def centroid(self, k):
        nodes = self.e2v[:, k]
        cx, cy = np.sum(self.VX[nodes]) / 3, np.sum(self.VY[nodes]) / 3
        return cx, cy

    def FindSource(self):
        mapS = []
        k = 0
        logger.info("Searching for element containing source...")
        while k < self.K:
            l1, l2, l3 = self.barycentric(k, self.src_position)
            # condition to calculate the element that contains the source
            if (
                (l1 >= 0.0 and l1 <= 1.0)
                and (l2 >= 0.0 and l2 <= 1.0)
                and (l3 >= 0.0 and l3 <= 1.0)
            ):
                mapS.append(k)
                logger.info("Element centroid located at %s[m]", self.centroid(k))
                break
            k = k + 1

        if len(mapS) == 0:
            logger.error("Source outside the computacional domain.")
            exit(-1)
        else:
            logger.info("Source found in triangle(s) %s.", mapS)

        return mapS

    def FindGather(self):
        mapG = [[], []]
        k = 0
        logger.info("Searching for element(s) containing gather(s)...")
        while k < self.K:
            for i, gth in enumerate(self.gather):
                l1, l2, l3 = self.barycentric(k, gth)
                # condition to calculate the element that contains the source
                if (
                    (l1 >= 0.0 and l1 <= 1.0)
                    and (l2 >= 0.0 and l2 <= 1.0)
                    and (l3 >= 0.0 and l3 <= 1.0)
                ):
                    mapG[0].append(i)
                    mapG[1].append(k)

            if len(mapG[1]) == len(self.gather):
                break

            k = k + 1

        result = np.array(sorted(np.array(mapG).T, key=operator.itemgetter(0), reverse=False)).T[1]

        if len(mapG[1]) == 0:
            logger.error("All gathers are outside the computacional domain.")
            exit(-1)
        else:
            logger.info("Gather found in triangle %s.", result)

        return result

    # Function to read the mesh
    def __neu(self):
        if os.path.exists(self.mesh_file):
            with open(self.mesh_file, "rt") as f:
                data = f.readlines()
                Nv, K = np.array(data[6].split()[:2]).astype(int)
                VX, VY = np.array(
                    list(map(lambda a: a.split()[1:], data[9 : Nv + 9])), dtype=float
                ).T
                e2v = (
                    np.array(
                        list(map(lambda a: a.split()[3:], data[Nv + 11 : K + Nv + 11])), dtype=int
                    ).T
                    - 1
                )

            logger.info(
                "Mesh file %s loadded successfully. (%s E, %s V)",
                self.filename,
                K,
                Nv,
            )

            return Nv, VX, VY, K, e2v
        else:
            logger.error("Especified file '%s' does not exists.", self.mesh_file)
            exit(-1)

    def __tri(self):
        if os.path.exists(self.ele_file):
            data = []
            with open(self.ele_file, "rt") as f:
                for line in f.readlines():
                    if not line.startswith("#"):
                        data.append(line)
                K = int(data[0].split()[0])
                e2v = np.array(list(map(lambda a: a.split()[1:], data[1 : K + 1])), dtype=int).T - 1
            logger.info("Mesh file %s loadded successfully.", self.ele_file)
            logger.info("There were found %s elements.", K)
        else:
            logger.error("Especified file '%s' does not exists.", self.mesh_file)
            exit(-1)

        if os.path.exists(self.node_file):
            data = []
            with open(self.node_file, "rt") as f:
                for line in f.readlines():
                    if not line.startswith("#"):
                        data.append(line)
                Nv = int(data[0].split()[0])
                VX, VY = np.array(
                    list(map(lambda a: a.split()[1:3], data[1 : Nv + 1])), dtype=float
                ).T
            logger.info("Mesh file %s loadded successfully.", self.node_file)
            logger.info("There were found %s vertices.", Nv)
        else:
            logger.error("Especified file '%s' does not exists.", self.mesh_file)
            exit(-1)

        return Nv, VX, VY, K, e2v

    def __msh(self):
        if os.path.exists(self.mesh_file):
            with open(self.mesh_file, "rt") as f:
                data = f.read().splitlines()
                # line identification for nodes
                nodeId = data.index("$Nodes") + 1
                # line identification for elements
                elementId = data.index("$Elements") + 1
                # Number of vertex
                Nv = int(data[nodeId])
                # number of elements
                K = int(data[elementId])

                # Nodes array
                nodes = np.array(
                    [
                        list(map(float, a.split(" ")[1:3]))
                        for a in data[nodeId + 1 : nodeId + Nv + 1]
                    ]
                )
                VX = nodes[:, 0]
                VY = nodes[:, 1]
                # Elements array (check if elements are of type 2)
                e2v = (
                    np.array(
                        [
                            list(map(int, a.split(" ")[1:]))[-3:]
                            for a in data[elementId + 1 : elementId + K + 1]
                            if list(map(int, a.split(" ")[1:]))[0] == 2
                        ]
                    ).T
                    - 1
                )
                # Update number of elements
                K = np.size(e2v, 1)

            logger.info(
                "Mesh file %s loadded successfully (%s E, %s V).",
                self.filename,
                K,
                Nv,
            )

            return Nv, VX, VY, K, e2v
        else:
            logger.error("Especified file '%s' does not exists.", self.mesh_file)
            exit(-1)
