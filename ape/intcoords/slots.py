import abc
import logging

from math import exp, sin
import numpy as np

from ape.intcoords.elem_data import COVALENT_RADII as CR
from ape.intcoords.derivatives import d2q_b, d2q_a, dq_lb, d2q_lb, dq_ld, d2q_ld, d2q_d, dq_oop, d2q_oop
from ape.intcoords.rotate import get_expmap, get_expmap_der, is_linear, calc_rot_vec_diff
from ape.intcoords import nifty, math_utils

class Primitive(metaclass=abc.ABCMeta):
    def __init__(self, indices, periodic=False, calc_kwargs=None):
        self.indices = list(indices)
        self.periodic = periodic
        if calc_kwargs is None:
            calc_kwargs = ()
        self.calc_kwargs = calc_kwargs
        self.logger = logging.getLogger("internal_coords")

    def log(self, msg, lvl=logging.DEBUG):
        self.logger.log(lvl, msg)

    @staticmethod
    def parallel(u, v, thresh=1e-6):
        dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return (1 - abs(dot)) < thresh

    @staticmethod
    def _get_cross_vec(coords3d, indices):
        m, o, n = indices
        # Select initial vector for cross product, similar to
        # geomeTRIC. It must NOT be parallel to u and/or v.
        x_dash = coords3d[n] - coords3d[m]
        x = x_dash / np.linalg.norm(x_dash)
        cross_vecs = np.eye(3)
        min_ind = np.argmin([np.dot(cv, x) ** 2 for cv in cross_vecs])
        return cross_vecs[min_ind]

    def set_cross_vec(self, coords3d, indices):
        self.cross_vec = self._get_cross_vec(coords3d, self.indices)
        self.log(f"Cross vector for {self} set to {self.cross_vec}")

    @abc.abstractmethod
    def _calculate(*, coords3d, indices, gradient, **kwargs):
        pass

    @abc.abstractmethod
    def _weight(self, atoms, coords3d, indices, f_damping):
        pass

    def weight(self, atoms, coords3d, f_damping=0.12):
        return self._weight(atoms, coords3d, self.indices, f_damping)

    @staticmethod
    def rho(atoms, coords3d, indices):
        i, j = indices
        distance = np.linalg.norm(coords3d[i] - coords3d[j])
        cov_rad_sum = CR[atoms[i].lower()] + CR[atoms[j].lower()]
        return exp(-(distance / cov_rad_sum - 1))

    def calculate(self, coords3d, indices=None, gradient=False):
        if indices is None:
            indices = self.indices

        # Gather calc_kwargs
        calc_kwargs = {key: getattr(self, key) for key in self.calc_kwargs}

        return self._calculate(
            coords3d=coords3d,
            indices=indices,
            gradient=gradient,
            **calc_kwargs,
        )

    def jacobian(self, coords3d, indices=None):
        if indices is None:
            indices = self.indices

        # Gather calc_kwargs
        calc_kwargs = {key: getattr(self, key) for key in self.calc_kwargs}

        return self._jacobian(
            coords3d=coords3d,
            indices=indices,
            **calc_kwargs,
        )

    def __str__(self):
        return f"{self.__class__.__name__}({self.indices})"

    def __repr__(self):
        return self.__str__()

class CartesianX(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, w=1.0):
        ind = indices[0]
        value = coords3d[ind][0]
        if gradient:
            row = np.zeros_like(coords3d)
            row[ind][0] = w
            row = row.flatten()
            return value, row
        return value

class CartesianY(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, w=1.0):
        ind = indices[0]
        value = coords3d[ind][1]
        if gradient:
            row = np.zeros_like(coords3d)
            row[ind][1] = w
            row = row.flatten()
            return value, row
        return value

class CartesianZ(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, w=1.0):
        ind = indices[0]
        value = coords3d[ind][2]
        if gradient:
            row = np.zeros_like(coords3d)
            row[ind][2] = w
            row = row.flatten()
            return value, row
        return value

class TranslationX(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass
    
    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        indices = np.array(indices)
        w = np.ones(len(indices))/len(indices)
        value = np.sum(coords3d[indices, 0] * w)
        if gradient:
            row = np.zeros_like(coords3d)
            for i, a in enumerate(indices):
                row[a][0] = w[i]
            row = row.flatten()
            return value, row
        return value

class TranslationY(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass
    
    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        indices = np.array(indices)
        w = np.ones(len(indices))/len(indices)
        value = np.sum(coords3d[indices, 1] * w)
        if gradient:
            row = np.zeros_like(coords3d)
            for i, a in enumerate(indices):
                row[a][1] = w[i]
            row = row.flatten()
            return value, row
        return value

class TranslationZ(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass
    
    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        indices = np.array(indices)
        w = np.ones(len(indices))/len(indices)
        value = np.sum(coords3d[indices, 2] * w)
        if gradient:
            row = np.zeros_like(coords3d)
            for i, a in enumerate(indices):
                row[a][2] = w[i]
            row = row.flatten()
            return value, row
        return value

class Rotator(object):

    def __init__(self, a, x0):
        self.a = list(tuple(sorted(a)))
        x0 = x0.reshape(-1, 3)
        self.x0 = x0.copy()
        self.stored_valxyz = np.zeros_like(x0)
        self.stored_value = None
        # A second set of xyz coordinates used only when computing
        # differences in rotation coordinates
        self.stored_valxyz2 = np.zeros_like(x0)
        self.stored_value2 = None
        self.stored_derxyz = np.zeros_like(x0)
        self.stored_deriv = None
        self.stored_deriv2xyz = np.zeros_like(x0)
        self.stored_deriv2 = None
        self.stored_norm = 0.0
        # Extra variables to account for the case of linear molecules
        # The reference axis used for computing dummy atom position
        self.e0 = None
        # Dot-squared measures alignment of molecule long axis with reference axis.
        # If molecule becomes parallel with reference axis, coordinates must be reset.
        self.stored_dot2 = 0.0
        # Flag that records linearity of molecule
        self.linear = False

    def reset(self, x0):
        x0 = x0.reshape(-1, 3)
        self.x0 = x0.copy()
        self.stored_valxyz = np.zeros_like(x0)
        self.stored_value = None
        self.stored_valxyz2 = np.zeros_like(x0)
        self.stored_value2 = None
        self.stored_derxyz = np.zeros_like(x0)
        self.stored_deriv = None
        self.stored_deriv2xyz = np.zeros_like(x0)
        self.stored_deriv2 = None
        self.stored_norm = 0.0
        self.e0 = None
        self.stored_dot2 = 0.0
        self.linear = False

    def __eq__(self, other):
        if type(self) is not type(other): return False
        eq = set(self.a) == set(other.a)
        if eq and np.sum((self.x0-other.x0)**2) > 1e-6:
            logger.warning("Warning: Rotator same atoms, different reference positions\n")
        return eq

    def __repr__(self):
        return "Rotator %s" % commadash(self.a)

    def __ne__(self, other):
        return not self.__eq__(other)
    
    @property
    def w(self):
        sel = self.x0[self.a,:] 
        sel -= np.mean(sel, axis=0)
        rg = np.sqrt(np.mean(np.sum(sel ** 2, axis=1)))
        return rg

    def calc_e0(self):
        """
        Compute the reference axis for adding dummy atoms. 
        Only used in the case of linear molecules.
        We first find the Cartesian axis that is "most perpendicular" to the molecular axis.
        Next we take the cross product with the molecular axis to create a perpendicular vector.
        Finally, this perpendicular vector is normalized to make a unit vector.
        """
        ysel = self.x0[self.a, :]
        vy = ysel[-1]-ysel[0]
        ev = vy / np.linalg.norm(vy)
        # Cartesian axes.
        ex = np.array([1.0,0.0,0.0])
        ey = np.array([0.0,1.0,0.0])
        ez = np.array([0.0,0.0,1.0])
        self.e0 = np.cross(vy, [ex, ey, ez][np.argmin([np.dot(i, ev)**2 for i in [ex, ey, ez]])])
        self.e0 /= np.linalg.norm(self.e0)

    def value(self, xyz, store=True):
        xyz = xyz.reshape(-1, 3)
        if np.max(np.abs(xyz-self.stored_valxyz)) < 1e-12:
            return self.stored_value
        else:
            xsel = xyz[self.a, :]
            ysel = self.x0[self.a, :]
            xmean = np.mean(xsel,axis=0)
            ymean = np.mean(ysel,axis=0)
            if not self.linear and is_linear(xsel, ysel):
                # print "Setting linear flag for", self
                self.linear = True
            if self.linear:
                # Handle linear molecules.
                vx = xsel[-1]-xsel[0]
                vy = ysel[-1]-ysel[0]
                # Calculate reference axis (if needed)
                if self.e0 is None: self.calc_e0()
                #log.debug(vx)
                ev = vx / np.linalg.norm(vx)
                # Measure alignment of molecular axis with reference axis
                self.stored_dot2 = np.dot(ev, self.e0)**2
                # Dummy atom is located one Bohr from the molecular center, direction
                # given by cross-product of the molecular axis with the reference axis
                xdum = np.cross(vx, self.e0)
                ydum = np.cross(vy, self.e0)
                exdum = xdum / np.linalg.norm(xdum)
                eydum = ydum / np.linalg.norm(ydum)
                xsel = np.vstack((xsel, exdum+xmean))
                ysel = np.vstack((ysel, eydum+ymean))
            answer = get_expmap(xsel, ysel)
            if store:
                self.stored_norm = np.linalg.norm(answer)
                self.stored_valxyz = xyz.copy()
                self.stored_value = answer.copy()
            return answer

    def calcDiff(self, xyz1, xyz2=None, val2=None):
        """
        Return the difference of the internal coordinate
        calculated for (xyz1 - xyz2).
        """
        if xyz2 is None and val2 is None:
            raise RuntimeError("Provide exactly one of xyz2 and val2")
        elif xyz2 is not None and val2 is not None:
            raise RuntimeError("Provide exactly one of xyz2 and val2")
        val1 = self.value(xyz1)
        if xyz2 is not None:
            # The "second" coordinate set is cached separately
            xyz2 = xyz2.reshape(-1, 3)
            if np.max(np.abs(xyz2-self.stored_valxyz2)) < 1e-12:
                val2 = self.stored_value2.copy()
            else:
                val2 = self.value(xyz2, store=False)
                self.stored_valxyz2 = xyz2.copy()
                self.stored_value2 = val2.copy()
        # Calculate difference in rotation vectors, modulo n*2pi displacement vectors
        return calc_rot_vec_diff(val1, val2)

    def derivative(self, xyz):
        xyz = xyz.reshape(-1, 3)
        if np.max(np.abs(xyz-self.stored_derxyz)) < 1e-12:
            return self.stored_deriv
        else:
            xsel = xyz[self.a, :]
            ysel = self.x0[self.a, :]
            xmean = np.mean(xsel,axis=0)
            ymean = np.mean(ysel,axis=0)
            if not self.linear and is_linear(xsel, ysel):
                # print "Setting linear flag for", self
                self.linear = True
            if self.linear:
                vx = xsel[-1]-xsel[0]
                vy = ysel[-1]-ysel[0]
                if self.e0 is None: self.calc_e0()
                xdum = np.cross(vx, self.e0)
                ydum = np.cross(vy, self.e0)
                exdum = xdum / np.linalg.norm(xdum)
                eydum = ydum / np.linalg.norm(ydum)
                xsel = np.vstack((xsel, exdum+xmean))
                ysel = np.vstack((ysel, eydum+ymean))
            deriv_raw = get_expmap_der(xsel, ysel)
            if self.linear:
                # Chain rule is applied to get terms from
                # dummy atom derivatives
                nxdum = np.linalg.norm(xdum)
                dxdum = math_utils.d_cross(vx, self.e0)
                dnxdum = math_utils.d_ncross(vx, self.e0)
                # Derivative of dummy atom position w/r.t. molecular axis vector
                dexdum = (dxdum*nxdum - np.outer(dnxdum,xdum))/nxdum**2
                # Here we may compute finite difference derivatives to check
                # h = 1e-6
                # fdxdum = np.zeros((3, 3), dtype=float)
                # for i in range(3):
                #     vx[i] += h
                #     dPlus = np.cross(vx, self.e0)
                #     dPlus /= np.linalg.norm(dPlus)
                #     vx[i] -= 2*h
                #     dMinus = np.cross(vx, self.e0)
                #     dMinus /= np.linalg.norm(dMinus)
                #     vx[i] += h
                #     fdxdum[i] = (dPlus-dMinus)/(2*h)
                # if np.linalg.norm(dexdum - fdxdum) > 1e-6:
                #     print dexdum - fdxdum
                #     raise Exception()
                # Apply terms from chain rule
                deriv_raw[0]  -= np.dot(dexdum, deriv_raw[-1])
                for i in range(len(self.a)):
                    deriv_raw[i]  += np.dot(np.eye(3), deriv_raw[-1])/len(self.a)
                deriv_raw[-2] += np.dot(dexdum, deriv_raw[-1])
                deriv_raw = deriv_raw[:-1]
            derivatives = np.zeros((xyz.shape[0], 3, 3), dtype=float)
            for i, a in enumerate(self.a):
                derivatives[a, :, :] = deriv_raw[i, :, :]
            self.stored_derxyz = xyz.copy()
            self.stored_deriv = derivatives.copy()
            return derivatives
        
    def second_derivative(self, xyz):
        xyz = xyz.reshape(-1, 3)
        if np.max(np.abs(xyz-self.stored_deriv2xyz)) < 1e-12:
            return self.stored_deriv2
        else:
            xsel = xyz[self.a, :]
            ysel = self.x0[self.a, :]
            xmean = np.mean(xsel,axis=0)
            ymean = np.mean(ysel,axis=0)
            if not self.linear and is_linear(xsel, ysel):
                # print "Setting linear flag for", self
                self.linear = True
            if self.linear:
                vx = xsel[-1]-xsel[0]
                vy = ysel[-1]-ysel[0]
                if self.e0 is None: self.calc_e0()
                xdum = np.cross(vx, self.e0)
                ydum = np.cross(vy, self.e0)
                exdum = xdum / np.linalg.norm(xdum)
                eydum = ydum / np.linalg.norm(ydum)
                xsel = np.vstack((xsel, exdum+xmean))
                ysel = np.vstack((ysel, eydum+ymean))
            deriv_raw, deriv2_raw = get_expmap_der(xsel, ysel, second=True)
            if self.linear:
                # Chain rule is applied to get terms from dummy atom derivatives
                def dexdum_(vx_):
                    xdum_ = np.cross(vx_, self.e0)
                    nxdum_ = np.linalg.norm(xdum_)
                    dxdum_ = math_utils.d_cross(vx_, self.e0)
                    dnxdum_ = math_utils.d_ncross(vx_, self.e0)
                    dexdum_ = (dxdum_*nxdum_ - np.outer(dnxdum_,xdum_))/nxdum_**2
                    return dexdum_.copy()
                # First indices: elements of vx that are being differentiated w/r.t.
                # Last index: elements of exdum itself
                dexdum = dexdum_(vx)
                dexdum2 = np.zeros((3, 3, 3), dtype=float)
                h = 1.0e-3
                for i in range(3):
                    vx[i] += h
                    dPlus = dexdum_(vx)
                    vx[i] -= 2*h
                    dMinus = dexdum_(vx)
                    vx[i] += h
                    dexdum2[i] = (dPlus-dMinus)/(2*h)
                # Build arrays that contain derivative of dummy atom position
                # w/r.t. real atom positions
                ddum1 = np.zeros((len(self.a), 3, 3), dtype=float)
                ddum1[0] = -dexdum
                ddum1[-1] = dexdum
                for i in range(len(self.a)):
                    ddum1[i] += np.eye(3)/len(self.a)
                ddum2 = np.zeros((len(self.a), 3, len(self.a), 3, 3), dtype=float)
                ddum2[ 0, : , 0, :] =  dexdum2
                ddum2[-1, : , 0, :] = -dexdum2
                ddum2[ 0, :, -1, :] = -dexdum2
                ddum2[-1, :, -1, :] =  dexdum2
                # =====
                # Do not delete - reference codes using loops for chain rule terms
                # for j in range(len(self.a)): # Loop over atom 1
                #     for m in range(3):       # Loop over xyz of atom 1
                #         for k in range(len(self.a)): # Loop over atom 2
                #             for n in range(3):       # Loop over xyz of atom 2
                #                 for i in range(3):   # Loop over elements of exponential map
                #                     for p in range(3): # Loop over xyz of dummy atom
                #                         deriv2_raw[j, m, k, n, i] += deriv2_raw[j, m, -1, p, i] * ddum1[k, n, p]
                #                         deriv2_raw[j, m, k, n, i] += deriv2_raw[-1, p, k, n, i] * ddum1[j, m, p]
                #                         deriv2_raw[j, m, k, n, i] += deriv_raw[-1, p, i] * ddum2[j, m, k, n, p]
                #                         for q in range(3):
                #                             deriv2_raw[j, m, k, n, i] += deriv2_raw[-1, p, -1, q, i] * ddum1[j, m, p] * ddum1[k, n, q]
                # =====
                deriv2_raw[:-1, :, :-1, :] += np.einsum('jmpi,knp->jmkni', deriv2_raw[:-1, :, -1, :, :], ddum1, optimize=True)
                deriv2_raw[:-1, :, :-1, :] += np.einsum('pkni,jmp->jmkni', deriv2_raw[-1, :, :-1, :, :], ddum1, optimize=True)
                deriv2_raw[:-1, :, :-1, :] += np.einsum('pi,jmknp->jmkni', deriv_raw[-1, :, :], ddum2, optimize=True)
                deriv2_raw[:-1, :, :-1, :] += np.einsum('pqi,jmp,knq->jmkni', deriv2_raw[-1, :, -1, :, :], ddum1, ddum1, optimize=True)
                deriv2_raw = deriv2_raw[:-1, :, :-1, :, :]
            second_derivatives = np.zeros((xyz.shape[0], 3, xyz.shape[0], 3, 3), dtype=float)
            for i, a in enumerate(self.a):
                for j, b in enumerate(self.a):
                    second_derivatives[a, :, b, :, :] = deriv2_raw[i, :, j, :, :]
            return second_derivatives

class RotationA(Primitive):

    def __init__(self, indices, coords3d):
        Primitive.__init__(self, indices)
        self.rotator =  Rotator(indices, coords3d)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass

    def _calculate(self, coords3d, indices, gradient=False):
        value = self.rotator.value(coords3d)[0] * self.rotator.w
        if gradient:
            row = self.rotator.derivative(coords3d)
            row = row[:, :, 0] * self.rotator.w
            row = row.flatten()
            return value, row
        return value

class RotationB(Primitive):

    def __init__(self, indices, coords3d):
        Primitive.__init__(self, indices)
        self.rotator =  Rotator(indices, coords3d)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass

    def _calculate(self, coords3d, indices, gradient=False):
        value = self.rotator.value(coords3d)[1] * self.rotator.w
        if gradient:
            row = self.rotator.derivative(coords3d)
            row = row[:, :, 1] * self.rotator.w
            row = row.flatten()
            return value, row
        return value

class RotationC(Primitive):

    def __init__(self, indices, coords3d):
        Primitive.__init__(self, indices)
        self.rotator =  Rotator(indices, coords3d)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        pass

    def _calculate(self, coords3d, indices, gradient=False):
        value = self.rotator.value(coords3d)[2] * self.rotator.w
        if gradient:
            row = self.rotator.derivative(coords3d)
            row = row[:, :, 2] * self.rotator.w
            row = row.flatten()
            return value, row
        return value

class Stretch(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return Stretch.rho(atoms, coords3d, indices)

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        n, m = indices
        bond = coords3d[m] - coords3d[n]
        bond_length = np.linalg.norm(bond)
        if gradient:
            bond_normed = bond / bond_length
            row = np.zeros_like(coords3d)
            # 1 / -1 correspond to the sign factor [1] Eq. 18
            row[m,:] =  bond_normed
            row[n,:] = -bond_normed
            row = row.flatten()
            return bond_length, row
        return bond_length

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_b(*coords3d[indices].flatten())

class Bend(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        m, o, n = indices
        rho_mo = Bend.rho(atoms, coords3d, (m, o))
        rho_on = Bend.rho(atoms, coords3d, (o, n))
        rad = Bend._calculate(coords3d, indices)
        return (rho_mo * rho_on)**0.5 * (f_damping + (1-f_damping)*sin(rad))

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm

        angle_rad = np.arccos(u.dot(v))

        if gradient:
            cross_vec1 = ( 1, -1, 1)
            cross_vec2 = (-1,  1, 1)

            # Determine second vector for the cross product, to get an
            # orthogonal direction. Eq. (24) in [1]
            uv_parallel = Bend.parallel(u, v)
            if not uv_parallel:
                cross_vec = v
            elif not Bend.parallel(u, cross_vec1):
                cross_vec = cross_vec1
            else:
                cross_vec = cross_vec2

            w_dash = np.cross(u, cross_vec)
            w = w_dash / np.linalg.norm(w_dash)

            uxw = np.cross(u, w)
            wxv = np.cross(w, v)

            row = np.zeros_like(coords3d)
            #                  |  m  |  n  |  o  |
            # -----------------------------------
            # sign_factor(amo) |  1  |  0  | -1  | first_term
            # sign_factor(ano) |  0  |  1  | -1  | second_term
            first_term = uxw / u_norm
            second_term = wxv / v_norm
            row[m,:] = first_term
            row[o,:] = -first_term - second_term
            row[n,:] = second_term
            row = row.flatten()
            return angle_rad, row
        return angle_rad

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_a(*coords3d[indices].flatten())

# [1] 10.1080/00268977200102361
#     Hoy, 1972
# [2] 10.1063/1.474377
#     Chuang, 1997
# [3] 10.1063/1.468630
#     Jackels, 1995
# Refs. [2] and [3] give a short discussion of the linear bends.

class LinearBend(Primitive):
    def __init__(self, *args, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", "cross_vec")
        super().__init__(*args, **kwargs)

        self.complement = complement
        self.cross_vec = None

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        m, o, n = indices
        rho_mo = LinearBend.rho(atoms, coords3d, (m, o))
        rho_on = LinearBend.rho(atoms, coords3d, (o, n))

        # Repeated code to avoid import of intcoords.Bend
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        rad = np.arccos(u.dot(v))

        return (rho_mo * rho_on) ** 0.5 * (f_damping + (1 - f_damping) * sin(rad))

    @staticmethod
    def _get_orthogonal_direction(coords3d, indices, complement=False, cross_vec=None):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        u = u_dash / u_norm

        if cross_vec is None:
            cross_vec = LinearBend._get_cross_vec(coords3d, indices)
        # Generate first orthogonal direction
        w_dash = np.cross(u, cross_vec)
        w = w_dash / np.linalg.norm(w_dash)

        # Generate second orthogonal direction
        if complement:
            w = np.cross(u, w)
        return w

    def calculate(self, coords3d, indices=None, gradient=False):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().calculate(coords3d, indices, gradient)

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False, cross_vec=None):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w = LinearBend._get_orthogonal_direction(
            coords3d, indices, complement, cross_vec
        )

        lb_rad = w.dot(np.cross(u_dash, v_dash)) / (u_norm * v_norm)

        if gradient:
            # Fourth argument is the orthogonal direction
            row = np.zeros_like(coords3d)
            row[indices] = dq_lb(*coords3d[indices].flatten(), *w).reshape(-1, 3)
            return lb_rad, row.flatten()
        return lb_rad

    def jacobian(self, coords3d, indices=None):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().jacobian(coords3d, indices)

    @staticmethod
    def _jacobian(coords3d, indices, complement=False, cross_vec=None):
        if cross_vec is None:
            cross_vec = LinearBend._get_cross_vec(coords3d, indices)

        w = LinearBend._get_orthogonal_direction(
            coords3d, indices, complement, cross_vec
        )
        return d2q_lb(*coords3d[indices].flatten(), *w)

    def __str__(self):
        return f"LinearBend({tuple(self.indices)}, complement={self.complement})"

class LinearDisplacement(Primitive):
    def __init__(self, *args, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", "cross_vec")
        super().__init__(*args, **kwargs)

        self.complement = complement
        self.cross_vec = None

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        raise Exception("Not yet implemented!")

    def calculate(self, coords3d, indices=None, gradient=False):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().calculate(coords3d, indices, gradient)

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False, cross_vec=None):
        m, o, n = indices
        w_dash = coords3d[n] - coords3d[m]
        w = w_dash / np.linalg.norm(w_dash)

        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u = u_dash / np.linalg.norm(u_dash)
        v = v_dash / np.linalg.norm(v_dash)

        # Vector for cross product to determine first orthogonal direction
        if cross_vec is None:
            cross_vec = LinearDisplacement._get_cross_vec(coords3d, indices)

        if complement:
            cross_vec = np.cross(w, cross_vec)
        cross_vec /= np.linalg.norm(cross_vec)

        # Orthogonal direction
        y = np.cross(w, cross_vec)
        y /= np.linalg.norm(y)

        lin_disp = y.dot(u) + y.dot(v)

        if gradient:
            row = np.zeros_like(coords3d)
            row[indices] = dq_ld(*coords3d[indices].flatten(), *cross_vec).reshape(-1, 3)
            return lin_disp, row.flatten()

        return lin_disp

    def jacobian(self, coords3d, indices=None):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().jacobian(coords3d, indices)

    @staticmethod
    def _jacobian(coords3d, indices, complement=False, cross_vec=None):
        if cross_vec is None:
            cross_vec = LinearDisplacement._get_cross_vec(coords3d, indices)

        if complement:
            m, _, n = indices
            w_dash = coords3d[n] - coords3d[m]
            w = w_dash / np.linalg.norm(w_dash)
            cross_vec = np.cross(w, cross_vec)

        return d2q_ld(*coords3d[indices].flatten(), *cross_vec)

    def __str__(self):
        return (
            f"LinearDisplacement({tuple(self.indices)}, complement={self.complement})"
        )

class Torsion(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        m, o, p, n = indices
        rho_mo = Torsion.rho(atoms, coords3d, (m, o))
        rho_op = Torsion.rho(atoms, coords3d, (o, p))
        rho_pn = Torsion.rho(atoms, coords3d, (p, n))
        rad_mop = Bend._calculate(coords3d, (m, o, p))
        rad_opn = Bend._calculate(coords3d, (o, p, n))
        return (
            (rho_mo * rho_op * rho_pn)**(1/3)
            * (f_damping + (1-f_damping)*sin(rad_mop))
            * (f_damping + (1-f_damping)*sin(rad_opn))
        )

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        m, o, p, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[p]
        w_dash = coords3d[p] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w_norm = np.linalg.norm(w_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        w = w_dash / w_norm
        phi_u = np.arccos(u.dot(w))
        phi_v = np.arccos(-w.dot(v))
        uxw = np.cross(u, w)
        vxw = np.cross(v, w)
        cos_dihed = uxw.dot(vxw)/(np.sin(phi_u)*np.sin(phi_v))
        # Restrict cos_dihed to the allowed interval for arccos [-1, 1]
        cos_dihed = min(1, max(cos_dihed, -1))

        dihedral_rad = np.arccos(cos_dihed)

        # Arccos only returns values between 0 and π, but dihedrals can
        # also be negative. This is corrected now.
        #
        # (v ⨯ w) · u will be < 0 when both vectors point in different directions.
        #
        #  M  --->  N
        #   \      /
        #    u    v    positive dihedral, M rotates into N clockwise
        #     \  /     (v ⨯ w) · u > 0, keep positive sign
        #      OwP
        #
        #  M
        #   \
        #  | u
        #  |  \
        #  |   OwP     negative dihedral, M rotates into N counter-clockwise
        #  v  /        (v ⨯ w) · u < 0, invert dihedral sign
        #    v
        #   /
        #  N
        #
        if (dihedral_rad != np.pi) and (vxw.dot(u) < 0):
            dihedral_rad *= -1

        if gradient:
            row = np.zeros_like(coords3d)
            #                  |  m  |  n  |  o  |  p  |
            # ------------------------------------------
            # sign_factor(amo) |  1  |  0  | -1  |  0  | 1st term
            # sign_factor(apn) |  0  | -1  |  0  |  1  | 2nd term
            # sign_factor(aop) |  0  |  0  |  1  | -1  | 3rd term
            # sign_factor(apo) |  0  |  0  | -1  |  1  | 4th term
            sin2_u = np.sin(phi_u)**2
            sin2_v = np.sin(phi_v)**2
            first_term  = uxw/(u_norm*sin2_u)
            second_term = vxw/(v_norm*sin2_v)
            third_term  = uxw*np.cos(phi_u)/(w_norm*sin2_u)
            fourth_term = -vxw*np.cos(phi_v)/(w_norm*sin2_v)
            row[m,:] = first_term
            row[n,:] = -second_term
            row[o,:] = -first_term + third_term - fourth_term
            row[p,:] = second_term - third_term + fourth_term
            row = row.flatten()
            return dihedral_rad, row
        return dihedral_rad

    @staticmethod
    def _jacobian(coords3d, indices):
        sign = np.sign(Torsion._calculate(coords3d, indices))
        return sign * d2q_d(*coords3d[indices].flatten())

class OutOfPlane(Primitive):
    """
    [1] https://doi.org/10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
        Lee, 1999
    """

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        raise Exception("Not yet implemented!")

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        """
              P
            / | \
           /  |  \
          u'  v'  w'
         /    |    \
        m     n     o
        """
        # p is apex
        m, n, o, p = indices

        u_dash = coords3d[m] - coords3d[p]
        v_dash = coords3d[n] - coords3d[p]
        w_dash = coords3d[o] - coords3d[p]

        u = u_dash / np.linalg.norm(u_dash)
        v = v_dash / np.linalg.norm(v_dash)
        w = w_dash / np.linalg.norm(w_dash)

        z_dash = np.cross(u, v) + np.cross(v, w) + np.cross(w, u)
        z = z_dash / np.linalg.norm(z_dash)

        oop_coord = z.dot(u)

        if gradient:
            grad = dq_oop(*coords3d[m], *coords3d[n], *coords3d[o], *coords3d[p])
            grad = grad.reshape(4, 3)
            row = np.zeros_like(coords3d)
            row[m, :] = grad[0]
            row[n, :] = grad[1]
            row[o, :] = grad[2]
            row[p, :] = grad[3]
            row = row.flatten()
            return oop_coord, row

        return oop_coord

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_oop(*coords3d[indices].flatten())