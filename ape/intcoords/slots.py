import abc
import logging

from math import exp, sin
import numpy as np

from ape.intcoords.elem_data import COVALENT_RADII as CR
from ape.intcoords.derivatives import d2q_b, d2q_a, dq_lb, d2q_lb, dq_ld, d2q_ld, d2q_d, dq_oop, d2q_oop

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