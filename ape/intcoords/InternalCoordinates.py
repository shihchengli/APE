#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844 , Pulay 1992

import math
import itertools as it
import logging
import copy

import numpy as np
from numpy.linalg.linalg import LinAlgError

from ape.intcoords.linalg import svd_inv
from ape.intcoords.slots import Stretch, Torsion
from ape.intcoords.eval import eval_primitives,check_primitives
from ape.intcoords.setup import setup_redundant, get_primitives, PrimTypes
from ape.intcoords.valid import check_typed_prims
from ape.intcoords.update import update_internals
from ape.intcoords.constants import BOHR2ANG
from ape.intcoords import nifty
from ape.exceptions import SamplingError

def getXYZ(atoms, cart_coords):
    """
    Return a string of the molecule in the XYZ file format.
    """
    natom = len(atoms)
    xyz = ''
    for i in range(natom):
        xyz += '{:s}       {:.10f}     {:.10f}     {:.10f}'.format(atoms[i],cart_coords[3*i],cart_coords[3*i+1],cart_coords[3*i+2])
        if i != natom-1: xyz += '\n'
    return xyz

def get_RedundantCoords(label, atoms, cart_coords, rotors_dict=None, nHcap=0, add_hrdrogen_bonds=True, addcart=False, addtr=False, add_interfragment_bonds=False):

    def set_typed_prims(internal, rotors_dict):

        pivots_list = [set([rotors_dict[i]['pivots'][0] - 1,
                        rotors_dict[i]['pivots'][1] - 1])
                        for i in rotors_dict]

        dihedral_indices = internal.dihedral_indices
        new_typed_prims = list()
        for typed_prim in internal.typed_prims:
            type_, *indices = typed_prim
            if len(indices) == 4:
                if set(indices[1:3]) in pivots_list:
                    print(indices)
                    continue
            new_typed_prims.append(typed_prim)
        
        for rotor in rotors_dict.values():
            scan = rotor['scan']
            typed_prim = (PrimTypes.PROPER_DIHEDRAL, scan[0] - 1, scan[1] - 1, scan[2] - 1, scan[3] - 1)
            new_typed_prims.append(typed_prim)

        return new_typed_prims
    
    internal = RedundantCoords(atoms, cart_coords, add_hrdrogen_bonds=add_hrdrogen_bonds, addcart=addcart, addtr=addtr, add_interfragment_bonds=add_interfragment_bonds)
    if rotors_dict != None and rotors_dict != []:
        typed_prims = set_typed_prims(internal, rotors_dict)
        internal = RedundantCoords(atoms, cart_coords, add_hrdrogen_bonds=add_hrdrogen_bonds, typed_prims=typed_prims, add_interfragment_bonds=True)

    internal.nHcap = nHcap

    return internal

class RedundantCoords:
    def __init__(
        self,
        atoms,
        coords3d,
        bond_factor=1.2,
        typed_prims=None,
        define_prims=None,
        bonds_only=False,
        check_bends=True,
        rebuild=True,
        bend_min_deg=15,
        dihed_max_deg=175.0,
        lb_min_deg=175.0,
        weighted=False,
        min_weight=0.3,
        # Corresponds to a threshold of 1e-7 for eigenvalues of G, as proposed by
        # Pulay in [5].
        svd_inv_thresh=3.16e-4,
        add_hrdrogen_bonds = True,
        addcart=False,
        addtr=False,
        add_interfragment_bonds=False,
    ):
        self.atoms = atoms
        self.coords3d = np.reshape(coords3d, (-1, 3)).copy()
        self.cart_coords = self.coords3d.reshape(-1,)
        self.bond_factor = bond_factor
        self.define_prims = define_prims
        self.bonds_only = bonds_only
        self.check_bends = check_bends
        self.rebuild = rebuild
        self.bend_min_deg = bend_min_deg
        self.dihed_max_deg = dihed_max_deg
        self.lb_min_deg = lb_min_deg
        self.weighted = weighted
        self.min_weight = float(min_weight)
        assert self.min_weight > 0.0, "min_weight must be a positive rational!"
        self.svd_inv_thresh = svd_inv_thresh
        self.add_hrdrogen_bonds = add_hrdrogen_bonds
        self.addcart = addcart
        self.addtr = addtr
        self.add_interfragment_bonds = add_interfragment_bonds

        self._B_prim = None
        # Lists for the other types of primitives will be created afterwards.
        # Linear bends may have been disabled, so we create the list here.
        self.linear_bend_indices = list()
        self.logger = None#logging.getLogger("internal_coords")

        if self.weighted:
            self.log(
                "Coordinate weighting requested, min_weight="
                f"{self.min_weight:.2f}. Calculating bond factor."
            )
            # Screening function is
            #   ρ(d) = exp(-(d/sum_cov_rad - 1)
            #
            # Swart proposed a min_weight of ρ(d) = 0.3. With this we can
            # calculate the appropriate factor for the bond detection.
            # d = (1 - ln(0.3)) * sum_cov_rad
            # bond_factor = (1 - ln(0.3)) ≈ 2.204
            #
            # The snippet below prints weights and corresponding bond_factors.
            # [f"{w:.2f}: {1-np.log(w):.4f}" for w in np.linspace(0.3, 1, 25)]
            self.bond_factor = -math.log(self.min_weight) + 1
        self.log(f"Using a factor of {self.bond_factor:.6f} for bond detection.")
        self.log(f"Using svd_inv_thresh={self.svd_inv_thresh:.4e} for inversions.")

        # Set up primitive coordinate indices
        if typed_prims is None:
            self.set_primitive_indices(
                self.atoms,
                self.coords3d,
            )
        # Use supplied typed_prims
        else:
            self.log(f"{len(typed_prims)} primitives were supplied. Checking them.")
            valid_typed_prims = check_typed_prims(
                self.coords3d,
                typed_prims,
                bend_min_deg=self.bend_min_deg,
                dihed_max_deg=self.dihed_max_deg,
                lb_min_deg=self.lb_min_deg,
                check_bends=self.check_bends,
            )
            self.log(
                f"{len(valid_typed_prims)} primitives are valid at the current Cartesians."
            )
            if len(valid_typed_prims) != len(typed_prims):
                self.log("Invalid primitives:")
                for i, invalid_prim in enumerate(set(typed_prims) - set(valid_typed_prims)):
                    self.log(f"\t{i:02d}: {invalid_prim}")
            self.typed_prims = valid_typed_prims
            self.set_inds_from_typed_prims(self.typed_prims)

        self.primitives = get_primitives(
            self.coords3d,
            self.typed_prims,
            logger=self.logger,
        )
        if self.bonds_only:
            self.bending_indices = list()
            self.dihedral_indices = list()
            self.linear_bend_indices = list()
            self.primitives = [
                prim for prim in self.primitives if isinstance(prim, Stretch)
            ]
        check_primitives(self.coords3d, self.primitives, logger=self.logger)

        self._prim_internals = self.eval(self.coords3d)
        self._prim_coords = np.array(
            [prim_int.val for prim_int in self._prim_internals]
        )

        carts = len(self.cart_indices)
        trans = len(self.translation_indices)
        rots = len(self.rotation_indices)
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices) + len(self.linear_bend_indices)
        dihedrals = len(self.dihedral_indices)
        assert carts + trans + rots + bonds + bends + dihedrals == len(self.primitives)
        self._carts_slice = slice(carts)
        self._trans_slice = slice(carts, carts + trans)
        self._rots_slice = slice(carts + trans, carts + trans + rots)
        self._bonds_slice = slice(carts + trans + rots, carts + trans + rots + bonds)
        self._bends_slice = slice(carts + trans + rots + bonds, carts + trans + rots + bonds + bends)
        self._dihedrals_slice = slice(carts + trans + rots + bonds + bends, carts + trans + rots + bonds + bends + dihedrals)

        self.nHcap = None

    def log(self, message):
        #logger.debug(message)
        pass

    @property
    def coords3d(self):
        return self._coords3d

    @coords3d.setter
    def coords3d(self, coords3d):
        self._coords3d = coords3d.reshape(-1, 3)
        self._B_prim = None
        self._prim_coords = None
        self._prim_internals = None

    @property
    def primitives(self):
        return self._primitives

    @primitives.setter
    def primitives(self, primitives):
        self._primitives = primitives

    @property
    def prim_indices(self):
        return [np.array(self.bond_indices), np.array(self.bending_indices), np.array(self.dihedral_indices)]

    @property
    def prim_indices_set(self):
        return set([tuple(prim_ind) for prim_ind in it.chain(*self.prim_indices)])

    @property
    def prim_internals(self):
        if self._prim_internals is None:
            self._prim_internals = self.eval(self.coords3d)
        return self._prim_internals

    @prim_internals.setter
    def prim_internals(self, prim_internals):
        self._prim_internals = prim_internals

    @property
    def prim_coords(self):
        return np.array([prim_int.val for prim_int in self.prim_internals])

    def return_inds(self, slice_):
        return np.array([prim_int.indices for prim_int in self.prim_internals[slice_]])

    @property
    def carts(self):
        return self.prim_internals[self._carts_slice]

    @property
    def bonds(self):
        return self.prim_internals[self._bonds_slice]

    @property
    def bends(self):
        return self.prim_internals[self._bends_slice]

    @property
    def dihedrals(self):
        return self.prim_internals[self._dihedrals_slice]

    @property
    def coords(self):
        return self.prim_coords

    @property
    def dihed_start(self):
        return len(self.bond_indices) + len(self.bending_indices)

    def get_index_of_prim_coord(self, prim_ind):
        """Index of primitive internal for the given atom indices."""
        prim_ind_set = set(prim_ind)
        for i, prim in enumerate(self.primitives):
            if set(prim.indices) == prim_ind_set:
                return i
        self.log(f"Primitive internal with indices {prim_ind} " "is not defined!")
        return None

    @property
    def B_prim(self):
        """Wilson B-Matrix"""
        if self._B_prim is None:
            self._B_prim = np.array([prim_int.grad for prim_int in self.prim_internals])

        return self._B_prim

    @property
    def B(self):
        """Wilson B-Matrix"""
        return self.B_prim

    def inv_B(self, B):
        try:
            return B.T.dot(np.linalg.pinv(B.dot(B.T)))
        except LinAlgError:
            logging.warning('LinAlgError appears. Using svd_inv_thresh={:.4e} for inversions.'.format(self.svd_inv_thresh))
            return B.T.dot(svd_inv(B.dot(B.T), thresh=self.svd_inv_thresh, hermitian=True))

    def inv_Bt(self, B):
        try:
            return np.linalg.pinv(B.dot(B.T)).dot(B)
        except LinAlgError:
            logging.warning('LinAlgError appears. Using svd_inv_thresh={:.4e} for inversions.'.format(self.svd_inv_thresh))
            return svd_inv(B.dot(B.T), thresh=self.svd_inv_thresh, hermitian=True).dot(B)

    @property
    def Bt_inv_prim(self):
        """Transposed generalized inverse of the primitive Wilson B-Matrix."""
        return self.inv_Bt(self.B_prim)

    @property
    def Bt_inv(self):
        """Transposed generalized inverse of the Wilson B-Matrix."""
        return self.inv_Bt(self.B)

    @property
    def B_inv_prim(self):
        """Generalized inverse of the primitive Wilson B-Matrix."""
        return self.inv_B(self.B_prim)

    @property
    def B_inv(self):
        """Generalized inverse of the Wilson B-Matrix."""
        return self.inv_B(self.B)

    @property
    def P(self):
        """Projection matrix onto B. See [1] Eq. (4)."""
        return self.B.dot(self.B_inv)

    def transform_forces(self, cart_forces):
        """Combination of Eq. (9) and (11) in [1]."""
        return self.Bt_inv.dot(cart_forces)

    def get_K_matrix(self, int_gradient=None):
        if int_gradient is not None:
            assert len(int_gradient) == len(self._primitives)

        size_ = self.coords3d.size
        if int_gradient is None:
            return np.zeros((size_, size_))

        K_flat = np.zeros(size_ * size_)
        coords3d = self.coords3d
        for primitive, int_grad_item in zip(self.primitives, int_gradient):
            # Contract with gradient
            val = np.rad2deg(primitive.calculate(coords3d))
            # self.log(f"K, {primitive}={val:.2f}°")
            # The generated code (d2q_d) seems unstable for these values...
            if isinstance(primitive, Torsion) and ((abs(val) < 1) or (abs(val) > 179)):
                self.log(f"Skipped 2nd derivative of {primitive} with val={val:.2f}°")
                continue
            # 2nd derivative of normal, but linear, bends is undefined.
            try:
                dg = int_grad_item * primitive.jacobian(coords3d)
            except (ValueError, ZeroDivisionError) as err:
                self.log(
                    "Error in calculation of 2nd derivative of primitive "
                    f"internal {primitive.indices}."
                )
                continue
            # Depending on the type of internal coordinate dg is a flat array
            # of size 36 (stretch), 81 (bend) or 144 (torsion).
            #
            # An internal coordinate contributes to an element K[j, k] of the
            # K matrix if the cartesian coordinate indices j and k belong to an
            # atom that contributes to the respective internal coordinate.
            #
            # As for now we build up the K matrix as flat array. To add the dg
            # entries at the appropriate places in K_flat we have to calculate
            # the corresponding flat indices of dg in K_flat.
            cart_inds = list(
                it.chain(*[range(3 * i, 3 * i + 3) for i in primitive.indices])
            )
            flat_inds = [
                row * size_ + col for row, col in it.product(cart_inds, cart_inds)
            ]
            K_flat[flat_inds] += dg
        K = K_flat.reshape(size_, size_)
        return K

    def log_int_grad_msg(self, int_gradient):
        if int_gradient is None:
            self.log(
                "Supplied 'int_gradient' is None. K matrix will be zero, "
                "so derivatives of the\nWilson-B-matrix are neglected in "
                "Hessian transformation."
            )

    def transform_hessian(self, cart_hessian, int_gradient=None):
        """Transform Cartesian Hessian to internal coordinates."""
        self.log_int_grad_msg(int_gradient)
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv_prim.dot(cart_hessian - K).dot(self.B_inv_prim)

    def backtransform_hessian(self, redund_hessian, int_gradient=None):
        """Transform Hessian in internal coordinates to Cartesians."""
        self.log_int_grad_msg(int_gradient)
        K = self.get_K_matrix(int_gradient)
        return self.B.T.dot(redund_hessian).dot(self.B) + K

    def project_hessian(self, H, shift=1000):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return P.dot(H).dot(P) + shift * (np.eye(P.shape[0]) - P)

    def project_vector(self, vector):
        """Project supplied vector onto range of B."""
        return self.P.dot(vector)

    def set_inds_from_typed_prims(self, typed_prims):
        linear_bend_types = (PrimTypes.LINEAR_BEND, PrimTypes.LINEAR_BEND_COMPLEMENT)
        translation_types = (PrimTypes.TRANSLATION_X, PrimTypes.TRANSLATION_Y, PrimTypes.TRANSLATION_Z)
        rotation_types = (PrimTypes.ROTATION_A, PrimTypes.ROTATION_B, PrimTypes.ROTATION_C)
        per_type = {
            1: list(),
            2: list(),
            3: list(),
            4: list(),
            "linear_bend": list(),
            "hydrogen_bond": list(),
            "translation": list(),
            "rotation": list(),
        }
        B_indices = list()
        for type_, *indices in typed_prims:
            key = len(indices)
            if type_ in (linear_bend_types):
                key = "linear_bend"
            elif type_ in (translation_types):
                key = "translation"
            elif type_ in (rotation_types):
                key = "rotation"
            per_type[key].append(indices)

            # Also keep hydrogen bonds
            if type_ == PrimTypes.HYDROGEN_BOND:
                per_type["hydrogen_bond"].append(indices)
            B_indices.append(indices)

        self.cart_indices = per_type[1]
        self.bond_indices = per_type[2]
        self.bending_indices = per_type[3]
        self.dihedral_indices = per_type[4]
        self.linear_bend_indices = per_type["linear_bend"]
        self.hydrogen_bond_indices = per_type["hydrogen_bond"]
        self.translation_indices = per_type["translation"]
        self.rotation_indices = per_type["rotation"]
        self.B_indices = B_indices

        # TODO
        # self.fragments = coord_info.fragments

    def set_primitive_indices(
        self,
        atoms,
        coords3d,
    ):
        coord_info = setup_redundant(
            atoms,
            coords3d,
            factor=self.bond_factor,
            define_prims=self.define_prims,
            min_deg=self.bend_min_deg,
            dihed_max_deg=self.dihed_max_deg,
            lb_min_deg=self.lb_min_deg,
            min_weight=self.min_weight if self.weighted else None,
            add_hrdrogen_bonds=self.add_hrdrogen_bonds,
            addcart=self.addcart,
            addtr=self.addtr,
            add_interfragment_bonds=self.add_interfragment_bonds,
            logger=self.logger,
        )

        self.typed_prims = coord_info.typed_prims
        self.set_inds_from_typed_prims(self.typed_prims)

        self.fragments = coord_info.fragments

    def eval(self, coords3d, attr=None):
        prim_internals = eval_primitives(coords3d, self.primitives)

        if attr is not None:
            return np.array(
                [getattr(prim_internal, attr) for prim_internal in prim_internals]
            )

        return prim_internals

    def transform_int_step(self, dq_in, ensure_convergence=False):
        """
        Transformation is done in primitive internals, so int_step must be given
        in primitive internals and not in DLC!
        """
        logging.info('\n\tBack-transformation to cartesian coordinates...')
        q_orig = self.prim_coords.copy()
        geom_orig = self.cart_coords.copy()
        q_target = q_orig + dq_in

        dq = dq_in.copy()
        conv = False  # is back-transformation converged?

        if ensure_convergence:
            cnt = -1

            while not conv:
                cnt += 1
                if cnt > 0:
                    logging.info("\tReducing step-size by a factor of {:d}.".format(2 * cnt))
                    dq[:] = dq_in / (2.0 * cnt)
                
                conv, dx = self.back_transformation(dq)

                if not conv:
                    self._prim_coords = q_orig
                    self.cart_coords = geom_orig
                    if cnt == 5:
                        logging.warning(
                            "\tUnable to back-transform even 1/10th of the desired step rigorously."
                            + "\tQuitting with previous geometry.")
                        conv, dx = self.back_transformation(dq)
                        break
            
                if conv and cnt > 0:  # We were able to take a modest step.  Try to complete it.
                    logging.info(
                        "\tAble to take a small step; trying another partial back-transformations.\n")

                    for j in range(1, 2 * cnt):
                        logging.info("\tMini-step {:d} of {:d}.\n".format(j + 1, 2 * cnt))
                        dq[:] = dq_in / (2 * cnt)

                        conv, mdx = self.back_transformation(dq)
                        dx += mdx

                        if not conv:
                            self._prim_coords = q_orig
                            self.cart_coords = geom_orig                            
                            if cnt == 5:
                                logging.warning(
                                    "\tCouldn't converge this mini-step; quitting with previous geometry.\n")
                                # raise SamplingError('Couldn\'t converge to targeted internal coordinate even with 1/10th of the desired step.')
                                dq = dq_in.copy()
                                conv, dx = self.back_transformation(dq)
                                conv = True
                            break
        
        else:  # try to back-transform, but continue even if desired dq is not achieved
            conv, dx = self.back_transformation(dq)
        
        intco_lbls, qShow_orig, qShow_target, dqShow, qShow_final = [], [], [], [], []
        i = 0
        for type_, *indices in self.typed_prims:
            indices = np.array(indices)
            if type_ in (PrimTypes.TRANSLATION_X, PrimTypes.TRANSLATION_Y, PrimTypes.TRANSLATION_Z):
                if type_ == PrimTypes.TRANSLATION_X:
                    symbol = 'TranX'
                elif type_ == PrimTypes.TRANSLATION_Y:
                    symbol = 'TranY'
                elif type_ == PrimTypes.TRANSLATION_Z:
                    symbol = 'TranZ'
                q = self.prim_coords[i] * BOHR2ANG
                intco_lbls.append(symbol + '({})'.format(nifty.commadash(indices)))
                qShow_orig.append(q_orig[i] * BOHR2ANG)
                qShow_target.append(q_target[i] * BOHR2ANG)
                dqShow.append(q - qShow_orig[i])
                qShow_final.append(q)
            elif type_ in (PrimTypes.ROTATION_A, PrimTypes.ROTATION_B, PrimTypes.ROTATION_C):
                if type_ == PrimTypes.ROTATION_A:
                    symbol = 'RotA'
                elif type_ == PrimTypes.ROTATION_B:
                    symbol = 'RotB'
                elif type_ == PrimTypes.ROTATION_C:
                    symbol = 'RotC'
                q = self.prim_coords[i] * 180 / np.pi
                intco_lbls.append(symbol + '({})'.format(nifty.commadash(indices)))
                qShow_orig.append(q_orig[i] * 180 / np.pi)
                qShow_target.append(q_target[i] * 180 / np.pi)
                dqShow.append(q - qShow_orig[i])
                qShow_final.append(q)
            elif len(indices) == 1:
                if type_ == PrimTypes.CARTESIAN_X:
                    symbol = 'X'
                elif type_ == PrimTypes.CARTESIAN_Y:
                    symbol = 'Y'
                elif type_ == PrimTypes.CARTESIAN_Z:
                    symbol = 'Z'
                q = self.prim_coords[i] * BOHR2ANG
                intco_lbls.append(symbol + '({})'.format(indices[0] + 1))
                qShow_orig.append(q_orig[i] * BOHR2ANG)
                qShow_target.append(q_target[i] * BOHR2ANG)
                dqShow.append(q - qShow_orig[i])
                qShow_final.append(q)
            elif len(indices) == 2:
                q = self.prim_coords[i] * BOHR2ANG
                intco_lbls.append('R' + str(tuple(indices + 1)).replace(" ", ""))
                qShow_orig.append(q_orig[i] * BOHR2ANG)
                qShow_target.append(q_target[i] * BOHR2ANG)
                dqShow.append(q - qShow_orig[i])
                qShow_final.append(q)
            elif len(indices) == 3:
                if type_ == PrimTypes.LINEAR_BEND:
                    symbol = 'L'
                elif type_ == PrimTypes.LINEAR_BEND_COMPLEMENT:
                    symbol = 'C'
                else:
                    symbol = 'B'
                q = self.prim_coords[i] * 180 / np.pi
                intco_lbls.append(symbol + str(tuple(indices + 1)).replace(" ", ""))
                qShow_orig.append(q_orig[i] * 180 / np.pi)
                qShow_target.append(q_target[i] * 180 / np.pi)
                dqShow.append(q - qShow_orig[i])
                qShow_final.append(q)
            elif len(indices) == 4:
                q = self.prim_coords[i] * 180 / np.pi
                intco_lbls.append('D' + str(tuple(indices + 1)).replace(" ", ""))
                qShow_orig.append(q_orig[i] * 180 / np.pi)
                qShow_target.append(q_target[i] * 180 / np.pi)
                dqShow.append(q - qShow_orig[i])
                qShow_final.append(q)
            i += 1
        # Make sure final Dq is actual change
        frag_report = "\tReport of back-transformation: (au)\n"
        frag_report += "\n\t int                       q_final      q_target         Error\n"
        frag_report += "\t -------------------------------------------------------------\n"
        for i in range(len(dq_in)):
            frag_report += ("\t %-16s=%16.6f%14.6f%14.6f\n"
                           % (intco_lbls[i], qShow_final[i], qShow_target[i], (qShow_final[i] - qShow_target[i])))
        frag_report += "\t -------------------------------------------------------------\n"
        logging.debug(frag_report)

        coordinate_change_report = (
            "\n\t---Internal Coordinate Step in ANG or DEG, aJ/ANG or AJ/DEG ---\n")
        coordinate_change_report += (
            "\t -------------------------------------------------------------\n")
        coordinate_change_report += (
            "\t Coordinate               Previous        Change           New\n")
        coordinate_change_report += (
            "\t ----------               --------        ------        ------\n")
        for i in range(len(dq_in)):
            coordinate_change_report += ("\t %-16s=%16.6f%14.6f%14.6f\n"
                                        % (intco_lbls[i], qShow_orig[i], dqShow[i], qShow_final[i]))
        coordinate_change_report += (
            "\t -------------------------------------------------------------\n")
        logging.info(coordinate_change_report)
        return dx

    def back_transformation(self, dq, bt_dx_conv=1.0e-6, bt_max_iter=100):

        dx_rms_last = -1

        q_orig = self.prim_coords.copy()
        q_target = q_orig + dq        

        prev_geom = self.cart_coords.copy() # cart geometry to start each iter
        geom = self.cart_coords.copy()

        B_prim = self.B_prim
        Bt_inv_prim = self.Bt_inv_prim

        prev_q = q_orig

        bt_iter_continue = True
        bt_converged = False
        bt_iter_cnt = 0

        while bt_iter_continue:

            dx = Bt_inv_prim.T.dot(dq)

            # Frozen the positions of dummy atoms and hydrogen caps of QMMM system
            if self.nHcap != 0:
                dx[-(self.nHcap * 3):] = 0 

            # Update cartesian coordinates
            geom += dx
            dx_rms = np.sqrt(np.mean(dx ** 2))
    
            # Met convergence thresholds
            if dx_rms < bt_dx_conv:
                bt_converged = True
                bt_iter_continue = False
            # No further progress toward convergence
            elif (np.absolute(dx_rms - dx_rms_last) < 1.0e-7
                  or bt_iter_cnt >= bt_max_iter or dx_rms > 100.0):
                bt_converged = False
                bt_iter_continue = False

            dx_rms_last = dx_rms

            # Determine new internal coordinates
            dihedral_inds = np.array(
                [i for i, primitive in enumerate(self.primitives) if isinstance(primitive, Torsion)]
            )
            new_prim_ints = update_internals(
                new_coords3d=geom.reshape(-1, 3),
                old_internals=prev_q,
                primitives=self.primitives,
                dihedral_inds=dihedral_inds,
                check_dihedrals=self.rebuild,
                logger=self.logger,
            )
            new_q = [prim.val for prim in new_prim_ints]
            dq[:] = q_target - new_q

            dq_rms = np.sqrt(np.mean(dq ** 2))
            if bt_iter_cnt == 0 or dq_rms < best_dq_rms:  # short circuit evaluation
                best_cycle = (copy.deepcopy(geom), copy.deepcopy(new_q))
                best_dq_rms = dq_rms
            
            bt_iter_cnt += 1
            prev_q = new_q

        if bt_converged:
            logging.info("\tSuccessfully converged to displaced geometry.")
        else:
            logging.warning("\tUnable to completely converge to displaced geometry.")
        
        if dq_rms > best_dq_rms:
            # logging.warning("\tPrevious geometry is closer to target in internal coordinates,"
            #                 + " so using that one.\n")
            # logging.warning("\tBest geometry has RMS(Delta(q)) = %8.2e\n" % best_dq_rms)
            geom, new_q = best_cycle

        self.prim_internals = self.eval(geom.reshape(-1,3))
        self.cart_coords = geom
        self.coords3d = np.reshape(geom, (-1, 3))
        
        dx = (geom - prev_geom)

        return bt_converged, dx

    def __str__(self):
        carts = len(self.cart_indices)
        trans = len(self.translation_indices)
        rots = len(self.rotation_indices)
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({carts} carts, {trans} trans, {rots} rots, {bonds} bonds, {bends} bends, {dihedrals} dihedrals)"

    def get_intco_log(self):
        log = "\t-------Internal Coordinate-------\n"
        log += "\t -------------------------------\n"
        log += "\t Coordinate                Value\n"
        log += "\t ----------                -----\n"
        i = 0
        for type_, *indices in self.typed_prims:
            indices = np.array(indices)
            if type_ in (PrimTypes.TRANSLATION_X, PrimTypes.TRANSLATION_Y, PrimTypes.TRANSLATION_Z):
                if type_ == PrimTypes.TRANSLATION_X:
                    symbol = 'TranX'
                elif type_ == PrimTypes.TRANSLATION_Y:
                    symbol = 'TranY'
                elif type_ == PrimTypes.TRANSLATION_Z:
                    symbol = 'TranZ'
                tran_string = '({})'.format(nifty.commadash(indices))
                value = self.prim_coords[i] * BOHR2ANG
                log += '\t {}{:11s}={:>14.6f}\n'.format(symbol, tran_string, value)
            elif type_ in (PrimTypes.ROTATION_A, PrimTypes.ROTATION_B, PrimTypes.ROTATION_C):
                if type_ == PrimTypes.ROTATION_A:
                    symbol = 'RotA'
                elif type_ == PrimTypes.ROTATION_B:
                    symbol = 'RotB'
                elif type_ == PrimTypes.ROTATION_C:
                    symbol = 'RotC'
                rot_string = '({})'.format(nifty.commadash(indices))
                value = self.prim_coords[i] / np.pi * 180
                log += '\t {}{:12s}={:>14.6f}\n'.format(symbol, rot_string, value)
            elif len(indices) == 1:
                if type_ == PrimTypes.CARTESIAN_X:
                    symbol = 'X'
                elif type_ == PrimTypes.CARTESIAN_Y:
                    symbol = 'Y'
                elif type_ == PrimTypes.CARTESIAN_Z:
                    symbol = 'Z'
                cart_string = '({})'.format(indices[0] + 1)
                value = self.prim_coords[i] * BOHR2ANG
                log += '\t {}{:15s}={:>14.6f}\n'.format(symbol, cart_string, value)
            elif len(indices) == 2:
                bond_string = str(tuple(indices + 1)).replace(" ", "")
                value = self.prim_coords[i] * BOHR2ANG
                log += '\t R{:15s}={:>14.6f}\n'.format(bond_string, value)
            elif len(indices) == 3:
                if type_ == PrimTypes.LINEAR_BEND:
                    symbol = 'L'
                elif type_ == PrimTypes.LINEAR_BEND_COMPLEMENT:
                    symbol = 'C'
                else:
                    symbol = 'B'
                bend_string = str(tuple(indices + 1)).replace(" ", "")
                value = self.prim_coords[i] / np.pi * 180
                log += '\t {}{:15s}={:>14.6f}\n'.format(symbol, bend_string, value)
            elif len(indices) == 4:
                dihedral_string = str(tuple(indices + 1)).replace(" ", "")
                value = self.prim_coords[i] / np.pi * 180
                log += '\t D{:15s}={:>14.6f}\n'.format(dihedral_string, value)
            i += 1
        log += "\t -------------------------------\n"
        return log
