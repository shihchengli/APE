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

from pysisyphus.linalg import svd_inv
from pysisyphus.intcoords import Stretch, Torsion
from pysisyphus.intcoords.update import transform_int_step
from pysisyphus.intcoords.eval import (
    eval_primitives,
    check_primitives,
)
from pysisyphus.intcoords.setup import setup_redundant, get_primitives, PrimTypes
from pysisyphus.intcoords.valid import 

import pybel

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

def geo_to_pybel_mol(atoms, cart_coords):
    xyz = getXYZ(atoms, cart_coords)
    natom = len(atoms)
    xyz = str(natom) + '\n\n' + xyz
    PYMol = pybel.readstring('xyz', xyz)
    return PYMol

def get_bond_indices(atoms, cart_coords, imaginary_bonds=None):
    PYMol = geo_to_pybel_mol(atoms, cart_coords)
    OBMol = PYMol.OBMol
    reactant_bond = sorted(
        [(bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, bond.GetBondOrder())
            for bond in pybel.ob.OBMolBondIter(OBMol)]
    )
    bond_indices = [sorted(np.array([bond[0],bond[1]])) for bond in reactant_bond]

    # Add bond indices of imaginary bonds
    if imaginary_bonds is not None:
        for bond_indice in imaginary_bonds:
            bond_indice = sorted([bond_indice[0] - 1, bond_indice[1] - 1])
            if bond_indice not in bond_indices:
                bond_indices.append(bond_indice)

    bond_indices = np.array(sorted(bond_indices))
    return bond_indices

def get_RedundantCoords(label, atoms, cart_coords, rotors_dict=None, nHcap=0, natoms_adsorbate=0, imaginary_bonds=None, save_log=True):

    def connect_fragments(atoms, cart_coords, internal, bond_indices, save_log=True):
        internal.bond_indices = get_bond_indices(atoms, cart_coords)
        coords3d = cart_coords.reshape(-1, 3)
        # Condensed distance matrix
        cdm = pdist(coords3d)
        # Merge bond index sets into fragments
        bond_ind_sets = [frozenset(bi) for bi in bond_indices]
        fragments = internal.merge_fragments(bond_ind_sets)

        # Look for unbonded single atoms and create fragments for them.
        bonded_set = set(tuple(bond_indices.flatten()))
        unbonded_set = set(range(len(internal.atoms))) - bonded_set
        fragments.extend(
            [frozenset((atom, )) for atom in unbonded_set]
        )
        internal.fragments = fragments

        # Check if there are any disconnected fragments. If there are some
        # create interfragment bonds between all of them.
        if len(fragments) != 1:
            interfragment_inds = internal.connect_fragments(cdm, fragments)
            if save_log:
                logging.info('Add interfragment bonds between {}'.format([(ind[0] + 1, ind[1] + 1) for ind in interfragment_inds]))
            bond_indices = np.concatenate((bond_indices, interfragment_inds))
        return bond_indices

    def set_primitive_indices(internal, bond_indices, define_prims=None):
        stretches, bends, dihedrals = internal.sort_by_prim_type(define_prims)
        internal.bond_indices = bond_indices
        internal.bending_indices = list()
        internal.set_bending_indices(bends)
        internal.dihedral_indices = list()
        internal.set_dihedral_indices(dihedrals)
        dihedral_indices = internal.dihedral_indices
        if rotors_dict is not None and rotors_dict != []:
            pivots_list = [set([rotors_dict[i]['pivots'][0] - 1,
                            rotors_dict[i]['pivots'][1] - 1])
                            for i in rotors_dict]
            
            scan_indices_set = set()
            for i in rotors_dict:
                scan = rotors_dict[i]['scan']
                scan_indices_set.add((scan[0] - 1, scan[1] - 1, scan[2] - 1, scan[3] - 1))
            
            new_dihedral_indices = []
            for ind in dihedral_indices:
                if set(ind[1:3]) not in pivots_list:
                    new_dihedral_indices.append(ind)
            
            new_dihedral_indices.extend(list(scan_indices_set))
            internal.dihedral_indices = np.array(new_dihedral_indices)
    
    # For QMMM system
    if natoms_adsorbate != 0:
        # Fragment 1: Adsorbate
        internal_1 = RedundantCoords(atoms[:natoms_adsorbate], cart_coords[:natoms_adsorbate * 3])
        bond_indices_1 = get_bond_indices(atoms[:natoms_adsorbate], cart_coords[:natoms_adsorbate * 3])
        bond_indices_1 = connect_fragments(atoms[:natoms_adsorbate], cart_coords[:natoms_adsorbate * 3], internal_1, bond_indices_1, save_log=False)

        # Fragment 2: Active site w/ hydrogen caps
        internal_2 = RedundantCoords(atoms[natoms_adsorbate:], cart_coords[natoms_adsorbate * 3:])
        bond_indices_2 = get_bond_indices(atoms[natoms_adsorbate:], cart_coords[natoms_adsorbate * 3:])
        bond_indices_2 = connect_fragments(atoms[natoms_adsorbate:], cart_coords[natoms_adsorbate * 3:], internal_2, bond_indices_2, save_log=False)

        # User defined imaginary bonds
        if imaginary_bonds is not None:
            imaginary_bonds_indice = [sorted([bond_indice[0] - 1, bond_indice[1] - 1]) for bond_indice in imaginary_bonds]
            bond_indices = np.concatenate((bond_indices_1, bond_indices_2 + natoms_adsorbate, imaginary_bonds_indice))
        else:
            bond_indices = np.concatenate((bond_indices_1, bond_indices_2 + natoms_adsorbate))

        # Concatecate to get the bond indices of QMMM system
        bond_indices = np.unique(bond_indices, axis=0)
    else:
        bond_indices = get_bond_indices(atoms, cart_coords, imaginary_bonds)
    
    # Setup RedundantCoords object
    internal = RedundantCoords(atoms, cart_coords)
    internal.nHcap = nHcap
    bond_indices = connect_fragments(atoms, cart_coords, internal, bond_indices, save_log=save_log)
    set_primitive_indices(internal, bond_indices)
    internal._prim_internals = internal.calculate(cart_coords)
    internal._prim_coords = np.array([pc.val for pc in internal._prim_internals])
    invalid_bends_list = internal.invalid_bends_list

    # Add dummy atoms to handle molecules with nearly linear bend
    if invalid_bends_list != []:
        if save_log:
            logging.info("Didn't create bend {0} for {1}".format([(bend[0] + 1, bend[1] + 1, bend[2] + 1) for bend in invalid_bends_list], label))
        addHcap = AddHcap(cart_coords, bond_indices, invalid_bends_list, save_log)
        cart_coords, new_primes, new_nHcap = addHcap.add_Hcap_xyzs()
        atoms = atoms + ['H'] * new_nHcap
        internal = RedundantCoords(atoms, cart_coords)
        internal.nHcap = nHcap + new_nHcap
        internal.number_of_dummy_atom = new_nHcap
        stretches, bends, dihedrals = internal.sort_by_prim_type(new_primes)
        bond_indices = np.concatenate((bond_indices,stretches))
        bond_indices = connect_fragments(atoms, cart_coords, internal, bond_indices, save_log=save_log)
        define_primes = bends + dihedrals
        set_primitive_indices(internal, bond_indices, define_prims=define_primes)
        internal._prim_internals = internal.calculate(cart_coords)
        internal._prim_coords = np.array([pc.val for pc in internal._prim_internals])

    return internal

class RedundantCoords:
    def __init__(
        self,
        atoms,
        coords3d,
        bond_factor=1.3,
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
    ):
        self.atoms = atoms
        self.coords3d = np.reshape(coords3d, (-1, 3)).copy()
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

        self._B_prim = None
        # Lists for the other types of primitives will be created afterwards.
        # Linear bends may have been disabled, so we create the list here.
        self.linear_bend_indices = list()
        self.logger = logging.getLogger("internal_coords")

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

        bonds = len(self.bond_indices)
        bends = len(self.bending_indices) + len(self.linear_bend_indices)
        dihedrals = len(self.dihedral_indices)
        assert bonds + bends + dihedrals == len(self.primitives)
        self._bonds_slice = slice(bonds)
        self._bends_slice = slice(bonds, bonds + bends)
        self._dihedrals_slice = slice(bonds + bends, bonds + bends + dihedrals)
        self.backtransform_counter = 0

        self.nHcap = None
        self.number_of_dummy_atom = None
        self.shift_pi = list()

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
        return [self.bond_indices, self.bending_indices, self.dihedral_indices]

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
    def B_indices(self):
        """Wilson B-Matrix indices"""
        return [c.inds.tolist() for c in self.calculate(self.cart_coords)]

    @property
    def B(self):
        """Wilson B-Matrix"""
        return self.B_prim

    def inv_B(self, B):
        return B.T.dot(svd_inv(B.dot(B.T), thresh=self.svd_inv_thresh, hermitian=True))
        # return B.T.dot(self.pinv(B.dot(B.T)))

    def inv_Bt(self, B):
        return svd_inv(B.dot(B.T), thresh=self.svd_inv_thresh, hermitian=True).dot(B)
        # return self.pinv(B.dot(B.T)).dot(B)

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
        per_type = {
            2: list(),
            3: list(),
            4: list(),
            "linear_bend": list(),
            "hydrogen_bond": list(),
        }
        for type_, *indices in typed_prims:
            key = len(indices)
            if type_ in (linear_bend_types):
                key = "linear_bend"
            per_type[key].append(indices)

            # Also keep hydrogen bonds
            if type_ == PrimTypes.HYDROGEN_BOND:
                per_type["hydrogen_bond"].append(indices)

        self.bond_indices = per_type[2]
        self.bending_indices = per_type[3]
        self.dihedral_indices = per_type[4]
        self.linear_bend_indices = per_type["linear_bend"]
        self.hydrogen_bond_indices = per_type["hydrogen_bond"]

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

    def update_internals(self, new_cartesians, prev_internals):
        new_internals = self.calculate(new_cartesians, attr="val")
        internal_diffs = np.array(new_internals - prev_internals)
        bond, bend, dihedrals = self.prim_indices

        dihedral_diffs = internal_diffs[-len(dihedrals):]
        # Find differences that are shifted by 2*pi
        shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
        org = dihedral_diffs.copy()
        new_dihedrals = new_internals[-len(dihedrals):]
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])
        new_internals[-len(dihedrals):] = new_dihedrals
        return new_internals

    def transform_int_step(self, dq_in, ensure_convergence=True):
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
        bonds, bends, dihedrals = self.prim_indices
        for i, bond in enumerate(bonds):
            q = self.prim_coords[i]
            intco_lbls.append('R' + str(tuple(bond + 1)).replace(" ", ""))
            qShow_orig.append(q_orig[i])
            qShow_target.append(q_target[i])
            dqShow.append(q - qShow_orig[i])
            qShow_final.append(q)
        for i, bend in enumerate(bends):
            q = self.prim_coords[len(bonds) + i] * 180 / np.pi
            intco_lbls.append('B' + str(tuple(bend + 1)).replace(" ", ""))
            qShow_orig.append(q_orig[len(bonds) + i] * 180 / np.pi)
            qShow_target.append(q_target[len(bonds) + i] * 180 / np.pi)
            dqShow.append(q - qShow_orig[len(bonds) + i])
            qShow_final.append(q)
        for i, dihedral in enumerate(dihedrals):
            q = self.prim_coords[len(bonds) + len(bends) + i] * 180 / np.pi
            intco_lbls.append('D' + str(tuple(dihedral + 1)).replace(" ", ""))
            qShow_orig.append(q_orig[len(bonds) + len(bends) + i] * 180 / np.pi)
            qShow_target.append(q_target[len(bonds) + len(bends) + i] * 180 / np.pi)
            dqShow.append(q - qShow_orig[len(bonds) + len(bends) + i])
            qShow_final.append(q)        
        
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

        bond, bend, dihedrals = self.prim_indices
        # for i in set(self.shift_pi):
        #        step[len(bond)+i] *= -1
        target_bends = q_target[len(bond):-(len(dihedrals))]
        for i, target_bend in enumerate(target_bends):
            bendi = tuple(bend[i] + 1)
            if target_bend > np.pi:
                # TODO solve target_bend > np.pi situation
                # target_bends[i] = 2*np.pi - target_bends[i]
                # self.shift_pi.append(i)
                raise Exception('A sampling bending angel of {} is over 180°.'.format(bendi))
            elif target_bend <= 0:
                raise Exception('A sampling bending angel of {} is below 0°.'.format(bendi))

        B_prim = self.B_prim
        Bt_inv_prim = np.linalg.pinv(B_prim.dot(B_prim.T)).dot(B_prim)

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
            new_q = self.update_internals(geom, prev_q)
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

        self._prim_coords = np.array(new_q)
        self.cart_coords = geom
        
        dx = (geom - prev_geom)
        if self.number_of_dummy_atom is not None:
            dx = dx[:-self.number_of_dummy_atom * 3]

        return bt_converged, dx

    def __str__(self):
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"

    def get_intco_log(self):
        log = "\t-------Internal Coordinate-------\n"
        log += "\t -------------------------------\n"
        log += "\t Coordinate                Value\n"
        log += "\t ----------                -----\n"
        bonds, bends, dihedrals = self.prim_indices
        for i, bond in enumerate(bonds):
            bond_string = str(tuple(bond + 1)).replace(" ", "")
            value = self.prim_coords[i]
            log += '\t R{:15s}={:>14.6f}\n'.format(bond_string, value)
        for i, bend in enumerate(bends):
            bend_string = str(tuple(bend + 1)).replace(" ", "")
            value = self.prim_coords[len(bonds) + i] / np.pi * 180
            log += '\t B{:15s}={:>14.6f}\n'.format(bend_string, value)
        for i, dihedral in enumerate(dihedrals):
            dihedral_string = str(tuple(dihedral + 1)).replace(" ", "")
            value = self.prim_coords[len(bonds) + len(bends) + i] / np.pi * 180
            log += '\t D{:15s}={:>14.6f}\n'.format(dihedral_string, value)
        log += "\t -------------------------------\n"
        return log

###############################################################################

class AddHcap(object):
    """
    Add dummy atoms to handle linear molecules or molecules with nearly linear bend.
    """
    def __init__(self, cart_coords, bond_indices, invalid_bends_list, save_log=True):
        self.cart_coords = cart_coords
        self.bond_indices = bond_indices
        self.invalid_bends_list = invalid_bends_list
        self.save_log = save_log

    def add_Hcap_xyzs(self):
        """
        Find the set of xyz of the dummy atoms.
        """
        if self.save_log:
            logging.info('Adding dummy atoms...')
        invalid_bends_list = self.invalid_bends_list
        nHcap = len(invalid_bends_list)
        self.new_cart_coords = self.cart_coords.copy()
        self.new_primes = list()
        for i, bend in enumerate(invalid_bends_list):
            terminal1, central, terminal2 = bend
            self.ind = central
            self.bend = bend
            Hxyz_guess = self.new_cart_coords[self.ind * 3:self.ind * 3 + 3] + np.array([1.09, 0, 0])
            result = minimize(self.objectiveFunction, Hxyz_guess, method='SLSQP',
                            constraints=[
                            {'type': 'eq', 'fun': self.constraintFunction1},
                            {'type': 'eq', 'fun': self.constraintFunction2}
                            ])
            Hxyzs = result.x
            self.new_cart_coords = np.concatenate((self.new_cart_coords, Hxyzs), axis=None)
            
            dummy_atom_ind = len(self.cart_coords) // 3 + i
            self.new_primes.extend([[central, dummy_atom_ind],
                                    [terminal1, central, dummy_atom_ind],
                                    [terminal2, central, dummy_atom_ind],
                                    [terminal1, central, dummy_atom_ind, terminal2]])
            if self.save_log:
                logging.info('Create a improper dihedral of ({0}, {1}, {2}, {3})'.format(terminal1 + 1, central + 1, dummy_atom_ind + 1, terminal2 + 1))
        return self.new_cart_coords, self.new_primes, nHcap
    
    def objectiveFunction(self, Hxyzs):
        """
        Sum of the distance between dummy atom and other atoms.
        """
        val = 0
        for i, xyz in enumerate(self.new_cart_coords.reshape(-1, 3)):
            val += np.sqrt(np.sum((xyz - Hxyzs[0:3]) ** 2))
        return -val
    
    def constraintFunction1(self, Hxyzs):
        """
        The distance between dummy atom and the central atom of the chosen bend is 1.09 Å.
        """
        Hxyz = Hxyzs[0:3]
        center = self.cart_coords[self.ind * 3:self.ind * 3 + 3]
        distance = np.sqrt(np.sum((center - Hxyz) ** 2))
        return distance - 1.09
    
    def constraintFunction2(self, Hxyzs):
        """
        Let the vector from the chosen atom to the dummy atom is perpendicular to the bond vector.
        """
        atomB_ind, atomA_ind, atomC_ind = self.bend
        bond_vector = self.cart_coords[atomB_ind * 3:atomB_ind * 3 + 3] - self.cart_coords[atomA_ind * 3:atomA_ind * 3 + 3]
        bond_vector /= np.linalg.norm(bond_vector)
        A2H_vector = Hxyzs[0:3] - self.cart_coords[atomA_ind * 3:atomA_ind * 3 + 3]
        A2H_vector /= np.linalg.norm(A2H_vector)
        val = bond_vector.dot(A2H_vector)
        return val

        
