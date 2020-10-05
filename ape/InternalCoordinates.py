#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844

from collections import namedtuple
from functools import reduce
import itertools as it
import logging
import typing
import copy

import attr
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize

from pysisyphus.constants import BOHR2ANG
from pysisyphus.elem_data import VDW_RADII, COVALENT_RADII as CR
from pysisyphus.intcoords.derivatives import d2q_b, d2q_a, d2q_d

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

def get_cov_radii_sum_array(atoms, coords):
    coords3d = coords.reshape(-1, 3)
    atom_indices = list(it.combinations(range(len(coords3d)), 2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = list()
    for i, j in atom_indices:
        atom1 = atoms[i].lower()
        atom2 = atoms[j].lower()
        cov_rad_sums.append(CR[atom1] + CR[atom2])
    cov_rad_sums = np.array(cov_rad_sums)
    return cov_rad_sums

@attr.s(auto_attribs=True)
class PrimitiveCoord:
    inds : typing.List[int]
    val : float
    grad : np.ndarray


class RedundantCoords:

    RAD_175 = 3.05432619
    RAD_5 = 0.08726646
    BEND_MIN_DEG = 15
    BEND_MAX_DEG = 170

    def __init__(self, atoms, cart_coords, bond_factor=1.3,
                 prim_indices=None, define_prims=None):
        self.atoms = atoms
        self.cart_coords = cart_coords
        self.bond_factor = bond_factor
        self.define_prims = define_prims

        self.bond_indices = list()
        self.bending_indices = list()
        self.dihedral_indices = list()
        self.hydrogen_bond_indices = list()

        if prim_indices is None:
            self.set_primitive_indices(self.define_prims)
        else:
            to_arr = lambda _: np.array(list(_), dtype=int)
            bonds, bends, dihedrals = prim_indices
            # We accept all bond indices. What could possibly go wrong?! :)
            self.bond_indices = to_arr(bonds)
            valid_bends = [inds for inds in bends
                           if self.is_valid_bend(inds)]
            self.bending_indices = to_arr(valid_bends)
            valid_dihedrals = [inds for inds in dihedrals if
                               self.is_valid_dihedral(inds)]
            self.dihedral_indices = to_arr(valid_dihedrals)
        self._prim_internals = self.calculate(self.cart_coords)
        self._prim_coords = np.array([pc.val for pc in self._prim_internals])

        self.nHcap = None
        self.number_of_dummy_atom = None
        self.shift_pi = list()

    def log(self, message):
        #logger = logging.getLogger("internal_coords")
        #logger.debug(message)
        pass

    @property
    def prim_indices(self):
        return [self.bond_indices, self.bending_indices, self.dihedral_indices]

    @property
    def prim_indices_set(self):
        return set([tuple(prim_ind) for prim_ind in it.chain(*self.prim_indices)])

    @property
    def prim_coords(self):
        if self._prim_coords is None:
           self._prim_coords = np.array(
                            [pc.val for pc in self.calculate(self.cart_coords)]
            )
        return self._prim_coords

    @property
    def coords(self):
        return self.prim_coords

    @property
    def coord_indices(self):
        ic_ind_tuples = [tuple(ic.inds) for ic in self._prim_internals]
        return {ic_inds: i for i, ic_inds in enumerate(ic_ind_tuples)}

    @property
    def dihed_start(self):
        return len(self.bond_indices) + len(self.bending_indices)

    def get_index_of_prim_coord(self, prim_ind):
        """Index of primitive internal for the given atom indices.

        TODO: simplify this so when we get a prim_ind of len 2
        (bond) we don't have to check the bending and dihedral indices."""
  
        prim_ind_set = set(prim_ind)
        indices = [i for i, pi in enumerate(it.chain(*self.prim_indices))
                 if set(pi) == prim_ind_set]
        index = None
        try:
            index = indices[0]
        except IndexError:
            self.log(f"Primitive internal with indices {prim_ind} "
                      "is not defined!")
        return index

    @property
    def c3d(self):
        return self.cart_coords.reshape(-1, 3)

    @property
    def B_prim(self):
        """Wilson B-Matrix"""
        return np.array([c.grad for c in self.calculate(self.cart_coords)])

    @property
    def B_indices(self):
        """Wilson B-Matrix indices"""
        return [c.inds.tolist() for c in self.calculate(self.cart_coords)]

    @property
    def B(self):
        """Wilson B-Matrix"""
        return self.B_prim

    @property
    def Bt_inv(self):
        """Transposed generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return np.linalg.pinv(B.dot(B.T)).dot(B)

    @property
    def B_inv(self):
        """Generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return B.T.dot(np.linalg.pinv(B.dot(B.T)))

    @property
    def P(self):
        """Projection matrix onto B. See [1] Eq. (4)."""
        return self.B.dot(self.B_inv)

    def transform_forces(self, cart_forces):
        """Combination of Eq. (9) and (11) in [1]."""
        return self.Bt_inv.dot(cart_forces)

    def get_K_matrix(self, int_gradient=None):
        assert len(int_gradient) == len(self._prim_internals)
        size_ = self.cart_coords.size
        if int_gradient is None:
            return np.zeros((size_, size_))

        dg_funcs = {
            2: d2q_b,
            3: d2q_a,
            4: d2q_d,
        }
        def grad_deriv_wrapper(inds):
            coords_flat = self.c3d[inds].flatten()
            dgrad = dg_funcs[len(inds)](*coords_flat)
            return dgrad

        K_flat = np.zeros(size_ * size_)
        for pc, int_grad_item in zip(self._prim_internals, int_gradient):
            # Contract with gradient
            try:
                dg = int_grad_item * grad_deriv_wrapper(pc.inds)
            except (ValueError, ZeroDivisionError) as err:
                self.log( "Error in calculation of 2nd derivative of primitive "
                         f"internal {pc.inds}."
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
            cart_inds = list(it.chain(*[range(3*i,3*i+3) for i in pc.inds]))
            flat_inds = [row*size_ + col for row, col in it.product(cart_inds, cart_inds)]
            K_flat[flat_inds] += dg
        K = K_flat.reshape(size_, size_)
        return K

    def transform_hessian(self, cart_hessian, int_gradient=None):
        if int_gradient is None:
            self.log("Supplied 'int_gradient' is None. K matrix will be zero, "
                     "so derivatives of the Wilson-B-matrix are neglected in "
                     "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv.dot(cart_hessian-K).dot(self.B_inv)

    def project_hessian(self, H, shift=1000):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return P.dot(H).dot(P) + shift*(np.eye(P.shape[0]) - P)

    def project_vector(self, vector):
        """Project supplied vector onto range of B."""
        P = self.P
        return self.P.dot(vector)

    def set_rho(self):
        # TODO: remove this as it is already in optimizers/guess_hessians
        atoms = [a.lower() for a in self.atoms]
        cov_radii = np.array([CR[a.lower()] for a in atoms])
        rref = np.array([r1+r2
                         for r1, r2 in it.combinations(cov_radii, 2)])
        coords3d = self.cart_coords.reshape(-1, 3)
        cdm = pdist(coords3d)
        self.rho = squareform(np.exp(-cdm/rref + 1))

    def get_initial_hessian(self):
        # TODO: remove this as it is already in optimizers/guess_hessians
        self.set_rho()
        k_dict = {
            2: 0.35,
            3: 0.15,
            4: 0.005,
        }
        k_diag = list()
        for primitive in self._prim_internals:
            rho_product = 1
            for i in range(primitive.inds.size-1):
                i1, i2 = primitive.inds[i:i+2]
                rho_product *= self.rho[i1, i2]
            k_diag.append(k_dict[len(primitive.inds)] * rho_product)
        return np.diagflat(k_diag)

    def merge_fragments(self, fragments):
        """Merge a list of sets."""
        # Hold the final fragments that can't be merged further, as they
        # contain distinct atoms.
        merged = list()
        while len(fragments) > 0:
            popped = fragments.pop(0)
            # Look for an intersection between the popped unmerged fragment
            # and the remaining unmerged fragments.
            for frag in fragments:
                if popped & frag:
                    fragments.remove(frag)
                    # If a intersecting unmerged fragment is found merge
                    # both fragments and append them at the end.
                    fragments.append(popped | frag)
                    break
            else:
                # Add the unmerged fragment into merged if it doesn't
                # intersect with any other unmerged fragment.
                merged.append(popped)
        return merged

    def connect_fragments(self, cdm, fragments):
        """Determine the smallest interfragment bond for a list
        of fragments and a condensed distance matrix."""
        dist_mat = squareform(cdm)
        interfragment_indices = list()
        for frag1, frag2 in it.combinations(fragments, 2):
            arr1 = np.array(list(frag1))[None,:]
            arr2 = np.array(list(frag2))[:,None]
            indices = [(i1, i2) for i1, i2 in it.product(frag1, frag2)]
            distances = np.array([dist_mat[ind] for ind in indices])
            min_index = indices[distances.argmin()]
            interfragment_indices.append(sorted(min_index))
        # Or as Philipp proposed: two loops over the fragments and only
        # generate interfragment distances. So we get a full matrix with
        # the original indices but only the required distances.
        return interfragment_indices

    def set_hydrogen_bond_indices(self, bond_indices):
        coords3d = self.cart_coords.reshape(-1, 3)
        tmp_sets = [frozenset(bi) for bi in bond_indices]
        # Check for hydrogen bonds as described in [1] A.1 .
        # Find hydrogens bonded to small electronegative atoms X = (N, O
        # F, P, S, Cl).
        hydrogen_inds = [i for i, a in enumerate(self.atoms)
                         if a.lower() == "h"]
        x_inds = [i for i, a in enumerate(self.atoms)
                  if a.lower() in "n o f p s cl".split()]
        hydrogen_bond_inds = list()
        for h_ind, x_ind in it.product(hydrogen_inds, x_inds):
            as_set = set((h_ind, x_ind))
            if not as_set in tmp_sets:
                continue
            # Check if distance of H to another electronegative atom Y is
            # greater than the sum of their covalent radii but smaller than
            # the 0.9 times the sum of their van der Waals radii. If the
            # angle X-H-Y is greater than 90° a hydrogen bond is asigned.
            y_inds = set(x_inds) - set((x_ind, ))
            for y_ind in y_inds:
                y_atom = self.atoms[y_ind].lower()
                cov_rad_sum = CR["h"] + CR[y_atom]
                distance = self.calc_stretch(coords3d, (h_ind, y_ind))
                vdw = 0.9 * (VDW_RADII["h"] + VDW_RADII[y_atom])
                angle = self.calc_bend(coords3d, (x_ind, h_ind, y_ind))
                if (cov_rad_sum < distance < vdw) and (angle > np.pi/2):
                    self.hydrogen_bond_indices.append((h_ind, y_ind))
                    self.log("Added hydrogen bond between {h_ind} and {y_ind}")
        self.hydrogen_bond_indices = np.array(self.hydrogen_bond_indices)

    def set_bond_indices(self, define_bonds=None, factor=None):
        """
        Default factor of 1.3 taken from [1] A.1.
        Gaussian uses somewhat less, like 1.2, or different radii than we do.
        """
        bond_factor = factor if factor else self.bond_factor
        coords3d = self.cart_coords.reshape(-1, 3)
        # Condensed distance matrix
        cdm = pdist(coords3d)
        # Generate indices corresponding to the atom pairs in the
        # condensed distance matrix cdm.
        atom_indices = list(it.combinations(range(len(coords3d)),2))
        atom_indices = np.array(atom_indices, dtype=int)
        cov_rad_sums = get_cov_radii_sum_array(self.atoms, self.cart_coords)
        cov_rad_sums *= bond_factor
        bond_flags = cdm <= cov_rad_sums
        bond_indices = atom_indices[bond_flags]


        if define_bonds:
            bond_indices = np.concatenate(((bond_indices, define_bonds)), axis=0)

        self.bare_bond_indices = bond_indices


        # Look for hydrogen bonds
        self.set_hydrogen_bond_indices(bond_indices)
        if self.hydrogen_bond_indices.size > 0:
            bond_indices = np.concatenate((bond_indices,
                                           self.hydrogen_bond_indices))

        # Merge bond index sets into fragments
        bond_ind_sets = [frozenset(bi) for bi in bond_indices]
        fragments = self.merge_fragments(bond_ind_sets)

        # Look for unbonded single atoms and create fragments for them.
        bonded_set = set(tuple(bond_indices.flatten()))
        unbonded_set = set(range(len(self.atoms))) - bonded_set
        fragments.extend(
            [frozenset((atom, )) for atom in unbonded_set]
        )
        self.fragments = fragments

        # Check if there are any disconnected fragments. If there are some
        # create interfragment bonds between all of them.
        if len(fragments) != 1:
            interfragment_inds = self.connect_fragments(cdm, fragments)
            bond_indices = np.concatenate((bond_indices, interfragment_inds))

        self.bond_indices = np.unique(bond_indices, axis=0)

    def are_parallel(self, vec1, vec2, angle_ind=None, thresh=1e-6):
        dot = max(min(vec1.dot(vec2), 1), -1)
        rad = np.arccos(dot)#vec1.dot(vec2))
        # angle > 175°
        if abs(rad) > self.RAD_175:
            # self.log(f"Nearly linear angle {angle_ind}: {np.rad2deg(rad)}")
            ind_str = f" ({angle_ind})" if (angle_ind is not None) else ""
            self.log(f"Nearly linear angle{ind_str}: {np.rad2deg(rad)}")
        return abs(rad) > (np.pi - thresh)

    def sort_by_central(self, set1, set2):
        """Determines a common index in two sets and returns a length 3
        tuple with the central index at the middle position and the two
        terminal indices as first and last indices."""
        central_set = set1 & set2
        union = set1 | set2
        assert len(central_set) == 1
        terminal1, terminal2 = union - central_set
        (central, ) = central_set
        return (terminal1, central, terminal2), central

    def is_valid_bend(self, bend_ind):
        val = self.calc_bend(self.c3d, bend_ind)
        deg = np.rad2deg(val)
        return self.BEND_MIN_DEG <= deg <= self.BEND_MAX_DEG

    def set_bending_indices(self, define_bends=None):
        bond_sets = {frozenset(bi) for bi in self.bond_indices}
        self.invalid_bends_list = list()
        for bond_set1, bond_set2 in it.combinations(bond_sets, 2):
            union = bond_set1 | bond_set2
            if len(union) == 3:
                as_tpl, _ = self.sort_by_central(bond_set1, bond_set2)
                if not self.is_valid_bend(as_tpl):
                    self.invalid_bends_list.append(as_tpl)
                    self.log(f"Didn't create bend ({as_tpl})")
                             # f" with value of {deg:.3f}°")
                    continue
                self.bending_indices.append(as_tpl)
        self.bending_indices = np.array(self.bending_indices, dtype=int)

        if define_bends:
            bis = np.concatenate(( (self.bending_indices, define_bends)), axis=0)
            self.bending_indices = np.unique(bis, axis=0)

    def is_valid_dihedral(self, dihedral_ind, thresh=1e-6):
        # Check for linear atoms
        first_angle = self.calc_bend(self.c3d, dihedral_ind[:3])
        second_angle = self.calc_bend(self.c3d, dihedral_ind[1:])
        pi_thresh = np.pi - thresh
        return ((abs(first_angle) < pi_thresh)
                and (abs(second_angle) < pi_thresh)
        )

    def set_dihedral_indices(self, define_dihedrals=None):
        dihedral_sets = list()
        self.dihedral_list = list()
        def set_dihedral_index(dihedral_ind):
            dihedral_set = set(dihedral_ind)
            self.dihedral_list.append(dihedral_ind)  #contain repeated dihedral indices
            # Check if this dihedral is already present
            if dihedral_set in dihedral_sets:
                return
            # Assure that the angles are below 175° (3.054326 rad)
            if not self.is_valid_dihedral(dihedral_ind, thresh=0.0873):
                self.log("Skipping generation of dihedral "
                               f"{dihedral_ind} as some of the the atoms "
                                "are linear."
                )
                return
            self.dihedral_indices.append(dihedral_ind)
            dihedral_sets.append(dihedral_set)

        improper_dihedrals = list()
        coords3d = self.cart_coords.reshape(-1, 3)
        for bond, bend in it.product(self.bond_indices, self.bending_indices):
            central = bend[1]
            bend_set = set(bend)
            bond_set = set(bond)
            # Check if the two sets share one common atom. If not continue.
            intersect = bend_set & bond_set
            if len(intersect) != 1:
                continue
            # When the common atom is a terminal atom of the bend, that is
            # it's not the central atom of the bend, we create a
            # proper dihedral. Before we create any improper dihedrals we
            # create these proper dihedrals.
            if central not in bond_set:
                # The new terminal atom in the dihedral is the one that
                # doesn' intersect.
                terminal = tuple(bond_set - intersect)[0]
                intersecting_atom = tuple(intersect)[0]
                if intersecting_atom == bend[0]:
                    dihedral_ind = [terminal] + bend.tolist()
                else:
                    dihedral_ind = bend.tolist() + [terminal]
                set_dihedral_index(dihedral_ind)
            # If the common atom is the central atom we try to form an out
            # of plane bend / improper torsion. They may be created later on.
            else:
                fourth_atom = list(bond_set - intersect)
                dihedral_ind = bend.tolist() + fourth_atom
                # This way dihedrals may be generated that contain linear
                # atoms and these would be undefinied. So we check for this.
                dihed = self.calc_dihedral(coords3d, dihedral_ind)
                if not np.isnan(dihed):
                    improper_dihedrals.append(dihedral_ind)
                else:
                    self.log("Dihedral {dihedral_ind} is undefinied. Skipping it!")

        # Now try to create the remaining improper dihedrals.
        if (len(self.atoms) >= 4) and (len(self.dihedral_indices) == 0):
            for improp in improper_dihedrals:
                set_dihedral_index(improp)
            self.log("Permutational symmetry not considerd in "
                            "generation of improper dihedrals.")

        self.dihedral_indices = np.array(self.dihedral_indices)

        if define_dihedrals:
            dis = np.concatenate(((self.dihedral_indices, define_dihedrals)), axis=0)
            self.dihedral_indices = np.unique(dis, axis=0)

    def sort_by_prim_type(self, to_sort):
        by_prim_type = [[], [], []]
        if to_sort is None:
            to_sort = list()
        for item in to_sort:
            len_ = len(item)
            by_prim_type[len_-2].append(item)
        return by_prim_type

    def set_primitive_indices(self, define_prims=None):
        stretches, bends, dihedrals = self.sort_by_prim_type(define_prims)
        self.set_bond_indices(stretches)
        self.set_bending_indices(bends)
        self.set_dihedral_indices(dihedrals)

    def calculate(self, coords, attr=None):
        coords3d = coords.reshape(-1, 3)
        def per_type(func, ind):
            val, grad = func(coords3d, ind, True)
            return PrimitiveCoord(ind, val, grad)
        self.bonds = list()
        self.bends = list()
        self.dihedrals = list()
        for ind in self.bond_indices:
            bonds = per_type(self.calc_stretch, ind)
            self.bonds.append(bonds)
        for ind in self.bending_indices:
            bend = per_type(self.calc_bend, ind)
            self.bends.append(bend)
        for ind in self.dihedral_indices:
            dihedral = per_type(self.calc_dihedral, ind)
            self.dihedrals.append(dihedral)
        int_coords = self.bonds + self.bends + self.dihedrals
        if attr:
            return np.array([getattr(ic,attr) for ic in int_coords])
        return int_coords

    def calculate_val_diffs(self, coords1, coords2):
        vals1 = np.array(self.calculate(coords1, attr="val"))
        vals2 = np.array(self.calculate(coords2, attr="val"))
        return vals1-vals2

    def calc_stretch(self, coords3d, bond_ind, grad=False):
        n, m = bond_ind
        bond = coords3d[m] - coords3d[n]
        bond_length = np.linalg.norm(bond)
        if grad:
            bond_normed = bond / bond_length
            row = np.zeros_like(coords3d)
            # 1 / -1 correspond to the sign factor [1] Eq. 18
            row[m,:] =  bond_normed
            row[n,:] = -bond_normed
            row = row.flatten()
            return bond_length, row
        return bond_length

    def calc_bend(self, coords3d, angle_ind, grad=False):
        m, o, n = angle_ind
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        angle_rad = np.arccos(u.dot(v))
        if grad:
            # Eq. (24) in [1]
            if self.are_parallel(u, v, angle_ind):
                tmp_vec = np.array((1, -1, 1))
                par = self.are_parallel(u, tmp_vec) and self.are_parallel(v, tmp_vec)
                tmp_vec = np.array((-1, 1, 1)) if par else tmp_vec
                w_dash = np.cross(u, tmp_vec)
            else:
                w_dash = np.cross(u, v)
            w_norm = np.linalg.norm(w_dash)
            w = w_dash / w_norm
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

    def calc_dihedral(self, coords3d, dihedral_ind, grad=False, cos_tol=1e-9):
        m, o, p, n = dihedral_ind
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

        # Restrict cos_dihed to [-1, 1]
        if cos_dihed >= 1 - cos_tol:
            dihedral_rad = 0
        elif cos_dihed <= -1 + cos_tol:
            dihedral_rad = np.arccos(-1)
        else:
            dihedral_rad = np.arccos(cos_dihed)

        if dihedral_rad != np.pi:
            # wxv = np.cross(w, v)
            # if wxv.dot(u) < 0:
            if vxw.dot(u) < 0:
                dihedral_rad *= -1
        if grad:
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

    def transform_int_step(self, dq_in, ensure_convergence=False):
        """
        This is always done in primitive internal coordinates so care
        has to be taken that the supplied step is given in primitive internal
        coordinates.
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
                    if cnt == 5:
                        logging.warning(
                            "\tUnable to back-transform even 1/10th of the desired step rigorously."
                            + "\tContinuing with best (small) step")
                    else:
                        self._prim_coords = q_orig
                        self.cart_coords = geom_orig
            
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
    
    def get_active_set(self, B, thresh=1e-6):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)

        nonzero_inds = np.abs(eigvals) > thresh
        active_eigvals = eigvals[nonzero_inds]
        return eigvectors[:,nonzero_inds]

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

        
