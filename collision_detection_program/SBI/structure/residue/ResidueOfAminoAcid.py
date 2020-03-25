"""
AminoAcid

author: jbonet
date:   02/2013

@oliva's lab
"""

from .          import Residue,        ResidueOfNucleotide
from SBI.data   import aminoacids3to1, aminoacids1to3, aminoacids_polarity_boolean, aminoacids_surface
from SBI        import SBIglobals

import numpy         as np
import scipy.spatial as sp


class ResidueOfAminoAcid(Residue):
    """
    A {AminoAcid} collects a series of {AminoAtoms}
    """
    def __init__(self, number, version, Rtype, mode):
        """
        @type  number: Integer
        @param number: Residue number

        @type  version: Char
        @param version: Optional char used on pdbs to change count

        @type  type: String
        @param type: Residue type

        @type  mode: String
        @param mode: Residue mode: ATOM or HETATM
        """
        super(ResidueOfAminoAcid, self).__init__(number = number, version = version, Rtype = Rtype, mode = mode)

        self._ca       = None
        self._cb       = None
        self._n        = None
        self._c        = None
        self._o        = None
        
        self._dssp     = None

    #
    # ATTRIBUTES
    #
    @property
    def ca(self):
    	"""
    	Returns the C-alpha {AminoAtom}
    	@rtype: {AminoAtom}
    	"""
    	return self._ca

    @property
    def cb(self):
    	"""
    	Returns the C-beta {AminoAtom}
    	@rtype: {AminoAtom}
    	"""
    	return self._cb

    @property
    def n(self):
        return self._n

    @property
    def c(self):
        return self._c

    @property
    def o(self):
        return self._o

    @property
    def backbone(self):
    	"""
    	Returns a list of the backbone {AminoAtom}s
    	@rtype: List of {AminoAtom}s
    	"""
    	return self._backbone_atoms

    @property
    def single_letter(self):
    	"""
    	Returns the AminoAcid identifier as a single letter code
    	@rtype: String
    	"""
    	return aminoacids3to1[self.type]

    @property
    def standard_type(self):
    	"""
    	For some AminoAcid are HETATM, we can return the standard AminoAcid
    	@rtype: String
    	"""
    	return aminoacids1to3[self.single_letter]

    @property
    def dssp(self):
        """
        @rtype: {DSSP}
        """
        return self._dssp

    @dssp.setter
    def dssp(self, value):
        """
        Sets {DSSP}
        """
        self._dssp = value

    @property
    def secondary_structure(self):
        if self._dssp is None:
            raise AttributeError("To call secondary structure DSSP needs to be calculated")

        return self._dssp.secondary_structure

    @property
    def accessibility(self):
        if self._dssp is None:
            raise AttributeError("To call accessibility DSSP needs to be calculated")

        return self._dssp.accessibility

    @property
    def accessibility10(self):
        if self._dssp is None:
            raise AttributeError("To call accessibility DSSP needs to be calculated")

        return self._dssp.accessibility10

    @property
    def accessibilitycoded(self):
        if self._dssp is None:
            raise AttributeError("To call accessibility DSSP needs to be calculated")

        return self._dssp.accesscode

    @property
    def exposed(self):
        if self._dssp is None:
            raise AttributeError("To call exposition DSSP needs to be calculated")

        return self._dssp.exposed

    @property
    def exposed_text(self):
        if self.exposed: return 'E'
        else:            return 'B'

    @property
    def surface_area(self):
        return aminoacids_surface[self.single_letter]

    #
    # BOOLEANS
    #
    @property
    def has_ca(self):
        """
        Checks if the AminoAcid has C-alpha
        @rtype: Boolean
        """
        return self._ca is not None

    @property
    def has_cb(self):
        """
        Checks if the AminoAcid has C-beta
        @rtype: Boolean
        """
        return self._cb is not None

    @property
    def has_full_backbone(self):
        return self._ca is not None and self._c is not None \
           and self._n  is not None and self._o is not None

    @property
    def is_only_ca(self):
    	"""
    	Some times the structure only contains the C-alpha
    	@rtype: Boolean
    	"""
    	return self.has_ca and len(self) == 1

    @property
    def is_only_backbone(self):
    	"""
    	The structure does not contain side chain
    	@rtype: Boolean
    	"""
    	return len(self.backbone) == len(self)

    @property
    def is_polar(self):
        return aminoacids_polarity_boolean[self.single_letter]

    #
    # METHODS
    #
    def normalize(self):
        self._type = self.standard_type
        self._mode = 'ATOM'
        
    def remove_side_chain(self):
        # self._atoms = self._backbone_atoms
        self._sidechain_atoms = []
        self._sidechain_coordinates = None

    def _distance_ca(self, aminoacid):
        """
        Calculates the ca-ca distance between two residues

        @type  aminoacid: {AminoAcid}
        @param aminoacid: AminoAcid to which we want to find the distance

        @rtype: {Contact}
        @rtype: None if one of them does not have ca
        """
        if not self.has_ca or not aminoacid.has_ca:
            # return self._distance_min(residue = aminoacid)
            return (None, None, -1)

        return (self.ca, aminoacid.ca, self.ca.distance(atom = aminoacid.ca))


    def _distance_cb(self, aminoacid):
        """
        Calculates the cb-cb distance between two residues.
        If one does not has cb, reverts to ca-cb or cb-ca.
        If none have cb, reverts to ca-ca.

        @type  aminoacid: {AminoAcid}
        @param aminoacid: AminoAcid to which we want to find the distance

        @rtype: {Contact}
        @rtype: None if one of them does not have ca
        """
        if (not self.has_ca and not self.has_cb) or (not aminoacid.has_ca and not aminoacid.has_cb):
            # return (self._distance_min(residue = aminoacid))
            return (None, None, -1)
        if not self.has_cb and not aminoacid.has_cb:
            # return self._distance_ca(aminoacid = aminoacid)
            return (None, None, -1)
        elif not self.has_cb:
            return (self.ca, aminoacid.cb, self.ca.distance(atom = aminoacid.cb))
        elif not aminoacid.has_cb:
            return (self.cb, aminoacid.ca, self.cb.distance(atom = aminoacid.ca))
        else:
            return (self.cb, aminoacid.cb, self.cb.distance(atom = aminoacid.cb))

    def _distance_cb_backbone(self, nucleotide):
        if self.has_cb:
            cb_coord = np.array(np.zeros(3))
            cb_coord = np.vstack((cb_coord, self.cb.coordinates))
            cb_coord = np.delete(cb_coord, 0, 0)
            atom     = self.cb
        elif self.has_ca:
            cb_coord = np.array(np.zeros(3))
            cb_coord = np.vstack((cb_coord, self.ca.coordinates))
            cb_coord = np.delete(cb_coord, 0, 0)
            atom     = self.ca
        else:
            return (None, None, -1)
        if nucleotide._backbone_coordinates is not None:
            n_coord = nucleotide._backbone_coordinates
            n_atoms = nucleotide._backbone_atoms
        else:
            n_coord = nucleotide._sidechain_coordinates
            n_atoms = nucleotide._sidechain_atoms
        SBIglobals.alert('deepdebug', self, '\tAminoAcid coordenate {0}'.format(cb_coord))
        SBIglobals.alert('deepdebug', self, '\tNucleotide coordenates {0}'.format(n_coord))
        distances = sp.distance.cdist(cb_coord, n_coord)
        index     = np.unravel_index(distances.argmin(), distances.shape)
        return (atom, n_atoms[index[1]], distances.min())

    #
    # OVERRIDE PARENT'S METHODS
    #
    def add_atom(self, atom):
        """
        Adds a new {AminoAtom} to the {AminoAcid}
        This includes:
        	filling the all_coordinates attribute (parent)
        	assign CA, CB and backbone
        @type  atom: {AminoAtom}
        @param atom: New {AminoAtom} added to the {AminoAcid}
        """
        super(ResidueOfAminoAcid, self).add_atom(atom = atom)

        if atom.is_Calpha:
        	self._ca = atom
        elif atom.is_Cbeta:
        	self._cb = atom
        elif atom.is_C:
            self._c  = atom
        elif atom.is_N:
            self._n  = atom
        elif atom.is_O:
            self._o  = atom

    def distance(self, residue, dist_type):
        """
        Calculates the distance between two residues

        @type  residue: {AminoAcid}
        @param residue: AminoAcid to which we want to find the distance

        @type  dist_type: String
        @param dist_type: Type of distance to evaluate
        @limit dist_type: Accepted: distance.types (ca, cb, min, geometric)

        @rtype: (AT1, AT2, DISTANCE)
        """
        if not isinstance(residue, ResidueOfAminoAcid):
            if isinstance(residue, ResidueOfNucleotide):
                if dist_type.strip().lower() != 'cbbackbone':
                    return super(ResidueOfAminoAcid, self).distance(residue = residue, dist_type = dist_type)
                else:
                    return self._distance_cb_backbone(nucleotide = residue)
            else:
                return super(ResidueOfAminoAcid, self).distance(residue = residue, dist_type = dist_type)
        if dist_type.strip().lower() == 'ca':
            return self._distance_ca(aminoacid = residue)
        if dist_type.strip().lower() == 'cb':
            return self._distance_cb(aminoacid = residue)
        else:
            return super(ResidueOfAminoAcid, self).distance(residue = residue, dist_type = dist_type)

    def follows(self, residue):
        return residue.is_followed(self)

    def is_followed(self, residue):
        c0, ca0, n0 = self.c,    self.ca,    self.n
        c1, ca1, n1 = residue.c, residue.ca, residue.n

        if c0  is not None and n1  is not None: return c0.distance(n1)   <= 1.5 # measured distance is 1.3 
        if ca0 is not None and ca1 is not None: return ca0.distance(ca1) <= 4   # measured distance around 3.8
        if c0  is not None and c1  is not None: return c0.distance(c1)   <= 4
        if n0  is not None and n1  is not None: return n0.distance(n1)   <= 4

    #
    # OVERRIDE DEFAULT OBJECT METHODS
    #
    def __repr__(self):
        repre = []
        if self._dssp is None:
            repre.append('<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(self))
        else:
            repre.append('<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}, ({0.secondary_structure}, {0.exposed_text})]>'.format(self))
        for atom in self.atoms:
            repre.append("\t" + repr(atom))
        return "\n".join(repre)