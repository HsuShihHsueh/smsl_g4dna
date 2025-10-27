from MDAnalysis import *
from MDAnalysis import coordinates
from typing import Union, Optional, List, Dict

## 1. u.trajectory.timeseries
## rewrite from:
## https://github.com/MDAnalysis/mdanalysis/blob/59f4e395178240d5e3f36088d7a4d98ddd0e3607/package/MDAnalysis/coordinates/base.py#L984
import numpy as np      
from tqdm import tqdm  

class _ProtoReader():
    def timeseries(self, asel: Optional['AtomGroup']=None,
                   start: Optional[int]=None, stop: Optional[int]=None,
                   step: Optional[int]=None,
                   order: Optional[str]='fac') -> np.ndarray:
        """Edit by Allen 2024/01/22
        Return a subset of coordinate data for an AtomGroup

        Parameters
        ----------
        asel : AtomGroup (optional)
            The :class:`~MDAnalysis.core.groups.AtomGroup` to read the
            coordinates from. Defaults to ``None``, in which case the full set
            of coordinate data is returned.
        start :  int (optional)
            Begin reading the trajectory at frame index `start` (where 0 is the
            index of the first frame in the trajectory); the default
            ``None`` starts at the beginning.
        stop : int (optional)
            End reading the trajectory at frame index `stop`-1, i.e, `stop` is
            excluded. The trajectory is read to the end with the default
            ``None``.
        step : int (optional)
            Step size for reading; the default ``None`` is equivalent to 1 and
            means to read every frame.
        order : str (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates)

        See Also
        --------
        :class:`MDAnalysis.coordinates.memory`


        .. versionadded:: 2.4.0
        """
        start, stop, step = self.check_slice_indices(start, stop, step)
        nframes = len(range(start, stop, step))

        if asel is not None:
            if len(asel) == 0:
                raise ValueError(
                    "Timeseries requires at least one atom to analyze")
            atom_numbers = asel.indices
            natoms = len(atom_numbers)
        else:
            natoms = self.n_atoms
            atom_numbers = np.arange(natoms)

        # allocate output array in 'fac' order
        coordinates = np.empty((nframes, natoms, 3), dtype=np.float32)
        for i, ts in enumerate(tqdm(self[start:stop:step])): ## only add this!!
            coordinates[i, :] = ts.positions[atom_numbers]

        # switch axes around
        default_order = 'fac'
        if order != default_order:
            try:
                newidx = [default_order.index(i) for i in order]
            except ValueError:
                raise ValueError(f"Unrecognized order key in {order}, "
                                 "must be permutation of 'fac'")

            try:
                coordinates = np.moveaxis(coordinates, newidx, [0, 1, 2])
            except ValueError:
                errmsg = ("Repeated or missing keys passed to argument "
                          f"`order`: {order}, each key must be used once")
                raise ValueError(errmsg)
        return coordinates
    
coordinates.XTC.XTCReader.timeseries = _ProtoReader.timeseries

## 2 u.sele_atoms('').write
## rewrite from:
## https://github.com/MDAnalysis/mdanalysis/blob/e776f124c18c41f8990b487e8557d3ad82fe7d1f/package/MDAnalysis/coordinates/CRD.py#L177
import itertools
import numpy as np
import warnings

class _CRDWriter():
    def write(self, selection, frame=None):
        from MDAnalysis.exceptions import NoDataError
        from MDAnalysis.lib import util

        """Write selection at current trajectory frame to file.
        Edit by Allen 2024/06/01

        Parameters
        ----------
        selection : AtomGroup
             group of atoms to be written
        frame : int (optional)
             Move the trajectory to frame `frame`; by default, write
             the current frame.

        """
        print('Edit by Allen 2024/06/01')
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single PDB (?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.positions  # can write from selection == Universe (Issue 49)

        n_atoms = len(atoms)
        # Detect which format string we're using to output (EXT or not)
        # *len refers to how to truncate various things,
        # depending on output format!
        if n_atoms > 0: ## edit 99999 -> 0
            at_fmt = self.fmt['ATOM_EXT']
            serial_len = 10
            resid_len = 8
            totres_len = 10
        else:
            at_fmt = self.fmt['ATOM']
            serial_len = 5
            resid_len = 4
            totres_len = 5

        # Check for attributes, use defaults for missing ones
        attrs = {}
        missing_topology = []
        for attr, default in (
                ('resnames', itertools.cycle(('UNK',))),
                # Resids *must* be an array because we index it later
                ('resids', np.ones(n_atoms, dtype=int)),
                ('names', itertools.cycle(('X',))),
                ('tempfactors', itertools.cycle((0.0,))),
        ):
            try:
                attrs[attr] = getattr(atoms, attr)
            except (NoDataError, AttributeError):
                attrs[attr] = default
                missing_topology.append(attr)
        # ChainIDs - Try ChainIDs first, fall back to Segids
        try:
            attrs['chainIDs'] = atoms.chainIDs
        except (NoDataError, AttributeError):
            # try looking for segids instead
            try:
                attrs['chainIDs'] = atoms.segids
            except (NoDataError, AttributeError):
                attrs['chainIDs'] = itertools.cycle(('',))
                missing_topology.append(attr)
        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "".format(miss=', '.join(missing_topology)))

        with util.openany(self.filename, 'w') as crd:
            # Write Title
            crd.write(self.fmt['TITLE'].format(
                frame=frame, where=u.trajectory.filename))
            crd.write("*\n")

            # Write NUMATOMS
            if n_atoms > 0: ## edit 99999 -> 0
                crd.write(self.fmt['NUMATOMS_EXT'].format(n_atoms))
            else:
                crd.write(self.fmt['NUMATOMS'].format(n_atoms))

            # Write all atoms

            current_resid = 1
            resids = attrs['resids']
            for i, pos, resname, name, chainID, resid, tempfactor in zip(
                    range(n_atoms), coor, attrs['resnames'], attrs['names'],
                    attrs['chainIDs'], attrs['resids'], attrs['tempfactors']):
                if not i == 0 and resids[i] != resids[i-1]:
                    current_resid += 1

                # Truncate numbers
                serial = util.ltruncate_int(i + 1, serial_len)
                resid = util.ltruncate_int(resid, resid_len)
                current_resid = util.ltruncate_int(current_resid, totres_len)

                crd.write(at_fmt.format(
                    serial=serial, totRes=current_resid, resname=resname,
                    name=name, pos=pos, chainID=chainID,
                    resSeq=resid, tempfactor=tempfactor))
                
coordinates.CRD.CRDWriter.write = _CRDWriter.write