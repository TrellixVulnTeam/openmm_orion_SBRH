#############################################################################
# Copyright (C) 2019 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np

from openeye import (oechem,
                     oedepict,
                     oegrapheme)
import mdtraj as md

from datarecord import OEField

import os, contextlib

import glob

from tempfile import TemporaryDirectory


from oeommtools import utils as oeommutils

def GetCardinalOrderOfProteinResNums(mol):
    # make map of protein res nums to the residue cardinal order index
    resmap = {}
    currRes = -10000
    currIdx = -1
    for atom in mol.GetAtoms(oechem.OEIsBackboneAtom()):
        thisRes = oechem.OEAtomGetResidue(atom)
        resnum = thisRes.GetResidueNumber()
        if resnum != currRes:
            currIdx += 1
            currRes = resnum
            resmap[currRes] = currIdx
    return resmap, currIdx


def ExtractProtLigActsiteResNums(mol, fromLigCutoff=5.0):
    '''Extracts the protein and ligand from a single OEMol containing a protein-ligand
    complex plus other components. A list of protein residues within a cutoff distance
    from the ligand is also returned.
    Inputs:
        mol: The OEMol containing protein, ligand, and all other components
        fromLigCutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
    Returns:
        protein: An OEMol containing only the protein
        ligand: An OEMol containing only the ligand
        actSiteResNums: A list of integers, one per residue number for protein residues
            with any atom within fromLigCutoff distance of the ligand.'''
    # perceive residue hierarchy of total system
    if not oechem.OEHasResidues(mol):
        oechem.OEPerceiveResidues(mol, oechem.OEPreserveResInfo_All)
    # split the total system into components

    protein, ligand, water, excipients = oeommutils.split(mol, ligand_res_name="LIG")

    # use residue-based distance cutoff (in angstroms) from ligand to define an active site residue
    nn = oechem.OENearestNbrs(protein, fromLigCutoff)

    actSiteResNums = set()

    for nbrs in nn.GetNbrs(ligand):
        residue = oechem.OEAtomGetResidue(nbrs.GetBgn())
        actSiteResNums.add(residue.GetResidueNumber())

    return protein, ligand, actSiteResNums

#
# def ExtractAlignedProtLigTraj(mol, traj_filename, fromLigCutoff=5.0, skip=0):
#     '''Extracts the aligned protein trajectory and aligned ligand trajectory from
#     a MD trajectory of a larger system that includes other components (eg water).
#     The passed in OEMol must have the topology that matches the trajectory, and its xyz
#     coordinates are the reference for the alignment. The alignment is done on the
#     alpha carbons (atom name CA) of the active site residues within fromLigCutoff
#     from the ligand. Once the alignment is done, the protein and ligand trajectories
#     are each placed into a separate OEMol, one conformer per trajectory frame.
#     Inputs:
#         mol: An OEMol giving the topology for the trajectory and the reference xyz
#             coordinates for the alignment.
#         traj_filename: The filename of the hdf5-format MD trajectory or Gromacs .xtc file format
#         fromLigCutoff: The cutoff distance in angstroms to include protein residues
#             close to the ligand.
#         skip: number of frames to skip at the beginning of the trajectory.
#     Outputs:
#         protTraj: A multiconformer OEMol for the protein, one conformer per frame.
#         ligTraj: A multiconformer OEMol for the ligand, one conformer per frame.'''
#
#     void, traj_ext = os.path.splitext(traj_filename)
#
#     traj_dir = os.path.dirname(traj_filename)
#
#     # get the topology from 1st frame of the traj file
#     if traj_ext == '.h5':
#         topologyTraj = md.load_hdf5(traj_filename, frame=1)
#
#     elif traj_ext == '.xtc':
#         pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
#         topologyTraj = md.load_xtc(traj_filename, top=pdb_fn, frame=1)
#     else:
#         raise ValueError("Trajectory file format {} not recognized in the trajectory {}".format(traj_ext, traj_filename))
#
#     # Put the reference mol xyz into the 1-frame topologyTraj to use as a reference in the fit
#     molXyz = oechem.OEDoubleArray(3*mol.GetMaxAtomIdx())
#     mol.GetCoords(molXyz)
#     molXyzArr = np.array(molXyz)
#     molXyzArr.shape = (-1, 3)
#
#     # convert from angstroms to nanometers and slice out the protein-ligand complex
#     topologyTraj.xyz[0] = molXyzArr/10.0
#
#     # extract protein and ligand molecules from the larger multicomponent system
#     # and identify residue numbers for residues within fromLigCutoff of the ligand.
#     protein, ligand, actSiteResNums = ExtractProtLigActsiteResNums(mol, fromLigCutoff)
#     protResMap, numProtRes = GetCardinalOrderOfProteinResNums(protein)
#     actSiteResIdxs = set()
#     for resnum in actSiteResNums:
#         actSiteResIdxs.add(protResMap[resnum])
#
#     # Extract protein atom indices: cannot trust mdtraj protein selection so
#     # assume they are contiguous and starting the atom list and just get the same
#     # number of atoms as in the OpenEye protein
#     protOEIdx = np.array([atom.GetIdx() for atom in protein.GetAtoms()])
#
#     # Extract ligand atom indices
#     #   Note: the ligand must have residue name 'MOL' or 'LIG' (bad, should change)
#     ligIdx = topologyTraj.topology.select('resname == MOL or resname == LIG')
#     protligIdx = np.append(protOEIdx, ligIdx)
#
#     # print( 'numAtoms prot, lig, protlig:', len(protOEIdx), len(ligIdx), len(protligIdx))
#     #protligIdx = topologyTraj.topology.select('protein or resname == MOL or resname == LIG')
#
#     # Read the protein-ligand subsystem of the trajectory file
#     if traj_ext == '.h5':  # OpenMM
#         trj_initial = md.load_hdf5(traj_filename, atom_indices=protligIdx)
#     elif traj_ext == '.xtc':  # Gromacs
#         trj_initial = md.load_xtc(traj_filename, top=pdb_fn, atom_indices=protligIdx)
#
#     if skip > 0 and len(trj_initial) > skip:
#         trj = trj_initial[skip:]
#     else:
#         trj = trj_initial
#
#     # Image the protein-ligand trajectory so the complex does not jump across box boundaries
#     protlig = topologyTraj.atom_slice(protligIdx)
#     protligAtoms = [atom for atom in protlig.topology.atoms]
#     inplace = False
#     trjImaged = trj.image_molecules(inplace, [protligAtoms])
#
#     # Make a list of the atom indices of the carbon-alphas of the active site residues;
#     # assume residue numbering matches the mol
#     actSiteCA = [atom.index for atom in topologyTraj.topology.atoms
#                  if ((atom.residue.resSeq in actSiteResIdxs) and (atom.name == 'CA'))]
#
#     # Fit the protein-ligand trajectory to the active site carbon-alphas of the reference
#     trjImaged.superpose(protlig, 0, actSiteCA)
#
#     # Generate a multiconformer representation of the ligand trajectory
#     ligTraj = oechem.OEMol(ligand)
#
#     ligTraj.DeleteConfs()
#     for frame in trjImaged.xyz:
#         xyzList = [10*frame[idx] for idx in ligIdx]
#         confxyz = oechem.OEFloatArray( np.array(xyzList).ravel())
#         conf = ligTraj.NewConf(confxyz)
#
#     # Generate a multiconformer representation of the protein trajectory
#     strNumProteinAtomsToSelect = 'index ' + str(protOEIdx[0]) + ' to ' + str(protOEIdx[-1])
#     protIdx = protlig.topology.select(strNumProteinAtomsToSelect)
#     protTraj = oechem.OEMol(protein)
#     protTraj.DeleteConfs()
#
#     for frame in trjImaged.xyz:
#         xyzList = [10*frame[idx] for idx in protIdx]
#         confxyz = oechem.OEFloatArray(np.array(xyzList).ravel())
#         conf = protTraj.NewConf(confxyz)
#
#     return protTraj, ligTraj


def extract_aligned_prot_lig_wat_traj(setup_mol, well, trj_fn, nmax, opt, cutoff=5.0):
    """Extracts the aligned protein trajectory and aligned ligand trajectory and aligned
    Water trajectory from a MD trajectory of a larger system that includes other
    components (eg water).
    The passed in setup mol must have the topology that matches the trajectory, and its xyz
    coordinates are the reference for the alignment. The alignment is done on the
    alpha carbons (atom name CA) of the active site residues within cutoff
    from the ligand. Once the alignment is done, the protein and ligand trajectories
    are each placed into a separate OEMol, one conformer per trajectory frame.
    Water trajectory is selecting the nmax waters from the ligand and protein CA
    within the cutoff distance for each trajectory snapshot

    Inputs:
        mol: An OEMol giving the topology for the trajectory and the reference xyz
            coordinates for the alignment.
        well: An OEMol of the system well

        trj_fn: The filename of the hdf5-format MD trajectory or Gromacs .xtc file format
        cutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
        nmax: max number of waters to select
    Outputs:
        multi_conf_protein: A multi conformer OEMol for the protein, one conformer per frame.
        multi_conf_ligand: A multi conformer OEMol for the ligand, one conformer per frame.
        multi_conf_water: A multi conformer OEMol for the waters, one conformer per frame."""

    def dist2(coord1, coord2, unit_cell_lengths):

        dx = coord1[0] - coord2[0]
        dy = coord1[1] - coord2[1]
        dz = coord1[2] - coord2[2]

        dx = dx - unit_cell_lengths[0] * int(round(dx / (2 * unit_cell_lengths[0])))
        dy = dy - unit_cell_lengths[1] * int(round(dy / (2 * unit_cell_lengths[1])))
        dz = dz - unit_cell_lengths[2] * int(round(dz / (2 * unit_cell_lengths[2])))

        distsq = dx ** 2 + dy ** 2 + dz ** 2

        return distsq

    void, traj_ext = os.path.splitext(trj_fn)

    traj_dir = os.path.dirname(trj_fn)

    # Determine the protein and ligand idx
    # Extract Ligand idx and ca protein idx within cutoff from the ligand
    protein, ligand, water, excipients = oeommutils.split(well, ligand_res_name="LIG")

    nn = oechem.OENearestNbrs(well, cutoff)
    pred = oechem.OEIsAlphaCarbon()
    ca_set = set()
    for nbrs in nn.GetNbrs(ligand):
        at = nbrs.GetBgn()
        if pred(at):
            ca_set.add(at.GetIdx())

    ca_idx = sorted(list(ca_set))

    lig_idx = []
    for at in well.GetAtoms():
        if oechem.OEAtomGetResidue(at).GetName() == 'LIG':
            lig_idx.append(at.GetIdx())

    prot_idx = []
    for at in protein.GetAtoms():
        prot_idx.append(at.GetIdx())

    lig_ca_idx = lig_idx + ca_idx

    # Load Trajectory
    if traj_ext == '.h5':
        trj = md.load_hdf5(trj_fn)

    elif traj_ext == '.xtc':
        pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
        trj = md.load_xtc(trj_fn, top=pdb_fn)
        trj = trj[1:]
    else:
        raise ValueError("Trajectory file format {} not recognized in the trajectory {}".format(traj_ext, trj_fn))

    # Image the protein-ligand trajectory so the complex does not jump across box boundaries
    protlig = trj[0].atom_slice(prot_idx + lig_idx)
    protligAtoms = [atom for atom in protlig.topology.atoms]

    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stderr(devnull):
            trjImaged = trj.image_molecules(inplace=False, anchor_molecules=[protligAtoms], make_whole=True)

    if nmax != 0:

        count = 0
        water_max_frames = []

        for frame in trjImaged:
            # print(count)
            opt['Logger'].info("count {}".format(count))
            # Extract frame coordinates
            xyz = frame.xyz * 10

            # Set Well Coordinates as the current frame
            well.SetCoords(xyz[0].flatten())

            atm_list_frame_idx = md.compute_neighbors(frame, cutoff/10.0, lig_ca_idx, periodic=True)

            atm_idx = atm_list_frame_idx[0].tolist()

            oe_water_atoms = set()

            pred = oechem.OEIsWater(checkHydrogens=False)

            for idx in atm_idx:
                atom = well.GetAtom((oechem.OEHasAtomIdx(idx)))
                if len(list(oechem.OEGetResidueAtoms(atom))) == 3:
                    # Check if the residue is a water molecule
                    for at in oechem.OEGetResidueAtoms(atom):
                        if pred(at):
                            oe_water_atoms.add(at)

            # Here we are assuming that all the water molecules
            # have been  imaged inside the box centered in (0,0,0)
            # and sides (Lx,Ly,Lz)
            # water_coords = waters.GetCoords()
            unit_cell_lengths = frame.unitcell_lengths[0] * 10

            well_coords = well.GetCoords()

            metric_distances = dict()

            pred = oechem.OEIsOxygen()
            for at_water in oe_water_atoms:
                if pred(at_water):

                    at_water_coords = well_coords[at_water.GetIdx()]
                    at_water_lig_distances2 = []
                    for idx in lig_idx:
                        at_lig_coords = well_coords[idx]
                        at_water_at_lig_distance2 = dist2(at_water_coords, at_lig_coords, unit_cell_lengths)
                        at_water_lig_distances2.append(at_water_at_lig_distance2)

                    min_distance2_at_water_lig = min(at_water_lig_distances2)

                    at_water_ca_distances2 = []
                    for idx in ca_idx:
                        at_ca_coords = well_coords[idx]
                        at_water_at_ca_distance2 = dist2(at_water_coords, at_ca_coords, unit_cell_lengths)
                        at_water_ca_distances2.append(at_water_at_ca_distance2)

                    min_distance2_at_water_cas = min(at_water_ca_distances2)

                    metric = min_distance2_at_water_lig + min_distance2_at_water_cas

                    metric_distances[at_water] = metric

            water_list_sorted_max = sorted(metric_distances.items(), key=lambda x: x[1])[:nmax]

            water_max_frames.append(water_list_sorted_max)
            count += 1

    # Put the reference mol xyz into the 1-frame topologyTraj to use as a reference in the fit
    setup_mol_array_coords = oechem.OEDoubleArray(3 * setup_mol.GetMaxAtomIdx())
    setup_mol.GetCoords(setup_mol_array_coords)

    setup_mol_xyzArr = np.array(setup_mol_array_coords)
    setup_mol_xyzArr.shape = (-1, 3)

    trj_reference = trjImaged[0]
    # convert from angstroms to nanometers
    trj_reference.xyz[0] = setup_mol_xyzArr / 10.0

    # Fitting
    trjImaged.superpose(trj_reference, 0, ca_idx)

    ligand_reference = oechem.OEMol(ligand)
    protein_reference = oechem.OEMol(protein)

    count = 0
    multi_conf_water = None
    for frame in trjImaged.xyz:
        # print(count)
        if nmax != 0:
            xyz = frame * 10
            # Set Well Coordinates as the current frame for the water extraction
            well.SetCoords(xyz.flatten())
            water_list_sorted_max = water_max_frames[count]

            bv = oechem.OEBitVector(nmax * 3)

            water_idx = []

            for pair in water_list_sorted_max:
                ow = pair[0]
                water_idx.append(ow.GetIdx())
                for atw in oechem.OEGetResidueAtoms(ow):
                    bv.SetBitOn(atw.GetIdx())
                    water_idx.append(atw.GetIdx())

            pred_vec = oechem.OEAtomIdxSelected(bv)
            water_nmax_reference = oechem.OEMol()
            oechem.OESubsetMol(water_nmax_reference, well, pred_vec)

        lig_xyz_list = [10 * frame[idx] for idx in lig_idx]
        lig_confxyz = oechem.OEFloatArray(np.array(lig_xyz_list).ravel())

        prot_xyz_list = [10 * frame[idx] for idx in prot_idx]
        prot_confxyz = oechem.OEFloatArray(np.array(prot_xyz_list).ravel())

        if count == 0:

            if nmax != 0:

                multi_conf_water = oechem.OEMol(water_nmax_reference)

                if multi_conf_water.NumAtoms() % 3 != 0:
                    raise ValueError("Number of Water atoms is not multiple of 3")

                # Clean ResNumber and Chain
                oechem.OEPerceiveResidues(multi_conf_water, oechem.OEPreserveResInfo_All)
                multi_conf_water.SetTitle("Water_" + str(nmax))

                res_num = 0
                i = 0
                for at in multi_conf_water.GetAtoms():
                    res = oechem.OEAtomGetResidue(at)
                    res.SetName("HOH")
                    res.SetChainID("A")
                    if i % 3 == 0:
                        res_num += 1
                    res.SetResidueNumber(res_num)
                    i += 1

            ligand_reference.SetCoords(lig_confxyz)
            protein_reference.SetCoords(prot_confxyz)
            multi_conf_ligand = oechem.OEMol(ligand_reference)
            multi_conf_protein = oechem.OEMol(protein_reference)
        else:
            if nmax != 0:
                water_confxyz = oechem.OEFloatArray(water_nmax_reference.NumAtoms() * 3)
                water_nmax_reference.GetCoords(water_confxyz)
                multi_conf_water.NewConf(water_confxyz)

            multi_conf_ligand.NewConf(lig_confxyz)
            multi_conf_protein.NewConf(prot_confxyz)

        count += 1

    return multi_conf_protein, multi_conf_ligand, multi_conf_water


def RequestOEField(record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))


def RequestOEFieldType( record, field):
    if not record.has_value(field):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field.get_name() ))
        raise ValueError('The record does not have field {}'.format( field.get_name() ))
    else:
        return record.get_value(field)


def ColorblindRGBMarkerColors( nColors=0):
    palette = [(0, 114, 178), (0, 158, 115), (213, 94, 0), (204, 121, 167),
        (240, 228, 66), (230, 159, 0), (86, 180, 233), (150, 150, 150)]
    if nColors<1:
        return palette
    elif nColors<9:
        return palette[:nColors]
    else:
        n = int(nColors/8)
        moreRGB = palette
        for i in range(n):
            moreRGB = moreRGB+palette
        return( moreRGB[:nColors])


def PoseInteractionsSVG(ligand, proteinOrig, width=400, height=300):
    """Generate a OEGrapheme interaction plot for a protein-ligand complex.
    The input protein may have other non-protein components as well so
    the input protein is first split into components to isolate the protein
    only for the plot. This may have to be changed if other components need
    to be included in the plot.
    """
    # perceive residue hierarchy of total system
    if not oechem.OEHasResidues(proteinOrig):
        oechem.OEPerceiveResidues(proteinOrig, oechem.OEPreserveResInfo_All)
        print('Perceiving residues')
    # split the total system into components
    ligandPsuedo = oechem.OEMol()
    protein = oechem.OEMol()
    water = oechem.OEMol()
    other = oechem.OEMol()
    #sopts = oechem.OESplitMolComplexOptions('MOL')
    sopts = oechem.OESplitMolComplexOptions()
    oechem.OESplitMolComplex(ligandPsuedo, protein, water, other, proteinOrig, sopts)
    #
    # make the OEHintInteractionContainer
    asite = oechem.OEInteractionHintContainer(protein, ligand)
    if not asite.IsValid():
        oechem.OEThrow.Fatal("Cannot initialize active site!")
    # do the perceiving
    oechem.OEPerceiveInteractionHints(asite)
    # set the depiction options
    opts = oegrapheme.OE2DActiveSiteDisplayOptions(width, height)
    opts.SetRenderInteractiveLegend(True)
    magnifyresidue = 1.0
    opts.SetSVGMagnifyResidueInHover(magnifyresidue)
    # make the depiction
    oegrapheme.OEPrepareActiveSiteDepiction(asite)
    adisp = oegrapheme.OE2DActiveSiteDisplay(asite, opts)
    # make the image
    image = oedepict.OEImage(width, height)
    oegrapheme.OERenderActiveSite(image, adisp)
    # Add a legend
    #iconscale = 0.5
    #oedepict.OEAddInteractiveIcon(image, oedepict.OEIconLocation_TopRight, iconscale)
    svgBytes = oedepict.OEWriteImageToString("svg", image)

    svgString = svgBytes.decode("utf-8")

    return svgString


def ligand_to_svg_stmd(ligand, ligand_name):

    class ColorLigandAtomByBFactor(oegrapheme.OEAtomGlyphBase):
        def __init__(self, colorg):
            oegrapheme.OEAtomGlyphBase.__init__(self)
            self.colorg = colorg

        def RenderGlyph(self, disp, atom):
            adisp = disp.GetAtomDisplay(atom)
            if adisp is None or not adisp.IsVisible():
                return False

            if not oechem.OEHasResidue(atom):
                return False

            res = oechem.OEAtomGetResidue(atom)
            bfactor = res.GetBFactor()
            color = self.colorg.GetColorAt(bfactor)

            pen = oedepict.OEPen(color, color, oedepict.OEFill_On, 1.0)
            radius = disp.GetScale() / 3.0

            layer = disp.GetLayer(oedepict.OELayerPosition_Below)
            circlestyle = oegrapheme.OECircleStyle_Default
            oegrapheme.OEDrawCircle(layer, circlestyle, adisp.GetCoords(), radius, pen)
            return True

        def CreateCopy(self):
            return ColorLigandAtomByBFactor(self.colorg).__disown__()

    with TemporaryDirectory() as output_directory:

        lig_copy = oechem.OEMol(ligand)

        if len(ligand_name) < 15:
            lig_copy.SetTitle(ligand_name)
        else:
            lig_copy.SetTitle(ligand_name[0:13] + '...')

        img_fn = os.path.join(output_directory, "img.svg")

        oegrapheme.OEPrepareDepictionFrom3D(lig_copy)

        colorg = oechem.OELinearColorGradient()
        colorg.AddStop(oechem.OEColorStop(0.0, oechem.OEDarkBlue))
        colorg.AddStop(oechem.OEColorStop(10.0, oechem.OELightBlue))
        colorg.AddStop(oechem.OEColorStop(25.0, oechem.OEYellowTint))
        colorg.AddStop(oechem.OEColorStop(50.0, oechem.OERed))
        colorg.AddStop(oechem.OEColorStop(100.0, oechem.OEDarkRose))

        color_bfactor = ColorLigandAtomByBFactor(colorg)

        width, height = 150, 150
        opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
        opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)
        disp = oedepict.OE2DMolDisplay(lig_copy, opts)

        oegrapheme.OEAddGlyph(disp, color_bfactor, oechem.OEIsTrueAtom())

        oedepict.OERenderMolecule(img_fn, disp)

        svg_lines = ""
        marker = False
        with open(img_fn, 'r') as file:
            for line in file:
                if marker:
                    svg_lines += line

                if line.startswith("<svg"):
                    marker = True
                    svg_lines += line
                    svg_lines += """<title>{}</title>\n""".format(ligand_name)

                if line.startswith("</svg>"):
                    marker = False

    return svg_lines
