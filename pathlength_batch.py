from dipy.tracking.utils import density_map, get_flexi_tvis_affine, path_length
from nipype.utils.filemanip import load_json
import os
import nibabel as nib
import numpy as np


def safesavenii(array, aff, path):
    saveim = nib.Nifti1Image(array, aff)
    saveim.set_qform(aff, 1)
    saveim.set_sform(aff, 1)
    saveim.to_filename(path)

basedirtxt = '/home/kjordan/python_code/mydata_4myscripts/radonc_basedir.txt'
basedir = open(basedirtxt).readline().strip()
ptjson = '/home/kjordan/python_code/mydata_4myscripts/radonc_ptlist.json'

aoi_pathfrag = 'tp1/GTV_diffsp.nii.gz'
trk_pathfrag = 'tracking/GTV.trk'

ptlist = load_json(ptjson)

for pt in ptlist:
    ptpath = os.path.join(basedir, pt)
    aoipath = os.path.join(ptpath, aoi_pathfrag)
    trkpath = os.path.join(ptpath, trk_pathfrag)
    putpath = aoipath.replace('.nii.gz', '_sls'+os.path.basename(trkpath).split('.')[0]+'_pmap.nii.gz')
    dmap_putpath = trkpath.replace('.trk', '_dmap.nii.gz')
    print putpath

    print ptpath
    print aoipath
    print trkpath

    trk, hdr = nib.trackvis.read(trkpath)
    sls = [item[0] for item in trk]
    aoiim = nib.load(aoipath)
    aoidata = aoiim.get_data()
    aoiaff = aoiim.get_affine()

    print hdr['dim']
    print hdr['voxel_order']
    grid2trk_aff = get_flexi_tvis_affine(hdr, aoiaff)
    pmap = path_length(sls, aoidata, grid2trk_aff)
    dmap = density_map(sls, aoidata.shape, affine=grid2trk_aff)

    safesavenii(pmap, aoiaff, putpath)
    safesavenii(dmap, aoiaff, dmap_putpath)
