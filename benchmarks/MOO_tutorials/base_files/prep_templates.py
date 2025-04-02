import shutil
import glob
import os
import pyemu

def prep_templates(tmpl_in, temp_d, nmax_inner):

    pst_files = glob.glob(os.path.join(tmpl_in, '*.pst'))
    pst = pyemu.Pst(pst_files[0])

    #prep outer iter pst files
    if os.path.exists(os.path.join(temp_d,'template_outer')):
        shutil.rmtree(os.path.join(temp_d, 'template_outer'))
    shutil.copytree(tmpl_in, os.path.join(temp_d, 'template_outer'))

    pst.control_data.noptmax = -1
    pst.write(os.path.join(temp_d, 'template_outer', os.path.basename(pst_files[0])))

    #prep outer repo update template
    pst = pyemu.Pst(pst_files[0])
    if os.path.exists(os.path.join(temp_d,'template_repo_update')):
        shutil.rmtree(os.path.join(temp_d,'template_repo_update'))
    shutil.copytree(tmpl_in, os.path.join(temp_d,'template_repo_update'))
    pst.pestpp_options['mou_dv_population_file'] = 'merged.dv_pop.csv'
    pst.pestpp_options['mou_obs_population_restart_file'] = 'merged.obs_pop.csv'
    pst.write(os.path.join(temp_d, 'template_repo_update', "outer_repo.pst"))

    #prep inner iter pst file
    pst = pyemu.Pst(pst_files[0])
    if os.path.exists(os.path.join(temp_d,'template_inner')):
        shutil.rmtree(os.path.join(temp_d,'template_inner'))
    shutil.copytree(tmpl_in, os.path.join(temp_d,'template_inner'))

    pst.control_data.noptmax = nmax_inner
    obs = pst.observation_data

