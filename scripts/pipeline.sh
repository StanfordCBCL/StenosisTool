# setup all models
simvascular --python -- scripts/dev_scripts/sv_dev_segmentation_gen.py
python scripts/dev_scripts/dev_configure_base_solver.py
python scripts/jc_model_construction.py -solver_dirs $BASE_0080_0001 $BASE_0082_0001 $BASE_0086_0001 $BASE_0118_1000 $BASE_SU0238 -o jc_solver_dir_0 -f

# optimize Healthy
python scripts/optimization_driver.py -solver_dirs $BASE_0080_0001 $BASE_0082_0001 $BASE_0086_0001 $JC0_0080_0001 $JC0_0082_0001 $JC0_0086_0001 -f -int
# add RCRS for stenosis
python scripts/add_solver_rcrs.py -solver_dirs $BASE_0118_1000 $BASE_SU0238 $JC0_0118_1000 $JC0_SU0238

# construct artificial stenosis
python scripts/stenosis_tool_scripts/artificial_stenosis_driver.py -solver_dirs $BASE_0080_0001 $BASE_0082_0001 $BASE_0086_0001 $JC0_0080_0001 $JC0_0082_0001 $JC0_0086_0001 -f 
# construct treated stenosis
python scripts/stenosis_tool_scripts/stenosis_detection_driver.py -solver_dirs $BASE_0118_1000 $BASE_SU0238 $JC0_0118_1000 $JC0_SU0238 -f

# get 3D simulation files
python scripts/three_d/correct_healthy_3D.py -solver_dirs $BASE_0080_0001 $BASE_0082_0001 $BASE_0086_0001 $JC0_0080_0001 $JC0_0082_0001 $JC0_0086_0001 -f

# run all
python scripts/run_zerod_sim.py -solver_dirs $BASE_0080_0001 $BASE_0082_0001 $BASE_0086_0001 $JC0_0080_0001 $JC0_0082_0001 $JC0_0086_0001 $BASE_0118_1000 $BASE_SU0238 $JC0_0118_1000 $JC0_SU0238 -b -s -f -r

simvascular --python -- scripts/viz_script/sv_zerod_to_geom.py -solver_dirs $BASE_0080_0001 $BASE_0082_0001 $BASE_0086_0001 $JC0_0080_0001 $JC0_0082_0001 $JC0_0086_0001 $BASE_0118_1000 $BASE_SU0238 $JC0_0118_1000 $JC0_SU0238 -s -r

python scripts
