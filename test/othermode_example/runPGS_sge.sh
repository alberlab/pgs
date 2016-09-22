python $PGS_SOURCE_DIR/pgs.py --input_config $WORKING_DIR/input_config.json --run_mode sge --nCores 300 --memMb 800000 --pyflow_dir $WORKING_DIR --schedulerArgList ['-q','YOUR_QUEUE_NAME','-l','time=77:00:00']

