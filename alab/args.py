####### utils.py

def add_pyflow_args(parser):
	"""
	Given an argument parser from the argparse library,
	adds a set of common pyflow parameters.

	Parameters added:
		run_mode
		nCores
		memMb
		mailTo
		dry
		isContinue
		pyflow_dir
		schedulerArgList
	"""
	parser.add_argument('--pyflow_dir', type=str, default='.', required=False)
	parser.add_argument('--nCores', type=int, default=128)
	parser.add_argument('--memMb', type=int, default=4096)
	parser.add_argument('--run_mode', type=str, default='sge')
	parser.add_argument('--dry', default=False, action="store_true")
	parser.add_argument('--isContinue', default='Auto', action="store_true")
	parser.add_argument('--forceContinue', default=False, action="store_true")
	parser.add_argument('--mailTo', type=str, default=None)
	parser.add_argument('--startFromTasks', type=str, default=None)
	parser.add_argument('--ignoreTasksAfter', type=str, default=None)
	parser.add_argument('--resetTasks', type=str, default=None)
	parser.add_argument('--schedulerArgList', type=str, default=None)
	
	return parser

def default_pyflow_args(args):
	"""
	Returns a dictionary of arguments commonly sent to
	a WorkflowRunner.run() method.
	"""
	startFromTasks = None
	ignoreTasksAfter = None
	resetTasks = None
	schedulerArgList = None
	
	if args.startFromTasks:
		startFromTasks = args.startFromTasks.split(',')
	if args.ignoreTasksAfter:
		ignoreTasksAfter = args.ignoreTasksAfter.split(',')
	if args.resetTasks:
		resetTasks = args.resetTasks.split(',')
	if args.schedulerArgList :
		schedulerArgList = args.schedulerArgList.strip("[|]")
		schedulerArgList = schedulerArgList.split(',')
	
	arg_dict = {
			'mode':       args.run_mode,
			'nCores':     args.nCores,
			'memMb':      args.memMb,
			'isDryRun':   args.dry,
			'isContinue': args.isContinue,
			'isForceContinue': args.forceContinue,
			'dataDirRoot':args.pyflow_dir,
			'mailTo':     args.mailTo,
			'startFromTasks': startFromTasks,
			'ignoreTasksAfter': ignoreTasksAfter,
			'resetTasks': resetTasks,
			'schedulerArgList' : schedulerArgList
			}
			
	return arg_dict

def extend_pyflow_docstring(docstring):
	"""
	Updates a doc string to include the usage instructions for default pyflow parameters
	"""
	default_docstring = """
	Pyflow arguments:
	Optional arguments:
		--run_mode   Valid options: (Default) 'local', 'sge'
		--nCores     Number of threads to use (Default: 1)
		--memMb      memory available in local mode in MB
		--mailTo     Email address to send updates on pipeline running progress.
		--pyflow_dir     Directory to store pyflow tracking data.
		--dry       Set this flag to run in dry mode (does not launch tasks).
		--isContinue    Modify the continuation method of pyflow. Default: auto
		--forceContinue    If isContinue is set, allows task definitions to change.
		--startFromTasks    A comma-delimited list of task labels as strings. Any tasks which are not in this set or descendants of this set will be marked as completed.
		--ignoreTasksAfter  A comma-delimited list of task labels as strings. All descendants of these task labels will be ignored.
		--resetTasks    A comma-delimited list of task labels as strings. These tasks and all of their descendants will be reset to the "waiting" state to be re-run. Note this option will only effect a workflow which has been continued from a previous run. This will not overide any nodes alterned by the startFromTasks setting in the case that both options are used together
		--schedulerArgList A list of additional Parameters related to scheduler
	"""
	return docstring + default_docstring