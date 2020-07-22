import sys, os
import ConfigParser 
import hashlib
import optparse
import subprocess
import pwd


def main():

    options = parse_options()
    run_SPServer_cluster(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("run_SPServer_cluster.py  -i INPUT_DIR -o OUTPUT_DIR [--dummy_dir=DUMMY_DIR] [--logs_dir=LOGS_DIR]")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Folder with all the input targets we want to analyze", metavar="INPUT_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Output folder where to store the results", metavar="OUTPUT_DIR")
    parser.add_option("--dummy_dir", default="dummy/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = ./)", metavar="DUMMY_DIR")
    parser.add_option("--logs_dir", default="logs/", action="store", type="string", dest="logs_dir", help="Logs directory (default = ./)", metavar="LOGS_DIR")


    (options, args) = parser.parse_args()

    if options.input_dir is None or options.output_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options




#################
#################
# MAIN FUNCTION #
#################
#################

def run_SPServer_cluster(options):
    """
    Runs SPSerer over a folder of targets with their corresponding models in the cluster
    """
    # Add "." to sys.path #
    src_path =  os.path.abspath(os.path.dirname(__file__))
    sys.path.append(src_path)

    # Read configuration file #     
    config = ConfigParser.ConfigParser()
    config_file = os.path.join(src_path, "config.ini")
    config.read(config_file)

    # Limit of jobs
    limit = config.get("Cluster", "max_jobs_in_queue")
    l = 1

    # Arguments & Options #
    options = parse_options()
    input_dir = options.input_dir
    output_dir = options.output_dir
    create_directory(output_dir)

    # Create the dummy directory (to store the commands)
    dummy_dir = os.path.abspath(options.dummy_dir)
    create_directory(dummy_dir)

    # Create log directory (to store the log files of the jobs)
    logs_dir = os.path.abspath(options.logs_dir)
    create_directory(logs_dir)

    # Define queue file
    #queue_file = os.path.join(src_path, 'command_queues_hydra.txt')
    # Define queue parameters
    max_mem = config.get("Cluster", "max_mem")
    queue = config.get("Cluster", "cluster_queue") # short, normal, long, bigmem
    modules = ['Python/2.7.11']
    queue_parameters = {'max_mem':max_mem, 'queue':queue, 'logs_dir':logs_dir, 'modules':modules}

    # Gather the target folders
    target_dirs = [f for f in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, f))]
    print(target_dirs)
    print(len(target_dirs))

    # Get the parameters
    pot_type = config.get("Parameters", "pot_type")
    spserverfold_path = os.path.join(config.get("Parameters", "spserverver_path"), 'SPServerFold.py')

    # Gather the results folders for each target
    for target in target_dirs:
        target_dir = os.path.join(input_dir, target)
        target_output_dir = os.path.join(output_dir, target)
        create_directory(target_output_dir)
        # Get the model files
        models = [model for model in os.listdir(target_dir) if fileExist(os.path.join(target_dir, model))]

        for model in models:

            if limit: # Break the loop if a limit of jobs is introduced
                if l > limit:
                    print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                    break

            # Output dir for the specific file
            file_output_dir = os.path.join(target_output_dir, model)
            create_directory(file_output_dir)
            output_xml = os.path.join(file_output_dir, '{}.xml'.format(model))
            if not fileExist(output_xml):

                # Define the command
                command = 'python {} -f {} -p {} -o {} -x {}'.format( spserverfold_path, os.path.join(target_dir, model), pot_type, file_output_dir, model )
                print(command)

                #To run in the cluster submitting files to queues
                script_name = '{}_{}.sh'.format(target, model)
                submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=script_name)

                l += 1

    return




#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return

#-------------#
# Cluster     #
#-------------#

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, queue_parameters={'max_mem':5000, 'queue':'short', 'logs_dir':'/tmp', 'modules':['Python/2.7.11']}, dummy_dir="/tmp", script_name=None):
    """
    This function submits any {command} to a cluster {queue}.

    @input:
    command {string}
    queue {string} by default it submits to any queue
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue

    """
    if max_jobs_in_queue is not None:
        while number_of_jobs_in_queue() >= max_jobs_in_queue: time.sleep(5)

    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    if not script_name:
        script_name = "submit_"+hashlib.sha224(command).hexdigest()+".sh"
    script= os.path.join(dummy_dir, script_name)
    if queue_file is not None:
        # Use a queue file as header of the bash file
        fd=open(script,"w")
        with open(queue_file,"r") as queue_standard:
            data=queue_standard.read()
            fd.write(data)
            fd.write("%s\n\n"%(command))
        fd.close()
        queue_standard.close()
        if queue is not None:
            os.system("sbatch -p %s %s" % (queue,script))
        else:
            os.system("sbatch %s"% (script))
    else:
        if queue_parameters is not None:
            # Write manually the bash file
            with open(script, "w") as fd:
                fd.write('#!/bin/bash\n') # bash file header
                fd.write('#SBATCH --mem={}\n'.format(queue_parameters['max_mem'])) # max memory
                fd.write('#SBATCH -p {}\n'.format(queue_parameters['queue'])) # queue name
                fd.write('#SBATCH -o {}.out\n'.format(os.path.join(queue_parameters['logs_dir'], script_name))) # standard output
                fd.write('#SBATCH -e {}.err\n'.format(os.path.join(queue_parameters['logs_dir'], script_name))) # standard error
                for module in queue_parameters['modules']:
                    fd.write('module load {}\n'.format(module)) # modules to load
                fd.write('{}\n'.format(command)) # command
            os.system("sbatch {}".format(script)) # execute bash file
        else:
            # Execute directly the command without bash file
            if queue is not None:
                os.system("echo \"%s\" | sbatch -p %s" % (command, queue))
            else:
                os.system("echo \"%s\" | sbatch" % command)


def number_of_jobs_in_queue():
    """
    This functions returns the number of jobs in queue for a given
    user.

    """

    # Initialize #
    user_name = get_username()

    process = subprocess.check_output(["squeue", "-u", user_name])

    return len([line for line in process.split("\n") if user_name in line])


def get_username():
    """
    This functions returns the user name.

    """

    return pwd.getpwuid(os.getuid())[0]


if  __name__ == "__main__":
    main()

