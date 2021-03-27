import subprocess
import sys
def giggle_sort(path):
    temp_path = path + '_sorted'
    # make sorted directory
    cmd_str = 'mkdir -p ' + temp_path
    proc = subprocess.check_call(cmd_str,
                                    shell=True,
                                    timeout=100)
    # sort files using giggle sort
    cmd_str = 'giggle/scripts/sort_bed \"{}/*.bed*\" {} {}'.format(path, temp_path, 12)
    print(cmd_str)
    proc = subprocess.check_output(cmd_str,
                                    shell=True,
                                    timeout=2*60)
    proc = subprocess.check_call("rm -R -f {}".format(path), shell=True)
giggle_sort(sys.argv[1])