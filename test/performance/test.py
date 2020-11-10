import os
import sys
import subprocess
import re
import time
import colorama
import argparse
from colorama import Fore, Back, Style

ME = os.path.realpath(os.path.dirname(__file__))
ROOTDIR = os.path.realpath(os.path.dirname(os.path.dirname(ME)))

class SingleTest:
    def __init__(self, rootdir):
        self.rootdir = rootdir
        self.git_checkout('master')

    def test_run (self, command, print_output=False):
        p = subprocess.Popen(command, cwd=self.rootdir,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        ret = 0
        while True:
            output = p.stdout.readline()
            if output.decode() == '' and p.poll() is not None:
                break
            if re.search('Execution time is', str(output)):
                et_str = str(output)
                ret = float(et_str.split()[3])
            p.poll()
        p.wait()

        return(ret)


    def cmd_run (self, command, print_output=False):
        p = subprocess.Popen(command, cwd=self.rootdir, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        if print_output:
            while True:
                output = p.stdout.readline()
                if output.decode() == '' and p.poll() is not None:
                    break
                if output:
                    print(output)
            p.poll()
        p.wait()

        return(p.returncode)


    def git_commits(self):
        p = subprocess.Popen('git log --pretty=format:"%H"', cwd=self.rootdir,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        commits = []

        while True:
            output = p.stdout.readline()
            if output.decode() == '' and p.poll() is not None:
                break
            if output:
                commits.append(output.decode().strip())
        p.poll()
        return(commits)


    def git_checkout(self, commit):
        return(self.cmd_run('git checkout ' + commit))


    def __del__(self):
        print("Restoring state to clear master...")
        self.git_checkout('master')
        self.cmd_run('make distclean')


def avg(values):
    count = 0
    sum = 0
    for i in values:
        count += 1
        sum += i

    return(sum/count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', type=int,
                        help='Amount of last commits (default 5)',
                        default=5)

    args = parser.parse_args()

    single_test = SingleTest(ROOTDIR)
    commits = single_test.git_commits()
    # prepare to testing
    sys.stdout.write("Preparing to test .")
    sys.stdout.flush()

    single_test.git_checkout('master')
    sys.stdout.write(".")
    sys.stdout.flush()

    single_test.cmd_run('./autogen.sh')
    sys.stdout.write(".")
    sys.stdout.flush()

    single_test.cmd_run('./configure')
    sys.stdout.write(". DONE\n")

    start_time = 0
    end_time = 0
    count = 1
    cmd_run_times = []
    test_run_times = []

    for c in commits[:args.l]:
        if count == 1:
            sys.stdout.write("Running 1st test for commit " + c + " ...")
        elif count == 2:
            sys.stdout.write("Running 2nd test for commit " + c + " ...")
        elif count == 3:
            sys.stdout.write("Running 3rd test for commit " + c + " ...")
        else:
            sys.stdout.write("Running " + str(count) + "th test for commit " + c + " ...")
        sys.stdout.flush()

        single_test.git_checkout(c)
        start_time = time.time()
        test_run_times.append(single_test.test_run('make test'))
        end_time = time.time()

        cmd_run_times.append(round(end_time - start_time, 2))

        sys.stdout.write(" DONE with times [command, app-run-only]: [" + str(cmd_run_times[-1]) + ", " + str(test_run_times[-1]) + "]\n")
        count=count+1
    single_test.git_checkout('master')

    avg_cmd = avg(cmd_run_times)
    avg_test = avg(test_run_times)

    print(Fore.YELLOW + "+--------+" + Style.RESET_ALL)
    print(Fore.YELLOW + "| REPORT |" + Style.RESET_ALL)
    print(Fore.YELLOW + "+--------+" + Style.RESET_ALL)
    if cmd_run_times[0] <= avg_cmd or cmd_run_times[0] / avg_cmd < 1.01:
        print(Fore.BLUE + "Command execution PASSED with " + str(cmd_run_times[0]) + ". Average is: " + str(round(avg_cmd, 2)) + Style.RESET_ALL)
    elif cmd_run_times[0] / avg_cmd < 1.1:
        print(Fore.YELLOW + "Command execution PASSED with warning " + str(cmd_run_times[0]) + ". Average is: " + str(round(avg_cmd, 2)) + Style.RESET_ALL)
    else:
        print(Fore.RED + "Command execution FAILED with " + str(cmd_run_times[0]) + ". Average is: " + str(round(avg_cmd, 2)) + Style.RESET_ALL)

    if test_run_times[0] <= avg_test or test_run_times[0] / avg_test < 1.01:
        print(Fore.BLUE + "App execution PASSED with " + str(test_run_times[0]) + ". Average is: " + str(round(avg_test, 2)) + Style.RESET_ALL)
    elif test_run_times[0] / avg_test < 1.1:
        print(Fore.YELLOW + "App execution PASSED with warning " + str(test_run_times[0]) + ". Average is: " + str(round(avg_test, 2)) + Style.RESET_ALL)
    else:
        print(Fore.RED + "App execution FAILED with " + str(test_run_times[0]) + ". Average is: " + str(round(avg_test, 2)) + Style.RESET_ALL)
