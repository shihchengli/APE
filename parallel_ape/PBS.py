# -*- coding: utf-8 -*-

"""
A module for giving commands in PBS job scheduler
"""

import os
import subprocess

def send_command_to_server(command, remote_path=''):
    """
    Send commands to server. `command` is either a sting or an array of string commands to send.
    If remote_path is not an empty string, the command will be executed in the directory path it points to.
    Returns lists of stdout, stderr corresponding to the commands sent.
    """
    wd = os.getcwd()

    if isinstance(command, list):
        command = ' '.join(command)

    if remote_path != '':
        os.chdir(remote_path)
    
    proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = [], []
    os.chdir(wd)

    if proc.stdout is not None: stdout = [line.decode().rstrip() for line in proc.stdout]
    if proc.stderr is not None: stderr = [line.decode().rstrip() for line in proc.stderr]
    return stdout, stderr

def check_job_status(job_id):
    """
    A helper function for checking job status.
    Possible status: `done`, `running`, `queue`, `errored`
    """
    username = get_username()
    cmd = 'qstat -u {un}'.format(un=username)
    stdout, stderr = send_command_to_server(cmd)
    for status_line in stdout:
        if str(job_id) in status_line:
            break
    else:
        return 'done'
    status = status_line.split()[9]
    if status.lower() in ['r']:
        return 'running'
    elif status.lower() in ['q']:
        return 'queue'
    else:
        return 'errored'

def get_username():
    proc = subprocess.Popen('whoami'.split(), stdout=subprocess.PIPE)
    username = (proc.communicate())[0].decode().strip()
    return username

def delete_job(job_id):
    """
    Deletes a running job
    """
    cmd = 'qdel' + ' ' + str(job_id)
    send_command_to_server(cmd)

def check_running_jobs_ids():
    """
    Return a list of ``int`` representing job IDs of all jobs submitted by the user on a server
    """
    running_jobs_ids = list()
    username = get_username()
    cmd = 'qstat -u {un}'.format(un=username)
    stdout, stderr = send_command_to_server(cmd)
    for i, status_line in enumerate(stdout):
        if i > 4:
            job_id = int(status_line.strip().split('.')[0])
            running_jobs_ids.append(job_id)
    return running_jobs_ids

def submit_job(submit_filename, remote_path=''):
    """
    Submit a job
    """
    job_status = ''
    job_id = 0
    cmd = 'qsub' + ' ' + submit_filename
    stdout, stderr = send_command_to_server(cmd, remote_path)
    if len(stderr) > 0 or len(stdout) == 0:
        job_status = 'errored'
        for line in stderr:
            if 'Requested node configuration is not available' in line:
                print('User may be requesting more resources than are available. \
                Please check server settings, such as cpus and memory')
            if 'No such file or directory' in line:
                print('User must check {filename} file is in the path of {path}'.format(filename=submit_filename,path=os.getcwd()))
    elif len(stdout[0]) > 0:
        job_status = 'queue'
        job_id = int(stdout[0].split('.')[0])
    return job_status, job_id