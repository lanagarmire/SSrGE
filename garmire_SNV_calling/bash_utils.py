from subprocess import call

from sys import stdout as STDOUT


def exec_cmd(cmd, stdout=STDOUT):
    """
    """
    if stdout is None:
        stdout = STDOUT

    try:
        answer = call(cmd.split(), stdout=stdout)
    except Exception:
        raise Exception('error when launching {0} \n cannot execute the command!'.format(cmd))

    try:
        assert(answer == 0)
    except Exception:
        raise Exception('{0} return a non 0 code!'.format(cmd))

    call('echo ### cmd: {0} succesfull ###\n'.format(cmd).split(), stdout=stdout)
