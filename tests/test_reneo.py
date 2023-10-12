import subprocess


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_cli():
    exec_command("reneo -v")
    exec_command("reneo -h")
    exec_command("reneo run -h")
    exec_command("reneo simulate -h")
    exec_command("reneo test -h")
    exec_command("reneo install -h")
    exec_command("reneo config -h")
    exec_command("reneo citation")


def test_reneo_pipelines():
    exec_command("reneo simulate")
    exec_command("reneo config")
