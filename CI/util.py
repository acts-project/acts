from fs.sshfs import SSHFS
import os
import contextlib
import sys
import gitlab as gitlab_api
import smtplib, ssl

try:
    from halo import Halo
except:
    Halo = None


def get_lxplus_fs(args):
    assert args.user is not None, "Need username for lxplus"
    assert args.pwd is not None, "Neew password for lxplus"
    with Spinner(text="Connecting to lxplus", persist=False):
        return SSHFS(
            host="lxplus.cern.ch",
            user=args.user,
            passwd=args.pwd,
            allow_agent=False,
            look_for_keys=False,
        )


def def_arguments(p, acc=False, gl=False):
    if acc:
        p.add_argument("--user", default="atsjenkins")
        p.add_argument("--pwd", default=os.getenv("ATSJENKINS_PASSWORD", None))
    if gl:
        p.add_argument(
            "--gitlab-access-token",
            "-t",
            default=os.getenv("ATSJENKINS_ACCESS_TOKEN", None),
            help="GitLab access token to authenticate with the GitLab API",
        )
    return p


def gitlab(args):
    assert args.gitlab_access_token is not None, "Need gitlab access token"
    gl = gitlab_api.Gitlab(
        "https://gitlab.cern.ch", private_token=args.gitlab_access_token
    )
    gl.auth()
    return gl


@contextlib.contextmanager
def smtp(args):
    server = "smtp.cern.ch"
    port = 587  # For SSL
    username = "atsjenkins@cern.ch"
    try:
        context = ssl.create_default_context()
        server = smtplib.SMTP(server, port)
        server.ehlo()  # Can be omitted
        server.starttls(context=context)  # Secure the connection
        server.ehlo()  # Can be omitted
        server.login(args.user, args.pwd)

        yield server
    finally:
        server.quit()


@contextlib.contextmanager
def Spinner(text, persist=True, *args, **kwargs):
    stream = kwargs.get("stream", sys.stdout)
    if stream.isatty() and Halo is not None:
        spinner = Halo(text, *args, **kwargs)
        spinner.start()
        try:
            yield
            if persist:
                spinner.succeed()
        except:
            if persist:
                spinner.fail()
            raise
        finally:
            if not persist:
                spinner.stop()
    else:
        sys.stdout.write(text + "\n")
        yield
