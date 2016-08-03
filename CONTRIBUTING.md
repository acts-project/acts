# Contributing to ACTS

Contributions to the ACTS project are very welcome and feedback on the documentation is greatly appreciated. In order to be able to contribute to the ACTS project, developers must have a valid CERN user account. Unfortunately, lightweight CERN accounts for external users do not have sufficient permissions to access certain CERN services used by the ACTS project.

## Mailing lists

1. [acts-users@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-users): Users of the ACTS project should subscribe to this list as it provides:
    - regular updates on the software,
    - access to the ACTS JIRA project for bug fixes/feature requests,
    - a common place for asking any kind of questions.
1. [acts-developers@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-developers): Developers are encouraged to also subscribe to this list as it provides you with:
    - a developer role in the ACTS JIRA project (allows you to handle tickets),
    - information about developer meetings,
    - a common place for technical discussions.

## Bug reports and feature requests

If you want to report or start a feature request, please open a ticket in the [ACTS JIRA](https://its.cern.ch/jira/projects/ACTS/) (**Note:** access is restricted to members of the mailing lists mentioned above). A comprehensive explanation will help the development team to respond in a timely manner. Therefore, the following details should be mentioned:

- bug reports
    - issue type: "Bug"
    - summary: short description of the problem
    - priority: will be set by the development team
    - components: if known, part of ACTS affected by this bug; leave empty otherwise
    - affects version: version of ACTS affected by this bug
    - a detailed description of the bug including a receipe on how to reproduce it and any hints which may help diagnosing the problem 
- feature requests
    - issue type: "Improvement" or "New Feature"
    - summary: short description of new feature
    - priority: will be set by the development team
    - a detailed description of the feature request including possible use cases and benefits for other users

## Make a contribution

The instructions below should help you getting started with development process in the ACTS project. If you have any questions, feel free to ask [acts-developers@cern](mailto:acts-developers@cern.ch) for help or guidance.

### Setting up your fork

The ACTS project uses a <tt>git</tt> repository which is hosted on the CERN GitLab server. In order to be able to create merge requests (and thus, contribute to the development of ACTS), you need to create a fork on the CERN GitLab server. A general introduction to the GitLab web interface can be found [here](https://gitlab.cern.ch/help/gitlab-basics/README.md). Very nice tutorials as well as explanations for concepts and workflows with <tt>git</tt> can be found on [Atlassian](https://www.atlassian.com/git/). For a shorter introduction and the full <tt>git</tt> documentation have a look [here](https://git-scm.com/docs/gittutorial).

##### Creating your fork

As a first step, you need to create your own fork of the ACTS project. For doing this, please go to the [ACTS GitLab page](https://gitlab.cern.ch/acts/a-common-tracking-sw), click on the fork button, and follow the instructions ([GitLab Help "How to fork a project"](https://gitlab.cern.ch/help/gitlab-basics/fork-project.md)).

##### Configuring your fork

**Important:** Due to some limitations in the GitLab JIRA integration, you need to fix the JIRA settings for your forked project. This can be done by starting from the GitLab project page of your fork and then going to "Settings -> Services -> JIRA". Remove anything in the field "Username" and leave it empty. If you fail to do so, the ACTS JIRA project will be spammed with hundreds of duplicates comments and you will likely receive an angry email from the development team ;-).

Once you have created your fork on the CERN GitLab server, you need to create a local copy to start coding. This is done by cloning your forked repository onto your local machine through the following command (if you are on a UNIX system):

`git clone <FORK_URL> <Destination>`

- &lt;FORK_URL&gt; is the web address of your forked repository which can be found on the project page in the GitLab web interface
- &lt;DESTINATION&gt; is optional and gives the location on your local machine where the clone will be created

You probably want to be able to pull in changes from the official ACTS repository to benefit from the latest and greatest improvements. This requires that you add the official ACTS repository as another remote to your local clone. 

`cd <DESTINATION>`<br />
`git remote add ACTS ssh://git@gitlab.cern.ch:7999/acts/a-common-tracking-sw.git`

You can check that everything went ok with

`git remote -v`

where the reference to the ACTS repository should appear (along with your forked repository on the CERN GitLab server). This procedure is also described [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/).

##### Keeping your fork up-to-date

At certain points you may want to sync your fork with the latest updates from the official ACTS repository. The following commands illustrate how to update the 'master' branch of fork. The same procedure can be used to sync any other branch, but you will rarely need this. Please mak sure to commit/stash all changes before proceeding to avoid any loss of data. The following commands must be run in the working directory of the local clone or your forked repository. 

`git fetch ACTS`<br />
`git checkout master`<br />
`git merge --ff-only ACTS/master`<br />
`git push origin master`

#### Workflow recommendations

#### Coding style and guidelines

### Submitting a merge request