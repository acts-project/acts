# Contributing to ACTS

Contributions to the ACTS project are very welcome and feedback on the documentation is greatly appreciated. In order to be able to contribute to the ACTS project, developers must have a valid CERN user account. Unfortunately, lightweight CERN accounts for external users do not have sufficient permissions to access certain CERN services used by the ACTS project.

1. [Mailing lists](#mailing-lists)
1. [Bug reports and feature requests](#bug-reports-and-feature-requests)
1. [Make a contribution](#make-a-contribution)
    1. [Setting up your fork](#setting-up-your-fork)
    1. [Workflow recommendations](#workflow-recommendations)
    1. [Coding style and guidelines](#coding-style-and-guidelines)
    1. [Creating a merge request](#creating-a-merge-request)
    1. [git tips and tricks](#git-tips-and-tricks)

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

The ACTS project uses a git repository which is hosted on the CERN GitLab server. In order to be able to create merge requests (and thus, contribute to the development of ACTS), you need to create a fork on the CERN GitLab server. A general introduction to the GitLab web interface can be found [here](https://gitlab.cern.ch/help/gitlab-basics/README.md). Very nice tutorials as well as explanations for concepts and workflows with git can be found on [Atlassian](https://www.atlassian.com/git/). For a shorter introduction and the full git documentation have a look at the [git tutorial](https://git-scm.com/docs/gittutorial).

##### Configuring git

Commits to repositories on the CERN GitLab server are only accepted from CERN users. Therefore, it is important that git is correctly configured on all machines you are working on. It is necessary to that the git user email address to the primary email address of your CERN account (usually: firstname.lastname@cern.ch). You can check the current values with:

    git config user.name
    git config user.email

You can change those settings by either editing the `.gitconfig` file in your home directory or by running:

    git config --global user.name "Donald Duck"
    git config --global user.email "donald.duck@cern.ch"

##### Creating your fork

As a first step, you need to create your own fork of the ACTS project. For doing this, please go to the [ACTS GitLab page](https://gitlab.cern.ch/acts/a-common-tracking-sw), click on the fork button, and follow the instructions ([GitLab Help "How to fork a project"](https://gitlab.cern.ch/help/gitlab-basics/fork-project.md)).

##### Configuring your fork

**Important:** Due to some limitations in the GitLab JIRA integration, you need to fix the JIRA settings for your forked project. This can be done by starting from the GitLab project page of your fork and then going to "Settings -> Services -> JIRA". Remove anything in the field "Username" and leave it empty. If you fail to do so, the ACTS JIRA project will be spammed with hundreds of duplicates comments and you will likely receive an angry email from the development team ;-).

Once you have created your fork on the CERN GitLab server, you need to create a local copy to start coding. This is done by cloning your forked repository onto your local machine through the following command (if you are on a UNIX system):

    git clone <FORK_URL> <Destination>

- &lt;FORK_URL&gt; is the web address of your forked repository which can be found on the project page in the GitLab web interface
- &lt;DESTINATION&gt; is optional and gives the location on your local machine where the clone will be created

You probably want to be able to pull in changes from the official ACTS repository to benefit from the latest and greatest improvements. This requires that you add the official ACTS repository as another remote to your local clone. 

    cd <DESTINATION>
    git remote add ACTS ssh://git@gitlab.cern.ch:7999/acts/a-common-tracking-sw.git

You can check that everything went ok with

    git remote -v

where the reference to the ACTS repository should appear (along with your forked repository on the CERN GitLab server). This procedure is also described on [github](https://help.github.com/articles/configuring-a-remote-for-a-fork/).

##### Keeping your fork up-to-date

At certain points you may want to sync your fork with the latest updates from the official ACTS repository. The following commands illustrate how to update the 'master' branch of fork. The same procedure can be used to sync any other branch, but you will rarely need this. Please mak sure to commit/stash all changes before proceeding to avoid any loss of data. The following commands must be run in the working directory of the local clone or your forked repository. 

    git fetch ACTS
    git checkout master
    git merge --ff-only ACTS/master
    git push origin master

#### Workflow recommendations

In the following a few recommendations are outlined which should help you to get familiar with development process in the ACTS project.

1. **Each development its own branch!**<br />
Branching in git is simple, it is fun and it helps you keep your working copy clean. Therefore, you should start a new branch for every development. All work which is logically/conceptually linked should happen in one branch. Keep your branches short. This helps immensly to understand the git history if you need to look at it in the future.<br />
If projects are complex (e.g. large code refactoring or complex new features), you may want to use _sub_-branches from the main development branch as illustrated in the picture below.<br />
<br />
<img src="doc/figures/sub_dev.png" alt="workflow for large feature">
1. **Never, ever directly work on any "official" branch!**<br />
Though not strictly necessary and in the end it is up to you, it is strongly recommended that you never commit directly on a branch which tracks an "official" branch. As all branches are equal in git, the definition of "official" branch is quite subjective. In the ACTS project you should not work directly on branches which are **protected** in the CERN GitLab repository. Usually, these are the _master_ and _release-X.Y.Z_ branches. The benefit of this strategy is that you will never have problems to update your fork. Any git merge in your local repository on such an "official" branch will always be a fast-forward merge.<br />
<br />
1. **Use atomic commits!**<br />
Similarly to the concept of branches, each commit should reflect a self-contained change. Try to avoid overly large commits (bad examples are for instance mixing logical change with code cleanup and typo fixes).<br />
<br />
1. **Write good commit messages!**<br />
Well-written commit messages are key to understand your changes. There are many guidelines available on how to write proper commit logs (e.g. [here](http://alistapart.com/article/the-art-of-the-commit), [here](http://chris.beams.io/posts/git-commit/), or [here](https://wiki.openstack.org/wiki/GitCommitMessages#Information_in_commit_messages)). As a short summary:
    - Structure your commit messages into short title (max 50 characters) and longer description (max width 72 characters)!<br />
      This is best achieved by avoiding the `commit -m` option. Instead write the commit message in an editor/git tool/IDE... 
    - Describe why you did the change (git diff already tells you what has changed)!
    - Mention any side effects/implications/consquences!<br /><br />

1. **Prefer git pull --rebase!**<br />
If you work with a colleague on a new development, you may want to include his latest changes. This is usually done by calling `git pull` which will synchronise your local working copy with the remote repository (which may have been updated by your colleague). By default, this action creates a merge commit if you have local commits which were not yet published to the remote repository. These merge commits are considered to contribute little information to the development process of the feature and they clutter the history (read more e.g.  [here](https://developer.atlassian.com/blog/2016/04/stop-foxtrots-now/) or [here](http://victorlin.me/posts/2013/09/30/keep-a-readable-git-history)). This problem can be avoided by using `git pull --rebase` which replays your local (un-pushed) commits on the tip of the remote branch. You can make this the default behaviour by running `git config pull.rebase true`. More about merging vs rebasing can be found [here](https://www.atlassian.com/git/tutorials/merging-vs-rebasing/).<br />
1. **Push your development branches as late as possible!**<br />
Unless required by other circumstances (e.g. collaboration with colleagues, code reviews etc) it is recommended to push your development branch once you are finished. This gives you more flexibility on what you can do with your local commits (e.g. rebase interactively) without affecting others. Thus, it minimises the risk for running into git rebase problems.

#### Coding style and guidelines

### Creating a merge request

### git tips and tricks

The following section gives some advise on how to solve certain situations you may encounter during your development process. Many of these commands have the potential to loose uncommitted data. So please make sure that you understand what you are doing before running the receipes below. Also, this collection is non-exhaustive and alternative approaches exist. If you want to contribute to this list, please drop an email to [acts-developers@cern.ch](mailto:acts-developers@cern.ch).

**Before doing anything**<br />
In the rare event that you end up in a situation you do not know how to solve, get to a clean state of working copy and create a (backup) branch, then switch back to the original branch. If anything goes wrong, you can always checkout the backup branch and you are back to where you started.<br /><br /> 
**Modify the author of a commit**<br />
If your git client is not correctly set up on the machine you are working on, it may derive the committer name and email address from login and hostname information. In this case your commits are likely rejected by the CERN GitLab server. As a first step, you should correctly configure git on this machine as described above so that this problems does not appear again.
-  You are lucky and only need to fix the author of the latest commit. You can use `git commit --amend`:<br />

    git commit --amend --no-edit --author "My Name <login@cern.ch>
    
- You need to fix (several) commit(s) which are not the current head. You can use `git rebase`:<br />
For the following it is assumed that all commits which need to be fixed are in the same branch &lt;BRANCH&gt;, and &lt;SHA&gt; is the hash of the earliest commit which needs to be corrected.<br />

    git checkout <BRANCH>
    git rebase -i -p <SHA>^
    
In the editor opened by the git rebase procedure, add the following line after each commit you want to fix:<br />

    exec git commit --amend --author="New Author Name <email@address.com>" -C HEAD
    
Then continue with the usual rebase procedure.

**Make a bugfix while working on a feature**<br />
    During the developmen of a new feature you discover a bug which needs to be fixed. In order to not mix bug fix and feature development, the bug fix should happen in a different branch. The recommended procedure for handling this situation is the following:
1. Get into a clean state of your working directory on your feature branche (either by commiting open changes or by stashing them).
1. Checkout the branch the bugfix should be merged into (either _master_ or _release-X.Y.Z_) and get the most recent version.
1. Create a new branch for the bugfix.
1. Fix the bug, write a test, update documentation etc.
1. Open a merge request for the bug fix.
1. Switch back to your feature branch.
1. Merge your local bugfix branch into the feature branch. Continue your feature development.
1. Eventually, the bugfix will be merged into _master_. Then, you can rebase your feature branch on master which will remove all duplicate commits related to the bugfix.    

In git commands this looks like:
1. `git stash`
1. `git checkout master && git pull`
1. `git checkout -b <bug_fix_branch>`
1. Implement the bug fix, add tests, commit.
1. Open a merge request.
1. `git checkout <feature_branch>`
1. `git merge <bug_fix_branch>`
1. Once the merge request for the bug fix is accepted in the upstream repository:
  `git fetch`
  `git rebase origin/master` 

This should give the following git history where the initial feature branch is blue, the bugfix branch is yellow and the feature branch after the rebase is red.
<img src="doc/figures/bugfix_while_on_feature.png" alt="fixing a bug while working on a feature">     
**Move most recent commit(s) to new branch**<br/>
Very enthusiastic about the cool feature you are going to implement, you started from master and made one (or more) commits for your development. That's when you realised that you did not create a new branch and committed directly to the master branch. As you know that you should not directly commit to any "official" branch, you are desperately looking for a solution. Assuming your current situation is: A -> B -> C -> D -> E where HEAD is pointing to E (= master) and the last "official" commit is B as shown below.
<img src="doc/figures/move_to_branch1.png" alt="moving commits to new branch"><br />
You can resolve this situation by running:<br />
<br />

    git checkout <new_branch_name>
    git reset --hard <hash of B>
    git checkout <new_branch_name>
    
<br />
which should give
<img src="doc/figures/move_to_branch2.png" alt="moving commits to new branch">
Now, master is pointing to B, HEAD and &lt;new\_branch\_name&gt; are pointing to E and you can happily continue with your work.