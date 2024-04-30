We plan to contribute our solver to open-source Tamaas project(https://tamaas.readthedocs.io/en/latest/) by creating merge requests(https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html#new-merge-request-from-a-fork).

We can start to learn how to play with git, especially how to merge branches, we have installed GitLab Workflow extension for VS Code (https://docs.gitlab.com/ee/editor_extensions/visual_studio_code/). 

We also installed blender, we will try to learn this rendering skill as well.

To initiate our work on GitLab, we follow the steps in [3]. We have been added in Tamaas project as developer, so we plan to play with our branch. Here we first clone from **master** branch of Tamaas[1].

Some useful commands are listed:

```bash
git remote set-url origin git@gitlab.com:tamaas/tamaas.git
git fetch --all
git remote --verbose
git checkout -b maxwell_viscoelastic-branch 
git branch -a
```

Before everything started, we should deal with SSH:

```bash
(base)  lizichen@lizichendeMacBook-Air  ~/Tamaas/tamaas   maxwell_viscoelastic-branch  cat ~/.ssh/id_rsa.pub 
```
And add this public key to our GitLab webpage.(Just like what we have done for Github)

Then we will be informed to set a passphrase for this key by typing:

```bash
(base)  lizichen@lizichendeMacBook-Air  ~/Tamaas/tamaas   maxwell_viscoelastic-branch  git fetch 
```

With the help of ssh-agent(integrated on macOS since Leopard, version 10.5 in 2007)[4], we don't need to type our passphrase every time when we push local to remote.

Start the SSH agent and add our private key to the SSH agent:
```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_rsa
```

List all keys loaded into ssh-agent:
```bash
ssh-add -l
```
If everything is set up correctly, we will see a welcome message stating that our SSH connection is now set up by typing:
```bash
ssh -T git@gitlab.com
```

**origin** is the default name of the remote repository. Establish a tracking relationship between the local branch and the remote branch, so that subsequent git push and git pull commands on this branch do not need to specify the remote repository name and branch name. The following command is usually used when we push a new branch to the remote repository for the first time to set up this tracking relationship:

```bash
git remote set-url origin git@gitlab.com:tamaas/tamaas.git
```

When we talk about git, we should be clear with local and remote repo, for GitLab, we can easily create a New branch by clicking **New branch** on GitLab webpage, or use command:

```bash
(base)  lizichen@lizichendeMacBook-Air  ~/Tamaas/tamaas   master  git checkout -b maxwell_viscoelastic-branch               
Switched to a new branch 'maxwell_viscoelastic-branch'
```
If we already link our local repo and remote repo with SSH, then a new branch will be automatic generated when we push our local repo to remote:

```bash
git push -u origin maxwell_viscoelastic-branch
```

##### Reference:

[1] https://tamaas.readthedocs.io/en/latest/

[2] https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html#new-merge-request-from-a-fork

[3] https://yang-xijie.github.io/LECTURE/Git/git/#_10

[4] https://en.wikipedia.org/wiki/Ssh-agent