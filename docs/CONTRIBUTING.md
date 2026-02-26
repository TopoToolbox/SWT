# Contribution Guidelines

See the TopoToolbox [Contribution Guidelines](https://topotoolbox.github.io/contributing.html) for general information about contributing to TopoToolbox. This document contains more detailed instructions for working with the SWT repository.

Development of SWT uses Git as a version control system and GitHub as a hosting platform for the git repository, [issue tracker](https://github.com/TopoToolbox/SWT/issues) and [pull request](https://github.com/TopoToolbox/SWT/pulls).

If you have any problems contributing to SWT, bug reports or any questions about changes you'd like to make to the code [open an issue](https://github.com/TopoToolbox/SWT/issues/new) and describe your problem and we can discuss it there. It is useful to document both challenges we run into and design decisions.

## Setting up Git and GitHub

To get started, you will need to have an account on GitHub. [Create an account](https://docs.github.com/en/get-started/start-your-journey/creating-an-account-on-github) if you have not already. The [GitHub Docs](https://docs.github.com) are very extensive, and I will refer to them frequently.

Next, you will need to [set up Git](https://docs.github.com/en/get-started/git-basics/set-up-git) on your own computer. There are many different ways to [install Git](https://git-scm.com/install/) depending on your operating system and level of comfort with command line interfaces: 

- If you don't want to use a command line, use [GitHub Desktop](https://github.com/apps/desktop)
- If you use Windows, use Git for Windows.
- If you use macOS, use the Git that comes with Xcode (`xcode-select --install` via the Terminal)
- If you use Linux or other Unix derivatives, install git via your package manager.

From now on, I will assume you are using a command line version of Git. If you are using GitHub Desktop, you should be able to map the ideas onto its interface fairly easily.

Once you have installed Git, you should set your [username](https://docs.github.com/en/get-started/git-basics/setting-your-username-in-git) and [email](https://docs.github.com/en/account-and-profile/how-tos/email-preferences/setting-your-commit-email-address). Git associates your name and email address with every change you make.

```
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

Your email should be added to your GitHub account to associate your commits to the repository with the GitHub account. You can have multiple emails added to your GitHub account, so even if you signed up using a personal email address, you can use a work email for contributing to SWT.

One last step before you get to coding. Every time you try to push code from your computer to GitHub, Git will ask you for your username and password. This will get annoying after a while, so you can [set up an SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) to authenticate automatically to GitHub.

## Forking and cloning the repository

The SWT code lives in a Git repository that is hosted on GitHub. Git tracks the changes to the code, and GitHub stores the repository and provides a nice interface for interacting with it. The two are separate interfaces, though, and they each have their own language for doing various things. To get a copy of the code to start working, you will "fork" the repository on GitHub and then "clone" the repository to your local machine.

The code repository is [TopoToolbox/SWT](https://github.com/TopoToolbox/SWT). This is protected so that the only way changes can be made is through approval by a TopoToolbox maintainer (i.e. Will). If you try to make changes directly to this repository, you will get an error and be prevented from continuing. This means that you cannot screw up the TopoToolbox/SWT repository, no matter how badly you mess up your own set up. 

It also means that to make changes you'll need to work indirectly on your own copy of the code. GitHub calls this a "fork," and creating one is as simple as going to [TopoToolbox/SWT](https://github.com/TopoToolbox/SWT) and clicking on the "Fork" button in the upper right corner. GitHub will ask you some questions on the next page, but the defaults should be fine. Click "Create fork" to copy the SWT repository into your own account on GitHub. You will make your changes to this fork and then use "Pull Requests" to get them merged into the TopoToolbox/SWT repository. 

Now you can technically make changes to the code on GitHub through the web interface, but you won't be able to compile or run the code. To get a copy of the repository on your local machine, you ["clone" the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) using Git. **Clone your fork rather than the TopoToolbox repository.** Open a command line, navigate to the directory you want to put the SWT repository in and run

```
git clone git@github.com:$GH_USERNAME/SWT.git
```

where you replace `$GH_USERNAME` with your GitHub user name. If you don't have the SSH key set up as described above, you can use an HTTPS URL:

```
git clone https://github.com/$GH_USERNAME/SWT.git
```

If you click on the green "<> Code" button on the front page of your fork, you'll get the options for URLs you can use.

`git clone` will download the complete history of the repository onto your machine. It should create a directory called `SWT` in the directory in which you ran `git clone`. Enter that directory and run `make check` to compile the code and run the tests. If you hit any errors, [open an issue](https://github.com/TopoToolbox/SWT/issues/new) to describe it and I will help you figure out what's gone wrong.

## Making changes

You should make small, incremental changes to the code and get them merged individually rather than making a bunch of unrelated changes at once. This makes it easier to make sure you don't mess anything up, and it is easier to review the changes before merging them into the TopoToolbox/SWT repository. 

Git is a very flexible tool, and there are many possible workflows for using it. Here I'll describe a basic one that works for most of the changes we make to TopoToolbox code.

First, create a new "branch" for the change you would like to make. Branches that track different versions of the code. The official current version of the SWT repository lives on a branch called `main`. If you make changes on the `main` branch, you'll easily run into conflicts when other changes have been merged into the `main` branch. Creating a branch for your change doesn't prevent the conflicts but makes them easier to resolve. Create a new branch with

```
git checkout -b <newbranchname>
```

where you replace `<newbranchname>` with the name of your branch. It doesn't really matter what you call it, but a name descriptive of the change you'd like to make is a good idea. This command will create the new branch and "check out" the new branch, so that your local copy is now tracking the new branch. If you want to switch back to the main branch, run `git checkout main`. Note the absence of the `-b` argument, which is what creates the new branch. Run `git checkout <newbranchname>` to get back onto your new branch.

Now use your favorite text editor to make changes to the code or other files in this repository. Exactly what this looks like is up to you. Save the files like you normally would. Run `make check` from the SWT/ directory to run the tests and make sure everything still works. 

Now you are ready to [record your changes](https://git-scm.com/book/en/v2/Git-Basics-Recording-Changes-to-the-Repository) in the Git repository by making a "commit." First run

```
git status
```

to get an overview of the current state of the repository. Files that you have made changes to should be listed here. Committing is a new step process. First, you "add" the changes you would like to record to a staging area, then you "commit" the changes listed in the staging area. You can add and commit changes to more than one file at one time. Assume that I've made some changes to this CONTRIBUTING.md file. I would first run


```
git add docs/CONTRIBUTING.md
```

to stage the file, then

```
git commit -m "Add instructions for committing to CONTRIBUTING.md"
```

The `-m` argument lets you add a "commit message" that should describe the changes you have made. If you run `git commit` without the `-m`, this will open a text editor that let's you write more extensive commit messages. The editor that launches depends on the value of your `EDITOR` environment variable or your [git configuration](https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup#_editor). Set this up if you wish, or use `-m` to add your commit messages for now.

Git will record your commit in the repository. Keep making changes and committing them. 

## Merging your changes

Once you are happy with the changes you've made, run the tests once more to make sure they continue to pass

```
make check
```

Now you are ready to "push" the changes back to your fork. This sends the changes from your local computer to GitHub, but does not merge them into the TopoToolbox/SWT repository yet. You can only push to a repository that you have access to. If you try to push to TopoToolbox/SWT, GitHub will refuse to allow it and Git will fail with an error message. A remote repository is called a "remote". Each remote has a name and a URL. If you cloned your GitHub fork as described above, Git will automatically create a remote with the name `origin` and a URL pointing to your fork. Push the changes on your feature branch using

```
git push origin <newbranchname>
```

where `<newbranchname>` is the name of the branch you created before. Use `git status` if you have forgotten what you called it.

This will create the branch on GitHub if it doesn't already exist and add your commits to the branch. You can keep making changes, committing and pushing after you've pushed a branch to GitHub. There are ways to delete commits that you have pushed, but they are somewhat advanced techniques. If you made a mistake and pushed a commit that you didn't want to, the easiest way to fix it is to make a new commit undoing the change and push it on top. If you really don't want that commit to be public (because it has sensitive information, for example), get in touch with Will who will help you fix it.

Now your fork has your changes. You'll propose that the changes be included in the TopoToolbox/SWT repository by making a "Pull Request" on GitHub. Go to [TopoToolbox/SWT](https://github.com/TopoToolbox/SWT). If you pushed your changes recently, GitHub will display a yellow banner that says your branch had recent pushes and give you the option to "Compare & pull request." Click on that to create a new pull request. If it has been a while since you pushed, you can go to the Pull requests tab of the repository and click "New pull request."

