
- When we are contributing project, we may run some unreleased version to test our results, such that we need:

```bash
pip install -e build-release/python
```

Then in our virtual environment, we can test our tamaas brunch.

- Merge or rebase

```bash
git config pull.rebase false
git pull --no-rebase origin maxwell_viscoelastic-branch
```


- Git tips: When we commit a change and we want to go back to the previous version and keep this change at the same time

```bash
git status #see status
git reset HEAD~ #HEAD~n, n can be preious commit
git status #see status
git add src/solvers/maxwell_viscoelastic.cpp #add change to staged zone
git commit -C ORIG_HEAD #we can see we keep the change and ready to push to remote repo
git log #see change info(Author, Date, commit message)
```
