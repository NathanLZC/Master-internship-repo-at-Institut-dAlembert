
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

