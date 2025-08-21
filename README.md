# localization_remapping_lipvipmst

## Gitlab steps

1. Create separate local branch (replace "$name" by own name):
```
git checkout -b $name
git push -u origin $name
```

2. To get the latest version of the shared repository: `git pull origin main`. 
While on your own branch, this will actually pull the origin repository to your own branch (confusingly, not to the main branch).

3. After significant changes to your own branch, push it to origin:
```
git add .
git commit -m "Description of changes"
git push origin $name
```

4. To merge own branch with shared repository:
    - Go to GitHub → Pull requests → New pull request
    - Set: Base: main, Compare: $name-features
    - Create pull request
    - Other person reviews and merges
