# Bash commands


## Permissions
Give lab mates (same group) read access to a folder and everything inside it:
```bash
chmod -R g+rX /path/to/out
```

- `-R`: apply to all subfolders and files
- `g+r`: group can read
- `X`: group can enter folders (X = only folders, not files)
- Re-run after adding new files, because new files don't inherit it


Check permissions:
```bash
ls -l
```
