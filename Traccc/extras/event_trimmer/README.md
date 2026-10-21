# Event trimmer

This tool takes an input event and extracts the data for a small number of particles. It is designed to be useful for debugging all kinds of reconstruction issues. For example, the invocation:

```
python main.py -p 1 -p 5 -i 2 input output
```

Will read event 2 (`-i`) from directory `input` and extract particles 1 and 5 (`-p`) from it. It will then write the trimmed event to the `output` directory.
