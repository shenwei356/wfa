# Changelog

- v0.3.0
    - faster speed and lower memory simply by changing the default value of `OFFSETS_BASE_SIZE` from 2048 to 128.
    - setting `WAVEFRONTS_BASE_SIZE` and `OFFSETS_BASE_SIZE` as global variables, not constants.
- v0.2.3
    - fix adaptive reduction.
- v0.2.2
    - fix a concurrency bug when using multiple aligners.
- v0.2.1
    - fix missing of bounds checking in adaptive reduction.
- v0.2.0
    - rewrite the whole package.
    - more accurate.
    - 2X faster than v0.1.0, lower memory.
- v0.1.0
    - first version.
    - The speed is no ideal as the data structures and backtrace algorithm are inefficient.
