# IWYU = Include what you use

This tool finds unused includes and suggestes additional includes and declarations.

It is not very stable at the moment and takes a few hours to complete within ACTS therefor it only runs once a week and can be triggered manually.

There is also not sufficient filtering offered by IWYU at the moment. For that reason there is a specific filtering script for now which tries to get rid of unwanted changes.

offline resources
- [GitHub Actions Workflow](../../.github/workflows/iwyu.yml)
- [Custom filter script](./filter.py)

online resources

- https://include-what-you-use.org/
- https://github.com/include-what-you-use/include-what-you-use
