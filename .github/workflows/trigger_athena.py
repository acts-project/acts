name: Trigger Athena build

on:
  push:
    branches:
      - main
jobs:
  trigger_athena_job:
    runs-on: ubuntu-latest
    steps:
