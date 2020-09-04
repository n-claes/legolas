---
title: Testing the Pylbo framework
layout: single
sidebar:
  nav: "leftcontents"
---

The Pylbo framework is also tested, we use [pytest](https://docs.pytest.org/en/stable/) for this.
To run these locally you'll need pytest installed, which is explained in their
[installation instructions](https://docs.pytest.org/en/stable/getting-started.html#install-pytest).
Once that's done, you navigate to `tests/pylbo_tests` and execute
```bash 
pytest
```
This will automatically search for all unit tests and execute them, and you'll see a detailed overview during runtime.

These tests are run automatically for each commit and pull request to the `master` and `develop` branches.