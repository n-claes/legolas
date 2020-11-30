Hi, thanks for the pull request! To make sure everything stays nice and consistent, here is a default template which you can use.
It also includes a checklist of small things that are usually forgotten or overlooked. You may remove parts that are not applicable to the current PR
(including these sentences).

Cheers,

Legolas dev team

## PR description

> fill or delete

## New features

> fill or delete

## Bugfixes

> fill or delete

## Checklist

**Source code stuff**
- [ ] everything is nicely formatted using [Black-like code style](https://black.readthedocs.io/en/stable/the_black_code_style.html)
- [ ] use-statements should use `use mod_xxx, only:` as much as possible
- [ ] things that can go in a separate module, should go in a separate module
- [ ] code additions|changes are also added|changed in the docstrings

**Testing stuff**
- [ ] all tests pass
- [ ] in case of new features, new tests are added

**Website stuff**
- [ ] relevant pages on the website are edited, if relevant
- [ ] on the edited pages, the `last_modified_at` frontmatter key is updated
- [ ] all links render as they should, check through a local `bundle exec jekyll serve` in the `docs` directory
- [ ] additions to the parfile are also added to `parameters_file.md`

**Release stuff**
(only if this is a merge-in-master release)
- [ ] bump version number
- [ ] changelog pages on the website are updated
