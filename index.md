---
title: Legolas
layout: splash
last_modified_at: 2020-08-28
header:
  image: /assets/images/logo_legolas_1280_trans.png
intro:
  - excerpt: "Legolas, short for **L**arge **E**igensystem **G**enerator for **O**ne-dimensional p**LAS**mas, is a novel
              finite element code tailored to MHD spectroscopy of 1D Cartesian/cylindrical
              equilibria with flow that balance pressure gradients, enriched with various non-adiabatic effects."
example:
  - image_path: /assets/images/example_suydam.png
    excerpt: Legolas can do a multitude of things, ranging from full spectrum calculations to eigenfunctions
             of specific modes, to full-on parametric studies of various equilibrium configurations in different
             geometries. Take a look at the examples to see what the code is capable of.
    title: "Examples"
    url: /content/getting-started/examples/
    btn_label: "Read More"
    btn_class: "btn--primary"
code_paper:
  - image_path: /assets/images/code_paper_fig.png
    excerpt: We have written a method paper, accepted in the Astrophysical Journal Supplement Series, showcasing
             Legolas to the scientific community. Here we explain the underlying mathematical formalism
             in great detail, including a plethora of examples ranging from p- and g-modes in gravitationally
             stratified atmospheres to modes relevant in coronal loop seismology and stability studies of astrophysical jets.
    title: "The Legolas paper"
    url: "https://arxiv.org/abs/2010.14148"
    btn_label: "Read More"
    btn_class: "btn--primary"
using_legolas:
  - image_path: /assets/images/console_output.png
    excerpt: Legolas is the result of months and months of developing, testing, fixing issues, testing again,
             thinking bugs are fixed, further development, discovering that bugs weren't fixed, headscratching, testing again, etc.
             In short, a typical development process of a brand new code. Since this took (and still takes) a lot of effort and time,
             we therefore kindly ask that _**the first published peer-reviewed paper from applying Legolas is done in
             co-authorship with at least one of the original authors**_.
             Since the code is brand new we would like to know how it is used and provide guidance if possible.
             Additionally, if you use Legolas in a publication, we kindly request that you cite our paper.
    title: "Using Legolas"
    url: "https://ui.adsabs.harvard.edu/abs/2020arXiv201014148C/exportcitation"
    btn_label: "BibTex citation"
    btn_class: "btn--primary"
funding:
  - image_path: /assets/images/prominent_logo.png
    excerpt: Legolas is supported by funding from the European Research Council (ERC) under the European
             Unions Horizon 2020 research and innovation programme,
             Grant agreement No. 833251 PROMINENT ERC-ADG 2018; from the VSC (Flemish Supercomputer Center),
             funded by the Research Foundation – Flanders (FWO) and the Flemish Government – department EWI;
             and from Internal Funds KU Leuven, project C14/19/089 TRACESpace.
    title: "Funding"
    url: "https://cordis.europa.eu/project/id/833251"
    btn_label: "Read more"
    btn_class: "btn--primary"
---

{% include feature_row id="intro" type="center" %}

{% include feature_row id="example" type="left" %}

{% include feature_row id="using_legolas" type="right" %}

{% include feature_row id="code_paper" type="left" %}

{% include feature_row id="funding" type="right" %}
