---
title: Legolas
layout: splash
header:
  image: /assets/logo_legolas_1280_trans.png
intro: 
  - excerpt: "Legolas, short for **L**arge **E**igensystem **G**enerator for **O**ne-dimensional p**LAS**mas, is a novel
              finite element code tailored to MHD spectroscopy of 1D Cartesian/cylindrical 
              equilibria with flow that balance pressure gradients, enriched with various non-adiabatic effects."
suydam_example:
  - image_path: /assets/example_suydam.png
    excerpt: Legolas can do a multitude of things, ranging from full spectrum calculations to eigenfunctions 
             of specific modes, to full-on parametric studies of various equilibrium configurations in different
             geometries. Take a look at the examples to see what the code is capable of.
    title: "Examples"
    url: /content/getting-started/examples/
    btn_label: "Read More"
    btn_class: "btn--primary"
code_paper:
  - image_path: /assets/code_paper_fig.png
    excerpt: We've written a method paper, published in the Astrophysical Journal Supplement Series,
             showcasing Legolas to the scientific community. If you use Legolas for a publication,
             we kindly ask you to cite it.
    title: "Our method paper"
    url: "https://arxiv.org"
    btn_label: "Paper pdf"
    btn_class: "btn--primary"
contact:
  - image_path: /assets/get_in_touch.png
    excerpt: If you run into issues installing/using the code, or encounter bugs, please get in touch with us!
             For bugs or ideas for improvements, log an issue or submit a pull request on the 
             [GitHub repository](https://github.com/n-claes/legolas). Simply chatting with the developers is
             also possible, Legolas has a dedicated Slack channel.
    title: "Get in touch"
    url: "https://join.slack.com/t/the-legolas-code/shared_invite/zt-h8whvb80-28yXwqafp8DCsPTw4OProQ"
    btn_label: "Join Legolas Slack channel"
    btn_class: "btn--primary"
funding:
  - image_path: /assets/prominent_logo.png
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

{% include feature_row id="suydam_example" type="left" %}

{% include feature_row id="code_paper" type="right" %}

{% include feature_row id="contact" type="left" %}

{% include feature_row id="funding" type="right" %}

