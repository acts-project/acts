## :bar_chart: Physics performance monitoring for {{ .commit }}
[Full report]({{ .url }}/)
CKF: [seeded]({{ .url }}/ckf_seeded.html), [truth smeared]({{ .url }}/ckf_truth_smeared.html), [truth estimated]({{ .url }}/ckf_truth_estimated.html)
IVF: [seeded]({{ .url }}/ivf_seeded.html), [truth smeared]({{ .url }}/ivf_truth_smeared.html), [truth estimated]({{ .url }}/ivf_truth_estimated.html) 
[Ambiguity resolution]({{ .url }}/ambi_seeded.html)
[Truth tracking]({{ .url }}/truth_tracking.html)

### Vertexing

<img src="{{ .url }}/vertexing_mu_scan.pdf?to_png=1" width="350"/>

<details>
  <summary><b>IVF seeded</b></summary>
  <img src="{{ .url }}/ivf_seeded_plots/covXX.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_seeded_plots/covYY.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ivf_seeded_plots/diffx.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_seeded_plots/diffy.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ivf_seeded_plots/diffz.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_seeded_plots/recoOverTrue.pdf?to_png=1" width="50%"/>
</details>

<details>
  <summary><b>IVF truth smeared</b></summary>
  <img src="{{ .url }}/ivf_truth_smeared_plots/covXX.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_truth_smeared_plots/covYY.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ivf_truth_smeared_plots/diffx.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_truth_smeared_plots/diffy.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ivf_truth_smeared_plots/diffz.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_truth_smeared_plots/recoOverTrue.pdf?to_png=1" width="50%"/>
</details>

<details>
  <summary><b>IVF truth estimated</b></summary>
  <img src="{{ .url }}/ivf_truth_estimated_plots/covXX.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_truth_estimated_plots/covYY.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ivf_truth_estimated_plots/diffx.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_truth_estimated_plots/diffy.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ivf_truth_estimated_plots/diffz.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ivf_truth_estimated_plots/recoOverTrue.pdf?to_png=1" width="50%"/>
</details>

### CKF

<details>
  <summary><b>seeded</b></summary>
  <img src="{{ .url }}/ckf_seeded_plots/trackeff_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ckf_seeded_plots/trackeff_vs_pT.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ckf_seeded_plots/nHoles_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ckf_seeded_plots/nMeasurements_vs_eta.pdf?to_png=1" width="50%"/>
</details>

<details>
  <summary><b>truth smeared</b></summary>
  <img src="{{ .url }}/ckf_truth_smeared_plots/trackeff_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ckf_truth_smeared_plots/trackeff_vs_pT.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ckf_truth_smeared_plots/nHoles_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ckf_truth_smeared_plots/nMeasurements_vs_eta.pdf?to_png=1" width="50%"/>
</details>

<details>
  <summary><b>truth estimated</b></summary>
  <img src="{{ .url }}/ckf_truth_estimated_plots/trackeff_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ckf_truth_estimated_plots/trackeff_vs_pT.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ckf_truth_estimated_plots/nHoles_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ckf_truth_estimated_plots/nMeasurements_vs_eta.pdf?to_png=1" width="50%"/>
</details>

### Ambiguity resolution

<details>
  <summary><b>seeded</b></summary>
  <img src="{{ .url }}/ambi_seeded_plots/trackeff_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ambi_seeded_plots/trackeff_vs_pT.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/ambi_seeded_plots/nHoles_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/ambi_seeded_plots/nMeasurements_vs_eta.pdf?to_png=1" width="50%"/>
</details>

### Truth tracking

<details>
  <summary><b>Truth tracking</b></summary>
  <img src="{{ .url }}/truth_tracking_plots/nHoles_vs_eta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/truth_tracking_plots/nMeasurements_vs_eta.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/truth_tracking_plots/pull_d0.pdf?to_png=1" width="50%"/><img src="{{ .url }}/truth_tracking_plots/pull_z0.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/truth_tracking_plots/pull_theta.pdf?to_png=1" width="50%"/><img src="{{ .url }}/truth_tracking_plots/pull_phi.pdf?to_png=1" width="50%"/>

  <img src="{{ .url }}/truth_tracking_plots/pull_qop.pdf?to_png=1" width="50%"/><img src="{{ .url }}/truth_tracking_plots/pull_t.pdf?to_png=1" width="50%"/>
</details>
